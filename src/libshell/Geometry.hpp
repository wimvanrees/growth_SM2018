//
//  Geometry.hpp
//  Elasticity
//
//  Created by Wim van Rees on 2/24/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef Geometry_h
#define Geometry_h

#include "common.hpp"

#ifdef USETRILIBRARY
#include "igl/triangle/triangulate.h"
#include "Triangulate_TRIANGLE.hpp"
#endif

#include <string>
#include <algorithm>

#include "ReadVTK.hpp"
#include <igl/readOFF.h>
#include <igl/readOBJ.h>

/*! \class Geometry
 * \brief Base class for different basic triangle mesh geometries.
 *
 * @details
 * This class provides one (pure virtual) method that fills the lists of vertices and faces. Every geometry that derives from here represents a different triangle geometry. Boundary conditions for now are not yet taken into account. Specifics for each geometry are initialized in the constructor.
 */

class Geometry
{
protected:
    /**
     * Whether or not to print stuff during running
     */
    mutable bool verbose;
    
    /**
     * A helper method for initializing the vertices_bc array to false
     */
    void initVertexBoundaries(const int nVertices, Eigen::MatrixXb & vertices_bc) const
    {
        vertices_bc.resize(nVertices,3);
        for(int i=0;i<nVertices;++i) vertices_bc.row(i) << false,false,false;
    }
    
public:
    Geometry():
    verbose(true)
    {}
    
    /**
     * The main method of this class: it fills/overwrites the input arrays vertices (Reals, nVertices x 3 entries -- vertex locations in 3D space); face2vertices(integers, nFaces x 3 entries -- each row contains the indices of the 3 vertices that make up that face .. ordering has to be consistent with the face normal);  and vertices_bc(booleans, nVertices x 3 entries -- each row says whether the vertex of that row is fixed or not in each of the 3 cartesian directions)
     */
    virtual void get(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const = 0;

    /**
     * whether or not we can provide an analytic expression of the normal over the entire mid-surface
     */
    virtual bool analyticNormals() const = 0;
    
    /**
     * returns analytic expression of the normal vector (in world-coordinates) at location pos on the mesh - only used in case analyticNormals returns true
     */
    virtual Eigen::Vector3d getNormal(const Eigen::Vector3d pos) const = 0;
    
    virtual ~Geometry()
    {}

    void setVerbose() const{verbose=true;}
    void setQuiet() const{verbose=false;}
};


/*! \class IOGeometry
 * \brief Geometry class that defines the topology and vertex locations from an input file
 *
 * @details
 * This class reads an external file from disk that contains the vertex location in three dimensions, as well as the vertx indices for each face. The external file can be an OBJ or OFF file (read using the libIGL methods), or a VTP file (read using the VTK library). Whichever method is used to read the file is chosen based on the file extension. The analyticNormals method returns false : this way the initial edge directors are aligned with the dihedral angles.
 */
class IOGeometry : public Geometry
{
protected:
    enum class MESHTYPE {OBJ, OFF, VTP};
    
    const std::string filename;
    MESHTYPE meshtype;
    
    std::string getFileExtension(const std::string& filename) const
    {
        if(filename.find_last_of(".") != std::string::npos)
            return filename.substr(filename.find_last_of(".")+1);
        return "";
    }
    
public:
    IOGeometry(const std::string filename):
    filename(filename)
    {
        std::string extension = getFileExtension(filename);
        std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower); // make everything lower case
        if(extension == "obj")
            meshtype = MESHTYPE::OBJ;
        else if(extension == "off")
            meshtype = MESHTYPE::OFF;
        else if(extension == "vtp")
            meshtype = MESHTYPE::VTP;
        else
        {
            std::cout << "invalid extension ( " << extension << " ) : no mesh reader available" << std::endl;
        }
    }
    
    virtual void get(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        // read the vertices and faces into the input arrays
        switch(meshtype)
        {
            case MESHTYPE::OBJ:
                igl::readOBJ(filename, vertices, face2vertices);
                break;
            case MESHTYPE::OFF:
                igl::readOFF(filename, vertices, face2vertices);
                break;
            case MESHTYPE::VTP:
                ReadVTK::read(filename, vertices, face2vertices);
                break;
            default:
                std::cout << "Invalid meshtype" << std::endl;
                break;
        }
                
        // set the boundary conditions
        initVertexBoundaries(vertices.rows(), vertices_bc);
    }
    
    /**
     * we have no analytic expression for the normals in this discrete mesh
     */
    bool analyticNormals() const override
    {
        return false;
    }
    
    /**
     * getNormal is not used - returns zero vector everywhere
     */
    virtual Eigen::Vector3d getNormal(const Eigen::Vector3d ) const override
    {
        return (Eigen::Vector3d() << 0,0,0).finished();
    }
};
                
/*! \class PlateGeometry
 * \brief Specialization of Geometry class for geometries that are flat and aligned with the xy-plane, so that the normal vector of each face points in the positive z direction
 */
class PlateGeometry : public Geometry
{
    virtual void getPlateGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const = 0;
public:
                                  
    virtual void get(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        getPlateGeometry(vertices, face2vertices, vertices_bc);
    }
    
    /**
     * for a plate, the normals are oriented vertically
     */
    virtual Eigen::Vector3d getNormal(const Eigen::Vector3d ) const override
    {
        return (Eigen::Vector3d() << 0,0,1).finished();
    }
    
    bool analyticNormals() const override
    {
        return true;
    }
};

/*! \class ShellGeometry
 * \brief Specialization of Geometry class for geometries that are curved. Sets analyticNormals true so that getNormal needs to be implemented - also provides a helper method
 */

class ShellGeometry : public Geometry
{
protected:
    virtual void getShellGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const  = 0;
    
    int orientNormalsOutwards(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices) const
    {
        const int nFaces = face2vertices.rows();
        int nFlippedNormals = 0;
        
        for(int i=0;i<nFaces;++i)
        {
            const int vIdx_0 = face2vertices(i,0);
            const int vIdx_1 = face2vertices(i,1);
            const int vIdx_2 = face2vertices(i,2);
            
            const Eigen::Vector3d v0 = vertices.row(vIdx_0);
            const Eigen::Vector3d v1 = vertices.row(vIdx_1);
            const Eigen::Vector3d v2 = vertices.row(vIdx_2);
            
            // get target normal
            const Eigen::Vector3d facepos = (v0 + v1 + v2)/3.0;
            const Eigen::Vector3d target_normal = getNormal(facepos);
            
            // compute actual normal
            const Eigen::Vector3d e1 = v1 - v0;
            const Eigen::Vector3d e2 = v2 - v0;
            const Eigen::Vector3d normal = e1.cross(e2);//.normalized();
            const Real r = normal.norm();
            if(r < 1e-9)
            {
                std::cout << "DEGENERATE FACE DETECTED : " << r << "\t" << i << "\t" << vIdx_0 << "\t" << vIdx_1 << "\t" << vIdx_2 << std::endl;
            }
            const Real dotProd = (normal.normalized()).dot(target_normal);
            if(dotProd < 0)
            {
                face2vertices(i,0) = vIdx_2;
                face2vertices(i,2) = vIdx_0;
                nFlippedNormals++;
            }
        }
        
        return nFlippedNormals;
    }
    
public:
    void get(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        getShellGeometry(vertices, face2vertices, vertices_bc);
    }
    
    bool analyticNormals() const override
    {
        return true;
    }
};

class Geometry_Dummy : public PlateGeometry
{
    Eigen::MatrixXd my_vertices;
    Eigen::MatrixXi my_face2vertices;
    
    void getPlateGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        vertices = my_vertices;
        face2vertices = my_face2vertices;
        initVertexBoundaries(vertices.rows(), vertices_bc);
    }
    
public:
    Geometry_Dummy(const Eigen::Ref<const Eigen::MatrixXd> my_vertices_in, const Eigen::Ref<const Eigen::MatrixXi> my_face2vertices_in)
    {
        my_vertices = my_vertices_in;
        my_face2vertices = my_face2vertices_in;
    }
};
                
                
/*! \class TwoTriangles
 * \brief Geometry consisting of two equilateral triangles  (L=1) sharing one edge
 *
 */

class TwoTriangles : public PlateGeometry
{
    /* the picture is this
     
         2
     0   |   3
         1

     triangle 0 -> 0,1,2
     triangle 1 -> 1,3,2
     edge 0 -> v0, v1
     edge 1 -> v0, v2
     edge 2 -> v1, v2
     edge 3 -> v1, v3
     edge 4 -> v2, v3
     */
    void getPlateGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        vertices.resize(4,3);
        face2vertices.resize(2,3);
        
        vertices.row(0) << 1 - 0.5*std::sqrt(3),0.5,0;
        vertices.row(1) << 1,0,0;
        vertices.row(2) << 1,1,0;
        vertices.row(3) << 1 + 0.5*std::sqrt(3),0.5,0;
        
        face2vertices.row(0) << 0,1,2;
        face2vertices.row(1) << 1,3,2;
        
        // set all free
        initVertexBoundaries(4, vertices_bc);
    }
    
public:

};

class RectangularPlate_RightAngle : public PlateGeometry
{
protected:
    const Real halfEdgeX;
    const Real halfEdgeY;
    const Real edgeLength;
    const bool closedX;
    const bool closedY;
    mutable int nX, nY;
    mutable Real dX, dY;

    virtual void getPlateGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        // compute number of vertices : want one vertex half-way so even number of closed, or odd number for open        
        const int nX_tmp = (int)std::ceil(2*halfEdgeX/edgeLength);
        nX = (((nX_tmp % 2 == 0) && closedX) or ((nX_tmp % 2 == 1) && not closedX)) ? nX_tmp : nX_tmp+1;
        dX = closedX ? 2*halfEdgeX / nX : 2*halfEdgeX / (nX-1);
        
        const int nY_tmp = (int)std::ceil(2*halfEdgeY/edgeLength);
        nY = (((nY_tmp % 2 == 0) && closedY) or ((nY_tmp % 2 == 1) && not closedY)) ? nY_tmp : nY_tmp+1;
        dY = closedY ? 2*halfEdgeY / nY : 2*halfEdgeY / (nY-1);
        
        const int nV = nX * nY;
        
        // compute the vertex positions
        vertices.resize(nV,3);
        for(int i=0;i<nX;++i)
            for(int j=0;j<nY;++j)
            {
                const int vIdx = i*nY + j;
                
                vertices(vIdx,0) = -halfEdgeX + i*dX;
                vertices(vIdx,1) = -halfEdgeY + j*dY;
                vertices(vIdx,2) = 0.0;
            }
        
        // compute the faces
        const int nF_x = closedX ? nX : (nX-1);
        const int nF_y = closedY ? nY : (nY-1);
        const int nF = 2 * nF_x * nF_y;
        
        face2vertices.resize(nF,3);
        for(int i=0;i<nF_x;++i)
            for(int j=0;j<nF_y;++j)
            {
                const int faceIdx = 2*nF_y*i + 2*j; // two faces per vertex (up to Nx-1)
                
                // get four indices : bottom left, bottom right, top right, top left : 0, 1, 2, 3 (CCW)
//                const int v0_idx = i*nY + j;
//                const int v1_idx = i*nY + j + 1;
//                const int v2_idx = (i+1)*nY + j + 1;
//                const int v3_idx = (i+1)*nY + j;
                
                // these exceptions are only active if nF_x and/or nF_y correspond to closed or not (if they go out of bounds)
                const int rightCol = (i==nX-1) ? 0 : (i+1)*nY; // instead if (i+1)*nY, shift to i=0
                const int topRow = (j==nY-1) ? 0 : j+1;  // instead of j+1, shift to j=0
                
                const int v0_idx = i*nY + j;
                const int v1_idx = rightCol + j;
                const int v2_idx = rightCol + topRow;
                const int v3_idx = i*nY + topRow;
                
                // construct the bottom face
                face2vertices(faceIdx,0) = v0_idx;
                face2vertices(faceIdx,1) = v1_idx;
                face2vertices(faceIdx,2) = v2_idx;
                
                // construct the top face
                face2vertices(faceIdx+1,0) = v0_idx;
                face2vertices(faceIdx+1,1) = v2_idx;
                face2vertices(faceIdx+1,2) = v3_idx;
            }
        
        initVertexBoundaries(nV, vertices_bc);
    }
    
public:
    RectangularPlate_RightAngle(const Real halfX, const Real halfY, const Real edgeLength, const bool closedX = false, const bool closedY = false):
    halfEdgeX(halfX),
    halfEdgeY(halfY),
    edgeLength(edgeLength),
    closedX(closedX),
    closedY(closedY)
    {}
    
    int getNx() const {return nX;}
    int getNy() const {return nY;}
    Real getDx() const {return dX;}
    Real getDy() const {return dY;}
};

class CylinderTubeCurved : public ShellGeometry
{
protected:
    std::function<Eigen::Vector3d(Real)> midline; // x(t) where x is vector, t between 0 and 1
    std::function<Eigen::Vector2d(Real)> radius_params; // a(t) and b(t) where t between 0 and 1, and a,b in normal/binormal dir respectively
    mutable Eigen::MatrixXd midline_normals;
    const Real edgeLength;
    const bool closedX;
    const bool closedY;

    Real approximateEllipsePerimeter(const Eigen::Vector2d a_and_b) const
    {
        const Real a = a_and_b(0);
        const Real b = a_and_b(1);
        
        const Real h = std::pow(a - b, 2) / std::pow(a + b, 2);
        const Real retval = M_PI * (a + b) * (1.0 + 3*h/(10 + std::sqrt(4.0 - 3*h))); // approximate elliptic integral of 2nd kind Ramanujan
        return retval;
    }
    
    Eigen::Vector3d getMidlineTangent(const Real t) const
    {
        const Real dt = 1e-4;
        Eigen::Vector3d retval;
        
        if(t==0)
            retval = midline(dt) - midline(0);
        else if(t==1)
            retval = midline(1) - midline(1-dt);
        else
            retval = midline(t + dt) - midline(t - dt);

//        std::cout << "\t\t" << t << "\t" << retval << "\t" << retval.normalized() << std::endl;
        return retval.normalized();
    }
    
    Eigen::Vector3d getMidlineNormal(const Real t, const Eigen::Vector3d last_norm) const
    {
        const Real dt = 1e-4;
        Eigen::Vector3d retval;
        if(t==0)
            retval = getMidlineTangent(dt) - getMidlineTangent(0);
        else if(t==1)
            retval = getMidlineTangent(1) - getMidlineTangent(1-dt);
        else
            retval = getMidlineTangent(t + dt) - getMidlineTangent(t - dt);
        
        if(retval.norm() < std::numeric_limits<Real>::epsilon())
        {
            // our curve is planar -- no curvature.
            // pick last normal as base curve and find orthogonal one based on that
            const Eigen::Vector3d tang = getMidlineTangent(t);
            if((tang-last_norm).norm() > std::numeric_limits<Real>::epsilon())
            {
                retval = last_norm - tang.dot(last_norm)*tang;
            }
            else
            {
                // we have a problem: last_norm and tang align
                // pick another last_norm
                Eigen::Vector3d new_last_norm = Eigen::Vector3d::Constant(1);
                if((last_norm - new_last_norm).norm() <  std::numeric_limits<Real>::epsilon())
                    new_last_norm << 1, -1, 0; // this one is orthogonal to 1,1,1
                retval = new_last_norm - tang.dot(new_last_norm)*tang;
            }
        }
        
        return retval.normalized();
    }
    
    virtual void getShellGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        // find the length of the section
        const int Ns = 1e4;
        const Real ds = closedX ? 1.0/Ns : 1.0/(Ns-1);
        Real length = 0.0;
        Eigen::Vector3d lastpos = midline(0);
        Real maxCircumference = approximateEllipsePerimeter(radius_params(0));
        for(int i=1;i<Ns;++i)
        {
            Eigen::Vector3d pos = midline(i*ds);
            length += (pos - lastpos).norm();
            lastpos = pos;
            maxCircumference = std::max(maxCircumference, approximateEllipsePerimeter(radius_params(i)));
        }
        
        const Real width = maxCircumference;
        std::cout << "length, width = " << length << "\t" << width << std::endl;
        
        // get a plate of those propertions
        RectangularPlate_RightAngle plate(0.5*length, 0.5*width, edgeLength, closedX, closedY);
        plate.get(vertices, face2vertices, vertices_bc);
        
        const int nX = plate.getNx();
        const int nY = plate.getNy();
        const Real dX = plate.getDx();
//        const Real dY = plate.getDy();
        
        std::cout << "nx, ny, nrows = " << nX << "\t" << nY << "\t" << vertices.rows() << std::endl;

        // then we fold it over : the x-axis is the length, the y-axis becomes the azimuthal component
        const Real dTheta = closedY ? 2.0*M_PI/nY : 2.0*M_PI/(nY-1);
        
        Eigen::Vector3d last_norm = Eigen::Vector3d::Constant(1);
        for(int i=0;i<nX;++i)
        {
            // compute center position of the curve - x-val of plate is along arclength (when normalized)
            const Real xx = i*dX;
            const Real tt = xx/length;
            
            // get curve parameters
            const Eigen::Vector3d posn = midline(tt);
            const Eigen::Vector3d tang = getMidlineTangent(tt);
            const Eigen::Vector3d norm = getMidlineNormal(tt,last_norm);
            const Eigen::Vector3d binorm = (tang.cross(norm)).normalized();
            const Eigen::Vector2d radii = radius_params(tt);
            const Real ellipse_a = radii(0);
            const Real ellipse_b = radii(1);
            for(int j=0;j<nY;++j)
            {
                const int vIdx = i*nY + j;
                
                // compute angle
                const Real theta = j*dTheta;
                
                const Eigen::Vector3d vertex_pos = posn + ellipse_a*std::cos(theta)*norm + ellipse_b*std::sin(theta)*binorm;
                
//                std::cout << i << "\t" << j << "\t" << vIdx << "\t\t" << xx << "\t" << tt << "\t" << theta << "\t\t" << posn(0) << "\t" << posn(1) << "\t" << posn(2) << "\t" << norm(0) << "\t" << norm(1) << "\t" << norm(2) << "\t" << binorm(0) << "\t" << binorm(1) << "\t" << binorm(2) << "\t\t " << vertex_pos(0) << "\t" << vertex_pos(1) << "\t" << vertex_pos(2) << std::endl;
                vertices(vIdx,0) = vertex_pos(0);
                vertices(vIdx,1) = vertex_pos(1);
                vertices(vIdx,2) = vertex_pos(2);
            }
            
            last_norm = norm;
        }

        
        // compute midline normals for the getNormal method later
        midline_normals = Eigen::MatrixXd(Ns, 3);
        last_norm = Eigen::Vector3d::Constant(1);
        for(int i=0;i<Ns;++i)
        {
            const Real tt = i*ds;
            const Eigen::Vector3d norm = getMidlineNormal(tt,last_norm);
            midline_normals.row(i) = norm;
            last_norm = norm;
        }
        
        
        // flip normals outwards
        //const int nFlippedNormals =
        orientNormalsOutwards(vertices, face2vertices);
        
        std::cout << "done with getShellGeometry" << std::endl;
    }
    
    std::pair<Real, Eigen::Vector3d> getClosestMidlinePointAndNormal(const Eigen::Vector3d pos) const
    {
        // have to find the closest point on the midline for which the normal intersects with this point. should be unique.
        // basically for all points along the midline: create a plane spanned by the normal and binormal to that point
        // project the vector from the midline point to our point (pos) onto that plane
        // find the projection with the smallest actual value -- this is the plane!
        
        const int Ns = (int)midline_normals.rows();
        const Real ds = closedX ? 1.0/Ns : 1.0/(Ns-1);
        //        Eigen::Vector3d lastpos = midline(0);
        //        Real maxCircumference = approximateEllipsePerimeter(radius_params(0));
        
        Real minval = 1e9;
        std::pair<Real, Eigen::Vector3d> minval_params;
        
        for(int i=0;i<Ns;++i)
        {
            const Real tt = i*ds;
            const Eigen::Vector3d posn = midline(tt);
            
            // ok now we have the two vectors
            const Eigen::Vector3d pt_to_midline = pos - posn; // not normalized : we also find min distance
            const Real deviation = pt_to_midline.norm();//std::pow(norm.dot(pt_to_midline),2) + std::pow(binorm.dot(pt_to_midline),2); // actual squared distance
            if(deviation < minval)
            {
                minval = deviation;
                minval_params.first = tt;
                minval_params.second = (Eigen::Vector3d)midline_normals.row(i);
            }
        }
        
        return minval_params;
    }
        
        
public:
    CylinderTubeCurved(std::function<Eigen::Vector3d(Real)> f, const Real rad, const Real el, const bool closedX, const bool closedY):
    ShellGeometry(),
    midline(f),
    radius_params([=](Real){Eigen::Vector2d retval; retval << rad, rad; return retval;}),
    edgeLength(el),
    closedX(closedX),
    closedY(closedY)
    {
    }
    
    CylinderTubeCurved(std::function<Eigen::Vector3d(Real)> f, std::function<Eigen::Vector2d(Real)> r, const Real el, const bool closedX, const bool closedY):
    ShellGeometry(),
    midline(f),
    radius_params(r),
    edgeLength(el),
    closedX(closedX),
    closedY(closedY)
    {
    }
    
    virtual Eigen::Vector3d getNormal(const Eigen::Vector3d pos) const override
    {
        const std::pair<Real, Eigen::Vector3d> minval_params = getClosestMidlinePointAndNormal(pos);
        
        // now we found the closest point on the midline : still have to compute the normal (we assume to be in the same cross-section)
        const Real tt = minval_params.first;
        const Eigen::Vector3d posn = midline(tt);
        const Eigen::Vector3d tang = getMidlineTangent(tt);
        const Eigen::Vector3d norm = minval_params.second;
        const Eigen::Vector3d binorm = (tang.cross(norm)).normalized();
        const Eigen::Vector2d radii = radius_params(tt);
        const Real ellipse_a = radii(0);
        const Real ellipse_b = radii(1);
        
        const Eigen::Vector3d relposn = pos - posn;
        
        const Real nX = relposn.dot(norm)/(ellipse_a*ellipse_a);
        const Real nY = relposn.dot(binorm)/(ellipse_b*ellipse_b);
//      const Eigen::Vector3d pos = posn + ellipse_a*std::cos(theta)*norm + ellipse_b*std::sin(theta)*binorm;
        
        const Eigen::Vector3d n = nX*norm + nY*binorm;
        return n.normalized();
    }
        
    Eigen::Vector3d getTangentOnSurface(const Eigen::Vector3d pos) const
    {
        const std::pair<Real, Eigen::Vector3d> minval_params = getClosestMidlinePointAndNormal(pos);
        const Real tt = minval_params.first;
        const Eigen::Vector3d tang = getMidlineTangent(tt);
        return tang;
    }
};
                
class CylinderTube : public ShellGeometry
{
protected:
    const Real radius;
    const Real length;
    const Real edgeLength;
    const bool fixBottomEdge, fixTopEdge;
    
    /* we triangulate a plane of height x perimeter and roll it up into a cylinder
       triangulation follows simple right-angle triangles
       __________
       | /| /| /|
       |/_|/_|/_|
   ^   | /| /| /|
   |   |/_|/_|/_|
   |
 perimeter
        --> height
       
       in this example we have 4 vertices along height and 3 along perimeter
       per row we have (4-1)*2 faces, so in total (4-1)*2*3 = 9 faces
       for one 'square' of two faces we have four vertices
     
       3--2
       | /|
       |/ |
       0--1
     
       the two faces for this row-column combination are then 0-1-2 and 0-2-3 (see code below)
     
       finally we impose an odd number of points along the cylinder height, and a %4 number of points along the perimeter
       this ensure that we have always a point in the middle of the cylinder, and symmetric points around the perimeter
     */
    
    virtual void getShellGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        const Real perimeter = 2.0*M_PI*radius;
        
        // in z-direction
        const int nZ_tmp = (int)std::ceil(length/edgeLength);
        const int nZ = (nZ_tmp % 2 == 0 ? nZ_tmp+1 : nZ_tmp); // make it an odd number
        const Real dz = length / (nZ-1);
        
        // in angular (hteta) direction
        const int nT_tmp = (int)std::ceil(perimeter/edgeLength);
        const int nT = (nT_tmp % 4 == 0 ? nT_tmp : nT_tmp+4-(nT_tmp%4)); // make it an even number
        const Real dt = perimeter / nT;
        
        // create the vertices
        const int nV = nT * nZ;
        
        // create the basic plane : inner loop in z (height) direction
        vertices.resize(nV,3);
        for(int j=0;j<nT;++j)
            for(int i=0;i<nZ;++i)
            {
                const int vIdx = j*nZ + i;
                
                // compute cylinder coordinates
                const Real theta = j*dt/radius;
                const Real x = radius*std::cos(theta);
                const Real y = radius*std::sin(theta);
                
                vertices(vIdx,0) = x;
                vertices(vIdx,1) = y;
                vertices(vIdx,2) = i*dz;
            }
        
        // compute number of faces
        const int nF = ((nZ-1)*2)*nT; // not Nt-1 because we are going to close the whole thing up
        face2vertices.resize(nF,3);
        for(int j=0;j<nT;++j)
            for(int i=0;i<nZ-1;++i)
            {
                const int faceIdx = 2*(nZ-1)*j + 2*i; // two faces per vertex (up to Nx-1)
                
                // get four indices : bottom left, bottom right, top right, top left : 0, 1, 2, 3 (CCW)
                const int v0_idx = j*nZ + i;
                const int v1_idx = j*nZ + i + 1;
                // wrap around for last point to close the cylinder
                const int v2_idx = (j==nT-1) ? i + 1 : (j+1)*nZ + i + 1;
                const int v3_idx = (j==nT-1) ? i     : (j+1)*nZ + i;
                
                // construct the bottom face
                face2vertices(faceIdx,0) = v0_idx;
                face2vertices(faceIdx,1) = v1_idx;
                face2vertices(faceIdx,2) = v2_idx;
                
                // construct the top face
                face2vertices(faceIdx+1,0) = v0_idx;
                face2vertices(faceIdx+1,1) = v2_idx;
                face2vertices(faceIdx+1,2) = v3_idx;
            }
        
        initVertexBoundaries(nV, vertices_bc);
        
        if(fixBottomEdge or fixTopEdge)
        {
            const Real tol = 1e-12;
            for(int i=0;i<nV;++i)
            {
                const bool isOnBottom = std::abs(vertices(i,2)) < tol;
                const bool isOnTop = std::abs(vertices(i,2) - length) < tol;
                if((fixBottomEdge and isOnBottom) or (fixTopEdge and isOnTop))
                    vertices_bc(i,0) = vertices_bc(i,1) = vertices_bc(i,2) = true;
            }
        }
    }
    
public:
    CylinderTube(const Real R, const Real L, const Real el, const bool fixBottomEdge = false, const bool fixTopEdge = false):
    ShellGeometry(),
    radius(R),
    length(L),
    edgeLength(el),
    fixBottomEdge(fixBottomEdge),
    fixTopEdge(fixTopEdge)
    {}
    
    virtual Eigen::Vector3d getNormal(const Eigen::Vector3d pos) const override
    {
        const Real x = pos(0);
        const Real y = pos(1);
        const Real rad = std::sqrt(x*x + y*y);
        const Real nX = -x/rad; // minus since face orientation is actually inwards
        const Real nY = -y/rad;
        
        Eigen::Vector3d n;
        n << nX, nY, 0.0;
        return n.normalized();
    }
};
        
        
        
class Triangulate : public PlateGeometry
{
    
protected:
    const Eigen::MatrixXd & polygon_vertices; // 2D vertex position array
    const Eigen::MatrixXi & polygon_boundary_edges; // list of (2) vertex ids forming unoriented edges
    const Eigen::MatrixXd polygon_holes;
    const Eigen::MatrixXd polygon_areas;
    const std::string & option_flags;
    
#ifdef USETRILIBRARY
    virtual void getPlateGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        Eigen::MatrixXd vertices_2D;
        TRIANGLE::triangulate(polygon_vertices, polygon_boundary_edges, polygon_holes, polygon_areas, option_flags, vertices_2D, face2vertices);
        //igl::triangle::triangulate(polygon_vertices, polygon_boundary_edges, polygon_holes, option_flags, vertices_2D, face2vertices);
        
        const int nVertices = vertices_2D.rows();
        vertices.resize(nVertices,3);
        for(int i=0;i<nVertices;++i)
        {
            vertices(i,0) = vertices_2D(i,0);
            vertices(i,1) = vertices_2D(i,1);
            vertices(i,2) = 0.0;
        }
        initVertexBoundaries(vertices.rows(), vertices_bc);
    }
#else
    virtual void getPlateGeometry(Eigen::MatrixXd & , Eigen::MatrixXi & , Eigen::MatrixXb & ) const override
    {
        std::cout << "Need USETRILIBRARY compilation flag to use the Triangulate geometry class" << std::endl;
        std::exit(1);
    }
#endif
        
public:
    
    Triangulate(const Eigen::MatrixXd & p_vertices, const Eigen::MatrixXi & p_boundaries, const std::string & flags):
    PlateGeometry(),
    polygon_vertices(p_vertices),
    polygon_boundary_edges(p_boundaries),
    option_flags(flags)
    {
    }

    Triangulate(const Eigen::MatrixXd & p_vertices, const Eigen::MatrixXi & p_boundaries, const Eigen::MatrixXd p_holes, const std::string & flags):
    PlateGeometry(),
    polygon_vertices(p_vertices),
    polygon_boundary_edges(p_boundaries),
    polygon_holes(p_holes),
    option_flags(flags)
    {
    }
    
    Triangulate(const Eigen::MatrixXd & p_vertices, const Eigen::MatrixXi & p_boundaries, const Eigen::MatrixXd p_holes, const Eigen::MatrixXd p_areas, const std::string & flags):
    PlateGeometry(),
    polygon_vertices(p_vertices),
    polygon_boundary_edges(p_boundaries),
    polygon_holes(p_holes),
    polygon_areas(p_areas),
    option_flags(flags)
    {
    }
};


class RectangularPlate : public PlateGeometry
{
protected:
    const Real halfX, halfY, maxRelArea;
    const std::array<bool,2> fixedEdgeX, fixedEdgeY;
    
    virtual void getPlateGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        Eigen::MatrixXd polygon_vertices(4,2);
        Eigen::MatrixXi polygon_edges(4,2);
        polygon_vertices(0,0) = -halfX;
        polygon_vertices(0,1) = -halfY;
        polygon_vertices(1,0) = +halfX;
        polygon_vertices(1,1) = -halfY;
        polygon_vertices(2,0) = +halfX;
        polygon_vertices(2,1) = +halfY;
        polygon_vertices(3,0) = -halfX;
        polygon_vertices(3,1) = +halfY;
        polygon_edges(0,0) = 0;
        polygon_edges(0,1) = 1;
        polygon_edges(1,0) = 1;
        polygon_edges(1,1) = 2;
        polygon_edges(2,0) = 2;
        polygon_edges(2,1) = 3;
        polygon_edges(3,0) = 3;
        polygon_edges(3,1) = 0;
        
        const Real totalArea = halfX*halfY;
        const std::string flags = "q20a"+std::to_string(maxRelArea*totalArea)+(verbose ? "" : "Q");
        Triangulate triangulate(polygon_vertices, polygon_edges, flags);
        triangulate.get(vertices, face2vertices, vertices_bc);
        
        const int nV = vertices.rows();
        for(int i=0;i<nV;++i)
        {
            // check if we are on the boundary
            if(std::abs(vertices(i,0) + halfX) < 1e-6*halfX and fixedEdgeX[0])
            {
                vertices_bc(i,0) = true;
                vertices_bc(i,1) = true;
                vertices_bc(i,2) = true;
            }
            
            if(std::abs(vertices(i,0) - halfX) < 1e-6*halfX and fixedEdgeX[1])
            {
                vertices_bc(i,0) = true;
                vertices_bc(i,1) = true;
                vertices_bc(i,2) = true;
            }
            
            if(std::abs(vertices(i,1) + halfY) < 1e-6*halfY and fixedEdgeY[0])
            {
                vertices_bc(i,0) = true;
                vertices_bc(i,1) = true;
                vertices_bc(i,2) = true;
            }
            
            if(std::abs(vertices(i,1) - halfY) < 1e-6*halfY and fixedEdgeY[1])
            {
                vertices_bc(i,0) = true;
                vertices_bc(i,1) = true;
                vertices_bc(i,2) = true;
            }
        }
    }
    
public:
    RectangularPlate(const Real halfX, const Real halfY, const Real maxRelArea, const std::array<bool,2> fixedEdgeX_, const std::array<bool,2> fixedEdgeY_):
    halfX(halfX),
    halfY(halfY),
    maxRelArea(maxRelArea),
    fixedEdgeX(fixedEdgeX_),
    fixedEdgeY(fixedEdgeY_)
    {}
};
        
        
class CylinderTube_unstructured : public ShellGeometry
{
protected:
    const Real radius;
    const Real length;
    const Real edgeLength;
    const bool fixBottomEdge, fixTopEdge;
    const bool cyl_verbose;
    
    virtual void getShellGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        const Real perimeter = 2.0*M_PI*radius;
        
        // in z-direction
        const int nZ_tmp = (int)std::ceil(length/edgeLength);
        const int nZ = (nZ_tmp % 2 == 0 ? nZ_tmp+1 : nZ_tmp); // make it an odd number
        const Real dz = length / (nZ-1);
        
        // in angular (hteta) direction (where we connect)
        const int nT_tmp = (int)std::ceil(perimeter/edgeLength);
        const int nT = (nT_tmp % 4 == 0 ? nT_tmp : nT_tmp+4-(nT_tmp%4)); // make it an even number
        const Real dt = perimeter / (nT-1);
        
        // create the vertices around the perimeter
        const int nV = 2*(nT-1) + 2*(nZ-1);
        
        Eigen::MatrixXd polygon_vertices(nV,2);
        Eigen::MatrixXi polygon_edges(nV,2);
        
        // bottom edge : height
        for(int i=0;i<nZ-1;++i)
        {
            polygon_vertices(i,0) = i*dz;
            polygon_vertices(i,1) = 0.0;
        }
        // right edge : perimeter
        for(int i=0;i<nT-1;++i)
        {
            polygon_vertices(i + nZ - 1, 0) = length;
            polygon_vertices(i + nZ - 1, 1) = i*dt;
        }
        // top edge : height (backwards)
        for(int i=0;i<nZ-1;++i)
        {
            polygon_vertices(i + nZ + nT - 2,0) = length - i*dz;
            polygon_vertices(i + nZ + nT - 2,1) = perimeter;
        }
        // left edge : perimeter (backwards)
        for(int i=0;i<nT-1;++i)
        {
            polygon_vertices(i + 2*nZ + nT - 3, 0) = 0.0;
            polygon_vertices(i + 2*nZ + nT - 3, 1) = perimeter - i*dt;
        }
        
        if(cyl_verbose) std::cout << polygon_vertices << std::endl;

        // do edges
        for(int i=0;i<nV;++i)
        {
            polygon_edges(i,0) = i;
            polygon_edges(i,1) = (i+1)%nV;
        }
        
        // triangulate
        const std::string flags = "q20a"+std::to_string(std::sqrt(3.0)/4.0*dt*dz)+(verbose ? "" : "Q");;
        
        Eigen::MatrixXd vertices_plane;
        Eigen::MatrixXb vertices_bc_plane;
        
        Triangulate triangulate(polygon_vertices, polygon_edges, flags);
        triangulate.get(vertices_plane, face2vertices, vertices_bc_plane);

        
        // now we need to identify the opposing vertices (they better be equal)
        const int nVertices_plane = vertices_plane.rows();
        std::vector<bool> botverts(nVertices_plane, false), topverts(nVertices_plane, false);
        int cnt_bot = 0;
        int cnt_top = 0;
        for(int i=0;i<nVertices_plane;++i)
        {
            if(std::abs(vertices_plane(i,1)) < std::numeric_limits<Real>::epsilon())
            {
                if(cyl_verbose) std::cout << "B : \t " << vertices_plane(i,0) << std::endl;
                botverts[i] = true;
                cnt_bot++;
            }
            
            if(std::abs(vertices_plane(i,1) - perimeter) < std::numeric_limits<Real>::epsilon())
            {
                if(cyl_verbose) std::cout << "T : \t " << vertices_plane(i,0) << std::endl;
                topverts[i] = true;
                cnt_top++;
            }
        }
        
        if(cnt_bot != cnt_top)
        {
            helpers::catastrophe("could not create unstructured cylinder mesh - botVerts and topVerts do not match",__FILE__,__LINE__);
        }
        
        // now we match all the bottom vertices to the top vertices
        std::vector<int> topverts_replacements(nVertices_plane, -1);
        const Real tol = 1e-6*radius;
        int cnt_topbot = 0;
        for(int i=0;i<nVertices_plane;++i)
        {
            if(not topverts[i]) continue;

            const Real topX = vertices_plane(i,0);
            
            // yay for N^2
            for(int j=0;j<nVertices_plane;++j)
            {
                if(not botverts[j]) continue;
                
                const Real botX = vertices_plane(j,0);
                if(std::abs(topX - botX) < tol)
                {
                    topverts_replacements[i] = j;
                    cnt_topbot++;
                    break;
                }
            }
        }
        
        if(cnt_topbot != cnt_top)
        {
            helpers::catastrophe("could not create unstructured cylinder mesh - x-locations between two mesh sides do not match",__FILE__,__LINE__);
        }
        
        // now we know which vertices we can cut
        const int nVertices = nVertices_plane - cnt_top;
        vertices.resize(nVertices,3);
        std::vector<int> planarverts_replacements(nVertices_plane, -1);
        
        int cnt_verts=0;
        for(int i=0;i<nVertices_plane;++i)
        {
            if(topverts[i]) continue;
            
            vertices(cnt_verts,0) = vertices_plane(i,0);
            vertices(cnt_verts,1) = vertices_plane(i,1);
            vertices(cnt_verts,2) = vertices_plane(i,2);
            planarverts_replacements[i] = cnt_verts;
            
            cnt_verts++;
        }
        
        // adjust faces by replacing all occurences of topverts with the corresponding botverts
        for(int i=0;i<face2vertices.rows();++i)
            for(int d=0;d<3;++d)
            {
                const int vidx = face2vertices(i,d);
                if(topverts[vidx])
                    face2vertices(i,d) = planarverts_replacements[topverts_replacements[vidx]];
                else
                    face2vertices(i,d) = planarverts_replacements[vidx];
            }
        
        // finally we roll up the vertices
        for(int i=0;i<nVertices;++i)
        {
            // z is x : height
            const Real x = vertices(i,0);
            const Real theta = vertices(i,1)/radius;
            
            vertices(i,0) = radius*std::cos(theta);
            vertices(i,1) = radius*std::sin(theta);
            vertices(i,2) = x; // height direction
        }
        
        initVertexBoundaries(nVertices, vertices_bc);
        
        if(fixBottomEdge or fixTopEdge)
        {
            const Real tol = 1e-12;
            for(int i=0;i<nV;++i)
            {
                const bool isOnBottom = std::abs(vertices(i,2)) < tol;
                const bool isOnTop = std::abs(vertices(i,2) - length) < tol;
                if((fixBottomEdge and isOnBottom) or (fixTopEdge and isOnTop))
                    vertices_bc(i,0) = vertices_bc(i,1) = vertices_bc(i,2) = true;
            }
        }
    }
    
public:
    CylinderTube_unstructured(const Real R, const Real L, const Real el, const bool fixBottomEdge = false, const bool fixTopEdge = false, const bool cyl_verbose = false):
    ShellGeometry(),
    radius(R),
    length(L),
    edgeLength(el),
    fixBottomEdge(fixBottomEdge),
    fixTopEdge(fixTopEdge),
    cyl_verbose(cyl_verbose)
    {}
    
    virtual Eigen::Vector3d getNormal(const Eigen::Vector3d pos) const override
    {
        const Real x = pos(0);
        const Real y = pos(1);
        const Real rad = std::sqrt(x*x + y*y);
        const Real nX = -x/rad; // minus since normal is actually pointed inwards
        const Real nY = -y/rad;
        
        Eigen::Vector3d n;
        n << nX, nY, 0.0;
        return n.normalized();
    }
};

class LilyPlate : public PlateGeometry
{
protected:
    const Real edgeLength;
    const bool trafo;    
    const Real uCutOff;
    
    virtual void getPlateGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        // create the parametrized domain (u,v) from the domain ( (-2, 184/25) x (0,1) )
        const auto trafo_uv0_to_uv = [](const Real uu, const Real vv0)
        {
            const Real tanarg = 0.75*std::pow(23./6 - uu/(1 + 0.125*uu), 0.75);
            const Real tanfac = 6.0 *std::pow(std::atan(tanarg), 2.25);
            
            const Real vext_1 = 0.5*(M_PI - tanfac);
            const Real vext_2 = 0.5*(M_PI + tanfac);
            const Real vv = vv0*(vext_2 - vext_1) + vext_1;
            
            return std::make_pair(uu, vv);
        };
        
        // create the reference domain (x,y) (planar shape of lily) from the domain ( (-2, 184/25) x (0,1) )
        const auto trafo_uv0_to_xy = [](const Real uu, const Real vv0)
        {
            const Real xx = (uu + 2.0) / (184./25 + 2.) * 13./2;
            const Real yy_ext = 0.6237263790955635 + xx * (-0.10180994569230434 + (0.3319727194151532 - 0.052416873356509235 * xx) * xx);
            // this is a fit to r E^(-(4/3)+z/2) (a Sin[\[Pi] (d(b-z))/(1+d(b-z))]^c)/.{a->4,b->6.5,c->2,d->.15,r->1/2}
            // or in latex
            // r e^{\frac{z}{2}-\frac{4}{3}} \left(a \sin ^c\left(\frac{\pi  (d (b-z))}{d (b-z)+1}\right)\right)\text{/.}\, \left\{a\to 4,b\to 6.5,c\to 2,d\to 0.15,r\to \frac{1}{2}\right\}
            
            const Real yy = yy_ext * (2 * vv0 - 1.0);
            return std::make_pair(xx,yy);
        };
        
        // discretize the edge for u between -2 and 184/25
        const Real umin = -2;
        const Real umax = 184/25. - uCutOff;
        
        const int Npts_u = (int)std::ceil((umax - umin)/edgeLength);
        const Real du = (umax - umin)/(Npts_u - 1.0); // Npts_u - 1 : in case uCutOff = 0, there will be a singularity and the two edge points of the top and bottom outline curves will overlap
        const int nVertices = 2*Npts_u;
        
        Eigen::MatrixXd polygon_vertices(nVertices, 2);
        
        // bottom border : v = 0
        for(int i=0;i<Npts_u;++i)
        {
            const Real uu = umin + i*du;
            const Real vv = 0.0;

            const auto vx_vy = trafo ? trafo_uv0_to_xy(uu, vv) : trafo_uv0_to_uv(uu, vv);
            polygon_vertices(i,0) = vx_vy.first;
            polygon_vertices(i,1) = vx_vy.second;
        }
        
        // top border : v = 1 and go backwards
        for(int i=0;i<Npts_u;++i)
        {
            const Real uu = umax - i*du;
            const Real vv = 1.0;

            const auto vx_vy = trafo ? trafo_uv0_to_xy(uu, vv) : trafo_uv0_to_uv(uu, vv);
            polygon_vertices(i + Npts_u,0) = vx_vy.first;
            polygon_vertices(i + Npts_u,1) = vx_vy.second;
        }


        // create edges (closed loop)
        Eigen::MatrixXi polygon_edges(nVertices,2); // closed loop
        for(int i=0;i<nVertices;++i)
        {
            polygon_edges(i,0) = i;
            polygon_edges(i,1) = (i+1)%nVertices;
        }
        
        const std::string flags = "q20a"+std::to_string(std::sqrt(3.0)/4.0*edgeLength*edgeLength)+(verbose ? "" : "Q");;
        Triangulate triangulate(polygon_vertices, polygon_edges, flags);
        triangulate.get(vertices, face2vertices, vertices_bc);
        // no boundary conditions
        
        
    }
    
public:
    LilyPlate(const Real edgeLength, const bool trafo, const Real uCutOff = 3./10):
    edgeLength(edgeLength),
    trafo(trafo),
    uCutOff(uCutOff)
    {}
};
        
        

class CurvedRectangularShell : public ShellGeometry
{
protected:
    const Real halfX;
    const Real halfY;
    const Real halfOpeningAngle;
    const Real maxRelArea;
    const std::array<bool,2> fixedEdgeX, fixedEdgeY;
    
    virtual void getShellGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        
        // get flat plate
        RectangularPlate plate(halfX, halfY, maxRelArea, fixedEdgeX, fixedEdgeY);
        plate.get(vertices, face2vertices, vertices_bc);
        
        // deform vertices
        if(std::abs(halfOpeningAngle) < std::numeric_limits<Real>::epsilon()) return;
        
        const int nVertices = vertices.rows();
        const Real radius = halfY / std::abs(halfOpeningAngle);
        const Real sign = halfOpeningAngle > 0 ? 1.0 : -1.0;
        for(int i=0;i<nVertices;++i)
        {
            const Eigen::Vector3d v_old = vertices.row(i);
            const Real newY =  radius*std::sin(v_old(1)/radius);
            const Real newZ = sign * radius * (1.0 - std::cos(v_old(1)/radius)); // 0 at the center
            
            vertices(i,0) = v_old(0);
            vertices(i,1) = newY;
            vertices(i,2) = newZ;
        }
    }
    
public:
    CurvedRectangularShell(const Real halfX, const Real halfY, const Real halfOpeningAngle, const Real maxRelArea, const std::array<bool,2> fixedEdgeX_, const std::array<bool,2> fixedEdgeY_):
    halfX(halfX),
    halfY(halfY),
    halfOpeningAngle(halfOpeningAngle),
    maxRelArea(maxRelArea),
    fixedEdgeX(fixedEdgeX_),
    fixedEdgeY(fixedEdgeY_)
    {}
    
    virtual Eigen::Vector3d getNormal(const Eigen::Vector3d pos) const override
    {
        Eigen::Vector3d n;
        if(std::abs(halfOpeningAngle) < std::numeric_limits<Real>::epsilon())
        {
            n << 0,0,1;
        }
        else
        {
            // each cross-section is a circle. but we need to know the center height
            const Real radius = halfY / std::abs(halfOpeningAngle);
            const Real sign = halfOpeningAngle > 0 ? 1.0 : -1.0;
            
            const Real y = -pos(1);
            const Real z = sign * (radius - pos(2));

            n << 0, y, z;
        }
        return n.normalized();
    }
};
        
class CircularPlate : public PlateGeometry
{
protected:
    const Real radius;
    const int nPointsAlongPerimeter;
    const bool fixedBoundary;
    
    virtual void getPlateGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        const Real dTheta = 2.0*M_PI/nPointsAlongPerimeter;
        const Real edgeLength = dTheta*radius;
        
        Eigen::MatrixXd polygon_vertices(nPointsAlongPerimeter,2);
        Eigen::MatrixXi polygon_edges(nPointsAlongPerimeter,2);
        
        polygon_edges(0,0) = 0;
        polygon_edges(nPointsAlongPerimeter-1,1) = 0;
        
        for(int i=0;i<nPointsAlongPerimeter;++i)
        {
            const Real theta = i*dTheta;
            const Real myX = std::cos(theta)*radius;
            const Real myY = std::sin(theta)*radius;
            polygon_vertices(i,0) = myX;
            polygon_vertices(i,1) = myY;
            if(i<nPointsAlongPerimeter-1)
            {
                polygon_edges(i,1) = i+1;
                polygon_edges(i+1,0) = i+1;
            }
        }
        
        const std::string flags = "q20a"+std::to_string(std::sqrt(3.0)/4.0*edgeLength*edgeLength)+(verbose ? "" : "Q");;
        
        Triangulate triangulate(polygon_vertices, polygon_edges, flags);
        triangulate.get(vertices, face2vertices, vertices_bc);
        
        if(fixedBoundary)
        {
            const int nV = vertices.rows();
            for(int i=0;i<nV;++i)
            {
                const Real radSq = vertices(i,0)*vertices(i,0) + vertices(i,1)*vertices(i,1);
                const Real diff = std::abs(radSq - radius*radius);
                if(diff < 1e-6*radius)
                {
                    vertices_bc(i,0) = true;
                    vertices_bc(i,1) = true;
                    vertices_bc(i,2) = true;
                }
            }
        }
    }
    
public:
    CircularPlate(const Real radius, const int nPointsAlongPerimeter, const bool fixedBoundary):
    PlateGeometry(),
    radius(radius),
    nPointsAlongPerimeter(nPointsAlongPerimeter),
    fixedBoundary(fixedBoundary)
    {}
};
                
class CircularPlateMultires : public PlateGeometry
{
protected:
    const Real outerRadius;
    const Real innerRadius;
    const Real outerRes;
    const Real innerRes;
    
    virtual void getPlateGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        const Real perimeter_inner = 2.0*M_PI*innerRadius;
        const Real perimeter_outer = 2.0*M_PI*outerRadius;
        
        const int Npoints_inner = (int)(perimeter_inner/innerRes);
        const int Npoints_outer = (int)(perimeter_outer/outerRes);
        
        const Real dTheta_inner = 2.0*M_PI/Npoints_inner;
        const Real dTheta_outer = 2.0*M_PI/Npoints_outer;
        
        Eigen::MatrixXd polygon_vertices(Npoints_inner + Npoints_outer,2);
        Eigen::MatrixXi polygon_edges(Npoints_inner + Npoints_outer ,2);
        
        // fill the inner radius
        for(int i=0;i<Npoints_inner;++i)
        {
            const Real theta = i*dTheta_inner;
            const Real myX = std::cos(theta)*innerRadius;
            const Real myY = std::sin(theta)*innerRadius;
            polygon_vertices(i,0) = myX;
            polygon_vertices(i,1) = myY;
        }
        
        // fill the outer radius
        for(int i=0;i<Npoints_outer;++i)
        {
            const int idx = i + Npoints_inner;
            const Real theta = i*dTheta_outer;
            const Real myX = std::cos(theta)*outerRadius;
            const Real myY = std::sin(theta)*outerRadius;
            polygon_vertices(idx,0) = myX;
            polygon_vertices(idx,1) = myY;
        }
        
        // edges for inner radius
        for(int i=0;i<Npoints_inner;++i)
        {
            polygon_edges(i,0) = i;
            polygon_edges(i,1) = (i+1)%Npoints_inner;
        }
        
        // edges for outer radius : go the other way around (not sure if I need to)
        for(int i=0;i<Npoints_outer;++i)
        {
            const int edge_idx = i + Npoints_inner;
            polygon_edges(edge_idx,0) = Npoints_inner + i;
            polygon_edges(edge_idx,1) = Npoints_inner + (i+1)%Npoints_outer;
        }
        
        // area constraints
        Eigen::MatrixXd polygon_areas(2,3);
        // region 1 : interior
        polygon_areas.row(0) << 0.0, 0.0, std::sqrt(3.0)/4.0*innerRes*innerRes;
        // region 2 : outer part
        polygon_areas.row(1) << 0.5*(innerRadius + outerRadius), 0.0, std::sqrt(3.0)/4.0*outerRes*outerRes;

        const std::string flags = "q20a"+std::string(verbose ? "" : "Q");
        
        Eigen::MatrixXd polygon_holes;

        Triangulate triangulate(polygon_vertices, polygon_edges, polygon_holes, polygon_areas, flags);
        triangulate.get(vertices, face2vertices, vertices_bc);
    }
    
public:
    CircularPlateMultires(const Real outerRadius, const Real innerRadius, const Real outerRes, const Real innerRes):
    PlateGeometry(),
    outerRadius(outerRadius),
    innerRadius(innerRadius),
    outerRes(outerRes),
    innerRes(innerRes)
    {}
};
                

        
class AnnulusPlate : public PlateGeometry
{
protected:
    const Real outerRadius;
    const Real innerRadius;
    const Real edgeLength;
    const bool fixedOuterBoundary;
    const bool fixedInnerBoundary;
    
    virtual void getPlateGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        const Real perimeter_inner = 2.0*M_PI*innerRadius;
        const Real perimeter_outer = 2.0*M_PI*outerRadius;
        
        const int Npoints_inner = (int)(perimeter_inner/edgeLength);
        const int Npoints_outer = (int)(perimeter_outer/edgeLength);
        
        const Real dTheta_inner = 2.0*M_PI/Npoints_inner;
        const Real dTheta_outer = 2.0*M_PI/Npoints_outer;
        
        Eigen::MatrixXd polygon_vertices(Npoints_inner + Npoints_outer,2);
        Eigen::MatrixXi polygon_edges(Npoints_inner + Npoints_outer ,2);
        
        // fill the inner radius
        for(int i=0;i<Npoints_inner;++i)
        {
            const Real theta = i*dTheta_inner;
            const Real myX = std::cos(theta)*innerRadius;
            const Real myY = std::sin(theta)*innerRadius;
            polygon_vertices(i,0) = myX;
            polygon_vertices(i,1) = myY;
        }
        
        // fill the outer radius
        for(int i=0;i<Npoints_outer;++i)
        {
            const int idx = i + Npoints_inner;
            const Real theta = i*dTheta_outer;
            const Real myX = std::cos(theta)*outerRadius;
            const Real myY = std::sin(theta)*outerRadius;
            polygon_vertices(idx,0) = myX;
            polygon_vertices(idx,1) = myY;
        }
        
        // edges for inner radius
        for(int i=0;i<Npoints_inner-1;++i)
        {
            polygon_edges(i,0) = i;
            polygon_edges(i,1) = i+1;
        }
        // last edge
        polygon_edges(Npoints_inner-1,0) = Npoints_inner-1;
        polygon_edges(Npoints_inner-1,1) = 0;
        
        // edges for outer radius : go the other way around (not sure if I need to)
        for(int i=0;i<Npoints_outer;++i)
        {
            const int edge_idx = i + Npoints_inner;
            polygon_edges(edge_idx,0) = edge_idx;
            polygon_edges(edge_idx,1) = edge_idx + 1;
        }
        // last connection
        polygon_edges(Npoints_inner + Npoints_outer - 1, 0) = Npoints_inner + Npoints_outer - 1;
        polygon_edges(Npoints_inner + Npoints_outer - 1, 1) = Npoints_inner;
        
        // hole
        Eigen::MatrixXd polygon_holes(1,3);
        polygon_holes.row(0) << 0,0,0; // center
        const std::string flags = "q20a"+std::to_string(std::sqrt(3.0)/4.0*edgeLength*edgeLength)+(verbose ? "" : "Q");;
        
        Triangulate triangulate(polygon_vertices, polygon_edges, polygon_holes, flags);
        triangulate.get(vertices, face2vertices, vertices_bc);
        
        // boundary conditions
        if(fixedInnerBoundary or fixedOuterBoundary)
        {
            const int nVertices = vertices.rows();
            for(int i=0;i<nVertices;++i)
            {
                const Real rSq = std::pow(vertices(i,0),2) + std::pow(vertices(i,1),2);
                
                const bool innerR = std::abs(rSq - innerRadius*innerRadius) < 1e-6;
                const bool outerR = std::abs(rSq - outerRadius*outerRadius) < 1e-6;
                
                if((innerR and fixedInnerBoundary) or (outerR and fixedOuterBoundary))
                {
                    vertices_bc(i,0) = true;
                    vertices_bc(i,1) = true;
                    vertices_bc(i,2) = true;
                }                
            }
        }
        
    }
    
public:
    AnnulusPlate(const Real outerRadius, const Real innerRadius, const Real edgeLength, const bool fixedOuterBoundary = false, const bool fixedInnerBoundary = false):
    PlateGeometry(),
    outerRadius(outerRadius),
    innerRadius(innerRadius),
    edgeLength(edgeLength),
    fixedOuterBoundary(fixedOuterBoundary),
    fixedInnerBoundary(fixedInnerBoundary)
    {}
};

class SlitAnnularPlate : public PlateGeometry
{
protected:
    const Real innerRadius, outerRadius;
    const Real edgeLength;
    
    void getPlateGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        const Real perimeter_inner = 2.0*M_PI*innerRadius;
        const Real perimeter_outer = 2.0*M_PI*outerRadius;
        
        const int Npoints_inner = (int)perimeter_inner/edgeLength;
        const int Npoints_outer = (int)perimeter_outer/edgeLength;
        
        const Real dTheta_inner = 2.0*M_PI/Npoints_inner;
        const Real dTheta_outer = 2.0*M_PI/Npoints_outer;
        
        Eigen::MatrixXd polygon_vertices(Npoints_inner + Npoints_outer + 2,2); // +2 to have the endpoints double
        Eigen::MatrixXi polygon_edges(Npoints_inner + Npoints_outer + 2 ,2); // + 2 to connect the inner and outer radii
        
        // fill the inner radius
        for(int i=0;i<Npoints_inner+1;++i)
        {
            const Real theta = i*dTheta_inner;
            const Real myX = std::cos(theta)*innerRadius;
            const Real myY = std::sin(theta)*innerRadius;
            polygon_vertices(i,0) = myX;
            polygon_vertices(i,1) = myY;
        }
        
        // fill the outer radius
        for(int i=0;i<Npoints_outer+1;++i)
        {
            const int idx = i + Npoints_inner + 1;
            const Real theta = i*dTheta_outer;
            const Real myX = std::cos(theta)*outerRadius;
            const Real myY = std::sin(theta)*outerRadius;
            polygon_vertices(idx,0) = myX;
            polygon_vertices(idx,1) = myY;
        }
        
        // edges for inner radius
        for(int i=0;i<Npoints_inner;++i)
        {
            polygon_edges(i,0) = i;
            polygon_edges(i,1) = i+1;
        }
        // connection!
        polygon_edges(Npoints_inner,0) = Npoints_inner;
        polygon_edges(Npoints_inner,1) = Npoints_inner + Npoints_outer + 1;
        
        // edges for outer radius : go the other way around (not sure if I need to)
        for(int i=0;i<Npoints_outer;++i)
        {
            const int edge_idx = i + Npoints_inner + 1;
            polygon_edges(edge_idx,0) = Npoints_inner + Npoints_outer + 1 - i;
            polygon_edges(edge_idx,1) = Npoints_inner + Npoints_outer - i;
        }
        // last connection
        polygon_edges(Npoints_inner + Npoints_outer + 1, 0) = Npoints_inner + 1;
        polygon_edges(Npoints_inner + Npoints_outer + 1, 1) = 0;
        
        // hole
        Eigen::MatrixXd polygon_holes(1,3);
        polygon_holes.row(0) << 0,0,0; // center
        const std::string flags = "q20a"+std::to_string(std::sqrt(3.0)/4.0*edgeLength*edgeLength)+(verbose ? "" : "Q");;
        
        Triangulate triangulate(polygon_vertices, polygon_edges, polygon_holes, flags);
        triangulate.get(vertices, face2vertices, vertices_bc);
    }
    
public:
    SlitAnnularPlate(const Real innerR, const Real outerR, const Real edgeLength):
    innerRadius(innerR),
    outerRadius(outerR),
    edgeLength(edgeLength)
    {
    }
};


class SphericalShell_RightAngles : public ShellGeometry
{
protected:
    const Real sphereRadius;
    const Real height; // should be smaller than radius
    const Real holeHeight; // should be smaller than height
    const Real edgeLength;
    
    void getShellGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        const Real perimeter = 2.0*M_PI*sphereRadius;
        const int nPhi_tmp = (int)std::ceil(perimeter/edgeLength);
        const int nPhi = (nPhi_tmp % 4 == 0 ? nPhi_tmp : nPhi_tmp+4-(nPhi_tmp%4));
        
        // theta is 0 means north pole, theta is pi/2 means equator, theta is pi means south pole
        const Real theta_start = std::acos((sphereRadius - holeHeight)/sphereRadius);
        const Real theta_end = std::acos((sphereRadius - height)/sphereRadius);
        
        const int nTheta = (int)std::ceil(nPhi * sphereRadius * (theta_end - theta_start) / perimeter);
        
        const Real dPhi = 2.0*M_PI/nPhi; // periodic
        const Real dTheta = (theta_end - theta_start)/(nTheta-1); // not periodic
        
        // create the vertices
        const int nVertices = nPhi * nTheta;
        vertices.resize(nVertices,3);
        
        for(int i=0;i<nTheta;++i)
            for(int j=0;j<nPhi;++j)
            {
                const Real theta = theta_end - i*dTheta;// go from bottom to top (equator to north pole)theta_start + i*dTheta;
                const Real phi = j*dPhi;
                
                const int vIdx = i*nPhi + j;
                vertices(vIdx,0) = sphereRadius*std::sin(theta)*std::cos(phi);
                vertices(vIdx,1) = sphereRadius*std::sin(theta)*std::sin(phi);
                vertices(vIdx,2) = sphereRadius*std::cos(theta);
            }

        // create the faces        
        const int nFaces = (nTheta-1)*nPhi*2;
        face2vertices.resize(nFaces,3);
        for(int i=0;i<nTheta-1;++i)
            for(int j=0;j<nPhi;++j)
            {
                // number the vertices : bottom left (0), bottom right (1), top right (2), top left (3) --> CCW
                
                // periodic in phi-direction (in j)
                const int nextJ = (j==nPhi-1) ? 0 : j+1;
                
                const int vIdx_0 = i*nPhi + j;
                const int vIdx_1 = i*nPhi + nextJ;
                const int vIdx_2 = (i+1)*nPhi + nextJ;
                const int vIdx_3 = (i+1)*nPhi + j;
                
                const int faceIdx_0 = 2*vIdx_0;
                const int faceIdx_1 = faceIdx_0 + 1;
                
                face2vertices(faceIdx_0,0) = vIdx_0;
                face2vertices(faceIdx_0,1) = vIdx_1;
                face2vertices(faceIdx_0,2) = vIdx_2;
                
                face2vertices(faceIdx_1,0) = vIdx_0;
                face2vertices(faceIdx_1,1) = vIdx_2;
                face2vertices(faceIdx_1,2) = vIdx_3;
            }
        
        // boundary conditions
        initVertexBoundaries(nVertices, vertices_bc);
    }
    
public:
    SphericalShell_RightAngles(const Real sphereRadius, const Real height, const Real holeHeight, const Real edgeLength):
    sphereRadius(sphereRadius),
    height(height),
    holeHeight(holeHeight),
    edgeLength(edgeLength)
    {
    }
    
    virtual Eigen::Vector3d getNormal(const Eigen::Vector3d pos) const override
    {
        Eigen::Vector3d n;
        
        // origin is at (x,y,z) = (0,0,0) so the normal is just the position vector
        n(0) = pos(0);
        n(1) = pos(1);
        n(2) = pos(2);
        
        return n.normalized();
    }
};
        
        
class SphericalShell_Lambert : public ShellGeometry
{
protected:
    const Real sphereRadius;
    const Real height; // should be smaller than radius
    const Real holeHeight; // should be smaller than height
    const Real edgeLength;
    const bool fixedPerimeter;
    
    void getShellGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        // sphere radius is always one (rescaling later)
        
        // compute the radius of the disk that will be mapped onto the sphere (should have the same area)
        const Real diskRadius = std::sqrt(2.0 * 1.0 * height / sphereRadius); // surface area of cap is 2 pi R h
        
        // create disk
        const bool withHole = holeHeight > std::numeric_limits<Real>::epsilon();
        Geometry * geometry = nullptr;
        if(withHole)
        {
            const Real innerRadius = std::sqrt(2.0 * 1.0 * holeHeight / sphereRadius);
            geometry = new AnnulusPlate (diskRadius, innerRadius, edgeLength / sphereRadius, fixedPerimeter, false);
        }
        else
        {
            // compute the points along the perimeter of the final cap
            const Real cap_radius = sphereRadius*std::sqrt(1.0 - std::min(1.0, std::pow((sphereRadius-height)/sphereRadius,2)));
            const int nPointsAlongPerimeter = std::ceil(2.0*M_PI*cap_radius / edgeLength);
            geometry = new CircularPlate(diskRadius, nPointsAlongPerimeter, fixedPerimeter);
        }
        if(not verbose) geometry->setQuiet();
        
        geometry->get(vertices, face2vertices, vertices_bc);
        delete geometry;
        
        // apply the projection and rescale
        const int nVertices = vertices.rows();
        for(int i=0;i<nVertices;++i)
        {
            Eigen::Vector3d v_old = vertices.row(i);
            const Real rSq = v_old(0)*v_old(0) + v_old(1)*v_old(1);
            
            const Real newX = std::sqrt(1 - rSq/4)*v_old(0);
            const Real newY = std::sqrt(1 - rSq/4)*v_old(1);
            const Real newZ = 0.5*rSq - 1;
            
            vertices(i,0) = newX*sphereRadius;
            vertices(i,1) = newY*sphereRadius;
            vertices(i,2) = -newZ*sphereRadius;
        }
    }
    
public:
    SphericalShell_Lambert(const Real sphereRadius, const Real height, const Real holeHeight, const Real edgeLength, const bool fixedPerimeter):
    sphereRadius(sphereRadius),
    height(height),
    holeHeight(holeHeight),
    edgeLength(edgeLength),
    fixedPerimeter(fixedPerimeter)
    {
    }
    
    virtual Eigen::Vector3d getNormal(const Eigen::Vector3d pos) const override
    {
        Eigen::Vector3d n;
        
        // move the position vector so that z=0 coincides with the sphere origin
        n(0) = pos(0);
        n(1) = pos(1);
        n(2) = pos(2);
        
        return n.normalized();
    }
};
        
        
class SphericalShell : public ShellGeometry
{
protected:
    const Real curvRadius;
    const Real halfOpeningAngle;
    const int nPointsAlongPerimeter;
    const bool fixedPerimeter;
    
    void getShellGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        // create a disk
        const Real diskRadius = curvRadius*std::sin(halfOpeningAngle);
        CircularPlate geometry(diskRadius, nPointsAlongPerimeter, fixedPerimeter);
        geometry.get(vertices, face2vertices, vertices_bc);
        
        // deform into spherical cap (the boundary/perimeter of the cap will be at z=0)
        const Real plane_height = curvRadius * std::cos(halfOpeningAngle); // distance from plane to origin
        
        const int nVertices = vertices.rows();
        for(int i=0;i<nVertices;++i)
        {
            Eigen::Vector3d v_old = vertices.row(i);
            
            const Real inplane_radius = std::sqrt(v_old(0)*v_old(0) + v_old(1)*v_old(1));
            const Real distance_to_center = std::sqrt( std::pow(inplane_radius,2) + std::pow(plane_height,2) );
            
            const Real deltaR = curvRadius  - distance_to_center;
            
            const Real x_over_z = v_old(0) / plane_height;
            const Real y_over_z = v_old(1) / plane_height;
            const Real deltaX = deltaR * x_over_z / (1.0 + std::pow(x_over_z,2)); // sin(atan(x/(Rcosalpha)))
            const Real deltaY = deltaR * y_over_z / (1.0 + std::pow(y_over_z,2)); // sin(atan(y/(Rcosalpha)))
            const Real deltaZ = std::sqrt(deltaR*deltaR - deltaX*deltaX - deltaY*deltaY);
            
            vertices(i,0) += deltaX;
            vertices(i,1) += deltaY;
            vertices(i,2) += deltaZ;

        }
    }
    
public:
    SphericalShell(const Real curvRadius, const Real halfOpeningAngle, const int nPointsAlongPerimeter, const bool fixedPerimeter):
    curvRadius(curvRadius),
    halfOpeningAngle(halfOpeningAngle),
    nPointsAlongPerimeter(nPointsAlongPerimeter),
    fixedPerimeter(fixedPerimeter)
    {
    }
    
    virtual Eigen::Vector3d getNormal(const Eigen::Vector3d pos) const override
    {
        Eigen::Vector3d n;
        
        // move the position vector so that z=0 coincides with the sphere origin
        const Real plane_height = curvRadius*std::cos(halfOpeningAngle);
        n(0) = pos(0);
        n(1) = pos(1);
        n(2) = pos(2) + plane_height;
        
        return n.normalized();
    }
};


class LoadFrom2DData : public PlateGeometry
{
protected:
    const std::string filename;
    const Real zval;
    
    void getPlateGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        std::ifstream ifs(filename);
        if(ifs.fail())
        {
            std::cout << "COULD NOT OPEN FILE " << filename << std::endl;
            std::cout << "exiting" << std::endl;
            std::exit(-1);
        }
        // read number of vertices, number of faces
        int nV, nF;
        ifs >> nV >> nF;
        
        // allocate data
        vertices.resize(nV, 3);
        face2vertices.resize(nF, 3);
        
        // read vertices
        for(int i=0;i<nV;++i)
        {
            ifs >> vertices(i,0) >> vertices(i,1);
            vertices(i,2) = zval;
        }
        
        // read face2vertices
        
        for(int i=0;i<nF;++i)
        {
            ifs >> face2vertices(i,0) >> face2vertices(i,1) >> face2vertices(i,2);
        }
        ifs.close();
        
        // boundary conditions? set by hand
        initVertexBoundaries(nV, vertices_bc);
    }

    
public:
    LoadFrom2DData(const std::string fname, const Real zval = 0.0):
    PlateGeometry(),
    filename(fname),
    zval(zval)
    {
    }
};
        

#endif /* Geometry_h */
