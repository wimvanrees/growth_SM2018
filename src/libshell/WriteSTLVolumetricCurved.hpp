//
//  WriteSTLVolumetricCurved.hpp
//  Elasticity
//
//  Created by Wim van Rees on 8/8/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef WriteSTLVolumetricCurved_h
#define WriteSTLVolumetricCurved_h

#include "common.hpp"
#include <igl/boundary_facets.h>
#include <igl/unique_edge_map.h>
#include <igl/writeSTL.h>

class WriteSTLVolumetricCurved
{
public:

    static void write(const TopologyData & topo, const Eigen::Ref<const Eigen::MatrixXd> vertices, const std::vector<ExtendedTriangleInfo> & vInfos, const Eigen::VectorXd & thickness, const std::string filename)
    {
        write(topo, vertices, vInfos, thickness, thickness, filename);
    }
        
    static void write(const TopologyData & topo, const Eigen::Ref<const Eigen::MatrixXd> vertices, const std::vector<ExtendedTriangleInfo> & vInfos, const Eigen::VectorXd & thickness_bot, const Eigen::VectorXd & thickness_top, const std::string filename)
    {
        // face2vertices gives the vertex indices for each face for the entire mesh (nFaces x 3)
        // vertices gives the vertex locations for the entire mesh (nVertices x 3)
        // thickness gives the thickness for each face
        
        
        std::cout << " === Writing volumetric STL file \t " << filename << std::endl;
        const int nVertices = vertices.rows();
        const int nFaces = topo.getNumberOfFaces();
        
        const auto face2vertices = topo.getFace2Vertices();
        const auto vertex2faces = topo.getVertex2Faces();
        
        // create the data structures for the 3D STL file : a face2vertices instance that contains all the faces, and a vertices instance that contains all the vertices
        // first we create all the vertices
        const int nVertices_STL = 2*nVertices;
        Eigen::MatrixXd vertices_STL(nVertices_STL,3);// twice the vertices because we have two layers
        for(int i=0;i<nVertices;++i)
        {
            // compute average thickness and average normal
            Real totalWeightedThickness_bot = 0.0;
            Real totalWeightedThickness_top = 0.0;
            Real totalArea = 0.0;
            Eigen::Vector3d totalWeightedNormal;
            totalWeightedNormal.setZero();
            for(size_t j=0;j<vertex2faces[i].size();++j)
            {
                const int faceidx = vertex2faces[i][j];
                const ExtendedTriangleInfo & info = vInfos[faceidx];
                
                const Real facearea = 0.5 * info.double_face_area;
                
                totalArea += facearea;
                totalWeightedThickness_bot += thickness_bot(faceidx)*facearea;
                totalWeightedThickness_top += thickness_top(faceidx)*facearea;
                totalWeightedNormal += info.face_normal;
            }
            const Real my_thickness_bot = totalWeightedThickness_bot / totalArea;
            const Real my_thickness_top = totalWeightedThickness_top / totalArea;
            
            totalWeightedNormal /= totalArea;
            const Eigen::Vector3d my_normal = totalWeightedNormal.normalized();
            
            // layer one
            for(int d=0;d<3;++d)
                vertices_STL(i,d) = vertices(i,d) - 0.5 * my_thickness_bot * my_normal(d);
            // layer two
            for(int d=0;d<3;++d)
                vertices_STL(i+nVertices,d) = vertices(i,d) + 0.5 * my_thickness_top * my_normal(d);
        }
        
        // get the boundaries of my structure
        Eigen::MatrixXi boundary_edges;
        {
            Eigen::MatrixXi face2vertices_copy(nFaces, 3);
            for(int i=0;i<nFaces;++i)
                for(int d=0;d<3;++d)
                    face2vertices_copy(i,d) = face2vertices(i,d);
            igl::boundary_facets(face2vertices_copy, boundary_edges);
        }
        const int nBoundaries = boundary_edges.rows();
        
        // get the total number of faces
        const int nFaces_STL = 2*nFaces + 2*nBoundaries;
        Eigen::MatrixXi face2vertices_STL(nFaces_STL, 3);

        // fill the interior faces
        for(int i=0;i<nFaces;++i)
        {
            for(int d=0;d<3;++d)
            {
                face2vertices_STL(i,d) = face2vertices(i,d);
                face2vertices_STL(i + nFaces, d) = face2vertices(i,d) + nVertices;
            }
        }
        
        // fill the boundary faces
        for(int i=0;i<nBoundaries;++i)
        {
            const int vidx_0 = boundary_edges(i,0);
            const int vidx_1 = boundary_edges(i,1);
            
            face2vertices_STL(2*nFaces + 2*i + 0, 0) = vidx_0;
            face2vertices_STL(2*nFaces + 2*i + 0, 1) = vidx_1;
            face2vertices_STL(2*nFaces + 2*i + 0, 2) = vidx_1 + nVertices;

            face2vertices_STL(2*nFaces + 2*i + 1, 0) = vidx_0;
            face2vertices_STL(2*nFaces + 2*i + 1, 1) = vidx_1 + nVertices;
            face2vertices_STL(2*nFaces + 2*i + 1, 2) = vidx_0 + nVertices;
        }
        
        
        // fill the normals array : not sure why I need this or why this should be correct - does not seem to matter
        Eigen::MatrixXd facenormals_STL(nFaces_STL,3);
        facenormals_STL.setZero();
        for(int i=0;i<nFaces_STL;++i)
            facenormals_STL(i,2) = 1.0;
        
        // write
        igl::writeSTL(filename+".stl", vertices_STL, face2vertices_STL, facenormals_STL, false); // no ascii (last argument)
    }
};


#endif /* WriteSTLVolumetricCurved_h */
