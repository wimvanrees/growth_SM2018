//
//  DCSConfigurations.hpp
//  Elasticity
//
//  Created by Wim van Rees on 9/24/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef DCSConfigurations_hpp
#define DCSConfigurations_hpp

//#include <igl/per_face_normals.h>
//#include <igl/edge_lengths.h>
//#include <igl/doublearea.h>

#include "common.hpp"
#include "TopologyData.hpp"
#include "TriangleInfo.hpp"
#include "ConfigurationData.hpp"
#include "DCSCurrentConfigurationData.hpp"

struct DCSBasicConfiguration
{
    // the raw data
    DCSConfigurationData dcsdata;
    
    // default constructor
    DCSBasicConfiguration()
    {}
    
    // assignment operator
    DCSBasicConfiguration & operator= (const DCSBasicConfiguration & other)
    {
        if (this != &other) // protect against invalid self-assignment
            dcsdata = other.dcsdata;
        // by convention, always return *this
        return *this;
    }
    
    // copy constructor
    DCSBasicConfiguration(const DCSBasicConfiguration &other):
    dcsdata(other.dcsdata)
    {
    }
    
    friend void swap(DCSBasicConfiguration& first, DCSBasicConfiguration& second)
    {
        // http://stackoverflow.com/questions/3279543/what-is-the-copy-and-swap-idiom
        using std::swap;
        
        swap(first.dcsdata, second.dcsdata);
    }
    
    virtual void clear()
    {
        dcsdata.clear();
    }
    
    Eigen::MatrixXd computeFaceNormals(const TopologyData & topology) const
    {
        const int nFaces = topology.getNumberOfFaces();
        Eigen::MatrixXd faceNormals(nFaces, 3);
        for(int i=0;i<nFaces;++i)
        {
            const Eigen::Vector3d v0 = dcsdata.vertexdata.row(topology.face2vertices(i,0));
            const Eigen::Vector3d v1 = dcsdata.vertexdata.row(topology.face2vertices(i,1));
            const Eigen::Vector3d v2 = dcsdata.vertexdata.row(topology.face2vertices(i,2));
            
            Eigen::Vector3d face_normal = (v0 - v2).cross(v1 - v0);
            const Real r = face_normal.norm();
            if(r==0)
                face_normal.setZero();
            else
                face_normal /= r;
            faceNormals.row(i) = face_normal;
        }
        return faceNormals;
    }
    
    void setEdgeDirectors(const TopologyData & topology, std::function<Eigen::Vector3d(const Eigen::Vector3d)> normal_function)
    {
        const Eigen::MatrixXd faceNormals = computeFaceNormals(topology);
        setEdgeDirectors(topology, faceNormals, normal_function);
    }
    
    void setEdgeDirectors_fromarray(const TopologyData & topology, const Eigen::Ref<const Eigen::MatrixXd> &normal_vectors)
    {
        const Eigen::MatrixXd faceNormals = computeFaceNormals(topology);
        setEdgeDirectors_fromarray(topology, faceNormals, normal_vectors);
    }
    
    void setEdgeDirectors(const TopologyData & topology, const Eigen::Ref<const Eigen::MatrixXd> &faceNormals, std::function<Eigen::Vector3d(const Eigen::Vector3d)> normal_function)
    {
        const int nEdges = topology.getNumberOfEdges();
        Eigen::MatrixXd normal_vectors(nEdges,3);
        for(int i=0;i<nEdges;++i)
        {
            const int idx_v0 = topology.edge2vertices(i,0);
            const int idx_v1 = topology.edge2vertices(i,1);
            
            const Eigen::Vector3d v0 = dcsdata.vertexdata.row(idx_v0);
            const Eigen::Vector3d v1 = dcsdata.vertexdata.row(idx_v1);
            
            const Eigen::Vector3d midedge_pos = 0.5*(v0 + v1);
            const Eigen::Vector3d midedge_normal = normal_function(midedge_pos);
            for(int d=0;d<3;++d)
                normal_vectors(i,d) = midedge_normal(d);
        }
        
        setEdgeDirectors_fromarray(topology, faceNormals, normal_vectors);
    }
    
    void setEdgeDirectors_fromarray(const TopologyData & topology, const Eigen::Ref<const Eigen::MatrixXd> &faceNormals, const Eigen::Ref<const Eigen::MatrixXd> &normal_vectors)
    {
        const int Nedges = topology.getNumberOfEdges();
        
        for(int i=0;i<Nedges;++i)
        {
            const int faceidx_0 = topology.edge2faces(i,0);
            const int faceidx_1 = topology.edge2faces(i,1);
            
            const bool isBoundary = (faceidx_0 < 0 or faceidx_1 < 0);
            
            // faceidx_0 has the edge vertices contained in order
            const int idx_v0 = topology.edge2vertices(i,0);
            const int idx_v1 = topology.edge2vertices(i,1);
            
            if(not isBoundary)
                assert( (idx_v0 == topology.face2vertices(faceidx_0,0) && idx_v1 == topology.face2vertices(faceidx_0,1)) or (idx_v0 == topology.face2vertices(faceidx_0,1) && idx_v1 == topology.face2vertices(faceidx_0,2)) or (idx_v0 == topology.face2vertices(faceidx_0,2) && idx_v1 == topology.face2vertices(faceidx_0,0)));
            
            
            // these are the counter clockwise vertices of the faceidx_0 plane (the one to the right of the edge if the edge is pointing out of the plane)
            const Eigen::Vector3d v0 = dcsdata.vertexdata.row(idx_v0);
            const Eigen::Vector3d v1 = dcsdata.vertexdata.row(idx_v1);
            const Eigen::Vector3d edge = v1 - v0; // this is the CCW edge viewed from the faceidx_0 triangle
            
            const Eigen::Vector3d midedge_normal_raw = normal_vectors.row(i);
            
            if(not isBoundary)
            {
                const Eigen::Vector3d face_normal_0 = faceNormals.row(faceidx_0);
                const Eigen::Vector3d face_normal_1 = faceNormals.row(faceidx_1);
                
                const Eigen::Vector3d edge_normal = (edge.cross(face_normal_0)).normalized();
                
                // project it onto the normal edge plane
                const Eigen::Vector3d midedge_normal = (midedge_normal_raw.dot(face_normal_0)*face_normal_0 + midedge_normal_raw.dot(edge_normal)*edge_normal).normalized();
                
                // get edge normal (average between face normals)
                const Eigen::Vector3d avg_edge_normal = (face_normal_0 + face_normal_1).normalized();
                
                assert(helpers::isclose(edge_normal.norm(), 1.0));
                assert(helpers::isclose(midedge_normal.norm(), 1.0));
                assert(helpers::isclose(avg_edge_normal.norm(), 1.0));
                
                // get angle between the two (should both be normalized)
                // put it between min/max to prevent nan for very small deviations from 1
                const Real cosphi = std::min(1.0, std::max(-1.0, avg_edge_normal.dot(midedge_normal)));
                
                // now we just need to figure out the sign of this (is the average face normal to the right (+1) or to the left (-1) of the mid edge normal if the edge vector points out of the plane).
                const Real sign = (avg_edge_normal.cross(midedge_normal)).dot(edge) > 0.0 ? 1.0 : -1.0;
                dcsdata.edgedata(i) = sign*std::acos(cosphi);
            }
            else
            {
                const Eigen::Vector3d face_normal = faceidx_0 >= 0 ? faceNormals.row(faceidx_0) : faceNormals.row(faceidx_1);
                const Eigen::Vector3d edge_normal = (edge.cross(face_normal)).normalized();
                
                const Eigen::Vector3d midedge_normal = (midedge_normal_raw.dot(face_normal)*face_normal + midedge_normal_raw.dot(edge_normal)*edge_normal).normalized();
                
                assert(helpers::isclose(edge_normal.norm(), 1.0));
                assert(helpers::isclose(midedge_normal.norm(), 1.0));
                
                // get angle between face normal and midedge normal
                // put it between min/max to prevent nan for very small deviations from 1
                const Real cosphi = std::min(1.0, std::max(-1.0, face_normal.dot(midedge_normal)));
                
                const Real sign = (face_normal.cross(midedge_normal)).dot(edge) > 0.0 ? 1.0 : -1.0;
                dcsdata.edgedata(i) = sign*std::acos(cosphi);
            }
        }
    }
    
    void changeVertices(std::function<Eigen::Vector3d(Eigen::Vector3d)> update_func)
    {
        auto currentvertices = getVertices();
        const int nVertices = getNumberOfVertices();
        for(int i=0;i<nVertices;++i)
        {
            const Eigen::Vector3d v_new = update_func(currentvertices.row(i));
            for(int j=0;j<3;++j)
                currentvertices(i,j) = v_new(j);
        }
    }
    
    void changeVertices(std::function<Eigen::Vector3d(Eigen::Vector3d)> update_func, const Eigen::Ref<const Eigen::MatrixXb> vertexBoundaryConditions)
    {
        auto currentvertices = getVertices();
        const int nVertices = getNumberOfVertices();
        for(int i=0;i<nVertices;++i)
        {
            const Eigen::Vector3d v_new = update_func(currentvertices.row(i));
            for(int j=0;j<3;++j)
                if(not vertexBoundaryConditions(i,j)) currentvertices(i,j) = v_new(j);
        }
    }

    void changeEdgeDirectors(std::function<Real(Real)> update_func)
    {
        auto currentdirectors = getEdgeDirectors();
        const int nEdges = getNumberOfEdges();
        for(int i=0;i<nEdges;++i)
            currentdirectors(i) = update_func(currentdirectors(i));
    }
    
    void changeEdgeDirectors(std::function<Real(Real)> update_func, const Eigen::Ref<const Eigen::VectorXb> edgeBoundaryConditions)
    {
        auto currentdirectors = getEdgeDirectors();
        const int nEdges = getNumberOfEdges();
        for(int i=0;i<nEdges;++i)
            if(not edgeBoundaryConditions(i)) currentdirectors(i) = update_func(currentdirectors(i));
    }
    
    void computeFaceNormalsFromDirectors(const TopologyData & topology, const BoundaryConditionsData & boundaryConditions, Eigen::Ref<Eigen::MatrixXd> computed_normals) const
    {
        // create extended configuration just for this case
        DCSCurrentConfigurationData dcs_ext_data;
        dcs_ext_data.update(dcsdata, topology, boundaryConditions);
        
        // call the method
        computeFaceNormalsFromDirectors(dcs_ext_data.ext_tinfo, computed_normals);
    }
    
    void computeFaceNormalsFromDirectors(const std::vector<ExtendedTriangleInfo> & vInfo, Eigen::Ref<Eigen::MatrixXd> computed_normals) const
    {
        const int nFaces = (int)vInfo.size();
        for(int i=0;i<nFaces;++i)
        {
            const ExtendedTriangleInfo & info = vInfo[i];
            computed_normals.row(i) = info.computeFaceNormalFromDirectors();
        }
    }
    
    Eigen::VectorXd computeDoubleFaceAreas(const TopologyData & topo) const
    {
        const int nFaces = topo.getNumberOfFaces();
        const auto vertexdata = dcsdata.getVertexData();
        
        Eigen::VectorXd doubleFaceAreas(nFaces);
        doubleFaceAreas.setZero();
        for(int f=0; f<nFaces; ++f)
        {
            for(int d=0; d<3; d++)
            {
                const int x = d;
                const int y = (d+1)%3;
                const auto rx = vertexdata(topo.face2vertices(f,0),x) - vertexdata(topo.face2vertices(f,2),x);
                const auto sx = vertexdata(topo.face2vertices(f,1),x) - vertexdata(topo.face2vertices(f,2),x);
                const auto ry = vertexdata(topo.face2vertices(f,0),y) - vertexdata(topo.face2vertices(f,2),y);
                const auto sy = vertexdata(topo.face2vertices(f,1),y) - vertexdata(topo.face2vertices(f,2),y);
                
                double dblAd = rx*sy - ry*sx;
                doubleFaceAreas(f) += dblAd*dblAd;
            }
        }
        doubleFaceAreas = doubleFaceAreas.array().sqrt().eval();
        return doubleFaceAreas;
    }
    
    void computeAreaVectors(const TopologyData & topo, Eigen::Ref<Eigen::VectorXd> areaVertices, Eigen::Ref<Eigen::VectorXd> areaEdges) const
    {
        const Eigen::VectorXd doubleFaceAreas = computeDoubleFaceAreas(topo);
        computeAreaVectors(topo, doubleFaceAreas, areaVertices, areaEdges);
    }
    
    void computeAreaVectors(const TopologyData & topo, const Eigen::Ref<const Eigen::VectorXd> doubleFaceAreas, Eigen::Ref<Eigen::VectorXd> areaVertices, Eigen::Ref<Eigen::VectorXd> areaEdges) const
    {
        const int nFaces = topo.getNumberOfFaces();
        
        areaVertices.setZero();
        areaEdges.setZero();
        
        for(int i=0;i<nFaces;++i)
        {
            const Real areaPart = doubleFaceAreas(i)/6.0;
            
            const int idx_v0 = topo.face2vertices(i,0);
            const int idx_v1 = topo.face2vertices(i,1);
            const int idx_v2 = topo.face2vertices(i,2);
            
            areaVertices(idx_v0) += areaPart;
            areaVertices(idx_v1) += areaPart;
            areaVertices(idx_v2) += areaPart;
            
            const int idx_e0 = topo.face2edges(i,0);
            const int idx_e1 = topo.face2edges(i,1);
            const int idx_e2 = topo.face2edges(i,2);
            
            areaEdges(idx_e0) += areaPart;
            areaEdges(idx_e1) += areaPart;
            areaEdges(idx_e2) += areaPart;
        }
    }
    
    virtual void init(const Eigen::MatrixXd & vertices_in, const int nEdges)
    {
        dcsdata.init(vertices_in, nEdges);
    }
    
    virtual void init(const Eigen::MatrixXd & vertices_in, const Eigen::VectorXd & edgeDirectors_in)
    {
        dcsdata.init(vertices_in, edgeDirectors_in);
    }
    
    int getNumberOfDataVariables() const
    {
        return dcsdata.getNumberOfDataVariables();
    }
    
    int getNumberOfVertices() const
    {
        return dcsdata.getNumberOfVertices();
    }
    
    int getNumberOfEdges() const
    {
        return dcsdata.getNumberOfEdges();
    }
    
    Real * getDataPointer() const
    {
        return dcsdata.getDataPointer();
    }
    
    Eigen::Ref<Eigen::VectorXd> getPositions()
    {
        return Eigen::Map<Eigen::VectorXd>(dcsdata.getDataPointer(), dcsdata.getNumberOfDataVariables());
    }
    
    Eigen::Ref<Eigen::MatrixXd> getVertices()
    {
        return dcsdata.getVertexData();
    }
    
    Eigen::Ref<Eigen::VectorXd> getEdgeDirectors()
    {
        return dcsdata.getEdgeData();
    }
    
    const Eigen::Ref<const Eigen::VectorXd> getPositions() const
    {
        return Eigen::Map<const Eigen::VectorXd>(dcsdata.getDataPointer(), dcsdata.getNumberOfDataVariables());
    }
    
    
    const Eigen::Ref<const Eigen::MatrixXd> getVertices() const
    {
        return dcsdata.getVertexData();
    }
    
    const Eigen::Ref<const Eigen::VectorXd> getEdgeDirectors() const
    {
        return dcsdata.getEdgeData();
    }
    
    TriangleInfo getTriangleInfoLite(const TopologyData & topo, const int face_idx) const
    {
        return dcsdata.getTriangleInfo(topo, face_idx);
    }
};

struct DCSCurrentConfiguration : DCSBasicConfiguration
{
    DCSCurrentConfigurationData dcs_ext_data;
    
    void clear() override
    {
        dcs_ext_data.clear();
        DCSBasicConfiguration::clear();
    }
    
    void update(const TopologyData & topology, const BoundaryConditionsData & boundaryConditions)
    {
        dcs_ext_data.update(dcsdata, topology, boundaryConditions);
    }
    
    // subset update
    void update(const TopologyData & topology, const BoundaryConditionsData & boundaryConditions, const std::vector<int> & face_indices)
    {
        dcs_ext_data.updateSubSet(dcsdata, topology, boundaryConditions, face_indices);
    }
    
    const std::vector<ExtendedTriangleInfo> & getTriangleInfos() const
    {
        return dcs_ext_data.ext_tinfo;
    }
    
    const ExtendedTriangleInfo & getTriangleInfo(const int face_idx) const
    {
        return dcs_ext_data.ext_tinfo[face_idx];
    }
};


struct DCSRestConfiguration : DCSBasicConfiguration
{
    tVecMat2d aform, bform;
    
    void clear() override
    {
        aform.clear();
        bform.clear();
        DCSBasicConfiguration::clear();
    }
    
    void setFormsFromVertices(const TopologyData & topology, const BoundaryConditionsData & boundaryConditions)
    {
        DCSCurrentConfigurationData dcs_ext_data;
        dcs_ext_data.update(dcsdata, topology, boundaryConditions);
        
        const int nFaces = topology.getNumberOfFaces();
        aform.resize(nFaces);
        bform.resize(nFaces);
        
        for(int i=0;i<nFaces;++i)
        {
            aform[i] = dcs_ext_data.ext_tinfo[i].computeFirstFundamentalForm();
            bform[i] = dcs_ext_data.ext_tinfo[i].computeSecondFundamentalForm();
        }
    }
    
    void update(const TopologyData & topology, const BoundaryConditionsData & boundaryConditions)
    {
        setFormsFromVertices(topology, boundaryConditions);
    }
    
    tVecMat2d & getFirstFundamentalForms()
    {
        return aform;
    }
    
    tVecMat2d & getSecondFundamentalForms()
    {
        return bform;
    }
    
    const tVecMat2d & getFirstFundamentalForms() const
    {
        return aform;
    }
    
    const tVecMat2d & getSecondFundamentalForms() const
    {
        return bform;
    }
    
    Eigen::Matrix2d getFirstFundamentalForm(const int face_idx) const
    {
        return aform[face_idx];
    }
    
    Eigen::Matrix2d getSecondFundamentalForm(const int face_idx) const
    {
        return bform[face_idx];
    }
    
    DCSRestConfiguration& operator=(const DCSCurrentConfiguration& rhs)
    {
        // assign the base data
        DCSBasicConfiguration::operator=(rhs);
        
        // need to call the assignFormsFromVertices by hand...
        return *this;
    }
};

struct DCSBilayerRestConfiguration : DCSBasicConfiguration
{
    tVecMat2d aform_bot, aform_top, bform;
    
    void clear() override
    {
        aform_bot.clear();
        aform_top.clear();
        bform.clear();
        DCSBasicConfiguration::clear();
    }
    
    void setFormsFromVertices(const TopologyData & topology, const BoundaryConditionsData & boundaryConditions)
    {
        DCSCurrentConfigurationData dcs_ext_data;
        dcs_ext_data.update(dcsdata, topology, boundaryConditions);
        
        const int nFaces = topology.getNumberOfFaces();
        aform_bot.resize(nFaces);
        aform_top.resize(nFaces);
        bform.resize(nFaces);
        
        for(int i=0;i<nFaces;++i)
        {
            aform_bot[i] = aform_top[i] = dcs_ext_data.ext_tinfo[i].computeFirstFundamentalForm();
            bform[i] = dcs_ext_data.ext_tinfo[i].computeSecondFundamentalForm();
        }
    }
    
    void update(const TopologyData & topology, const BoundaryConditionsData & boundaryConditions)
    {
        setFormsFromVertices(topology, boundaryConditions);
    }
    
    template<MeshLayer layer>
    tVecMat2d & getFirstFundamentalForms()
    {
        static_assert(layer != single, "Can not use DCSBilayerRestConfiguration with single layer type");
        return (layer == bottom ? aform_bot : aform_top);
    }
    
    tVecMat2d & getSecondFundamentalForms()
    {
        return bform;
    }
    
    template<MeshLayer layer>
    const tVecMat2d & getFirstFundamentalForms() const
    {
        static_assert(layer != single, "Can not use DCSBilayerRestConfiguration with single layer type");
        return (layer == bottom ? aform_bot : aform_top);
    }
    
    const tVecMat2d & getSecondFundamentalForms() const
    {
        return bform;
    }
    
    template<MeshLayer layer>
    Eigen::Matrix2d getFirstFundamentalForm(const int face_idx) const
    {
        static_assert(layer != single, "Can not use DCSBilayerRestConfiguration with single layer type");
        return (layer == bottom ? aform_bot[face_idx] : aform_top[face_idx]); // compile-time conditional
    }
    
    Eigen::Matrix2d getSecondFundamentalForm(const int face_idx) const
    {
        return bform[face_idx];
    }
};


struct DCSDynamicCurrentConfiguration : DCSCurrentConfiguration
{
    // add the raw data for the dynamics
    DCSConfigurationData dyndcsdata;
    
    // default constructor
    DCSDynamicCurrentConfiguration():
    DCSCurrentConfiguration()
    {}
    
    // assignment operator
    DCSDynamicCurrentConfiguration & operator= (const DCSDynamicCurrentConfiguration & other)
    {
        if (this != &other) // protect against invalid self-assignment
        {
            DCSCurrentConfiguration::operator=(other);
            dyndcsdata = other.dyndcsdata;
        }
        // by convention, always return *this
        return *this;
    }
    
    // copy constructor
    DCSDynamicCurrentConfiguration(const DCSDynamicCurrentConfiguration &other):
    DCSCurrentConfiguration(other),
    dyndcsdata(other.dyndcsdata)
    {
    }
    
    friend void swap(DCSDynamicCurrentConfiguration& first, DCSDynamicCurrentConfiguration& second)
    {
        // http://stackoverflow.com/questions/3279543/what-is-the-copy-and-swap-idiom
        using std::swap;
        
        swap(first.dyndcsdata, second.dyndcsdata);
        
        swap(static_cast<DCSCurrentConfiguration&>(first), static_cast<DCSCurrentConfiguration&>(second));
    }
    
    void clear() override
    {
        dyndcsdata.clear();
        DCSCurrentConfiguration::clear();
    }
    
    virtual void init(const Eigen::MatrixXd & vertices_in, const int nEdges) override
    {
        DCSCurrentConfiguration::init(vertices_in, nEdges);
        dyndcsdata.init(vertices_in.rows(), nEdges);
    }
    
    virtual void init(const Eigen::MatrixXd & vertices_in, const Eigen::VectorXd & edgeDirectors_in) override
    {
        DCSCurrentConfiguration::init(vertices_in, edgeDirectors_in);
        dyndcsdata.init(vertices_in.rows(), edgeDirectors_in.rows());
    }
    
    Real * getDynamicDataPointer() const
    {
        return dyndcsdata.getDataPointer();
    }
    
    Eigen::Ref<Eigen::VectorXd> getVelocities()
    {
        return Eigen::Map<Eigen::VectorXd>(dyndcsdata.getDataPointer(), dyndcsdata.getNumberOfDataVariables());
    }
    
    Eigen::Ref<Eigen::MatrixXd> getvVertices()
    {
        return dyndcsdata.getVertexData();
    }
    
    Eigen::Ref<Eigen::VectorXd> getvEdgeDirectors()
    {
        return dyndcsdata.getEdgeData();
    }
    
    const Eigen::Ref<const Eigen::VectorXd> getVelocities() const
    {
        return Eigen::Map<const Eigen::VectorXd>(dyndcsdata.getDataPointer(), dyndcsdata.getNumberOfDataVariables());
    }
    
    const Eigen::Ref<const Eigen::MatrixXd> getvVertices() const
    {
        return dyndcsdata.getVertexData();
    }
    
    const Eigen::Ref<const Eigen::VectorXd> getvEdgeDirectors() const
    {
        return dyndcsdata.getEdgeData();
    }
    
    TriangleInfo getvTriangleInfo(const TopologyData & topo, const int face_idx) const
    {
        return dyndcsdata.getTriangleInfo(topo, face_idx);
    }
};

#endif /* DCSConfigurations_hpp */
