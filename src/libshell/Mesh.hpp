//
//  Mesh.hpp
//  Elasticity
//
//  Created by Wim van Rees on 2/18/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//
//! \file Mesh.hpp


#ifndef Mesh_hpp
#define Mesh_hpp

#include "common.hpp"
#include "TopologyData.hpp"
#include "DCSConfigurations.hpp"
#include "BoundaryConditionsData.hpp"
#include "Geometry.hpp"

#ifdef USELIBLBFGS
#include "lbfgs.h"
#endif

/*! \class Mesh
 * \brief To store the geometry of helper arrays related to the current triangle mesh.
 *
 * More details here
 */
template<typename tCCD, typename tRCD>
class BaseMesh
{
public:
    typedef tCCD tCurrentConfigData;
    typedef tRCD tReferenceConfigData;
    
protected:
    tCurrentConfigData currentState;
    tReferenceConfigData restState;
    
    TopologyData topology;
    BoundaryConditionsData boundaryConditions;
    
public:
    BaseMesh()
    {
        // no initialization: happens on init instead
        // maybe need a bool isInitialized flag or something
    }
    
    /*! initialize the grid with the these vertex and faces list, as well as BCs. Will set the current state as rest state */
    virtual void init(const Geometry & geometry, const bool bClampFixedEdges = false)
    {
        currentState.clear();
        restState.clear();
        topology.clear();
        boundaryConditions.clear();
        
        Eigen::MatrixXd vertices_in;
        Eigen::MatrixXi faces_in;
        Eigen::MatrixXb vertices_bc_in;
        
        geometry.get(vertices_in, faces_in, vertices_bc_in);
        topology.init(vertices_in, faces_in);
        
        const int nEdges = getNumberOfEdges();
        
        // now we can initialize our variables: assume current state is rest state
        restState.init(vertices_in, nEdges);
        currentState.init(vertices_in, nEdges);
        
        // set the rest and edge directors using the geometry function
        const bool bAnalyticNormals = geometry.analyticNormals();
        if(not bAnalyticNormals && bClampFixedEdges)
        {
            helpers::catastrophe("When specifying clamped edges - need to specify analytic normals as well!", __FILE__, __LINE__);
        }
        
        std::function<Eigen::Vector3d(const Eigen::Vector3d)> normal_function = std::bind(&Geometry::getNormal, &geometry, std::placeholders::_1);
        
        if(bAnalyticNormals)
        {
            restState.setEdgeDirectors(topology, normal_function);
            currentState.setEdgeDirectors(topology, normal_function);
        }
        // else : by default edge directors are set to zero

        // boundary conditions
        boundaryConditions.init(vertices_bc_in, nEdges);
        if(bClampFixedEdges) boundaryConditions.clampFixedEdges(topology, currentState.getVertices(), normal_function); // fix wrt current configuration : this could also be changed
        
        restState.update(topology, boundaryConditions);
        currentState.update(topology, boundaryConditions);
        

    }
    
    void updateDeformedConfiguration()
    {
        currentState.update(topology, boundaryConditions);
    }
    
    void updateDeformedConfiguration(const std::vector<int> & face_indices)
    {
        currentState.update(topology, boundaryConditions, face_indices);
    }
    
    virtual void resetToRestState()
    {
        currentState.dcsdata.getData() = restState.dcsdata.getData();
        updateDeformedConfiguration();
    }
    
    // boundary conditions according to a function
    void setVertexBoundaryConditions(std::function<Eigen::Vector3b(Real,Real,Real)> vertex_bc)
    {
        const int nVertices = getNumberOfVertices();
        const auto restvertices = restState.getVertices();
        
        for(int i=0;i<nVertices;++i)
        {
            const Real vx = restvertices(i,0);
            const Real vy = restvertices(i,1);
            const Real vz = restvertices(i,2);
            const Eigen::Vector3b isFixed = vertex_bc(vx,vy,vz);
            boundaryConditions.vertices_bc(i,0) = isFixed(0);
            boundaryConditions.vertices_bc(i,1) = isFixed(1);
            boundaryConditions.vertices_bc(i,2) = isFixed(2);
        }
    }
    
    void setEdgeBoundaryConditions(std::function<bool(Real,Real,Real)> edge_bc)
    {
        const int nEdges = getNumberOfEdges();
        const auto restvertices = restState.getVertices();
        
        int count = 0;
        for(int i=0;i<nEdges;++i)
        {
            const int idx_v0 = topology.edge2vertices(i,0);
            const int idx_v1 = topology.edge2vertices(i,1);
            const Real ex = 0.5*(restvertices(idx_v0,0) + restvertices(idx_v1,0));
            const Real ey = 0.5*(restvertices(idx_v0,1) + restvertices(idx_v1,1));
            const Real ez = 0.5*(restvertices(idx_v0,2) + restvertices(idx_v1,2));
            
            const bool isFixed = edge_bc(ex, ey, ez);
            boundaryConditions.edges_bc(i) = isFixed;
            
            if(isFixed) count++;
        }
        std::cout << "number of fixed edges = " << count << std::endl;
    }
    
    int clampFixedEdges()
    {
        auto vertical_normals = [](const Eigen::Vector3d) -> Eigen::Vector3d
        {
            Eigen::Vector3d retval;
            retval << 0,0,1;
            return retval;
        };
        
        return boundaryConditions.clampFixedEdges(topology, currentState.getVertices(), vertical_normals);
    }
    
    int clampFixedEdges(std::function<Eigen::Vector3d(const Eigen::Vector3d)> midedge_normals)
    {
        return boundaryConditions.clampFixedEdges(topology, currentState.getVertices(), midedge_normals);
    }
    
    int clampFixedEdgesDiscrete(const Eigen::Ref<const Eigen::MatrixXd> & midedge_normals)
    {
        return boundaryConditions.clampFixedEdgesDiscrete(topology, midedge_normals);
    }
    
    void changeVertices(std::function<Eigen::Vector3d(Eigen::Vector3d)> update_func, const bool keepBoundaryConditions=true)
    {
        if(keepBoundaryConditions)
        { 
            currentState.changeVertices(update_func, boundaryConditions.getVertexBoundaryConditions());
        }
        else
        {
            currentState.changeVertices(update_func);
        }
    }
    
    void changeEdgeDirectors(std::function<Real(Real)> update_func, const bool keepBoundaryConditions=true)
    {
        if(keepBoundaryConditions)
        {
            currentState.changeEdgeDirectors(update_func, boundaryConditions.getEdgeBoundaryConditions());
        }
        else
        {
            currentState.changeEdgeDirectors(update_func);
        }
    }
    
    // getter for data pointer (optimization)
    virtual Real * getDataPointer()
    {
        return currentState.getDataPointer();
    }
    
    // non-const getters
    TopologyData & getTopology()
    {
        return topology;
    }
    
    BoundaryConditionsData & getBoundaryConditions()
    {
        return boundaryConditions;
    }
    
    tCurrentConfigData & getCurrentConfiguration()
    {
        return currentState;
    }
    
    tReferenceConfigData & getRestConfiguration()
    {
        return restState;
    }
    
    
    // const getters
    int getNumberOfVertices() const
    {
        return currentState.getNumberOfVertices();
    }
    
    int getNumberOfFaces() const
    {
        return topology.getNumberOfFaces();
    }
    
    int getNumberOfEdges() const
    {
        return topology.getNumberOfEdges();
    }
    
    const TopologyData & getTopology() const
    {
        return topology;
    }
    
    const BoundaryConditionsData & getBoundaryConditions() const
    {
        return boundaryConditions;
    }
    
    const tCurrentConfigData & getCurrentConfiguration() const
    {
        return currentState;
    }
    
    const tReferenceConfigData & getRestConfiguration() const
    {
        return restState;
    }

    virtual ~BaseMesh()
    {
    }
};

class Mesh : public BaseMesh<DCSCurrentConfiguration, DCSRestConfiguration>
{
};

class BilayerMesh : public BaseMesh<DCSCurrentConfiguration, DCSBilayerRestConfiguration>
{
};



#endif /* Mesh_hpp */
