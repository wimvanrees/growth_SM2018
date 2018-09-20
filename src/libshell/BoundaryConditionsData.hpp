//
//  BoundaryConditionsData.hpp
//  Elasticity
//
//  Created by Wim van Rees on 5/4/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef BoundaryConditionsData_hpp
#define BoundaryConditionsData_hpp

#include "common.hpp"
#include "TopologyData.hpp"

struct BoundaryConditionsData
{
    Eigen::MatrixXb vertices_bc;
    Eigen::VectorXb edges_bc;
    
    Eigen::MatrixXd clampedEdgeNormals;
    std::map<int, int> clampedEdgeNormalsIndices;
    
    void init(const Eigen::MatrixXb & vertices_bc_in, const int nEdges)
    {
        vertices_bc = vertices_bc_in;
        
        edges_bc.resize(nEdges);
        for(int i=0;i<nEdges;++i)
            edges_bc(i) = false;
        
    }
    
    void clear()
    {
        vertices_bc.resize(0,0);
        edges_bc.resize(0);
        clampedEdgeNormals.resize(0,0);
        clampedEdgeNormalsIndices.clear();
    }
    
    int clampFixedEdges(const TopologyData & topology, const Eigen::Ref<const Eigen::MatrixXd> & vertices, std::function<Eigen::Vector3d(const Eigen::Vector3d)> midedge_normals)
    {
        const int nEdges = topology.getNumberOfEdges();
        
        // count number of fixed edges
        int cnt=0;
        for(int i=0;i<nEdges;++i)
        {
            const int idx_v0 = topology.edge2vertices(i,0);
            const int idx_v1 = topology.edge2vertices(i,1);
            
            const bool fixed_v0 = (vertices_bc(idx_v0,0) and vertices_bc(idx_v0,1) and vertices_bc(idx_v0,2));
            const bool fixed_v1 = (vertices_bc(idx_v1,0) and vertices_bc(idx_v1,1) and vertices_bc(idx_v1,2));
            
            if(fixed_v0 and fixed_v1) cnt++;
        }
        
        // allocate memory
        clampedEdgeNormals.resize(cnt,3);
        clampedEdgeNormals.setZero();

        // fill arrays
        cnt=0;
        for(int i=0;i<nEdges;++i)
        {
            const int idx_v0 = topology.edge2vertices(i,0);
            const int idx_v1 = topology.edge2vertices(i,1);
            
            const bool fixed_v0 = (vertices_bc(idx_v0,0) and vertices_bc(idx_v0,1) and vertices_bc(idx_v0,2));
            const bool fixed_v1 = (vertices_bc(idx_v1,0) and vertices_bc(idx_v1,1) and vertices_bc(idx_v1,2));
            
            if(fixed_v0 and fixed_v1)
            {
                clampedEdgeNormalsIndices.insert( std::pair<int,int>(i, cnt));
                const Eigen::Vector3d pos_v0 = vertices.row(idx_v0);
                const Eigen::Vector3d pos_v1 = vertices.row(idx_v1);
                const Eigen::Vector3d pos_edge = 0.5*(pos_v0 + pos_v1);
                const Eigen::Vector3d normal = midedge_normals(pos_edge);
                clampedEdgeNormals(cnt,0) = normal(0);
                clampedEdgeNormals(cnt,1) = normal(1);
                clampedEdgeNormals(cnt,2) = normal(2);
                
                edges_bc(i) = true;
                
                cnt++;
            }
            else
                edges_bc(i) = false;
        }
        
        return cnt;
    }
    
    int clampFixedEdgesDiscrete(const TopologyData & topology, const Eigen::Ref<const Eigen::MatrixXd> & midedge_normals)
    {
        const int nEdges = topology.getNumberOfEdges();
        
        assert(midedge_normals.rows() == nEdges);
        assert(midedge_normals.cols() == 3);
        
        // count number of fixed edges
        int cnt=0;
        for(int i=0;i<nEdges;++i)
        {
            const int idx_v0 = topology.edge2vertices(i,0);
            const int idx_v1 = topology.edge2vertices(i,1);
            
            const bool fixed_v0 = (vertices_bc(idx_v0,0) and vertices_bc(idx_v0,1) and vertices_bc(idx_v0,2));
            const bool fixed_v1 = (vertices_bc(idx_v1,0) and vertices_bc(idx_v1,1) and vertices_bc(idx_v1,2));
            
            if(fixed_v0 and fixed_v1) cnt++;
        }
        
        // allocate memory
        clampedEdgeNormals.resize(cnt,3);
        clampedEdgeNormals.setZero();
        
        // fill arrays
        cnt=0;
        for(int i=0;i<nEdges;++i)
        {
            const int idx_v0 = topology.edge2vertices(i,0);
            const int idx_v1 = topology.edge2vertices(i,1);
            
            const bool fixed_v0 = (vertices_bc(idx_v0,0) and vertices_bc(idx_v0,1) and vertices_bc(idx_v0,2));
            const bool fixed_v1 = (vertices_bc(idx_v1,0) and vertices_bc(idx_v1,1) and vertices_bc(idx_v1,2));
            
            if(fixed_v0 and fixed_v1)
            {
                clampedEdgeNormalsIndices.insert( std::pair<int,int>(i, cnt));
                const Eigen::Vector3d normal = midedge_normals.row(i);
                clampedEdgeNormals(cnt,0) = normal(0);
                clampedEdgeNormals(cnt,1) = normal(1);
                clampedEdgeNormals(cnt,2) = normal(2);
                
                edges_bc(i) = true;
                
                cnt++;
            }
            else
            edges_bc(i) = false;
        }
        
        return cnt;
    }
    
    Eigen::Ref<Eigen::MatrixXb> getVertexBoundaryConditions()
    {
        return vertices_bc;
    }
    
    const Eigen::Ref<const Eigen::MatrixXb> getVertexBoundaryConditions() const
    {
        return vertices_bc;
    }
    
    const Eigen::Ref<const Eigen::VectorXb> getEdgeBoundaryConditions() const
    {
        return edges_bc;
    }
    
    const Eigen::Ref<const Eigen::MatrixXd> getClampedEdgeNormals() const
    {
        return clampedEdgeNormals;
    }
    
    const std::map<int, int> & getClampedEdgeNormalsIndices() const
    {
        return clampedEdgeNormalsIndices;
    }
    
};

#endif /* BoundaryConditionsData_hpp */
