//
//  CreateExtendedTriangleInfos.hpp
//  Elasticity
//
//  Created by Wim van Rees on 1/18/17.
//  Copyright © 2017 Wim van Rees. All rights reserved.
//

#ifndef CreateExtendedTriangleInfos_hpp
#define CreateExtendedTriangleInfos_hpp

#include "common.hpp"
#include "ConfigurationData.hpp"
#include "ExtendedTriangleInfo.hpp"
#include "TopologyData.hpp"
#include "BoundaryConditionsData.hpp"

struct CreateExtendedTriangleInfos_tbb_base
{
    const Eigen::Ref<const Eigen::MatrixXi> face2vertices;
    const Eigen::Ref<const Eigen::MatrixXi> face2edges;
    const Eigen::Ref<const Eigen::MatrixXi> edge2faces;
    const Eigen::Ref<const Eigen::MatrixXi> edge2vertices;
    
    const DCSConfigurationData & data;
    
    const Eigen::Ref<const Eigen::MatrixXb> vertices_bc;
    const Eigen::Ref<const Eigen::VectorXb> edges_bc;
    const Eigen::Ref<const Eigen::MatrixXd>  clampedEdgeNormals;
    const std::map<int, int> & clampedEdgeNormalsIndices;
    
    std::vector<ExtendedTriangleInfo> & vInfos;
    
    CreateExtendedTriangleInfos_tbb_base(const TopologyData & topo, const DCSConfigurationData & data, const BoundaryConditionsData & boundaryConditions, std::vector<ExtendedTriangleInfo> & vInfos_in):
    face2vertices(topo.getFace2Vertices()),
    face2edges(topo.getFace2Edges()),
    edge2faces(topo.getEdge2Faces()),
    edge2vertices(topo.getEdge2Vertices()),
    data(data),
    vertices_bc(boundaryConditions.getVertexBoundaryConditions()),
    edges_bc(boundaryConditions.getEdgeBoundaryConditions()),
    clampedEdgeNormals(boundaryConditions.getClampedEdgeNormals()),
    clampedEdgeNormalsIndices(boundaryConditions.getClampedEdgeNormalsIndices()),
    vInfos(vInfos_in)
    {}
    
    // split constructor (dont need copy constructor)
    CreateExtendedTriangleInfos_tbb_base(const CreateExtendedTriangleInfos_tbb_base& c):
    face2vertices(c.face2vertices),
    face2edges(c.face2edges),
    edge2faces(c.edge2faces),
    edge2vertices(c.edge2vertices),
    data(c.data),
    vertices_bc(c.vertices_bc),
    edges_bc(c.edges_bc),
    clampedEdgeNormals(c.clampedEdgeNormals),
    clampedEdgeNormalsIndices(c.clampedEdgeNormalsIndices),
    vInfos(c.vInfos)
    {
    }
    
    void computeVertexData(const int face_idx) const
    {
        ExtendedTriangleInfo & info = vInfos[face_idx];
        
        // =============== INDEX INFORMATION ===============//
        info.face_idx = face_idx;
        
        info.idx_v0 = face2vertices(info.face_idx,0);
        info.idx_v1 = face2vertices(info.face_idx,1);
        info.idx_v2 = face2vertices(info.face_idx,2);
        
        info.idx_e0 = face2edges(info.face_idx, 0);
        info.idx_e1 = face2edges(info.face_idx, 1);
        info.idx_e2 = face2edges(info.face_idx, 2);
        
        // =============== RAW DATA ===============//
        info.v0 = data.getVertexData().row(info.idx_v0);
        info.v1 = data.getVertexData().row(info.idx_v1);
        info.v2 = data.getVertexData().row(info.idx_v2);
        
        info.phi(0) = data.getEdgeData()(info.idx_e0);
        info.phi(1) = data.getEdgeData()(info.idx_e1);
        info.phi(2) = data.getEdgeData()(info.idx_e2);
        
        info.e0 = info.v1 - info.v0;
        info.e1 = info.v2 - info.v1;
        info.e2 = info.v0 - info.v2;
        
        // =============== TOPOLOGY INFORMATION ===============//
        // get the sign (does not matter as long as its consistent for each of the two faces around an edge)
        info.sign(0) = (edge2faces(info.idx_e0,0) == face_idx ? +1 : -1);
        info.sign(1) = (edge2faces(info.idx_e1,0) == face_idx ? +1 : -1);
        info.sign(2) = (edge2faces(info.idx_e2,0) == face_idx ? +1 : -1);
        
        // precompute the other (opposite) face indices to each edge (no need to store them : we store other_faces which has the face_idx in there)
        const int other_face_idx_e0 = edge2faces(info.idx_e0,0) == face_idx ? edge2faces(info.idx_e0,1) : edge2faces(info.idx_e0,0);
        const int other_face_idx_e1 = edge2faces(info.idx_e1,0) == face_idx ? edge2faces(info.idx_e1,1) : edge2faces(info.idx_e1,0);
        const int other_face_idx_e2 = edge2faces(info.idx_e2,0) == face_idx ? edge2faces(info.idx_e2,1) : edge2faces(info.idx_e2,0);
        
        // also precompute the relative index of each edge from the perspective of the opposite face (we need that when computing gradients of bform)
        if(other_face_idx_e0 >= 0)
        {
            info.other_face_edge_idx(0) = (face2edges(other_face_idx_e0,0) == info.idx_e0 ? 0 : (face2edges(other_face_idx_e0,1) == info.idx_e0 ? 1 : 2) );
            info.other_faces[0] = &vInfos[other_face_idx_e0];
        }
        if(other_face_idx_e1 >= 0)
        {
            info.other_face_edge_idx(1) = (face2edges(other_face_idx_e1,0) == info.idx_e1 ? 0 : (face2edges(other_face_idx_e1,1) == info.idx_e1 ? 1 : 2) );
            info.other_faces[1] = &vInfos[other_face_idx_e1];
        }
        if(other_face_idx_e2 >= 0)
        {
            info.other_face_edge_idx(2) = (face2edges(other_face_idx_e2,0) == info.idx_e2 ? 0 : (face2edges(other_face_idx_e2,1) == info.idx_e2 ? 1 : 2) );
            info.other_faces[2] = &vInfos[other_face_idx_e2];
        }
        
        // =============== DERIVED INFORMATION ===============//
        
        // face area
        {
            info.double_face_area = 0.0;
            for(int d=0;d<3;++d)
            {
                const int x = d;
                const int y = (d+1)%3;
                const auto rx = info.v0(x) - info.v2(x);
                const auto sx = info.v1(x) - info.v2(x);
                const auto ry = info.v0(y) - info.v2(y);
                const auto sy = info.v1(y) - info.v2(y);
                
                double dblAd = rx*sy - ry*sx;
                info.double_face_area += dblAd*dblAd;
            }
            info.double_face_area = std::sqrt(info.double_face_area);
        }
        
        // edge lengths
        {
            info.edgelength(0) = info.e0.norm();
            info.edgelength(1) = info.e1.norm();
            info.edgelength(2) = info.e2.norm();
        }
        
        // normal vectors
        {
            info.face_normal = (info.e2).cross(info.e0);
            const Real r = info.face_normal.norm();
            if(r==0)
                info.face_normal.setZero();
            else
                info.face_normal /= r;
        }
        
        // triangle heights
        {
            info.height(0) = info.double_face_area / info.edgelength(0);
            info.height(1) = info.double_face_area / info.edgelength(1);
            info.height(2) = info.double_face_area / info.edgelength(2);
        }
        
        // cos angles
        {
            info.cosGamma(0) = -info.e1.dot(info.e2)/(info.edgelength(1)*info.edgelength(2));
            info.cosGamma(1) = -info.e0.dot(info.e2)/(info.edgelength(0)*info.edgelength(2));
            info.cosGamma(2) = -info.e0.dot(info.e1)/(info.edgelength(0)*info.edgelength(1));
        }
        
        
        // edge normals
        {
            info.edge_normals.col(0) = (info.e0.cross(info.face_normal)).normalized();
            info.edge_normals.col(1) = (info.e1.cross(info.face_normal)).normalized();
            info.edge_normals.col(2) = (info.e2.cross(info.face_normal)).normalized();
        }
        
        // dihedral angles : do it in a second loop
        info.theta.setZero();
        
        // alpha angles : add the 1/2*theta in a second loop (this way, if theta=0, we still have the correct value)
        info.alpha = (info.sign.array() * info.phi.array()).matrix();
        //            info.alpha(0) = 0.5*info.theta(0) + info.sign(0)*info.phi(0);
        //            info.alpha(1) = 0.5*info.theta(1) + info.sign(1)*info.phi(1);
        //            info.alpha(2) = 0.5*info.theta(2) + info.sign(2)*info.phi(2);
        
        
        // =============== BOUNDARY CONDITIONS ===============//
        for(int d=0;d<3;++d)
        {
            info.v0_fixed(d) = vertices_bc(info.idx_v0, d);
            info.v1_fixed(d) = vertices_bc(info.idx_v1, d);
            info.v2_fixed(d) = vertices_bc(info.idx_v2, d);
        }
        
        info.e0_clamped = edges_bc(info.idx_e0);
        info.e1_clamped = edges_bc(info.idx_e1);
        info.e2_clamped = edges_bc(info.idx_e2);
        
        info.clamped_normal_e0 = (info.e0_clamped ? (Eigen::Vector3d)clampedEdgeNormals.row(clampedEdgeNormalsIndices.find(info.idx_e0)->second) : Eigen::Vector3d::Zero());
        info.clamped_normal_e1 = (info.e1_clamped ? (Eigen::Vector3d)clampedEdgeNormals.row(clampedEdgeNormalsIndices.find(info.idx_e1)->second) : Eigen::Vector3d::Zero());
        info.clamped_normal_e2 = (info.e2_clamped ? (Eigen::Vector3d)clampedEdgeNormals.row(clampedEdgeNormalsIndices.find(info.idx_e2)->second) : Eigen::Vector3d::Zero());
    }
    
    void computeEdgeData(const int edge_idx) const
    {
        const int face_idx_0 = edge2faces(edge_idx,0);
        const int face_idx_1 = edge2faces(edge_idx,1);
        
        if(face_idx_0 < 0 or face_idx_1 < 0) return; // this edge is on the boundary
        
        // get the normal vectors
        const Eigen::Vector3d normal_0 = vInfos[face_idx_0].face_normal;
        const Eigen::Vector3d normal_1 = vInfos[face_idx_1].face_normal;
        
        // face 0 lies to the left of the edge and has the edge in its counter clock-wise orientation: so we assume the perspective of face 0 for this
        const int start = face2edges(face_idx_0,0) == edge_idx ? 0 : (face2edges(face_idx_0,1) == edge_idx ? 1 : (face2edges(face_idx_0,2) == edge_idx ? 2 : -1));
        assert(start>=0);
        const int v0_idx = face2vertices(face_idx_0,start);
        const int v1_idx = face2vertices(face_idx_0,(start+1)%3);
        
        assert((v0_idx == edge2vertices(edge_idx,0) && v1_idx == edge2vertices(edge_idx,1)) or (v0_idx == edge2vertices(edge_idx,1) && v1_idx == edge2vertices(edge_idx,0)));
        
        // compute edge
        const Eigen::Vector3d v0 = data.getVertexData().row(v0_idx);
        const Eigen::Vector3d v1 = data.getVertexData().row(v1_idx);
        const Eigen::Vector3d edge = v1 - v0;
        
        const Real signTheta = ( normal_0.cross(normal_1) ).dot(edge) > 0.0 ? 1.0 : -1.0; // theta has same sign as (n1 cross n2) dot e0 --> positive if the normals points away from each other when the edge points out of the plane
        
        const Real denum = (normal_0 + normal_1).norm();
        const Real theta = 2.0 * signTheta * std::atan2((normal_0 - normal_1).norm(), denum); // use atan2 because we are doing numerics
        
        // store them in the triangle infos
        vInfos[face_idx_0].setTheta(start, theta); // is this thread-safe?
        const int start_1 = face2edges(face_idx_1,0) == edge_idx ? 0 : (face2edges(face_idx_1,1) == edge_idx ? 1 : (face2edges(face_idx_1,2) == edge_idx ? 2 : -1));
        vInfos[face_idx_1].setTheta(start_1, theta); // is this thread-safe?
    }
};


template<int stage>
struct CreateExtendedTriangleInfos_tbb : CreateExtendedTriangleInfos_tbb_base
{
    CreateExtendedTriangleInfos_tbb(const TopologyData & topo, const DCSConfigurationData & data, const BoundaryConditionsData & boundaryConditions, std::vector<ExtendedTriangleInfo> & vInfos_in):
    CreateExtendedTriangleInfos_tbb_base(topo, data, boundaryConditions, vInfos_in)
    {}
    
    // split constructor (dont need copy constructor)
    CreateExtendedTriangleInfos_tbb(const CreateExtendedTriangleInfos_tbb& c):
    CreateExtendedTriangleInfos_tbb_base(c)
    {
    }
    
    void operator()(const tbb::blocked_range<int> & idx_range) const
    {
        for (int i=idx_range.begin(); i != idx_range.end(); ++i)
        {
            if(stage == 0) this->computeVertexData(i);
            else this->computeEdgeData(i);
        }
    }
};

template<int stage>
struct CreateExtendedTriangleInfos_tbb_subset : CreateExtendedTriangleInfos_tbb_base
{
    const std::vector<int> & face_indices;
    const std::vector<int> & edge_indices;
    
    CreateExtendedTriangleInfos_tbb_subset(const TopologyData & topo, const DCSConfigurationData & data, const BoundaryConditionsData & boundaryConditions, std::vector<ExtendedTriangleInfo> & vInfos_in, const std::vector<int> & face_indices_in, const std::vector<int> & edge_indices_in):
    CreateExtendedTriangleInfos_tbb_base(topo, data, boundaryConditions, vInfos_in),
    face_indices(face_indices_in),
    edge_indices(edge_indices_in)
    {}
    
    // split constructor (dont need copy constructor)
    CreateExtendedTriangleInfos_tbb_subset(const CreateExtendedTriangleInfos_tbb_subset& c):
    CreateExtendedTriangleInfos_tbb_base(c),
    face_indices(c.face_indices),
    edge_indices(c.edge_indices)
    {
    }
    
    void operator()(const tbb::blocked_range<int> & idx_range) const
    {
        for (int i=idx_range.begin(); i != idx_range.end(); ++i)
        {
            if(stage == 0) this->computeVertexData(face_indices[i]);
            else this->computeEdgeData(edge_indices[i]);
        }
    }
};

#endif /* CreateExtendedTriangleInfos_hpp */
