//
//  MergePerFaceQuantities.hpp
//  Elasticity
//
//  Created by Wim van Rees on 26/03/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef MergePerFaceQuantities_hpp
#define MergePerFaceQuantities_hpp

#include "common.hpp"
#include "TopologyData.hpp"
#include "BoundaryConditionsData.hpp"
#include "ExtendedTriangleInfo.hpp"

#include "tbb/blocked_range.h"

struct MergeGradVertices
{
    const TopologyData & topo;
    const BoundaryConditionsData & boundaryConditions;
    const Eigen::Ref<const Eigen::MatrixXd> per_face_gradients;
    Eigen::Ref<Eigen::MatrixXd> gradient_vertices;
    
    MergeGradVertices(const TopologyData & topo, const BoundaryConditionsData & boundaryConditions, const Eigen::Ref<const Eigen::MatrixXd> pfg, Eigen::Ref<Eigen::MatrixXd> g):
    topo(topo),
    boundaryConditions(boundaryConditions),
    per_face_gradients(pfg),
    gradient_vertices(g)
    {}
    
    MergeGradVertices(const MergeGradVertices & c, tbb::split):
    topo(c.topo),
    boundaryConditions(c.boundaryConditions),
    per_face_gradients(c.per_face_gradients),
    gradient_vertices(c.gradient_vertices)
    {}
    
    void join(const MergeGradVertices & )
    {
//        nothing here : but it has to be here to allow parallel_reduce, which in turn has to be to allow a non-const operator() (we technically are const but we change a const & in the form of const Eigen::Ref<const Eigen::...> ). Same for MergeGradEdges below.
    }
    
    void operator () (const tbb::blocked_range<int>& vertex_range)
    {
        const auto face2edges = topo.getFace2Edges();
        const auto edge2faces = topo.getEdge2Faces();
        const auto face2vertices = topo.getFace2Vertices();
        const auto vertices_bc = boundaryConditions.getVertexBoundaryConditions();
        const auto & vertex2faces = topo.getVertex2Faces();
        
        for (int i=vertex_range.begin(); i != vertex_range.end(); ++i)
        {
            const std::vector<int> & faces = vertex2faces[i];
            
            if(vertices_bc(i,0) && vertices_bc(i,1) && vertices_bc(i,2)) continue;
            
            for(size_t f=0;f<faces.size();++f)
            {
                const int face_idx = faces[f];
                
                // figure out which vertex I am
                const int idx_v_rel = face2vertices(face_idx,0) == i ? 0 : (face2vertices(face_idx,1) == i ? 1 : 2);
                
                for(int j=0;j<3;++j)
                    if(not vertices_bc(i,j)) gradient_vertices(i,j) += per_face_gradients(face_idx, 3*idx_v_rel+j);
                
                // deal with the opposite faces
                // 1. get the opposite edge to my vertex
                const int idx_opp_edge_rel = (idx_v_rel+1)%3; // v0 --> e1, v1 --> e2, v2 --> e0
                const int idx_opp_edge = face2edges(face_idx, idx_opp_edge_rel);
                
                // 2. get the opposite face to that edge
                const int idx_opp_face_rel = (edge2faces(idx_opp_edge,0) == face_idx ? 1 : 0);
                const int idx_opp_face = edge2faces(idx_opp_edge, idx_opp_face_rel);
                
                // add it to our gradient
                if(idx_opp_face >= 0)
                {
                    // 3. get the index of the neighboring edge relative to the opposite face
                    const int idx_opp_edge_rel_other = (face2edges(idx_opp_face,0) == idx_opp_edge ? 0 : (face2edges(idx_opp_face,1) == idx_opp_edge ? 1 : 2) );
                    
                    for(int j=0;j<3;++j)
                        if(not vertices_bc(i,j)) gradient_vertices(i,j) += per_face_gradients(idx_opp_face, 3*(idx_opp_edge_rel_other+3)+j);
                }
            }
        }
    }
};

struct MergeGradEdges
{
    const TopologyData & topo;
    const BoundaryConditionsData & boundaryConditions;
    const Eigen::Ref<const Eigen::MatrixXd> per_face_gradients;
    Eigen::Ref<Eigen::VectorXd> gradient_edges;
    
    MergeGradEdges(const TopologyData & topo, const BoundaryConditionsData & boundaryConditions, const Eigen::Ref<const Eigen::MatrixXd> pfg, Eigen::Ref<Eigen::VectorXd> g):
    topo(topo),
    boundaryConditions(boundaryConditions),
    per_face_gradients(pfg),
    gradient_edges(g)
    {}
    
    MergeGradEdges(const MergeGradEdges & c, tbb::split):
    topo(c.topo),
    boundaryConditions(c.boundaryConditions),
    per_face_gradients(c.per_face_gradients),
    gradient_edges(c.gradient_edges)
    {}
    
    void join(const MergeGradEdges & )
    {
        //        nothing here
    }
    
    void operator () (const tbb::blocked_range<int>& edge_range)
    {
        const auto face2edges = topo.getFace2Edges();
        const auto edge2faces = topo.getEdge2Faces();
        const auto edges_bc = boundaryConditions.getEdgeBoundaryConditions();
        
        for (int i=edge_range.begin(); i != edge_range.end(); ++i)
        {
            const int face_idx_0 = edge2faces(i,0);
            const int face_idx_1 = edge2faces(i,1);
            
            if(edges_bc(i)) continue;
            
            if(face_idx_0 >= 0)
            {
                const int edge_idx_rel = face2edges(face_idx_0,0) == i ? 0 : (face2edges(face_idx_0,1) == i ? 1 : 2);
                gradient_edges(i) += per_face_gradients(face_idx_0, 3*6 + edge_idx_rel);
            }
            if(face_idx_1 >= 0)
            {
                const int edge_idx_rel = face2edges(face_idx_1,0) == i ? 0 : (face2edges(face_idx_1,1) == i ? 1 : 2);
                gradient_edges(i) += per_face_gradients(face_idx_1, 3*6 + edge_idx_rel);
            }
        }
    }
};


struct ScalarGradientHelper
{
    const bool withGradient;
    const Real avgVal;
    Eigen::Map<Eigen::VectorXd> gradient;
    
    ScalarGradientHelper(const bool withGradient_in, const Real avgVal_in):
    withGradient(withGradient_in),
    avgVal(avgVal_in),
    gradient(nullptr,0)
    {
    }
};


struct MergeHessian
{
    const TopologyData & topo;
    const BoundaryConditionsData & boundaryConditions;
    const std::vector<ExtendedTriangleInfo> & vInfo;
    const tVecMat2d & aforms;
    const Eigen::Ref<const Eigen::MatrixXd> gradient_vertices; // the energy gradient per-vertex (all faces contributions summed together)
    const Eigen::Ref<const Eigen::VectorXd> gradient_edges; // the energy gradient per-edge (all faces contributions summed together)
    const Eigen::Ref<const Eigen::MatrixXd> per_face_hessian; // the energy hessian per-face, separated between all individual components (all vertices, edges, opposite vertices, etc)
    const Eigen::Ref<const Eigen::VectorXd> per_vertex_areas;
    const Eigen::Ref<const Eigen::VectorXd> per_edge_areas;
    
    Eigen::Ref<Eigen::MatrixXd> hessian_faces_abar; // the accumulation of all gradient contributions wrt abar
    ScalarGradientHelper & grad_h;
    ScalarGradientHelper & grad_Y;
    const Real energy_before_avg;
    
    MergeHessian(const TopologyData & topo, const BoundaryConditionsData & boundaryConditions, const std::vector<ExtendedTriangleInfo> & vInfo, const tVecMat2d & aforms, const Eigen::Ref<const Eigen::MatrixXd> gv, const Eigen::Ref<const Eigen::VectorXd> ge, const Eigen::Ref<const Eigen::MatrixXd> pfh, const Eigen::Ref<const Eigen::VectorXd> pva, const Eigen::Ref<const Eigen::VectorXd> pea, Eigen::Ref<Eigen::MatrixXd> hf_abar, ScalarGradientHelper & grad_h, ScalarGradientHelper & grad_Y, const Real energy_before_avg):
    topo(topo),
    boundaryConditions(boundaryConditions),
    vInfo(vInfo),
    aforms(aforms),
    gradient_vertices(gv),
    gradient_edges(ge),
    per_face_hessian(pfh),
    per_vertex_areas(pva),
    per_edge_areas(pea),
    hessian_faces_abar(hf_abar),
    grad_h(grad_h),
    grad_Y(grad_Y),
    energy_before_avg(energy_before_avg)
    {}
    
    MergeHessian(const MergeHessian & c, tbb::split):
    topo(c.topo),
    boundaryConditions(c.boundaryConditions),
    vInfo(c.vInfo),
    aforms(c.aforms),
    gradient_vertices(c.gradient_vertices),
    gradient_edges(c.gradient_edges),
    per_face_hessian(c.per_face_hessian),
    per_vertex_areas(c.per_vertex_areas),
    per_edge_areas(c.per_edge_areas),
    hessian_faces_abar(c.hessian_faces_abar),
    grad_h(c.grad_h),
    grad_Y(c.grad_Y),
    energy_before_avg(c.energy_before_avg)
    {}
    
    void join(const MergeHessian & )
    {
        //        nothing here : but it has to be here to allow parallel_reduce, which in turn has to be to allow a non-const operator() (we technically are const but we change a const & in the form of const Eigen::Ref<const Eigen::...> ).
    }
    
    void operator () (const tbb::blocked_range<int>& face_range)
    {
        // back to the per-face approach
        const auto face2edges = topo.getFace2Edges();
        const auto edge2faces = topo.getEdge2Faces();
        const auto face2vertices = topo.getFace2Vertices();
        
        const int nDims = 21;

        // one or the other for now --> note: grad_h not implemented yet!!
        const Real hY_prefac = grad_h.withGradient ? 1.0 / grad_h.avgVal : (grad_Y.withGradient ? 1.0 / grad_Y.avgVal : 1.0);
        
        for (int i=face_range.begin(); i != face_range.end(); ++i)
        {
            const ExtendedTriangleInfo & info = vInfo[i];
            
            int idx_v_other_e0 = -1;
            int idx_v_other_e1 = -1;
            int idx_v_other_e2 = -1;
            
            if(info.other_faces[0] != nullptr)
            {
                const int other_face_idx = info.other_faces[0]->face_idx;
                
                const int startidx = (face2edges(other_face_idx,0) == info.idx_e0 ? 0 : (face2edges(other_face_idx,1) == info.idx_e0 ? 1 : 2) );
                const int idx_other_v2 = face2vertices(other_face_idx,(startidx+2)%3);
                idx_v_other_e0 = idx_other_v2;
            }
            
            if(info.other_faces[1] != nullptr)
            {
                const int other_face_idx = info.other_faces[1]->face_idx;
                
                const int startidx = (face2edges(other_face_idx,0) == info.idx_e1 ? 0 : (face2edges(other_face_idx,1) == info.idx_e1 ? 1 : 2) );
                const int idx_other_v0 = face2vertices(other_face_idx,(startidx+2)%3);
                idx_v_other_e1 = idx_other_v0;
                
            }
            
            if(info.other_faces[2] != nullptr)
            {
                const int other_face_idx = info.other_faces[2]->face_idx;
                
                const int startidx = (face2edges(other_face_idx,0) == info.idx_e2 ? 0 : (face2edges(other_face_idx,1) == info.idx_e2 ? 1 : 2) );
                const int idx_other_v1 = face2vertices(other_face_idx,(startidx+2)%3);
                idx_v_other_e2 = idx_other_v1;
            }
            
            // compute the gradient of the area
            const Real det_abar = std::sqrt(aforms[i].determinant()); // twice the area
            const Real grad_det_11 = aforms[i](1,1);
            const Real grad_det_12 = - aforms[i](0,1) - aforms[i](1,0);
            const Real grad_det_22 = aforms[i](0,0);
            const Eigen::Vector3d grad_area = 0.5 * (Eigen::Vector3d() << 0.5/det_abar * grad_det_11, 0.5/det_abar * grad_det_12, 0.5/det_abar * grad_det_22).finished();
            
            // now we can accumulate the different terms
            for(int j=0;j<3;++j) // loop over a11, a12, a22
            {
                for(int d=0;d<3;++d)
                {
                    // v0, v1, v2
                    hessian_faces_abar(i,j) += hY_prefac * (2.0 * gradient_vertices(info.idx_v0, d) * per_face_hessian(i, nDims*j + 3*0 + d) / per_vertex_areas(info.idx_v0) - std::pow(gradient_vertices(info.idx_v0, d) / per_vertex_areas(info.idx_v0),2) * grad_area(j) / 3.0); // onethird area
                    hessian_faces_abar(i,j) += hY_prefac * (2.0 * gradient_vertices(info.idx_v1, d) * per_face_hessian(i, nDims*j + 3*1 + d) / per_vertex_areas(info.idx_v1) - std::pow(gradient_vertices(info.idx_v1, d) / per_vertex_areas(info.idx_v1),2) * grad_area(j) / 3.0);
                    hessian_faces_abar(i,j) += hY_prefac * (2.0 * gradient_vertices(info.idx_v2, d) * per_face_hessian(i, nDims*j + 3*2 + d) / per_vertex_areas(info.idx_v2) - std::pow(gradient_vertices(info.idx_v2, d) / per_vertex_areas(info.idx_v2),2) * grad_area(j) / 3.0);
                    
                    // v_other_e0, v_other_e1, v_other_e2 (those areas are independent of our aform - no second term there)
                    if(idx_v_other_e0 >= 0) hessian_faces_abar(i,j) += hY_prefac * (2.0 * gradient_vertices(idx_v_other_e0, d) * per_face_hessian(i, nDims*j + 3*3 + d) / per_vertex_areas(idx_v_other_e0));
                    if(idx_v_other_e1 >= 0) hessian_faces_abar(i,j) += hY_prefac * (2.0 * gradient_vertices(idx_v_other_e1, d) * per_face_hessian(i, nDims*j + 3*4 + d) / per_vertex_areas(idx_v_other_e1));
                    if(idx_v_other_e2 >= 0) hessian_faces_abar(i,j) += hY_prefac * (2.0 * gradient_vertices(idx_v_other_e2, d) * per_face_hessian(i, nDims*j + 3*5 + d) / per_vertex_areas(idx_v_other_e2));
                }
                
                // edge e0, e1, e2
                hessian_faces_abar(i,j) += hY_prefac * (2.0 * gradient_edges(info.idx_e0) * per_face_hessian(i, nDims*j + 3*6 + 0) / per_edge_areas(info.idx_e0) - std::pow(gradient_edges(info.idx_e0) / per_edge_areas(info.idx_e0), 2) * grad_area(j) / 2.0); // onehalf area
                hessian_faces_abar(i,j) += hY_prefac * (2.0 * gradient_edges(info.idx_e1) * per_face_hessian(i, nDims*j + 3*6 + 1) / per_edge_areas(info.idx_e1) - std::pow(gradient_edges(info.idx_e1) / per_edge_areas(info.idx_e1), 2) * grad_area(j) / 2.0);
                hessian_faces_abar(i,j) += hY_prefac * (2.0 * gradient_edges(info.idx_e2) * per_face_hessian(i, nDims*j + 3*6 + 2) / per_edge_areas(info.idx_e2) - std::pow(gradient_edges(info.idx_e2) / per_edge_areas(info.idx_e2), 2) * grad_area(j) / 2.0);
            }
            
            if( (not grad_h.withGradient) and (not grad_Y.withGradient) ) continue;
            
            // NOT grad_h.withGradient implemented yet!!!!
            // only grad_Y.withGradient
            
            Real tmpSum = 0.0;
            
            for(int d=0;d<3;++d)
            {
                // v0, v1, v2
                tmpSum += gradient_vertices(info.idx_v0, d) * per_face_hessian(i, nDims*3 + 3*0 + d) / per_vertex_areas(info.idx_v0);
                tmpSum += gradient_vertices(info.idx_v1, d) * per_face_hessian(i, nDims*3 + 3*1 + d) / per_vertex_areas(info.idx_v1);
                tmpSum += gradient_vertices(info.idx_v2, d) * per_face_hessian(i, nDims*3 + 3*2 + d) / per_vertex_areas(info.idx_v2);
                
                // v_other_e0, v_other_e1, v_other_e2 (those areas are independent of our aform - no second term there)
                if(idx_v_other_e0 >= 0) tmpSum += gradient_vertices(idx_v_other_e0, d) * per_face_hessian(i, nDims*3 + 3*3 + d) / per_vertex_areas(idx_v_other_e0);
                if(idx_v_other_e1 >= 0) tmpSum += gradient_vertices(idx_v_other_e1, d) * per_face_hessian(i, nDims*3 + 3*4 + d) / per_vertex_areas(idx_v_other_e1);
                if(idx_v_other_e2 >= 0) tmpSum += gradient_vertices(idx_v_other_e2, d) * per_face_hessian(i, nDims*3 + 3*5 + d) / per_vertex_areas(idx_v_other_e2);
            }
            
            // edge e0, e1, e2
            tmpSum += gradient_edges(info.idx_e0) * per_face_hessian(i, nDims*3 + 3*6 + 0) / per_edge_areas(info.idx_e0);
            tmpSum += gradient_edges(info.idx_e1) * per_face_hessian(i, nDims*3 + 3*6 + 1) / per_edge_areas(info.idx_e1);
            tmpSum += gradient_edges(info.idx_e2) * per_face_hessian(i, nDims*3 + 3*6 + 2) / per_edge_areas(info.idx_e2);

            grad_Y.gradient(i) += 2.0 * hY_prefac * tmpSum;
            
            // global term
            const Real h_postfac = -1.0 / (grad_Y.avgVal*grad_Y.avgVal) * (0.5*vInfo[i].double_face_area);
            grad_Y.gradient(i) += energy_before_avg * h_postfac;
            
            // derivative of h_prefac = 1/avgThickness wrt current face thickness :
            //
            //            for(int d=0;d<3;++d)
            //            {
            //                // v0, v1, v2
            //                hessian_faces_h(i) += h_prefac * (2.0 * gradient_vertices(idx_v0, d) * per_face_hessian(i, nDims*3 + 3*0 + d) / per_vertex_areas(idx_v0));
            //                hessian_faces_h(i) += h_prefac * (2.0 * gradient_vertices(idx_v1, d) * per_face_hessian(i, nDims*3 + 3*1 + d) / per_vertex_areas(idx_v1));
            //                hessian_faces_h(i) += h_prefac * (2.0 * gradient_vertices(idx_v2, d) * per_face_hessian(i, nDims*3 + 3*2 + d) / per_vertex_areas(idx_v2));
            //
            //                // v_other_e0, v_other_e1, v_other_e2
            //                if(idx_v_other_e0 >= 0) hessian_faces_h(i) += h_prefac * (2.0 * gradient_vertices(idx_v_other_e0, d) * per_face_hessian(i, nDims*3 + 3*3 + d) / per_vertex_areas(idx_v_other_e0));
            //                if(idx_v_other_e1 >= 0) hessian_faces_h(i) += h_prefac * (2.0 * gradient_vertices(idx_v_other_e1, d) * per_face_hessian(i, nDims*3 + 3*4 + d) / per_vertex_areas(idx_v_other_e1));
            //                if(idx_v_other_e2 >= 0) hessian_faces_h(i) += h_prefac * (2.0 * gradient_vertices(idx_v_other_e2, d) * per_face_hessian(i, nDims*3 + 3*5 + d) / per_vertex_areas(idx_v_other_e2));
            //            }
            //
            //            // edge e0, e1, e2
            //            hessian_faces_h(i) += h_prefac * (2.0 * gradient_edges(idx_e0) * per_face_hessian(i, nDims*3 + 3*6 + 0) / per_edge_areas(idx_e0));
            //            hessian_faces_h(i) += h_prefac * (2.0 * gradient_edges(idx_e1) * per_face_hessian(i, nDims*3 + 3*6 + 1) / per_edge_areas(idx_e1));
            //            hessian_faces_h(i) += h_prefac * (2.0 * gradient_edges(idx_e2) * per_face_hessian(i, nDims*3 + 3*6 + 2) / per_edge_areas(idx_e2));
            //
            
            // global term
            //const Real h_postfac = -1.0 / (avgThickness*avgThickness) * (0.5*double_face_areas(i));
            //hessian_faces_h(i) += energy_before_avg * h_postfac;
        }
    }
};

#endif /* MergePerFaceQuantities_hpp */
