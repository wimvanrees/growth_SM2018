//
//  Derivative_FD.hpp
//  Elasticity
//
//  Created by Wim van Rees on 4/28/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef Derivative_FD_hpp
#define Derivative_FD_hpp

#include "Mesh.hpp"
#include "MaterialProperties.hpp"
#include "common.hpp"

struct Derivative_FD_Base
{
    const int order;
    
    std::vector< std::vector <Real>> FD1_prefac; // coefficients for first derivative
    std::vector< std::vector <Real>> FD1_offset; // offsets for first derivative
    std::vector<Real> FD1_denum; // denumerator for first derivative
    
    std::vector< std::vector <Real>> FD2_prefac; // coefficients for second derivative
    std::vector< std::vector <Real>> FD2_offset; // offsets for second derivative
    std::vector<Real> FD2_denum; // denumerator for second derivative
    
    Derivative_FD_Base(const int order):
    order(order),
    FD1_prefac({ {1, -1}, {1, -8, 8, -1}, {-1, 9, -45, 45, -9, 1}, {3, -32, 168, -672, 672, -168, 32, -3} }),
    FD1_offset({ {-1, 1}, {-2, -1, 1, 2}, {-3, -2, -1, 1, 2, 3}, {-4, -3, -2, -1, 1, 2, 3, 4} }),
    FD1_denum({2, 12, 60, 840}),
    FD2_prefac({ {1,-2, 1}, {-1, 16, -30, 16, -1}, {2, -27, 270, -490, 270, -27, 2}, {-9, 128, -1008, 8064, -14350, 8064, -1008, 128, -9} }),
    FD2_offset({ {-1,0,1}, {-2, -1, 0, 1, 2}, {-3, -2, -1, 0, 1, 2, 3}, {-4, -3, -2, -1, 0, 1, 2, 3, 4} }),
    FD2_denum({1, 12, 180, 5040})
    {}
};


template<typename tMesh, typename tMeshOperator>
struct Derivative_FD : public Derivative_FD_Base
{
    tMesh & mesh;
    const tMeshOperator & meshOperator;
    
    const Real eps_vert;
    const Real eps_edge;
    
    Derivative_FD(const int order, tMesh & mesh, const tMeshOperator & meshOperator, const Real eps_vert, const Real eps_edge):
    Derivative_FD_Base(order),
    mesh(mesh),
    meshOperator(meshOperator),
    eps_vert(eps_vert),
    eps_edge(eps_edge)
    {}

    Derivative_FD(const int order, tMesh & mesh, const tMeshOperator & meshOperator, const Real eps):
    Derivative_FD(order, mesh, meshOperator, eps, eps)
    {}

    
    void computeGradient(Eigen::MatrixXd & gradient_vertices_approx, Eigen::VectorXd & gradient_edges_approx)
    {
        const int nPoints = FD1_prefac[order].size();
        assert(nPoints%2==0);
        
        const int nVertices = mesh.getNumberOfVertices();
        const int nEdges = mesh.getNumberOfEdges();

        gradient_vertices_approx.resize(nVertices,3);
        gradient_edges_approx.resize(nEdges);
        gradient_vertices_approx.setZero();
        gradient_edges_approx.setZero();
        
        const auto vertices_bc = mesh.getBoundaryConditions().getVertexBoundaryConditions();
        const auto edges_bc = mesh.getBoundaryConditions().getEdgeBoundaryConditions();

        auto vertices = mesh.getCurrentConfiguration().getVertices();
        auto edgeDirectors = mesh.getCurrentConfiguration().getEdgeDirectors();
        
        // wrt vertices
        for(int v=0;v<nVertices;++v)
            for(int d=0;d<3;++d)
            {
                gradient_vertices_approx(v,d) = 0.0;
                if(vertices_bc(v,d)) continue;
                
                for(int i=0;i<nPoints;++i)
                {
                    // perturb
                    vertices(v,d) += FD1_offset[order][i]*eps_vert;
                    mesh.updateDeformedConfiguration();
                    
                    gradient_vertices_approx(v,d) += FD1_prefac[order][i] * meshOperator.compute(mesh);
                    
                    // reset before next point
                    vertices(v,d) -= FD1_offset[order][i]*eps_vert;
                }
                
                gradient_vertices_approx(v,d) /= (eps_vert*FD1_denum[order]);
            }
        
        // wrt edge directors
        for(int e=0;e<nEdges;++e)
        {
            gradient_edges_approx(e) = 0.0;
            if(edges_bc(e)) continue;
            
            for(int i=0;i<nPoints;++i)
            {
                // perturb
                edgeDirectors(e) += FD1_offset[order][i]*eps_edge;
                mesh.updateDeformedConfiguration();
                
                gradient_edges_approx(e) += FD1_prefac[order][i] * meshOperator.compute(mesh);
                
                // reset before next point
                edgeDirectors(e) -= FD1_offset[order][i]*eps_edge;
            }
            // normalize
            gradient_edges_approx(e) /= (eps_edge*FD1_denum[order]);
        }
    }
};

template<typename tMesh, typename tMeshOperator>
struct Derivative_FD_Inverse : public Derivative_FD_Base
{
    tMesh & mesh;
    const tMeshOperator & meshOperator;
    tVecMat2d & vAbar;
    
    const Real eps;
    
    Derivative_FD_Inverse(const int order, tMesh & mesh, const tMeshOperator & meshOperator, tVecMat2d & vAbar, const Real eps):
    Derivative_FD_Base(order),
    mesh(mesh),
    meshOperator(meshOperator),
    vAbar(vAbar),
    eps(eps)
    {}

    void computeGradient(Eigen::VectorXd & gradient_abar_approx_vec)
    {
        const int nPoints = FD1_prefac[order].size();
        assert(nPoints%2==0);
        
        const int nFaces = mesh.getNumberOfFaces();
        
        gradient_abar_approx_vec.resize(3*nFaces);
        gradient_abar_approx_vec.setZero();
        
        Eigen::Map<Eigen::MatrixXd> gradient_abar_approx(gradient_abar_approx_vec.data(), nFaces, 3);
        
        // wrt abar
        for(int f=0;f<nFaces;++f)
            for(int d=0;d<3;++d) // a11, a12, a22
            {
                const int idx_1 = (d==2 ? 1 : 0);
                const int idx_2 = (d==0 ? 0 : 1);
                
                gradient_abar_approx(f,d) = 0.0;
                
                for(int i=0;i<nPoints;++i)
                {
                    // perturb
                    vAbar[f](idx_1, idx_2) += FD1_offset[order][i]*eps;
                    vAbar[f](1,0) = vAbar[f](0,1); // to be sure
                    
                    gradient_abar_approx(f,d) += FD1_prefac[order][i] * meshOperator.compute(mesh);
                    
                    // reset before next point
                    vAbar[f](idx_1, idx_2) -= FD1_offset[order][i]*eps;
                    vAbar[f](1,0) = vAbar[f](0,1); // to be sure
                }
                
                gradient_abar_approx(f,d) /= (eps*FD1_denum[order]);
            }
    }
    
    
    void computeGradient_thickness(MaterialProperties_Iso_Array & matprop, Eigen::VectorXd & gradient_h_approx)
    {
        const int nPoints = FD1_prefac[order].size();
        assert(nPoints%2==0);
        
        const int nFaces = mesh.getNumberOfFaces();
                
        // wrt abar
        for(int f=0;f<nFaces;++f)
        {
            gradient_h_approx(f) = 0.0;
            
            for(int i=0;i<nPoints;++i)
            {
                // perturb
                matprop.materials[f].thickness += FD1_offset[order][i]*eps;
                
                gradient_h_approx(f) += FD1_prefac[order][i] * meshOperator.compute(mesh);
                
                // reset before next point
                matprop.materials[f].thickness -= FD1_offset[order][i]*eps;
            }
            
            gradient_h_approx(f) /= (eps*FD1_denum[order]);
        }
    }
    
    void computeGradient_Young(MaterialProperties_Iso_Array & matprop, Eigen::VectorXd & gradient_Y_approx)
    {
        const int nPoints = FD1_prefac[order].size();
        assert(nPoints%2==0);
        
        const int nFaces = mesh.getNumberOfFaces();
        
        // wrt abar
        for(int f=0;f<nFaces;++f)
        {
            gradient_Y_approx(f) = 0.0;
            
            for(int i=0;i<nPoints;++i)
            {
                // perturb
                matprop.materials[f].Young += FD1_offset[order][i]*eps;
                
                gradient_Y_approx(f) += FD1_prefac[order][i] * meshOperator.compute(mesh);
                
                // reset before next point
                matprop.materials[f].Young -= FD1_offset[order][i]*eps;
            }
            
            gradient_Y_approx(f) /= (eps*FD1_denum[order]);
        }
    }
    
};

#endif /* Derivative_FD_hpp */
