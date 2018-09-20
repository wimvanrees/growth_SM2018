//
//  EnergyOperator.hpp
//  Elasticity
//
//  Created by Wim van Rees on 2/24/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef EnergyOperator_h
#define EnergyOperator_h

#include "common.hpp"
#include "Profiler.hpp"

/*! \class EnergyOperator
 * \brief Base class for energy-based operators.
 *
 * This class should implement here the energy, energy-gradient (-force vector), and energy hessian (-force jacobian) for each of the energy metrics that derives from here.
 */

template<typename tMesh>
class EnergyOperator
{
protected:
    mutable Profiler profiler;
    
public:
    
    virtual Real compute(const tMesh & mesh) const = 0;
    virtual Real compute(const tMesh & mesh, Eigen::Ref<Eigen::VectorXd> gradient) const = 0;
    virtual Real computeSubsetEnergies(const tMesh & /*mesh*/, const std::vector<int> & /*face_indices*/) const
    {return 0;}
    virtual bool isSubsetImplemented() const {return false;}
    
    virtual void printProfilerSummary() const
    {
        profiler.printSummary();
    }
    
    virtual void addEnergy(std::vector<std::pair<std::string, Real>> & ) const
    {        
    }
    
    virtual int getNumberOfVariables(const tMesh & mesh) const = 0;
    
    virtual ~EnergyOperator(){}
};

template<typename tMesh>
class EnergyOperator_DKS : public EnergyOperator<tMesh> // vertex only
{
protected:
    
    virtual Real computeAll(const tMesh & mesh, Eigen::Ref<Eigen::MatrixXd> gradient_vertices, const bool computeGradient) const = 0;
    
public:
    
    virtual Real compute(const tMesh & mesh, Eigen::Ref<Eigen::VectorXd> gradient) const override
    {
        assert(mesh.getNumberOfVertices() > 0);
        assert(gradient.rows() >= 3*mesh.getNumberOfVertices());
     
        const int nVertices = mesh.getNumberOfVertices();
        
        // map the gradient vector into a column-major matrix
        // memory layout is [v0x, v1x, .. ,vNx, v0y, v1y, .. ,vNy, v0z, v1z, .. ,vNz]^T
        Eigen::Map<Eigen::MatrixXd> grad_vertices(gradient.data(), nVertices, 3);
        const Real energy = computeAll(mesh, grad_vertices, true);
        return energy;
    }
    
    Real compute(const tMesh & mesh) const override
    {
        assert(mesh.getNumberOfVertices() > 0);
        
        Eigen::MatrixXd dummy1;
        const Real energy = computeAll(mesh, dummy1, false);
        return energy;
    }
    
    virtual Real compute(const tMesh & mesh, Eigen::Ref<Eigen::MatrixXd> gradient_vertices) const
    {
        assert(gradient_vertices.rows() >= mesh.getNumberOfVertices());
        assert(gradient_vertices.cols() == 3);
        assert(mesh.getNumberOfVertices() > 0);
        
        const Real energy = computeAll(mesh, gradient_vertices, true);
        return energy;
    }
    
    virtual int getNumberOfVariables(const tMesh & mesh) const override
    {
        return 3*mesh.getNumberOfVertices();
    }
};


template<typename tMesh>
class EnergyOperator_DCS : public EnergyOperator<tMesh> //vertex and edge directors
{
protected:
    
    virtual Real computeAll(const tMesh & mesh, Eigen::Ref<Eigen::MatrixXd> gradient_vertices, Eigen::Ref<Eigen::VectorXd> gradient_edges, const bool computeGradient) const = 0;
    
public:
    
    virtual Real compute(const tMesh & mesh, Eigen::Ref<Eigen::VectorXd> gradient) const override
    {
        assert(mesh.getNumberOfVertices() > 0);
        assert(gradient.rows() >= 3*mesh.getNumberOfVertices() + mesh.getNumberOfEdges());
        
        const int nVertices = mesh.getNumberOfVertices();
        const int nEdges = mesh.getNumberOfEdges();
        
        // map the gradient vector into two separate parts: a column-major matrix and a vector
        // memory layout is [v0x, v1x, .. ,vNx, v0y, v1y, .. ,vNy, v0z, v1z, .. ,vNz, e0, e1, .. ,eN]^T
        Eigen::Map<Eigen::MatrixXd> grad_vertices(gradient.data(), nVertices, 3);
        Eigen::Map<Eigen::VectorXd> grad_edges(gradient.data() + 3*nVertices, nEdges);
        const Real energy = computeAll(mesh, grad_vertices, grad_edges, true);
        return energy;
    }
    
    Real compute(const tMesh & mesh) const override
    {
        assert(mesh.getNumberOfVertices() > 0);
        
        Eigen::MatrixXd dummy1;
        Eigen::VectorXd dummy2;
        const Real energy = computeAll(mesh, dummy1, dummy2, false);
        return energy;
    }
    
    virtual Real compute(const tMesh & mesh, Eigen::Ref<Eigen::MatrixXd> gradient_vertices, Eigen::Ref<Eigen::VectorXd> gradient_edges) const
    {
        assert(mesh.getNumberOfVertices() > 0);
        assert(gradient_vertices.rows() >= mesh.getNumberOfVertices());
        assert(gradient_vertices.cols() == 3);
        assert(gradient_edges.rows() >= mesh.getNumberOfEdges());        
        
        const Real energy = computeAll(mesh, gradient_vertices, gradient_edges, true);
        return energy;
    }
    
    virtual int getNumberOfVariables(const tMesh & mesh) const override
    {
        return 3*mesh.getNumberOfVertices() + mesh.getNumberOfEdges();
    }
};
#endif /* EnergyOperator_h */
