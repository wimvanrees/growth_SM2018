//
//  Parametrizer.hpp
//  Elasticity
//
//  Created by Wim van Rees on 8/17/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef Parametrizer_hpp
#define Parametrizer_hpp

#include "common.hpp"

/**
 * @brief Abstract class that wraps around a set of data which is the input for a minimization routine
 * @tparam tMesh Type of the mesh
 *
 * @details
 * Base class for parametrization of a set of variables that will be used to minimize a cost function. 
 *
 * @see HLBFGS_Energy_Parametrized
 */
template<typename tMesh>
class Parametrizer
{
protected:
    tMesh & mesh;
    Real * data;
public:
    Parametrizer(tMesh & mesh_in):
    mesh(mesh_in),
    data(nullptr)
    {}

    /**
     * Initialize the solution / unknown vector and allocate if necessary
     *
     * Here we will initialize the data pointer by allocating the corresponding data and filling it with an initial guess. 
     *
     * @param [in] initVals a vector with the initial values for the unknowns
     * @param [in] reAllocate whether we need to reallocate the array or not (if we reuse the parametrizer between different optimizations we might not want to reallocate)
     */
    virtual void initSolution(const Eigen::Ref<const Eigen::VectorXd> initVals, const bool reAllocate = false)
    {
        if((data==nullptr) or reAllocate)
        {
            if(data!=nullptr)
            {
                delete [] data;
                data = nullptr;
            }
            const int nVariables = getNumberOfVariables();
            data = new Real[nVariables];
            assert(initVals.rows() == nVariables);
            for(int i=0;i<nVariables;++i)
                data[i] = initVals(i);
        }
        // also update the mesh
        updateSolution();
    }
    
    virtual int getNumberOfVariables() const = 0;
    virtual Real * getDataPointer() const
    {
        return data;
    }

    
    virtual void updateSolution() = 0;
    virtual void updateGradient(const int nVars, const Eigen::Ref<const Eigen::VectorXd> energyGradient, Real * const grad_ptr) = 0;
    virtual Real computeEnergyContribution() {return 0;}
    
    virtual ~Parametrizer()
    {
        if(data!=nullptr)
            delete [] data;
        data = nullptr;
    }
};



#endif /* Parametrizer_hpp */
