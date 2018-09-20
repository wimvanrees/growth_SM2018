//
//  EnergyOperatorList.hpp
//  Elasticity
//
//  Created by Wim van Rees on 27/02/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef EnergyOperatorList_h
#define EnergyOperatorList_h

#include "EnergyOperator.hpp"
#include "common.hpp"

/*! \class EnergyOperatorList
 * \brief List of energy operators
 *
 * In case multiple energy operators are used, this class can take them together and create a single method to compute all energies, gradients and/or hessians. Based on the composite design pattern
 */

template<typename tMesh>
class EnergyOperatorList : public EnergyOperator<tMesh>
{
protected:
    std::vector<EnergyOperator<tMesh>*> operators;
    
public:
    
    EnergyOperatorList(std::initializer_list<EnergyOperator<tMesh>*> init_list):
    EnergyOperator<tMesh>()
    {
        operators.clear();
        for( auto elem : init_list )
            operators.push_back(elem);
    }
    
    
    Real compute(const tMesh & mesh) const override
    {
        Real energy = 0.0;
        for(size_t i=0;i<operators.size();++i)
            energy += operators[i]->compute(mesh);
        
        return energy;
    }
    
    Real compute(const tMesh & mesh, Eigen::Ref<Eigen::VectorXd> gradient) const override
    {
        Real energy = 0.0;
        
        for(size_t i=0;i<operators.size();++i)
            energy += operators[i]->compute(mesh, gradient);
        
        return energy;
    }
    
    Real computeSubsetEnergies(const tMesh & mesh, const std::vector<int> & face_indices) const override
    {
        Real energy = 0.0;
        
        for(size_t i=0;i<operators.size();++i)
            energy += operators[i]->computeSubsetEnergies(mesh, face_indices);
        
        return energy;        
    }
    
    bool isSubsetImplemented() const override
    {
        bool retval = true;
        for(size_t i=0;i<operators.size();++i)
            retval = (retval && operators[i]->isSubsetImplemented());
        return retval;
    }
    
    virtual void printProfilerSummary() const override
    {
        for(size_t i=0;i<operators.size();++i)
            operators[i]->printProfilerSummary();
    }
        
    virtual void addEnergy(std::vector<std::pair<std::string, Real>> & out) const override
    {
        for(size_t i=0;i<operators.size();++i)
            operators[i]->addEnergy(out);
    }
    
    int getNumberOfVariables(const tMesh & mesh) const override
    {
        int retval = operators[0]->getNumberOfVariables(mesh);
        for(size_t i=1;i<operators.size();++i)
            retval = std::max(retval, operators[i]->getNumberOfVariables(mesh));
        return retval;
    }
};


#endif /* EnergyOperatorList_h */
