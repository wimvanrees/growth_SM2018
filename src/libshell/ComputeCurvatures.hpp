//
//  ComputeCurvatures.hpp
//  Elasticity
//
//  Created by Wim van Rees on 5/20/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef ComputeCurvatures_hpp
#define ComputeCurvatures_hpp

#include "common.hpp"

/**
 * @brief Class to compute the gauss and mean curvatures of the current configuration of a mesh
 * @tparam tMesh Type of the mesh 
 *
 * @details
 * This class is stateless. It provides a compute() method that computes the gaussian and mean curvature, as the determinant and 1/2 times the trace of the shape operator
 */

template<typename tMesh>
class ComputeCurvatures
{
protected:
    
public:
        
    ComputeCurvatures()
    {}
    
    /**
     * Computes the curvature of the current configuration in the mesh
     *
     * We compute the shape operator as \f$S = a^{-1} b\f$. Then we can
     * compute the gauss curvature as \f$K = \det S\f$ and the mean curvature
     * as \f$H = \frac{1}{2} \mathrm{Tr} S\f$.
     
     * @param [in] mesh the mesh from which to take the current configuration
     * @param [out] gauss reference to a vector in which the gauss curvature will be stored (resized if necessary)
     * @param [out] mean reference to a vector in which the mean curvature will be stored (resized if necessary)
     */
    void compute(const tMesh & mesh, Eigen::VectorXd & gauss, Eigen::VectorXd & mean) const;
};


#endif /* ComputeCurvatures_hpp */
