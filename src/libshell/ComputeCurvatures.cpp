//
//  ComputeCurvatures.cpp
//  Elasticity
//
//  Created by Wim van Rees on 5/20/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#include "ComputeCurvatures.hpp"
#include "ExtendedTriangleInfo.hpp"

template<typename tMesh>
void ComputeCurvatures<tMesh>::compute(const tMesh & mesh, Eigen::VectorXd & gauss, Eigen::VectorXd & mean) const
{
    // get current state quantities
    const int nFaces = mesh.getNumberOfFaces();
    const auto & currentState = mesh.getCurrentConfiguration();
    
    for(int i=0;i<nFaces;++i)
    {
        const ExtendedTriangleInfo & info = currentState.getTriangleInfo(i);
        
        const Eigen::Matrix2d aform = info.computeFirstFundamentalForm();
        const Eigen::Matrix2d bform = info.computeSecondFundamentalForm();
        
        const Eigen::Matrix2d shapeOp = aform.inverse() * bform;
        
        const Real gauss_curv = shapeOp.determinant();
        const Real mean_curv = 0.5*shapeOp.trace();
        
        gauss(i) = gauss_curv;
        mean(i) = mean_curv;
    }
}


#include "Mesh.hpp"
template class ComputeCurvatures<Mesh>;
template class ComputeCurvatures<BilayerMesh>;

