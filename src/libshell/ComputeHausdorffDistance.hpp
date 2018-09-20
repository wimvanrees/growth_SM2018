//
//  ComputeHausdorffDistance.hpp
//  Elasticity
//
//  Created by Wim van Rees on 12/28/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef ComputeHausdorffDistance_hpp
#define ComputeHausdorffDistance_hpp

#include "common.hpp"

class ComputeHausdorffDistance
{
protected:
    const bool bDump;
    
public:
        
    ComputeHausdorffDistance(const bool bDump = false):
    bDump(bDump)
    {}

    Real compute(const Eigen::Ref<const Eigen::MatrixXd> vertices_A, const Eigen::Ref<const Eigen::MatrixXd> vertices_B, const Eigen::Ref<const Eigen::MatrixXi> faces, const Real rescale = false) const;
};

#endif /* ComputeHausdorffDistance_hpp */
