//
//  Test_WriteBilayerSTL.hpp
//  Elasticity
//
//  Created by Wim van Rees on 7/29/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef Test_WriteBilayerSTL_hpp
#define Test_WriteBilayerSTL_hpp

#include "Sim.hpp"
#include "Mesh.hpp"

class Test_WriteBilayerSTL : public Sim<Mesh>
{
protected:
    
    void test_annulus();
    void test_twotriangles();
    void test_writeSingleLayerCurved();
    
public:
    
    void init() override;
    void run() override;
};
#endif /* Test_WriteBilayerSTL_hpp */
