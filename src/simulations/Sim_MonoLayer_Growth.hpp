//
//  Sim_MonoLayer_Growth.hpp
//  Elasticity
//
//  Created by Wim van Rees on 8/17/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef Sim_MonoLayer_Growth_hpp
#define Sim_MonoLayer_Growth_hpp

#include "Sim.hpp"
#include "Mesh.hpp"
#include "GrowthHelper.hpp"

class Sim_MonoLayer_Growth : public Sim<Mesh>
{
    typedef Mesh tMesh;
protected:
    
    void runCone();
    void runEdgeGrowth();
    void run_basic_disk();
    
    void dumpWithGrowthRates(const Eigen::Ref<const Eigen::VectorXd> growthRates, const std::string filename, const bool addCurvature = false);
    
    void dumpOrthoNew(const std::vector<GrowthState> & growth, const std::string filename, const bool restConfig = false);
    void assignGrowthToMetric(const std::vector<GrowthState> & growth, const Real t, const bool interp_logeucl);

public:
    
    Sim_MonoLayer_Growth(ArgumentParser & parser):
    Sim<tMesh>(parser)
    {}
    
    void init() override;
    void run() override;


};

#endif /* Sim_MonoLayer_Growth_hpp */
