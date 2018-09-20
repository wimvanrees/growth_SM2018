//
//  Sim_IOops.hpp
//  Elasticity
//
//  Created by Wim van Rees on 12/19/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef Sim_IOops_hpp
#define Sim_IOops_hpp

#include "Sim.hpp"
#include "Mesh.hpp"


class Sim_IOops: public Sim<Mesh>
{
    typedef Mesh tMesh;
    
protected:
    
    void convertSTL(const std::string & fname_out);
    void convertOBJ(const std::string & fname_out);
    void convertCSV(const std::string & fname_in, const std::string & fname_out);
    void convertINP(const std::string & fname_out);

public:
    
    Sim_IOops(ArgumentParser & parser):
    Sim<tMesh>(parser)
    {}
    
    void init() override;
    void run() override;
};

#endif /* Sim_IOops_hpp */
