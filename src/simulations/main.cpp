//
//  main.cpp
//  Elasticity
//
//  Created by Wim van Rees on 2/15/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#include "common.hpp"
#include "ArgumentParser.hpp"
#include "Sim.hpp"

#include "Sim_IOops.hpp"
#include "Sim_Bilayer_4DFilaments.hpp"
#include "Sim_MonoLayer_Growth.hpp"


#ifdef USETBB
#include "tbb/task_scheduler_init.h"
#endif

int main(int argc,  const char ** argv)
{
#ifdef USETBB
    const int num_threads = tbb::task_scheduler_init::default_num_threads();
    tbb::task_scheduler_init init(num_threads);
    std::cout << "Starting TBB with " << num_threads << " threads " << std::endl;
#endif
    
    BaseSim * sim = nullptr;
    
    ArgumentParser parser(argc,argv);

    const std::string simCase = parser.parse<std::string>("-sim", "");

    if(simCase == "IO")
        sim = new Sim_IOops(parser);
    else if(simCase == "bilayer_4dfilaments")
        sim = new Sim_Bilayer_4DFilaments(parser);
    else if(simCase == "monolayer_growth")
        sim = new Sim_MonoLayer_Growth(parser);
    else
    {
        std::cout << "No valid sim case defined. Options are \n";
        std::cout << "\t -sim IO\n";
        std::cout << "\t -sim bilayer_4dfilaments\n";
        std::cout << "\t -sim monolayer_growth\n";

        helpers::catastrophe("sim case does not exist",__FILE__,__LINE__);
    }
    
    sim->init();
    sim->run();
    
    delete sim;

    parser.finalize();
    
    return 0;
}
