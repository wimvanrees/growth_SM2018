//
//  WriteINP.hpp
//  Elasticity
//
//  Created by Wim M. van Rees on 7/18/18.
//  Copyright © 2018 Wim van Rees. All rights reserved.
//

#ifndef WriteINP_hpp
#define WriteINP_hpp

#include "TopologyData.hpp"
#include "DCSConfigurations.hpp"
#include <ctime>


class WriteINP
{
    
public:
    
    static void write(const TopologyData & topo, const DCSBasicConfiguration & config, const std::string filename, const std::string elType = "STRI3")
    {
        // include file for abaqus .inp format
        const int nVertices = config.getNumberOfVertices();
        const int nFaces = topo.getNumberOfFaces();
        const auto vertices = config.getVertices();
        const auto face2vertices = topo.getFace2Vertices();
        
        // include time of day
        const auto curtime = std::chrono::system_clock::now();
        const std::time_t curtime_t = std::chrono::system_clock::to_time_t(curtime);
        std::string curtime_str = std::ctime(&curtime_t);
        curtime_str.erase(std::remove(curtime_str.begin(), curtime_str.end(), '\n'), curtime_str.end()); // remove new line
        
        FILE * f = fopen((filename+".inp").c_str(), "w");
        fprintf(f, "**\n");
        fprintf(f, "** part automatically generated on %s\n", curtime_str.c_str());
        fprintf(f, "**\n");
        fprintf(f, "*Node\n");
        for(int i=0;i<nVertices;++i)
        {
            fprintf(f, "%d, %5.5f, %5.5f, %5.5f\n", i+1, vertices(i,0), vertices(i,1), vertices(i,2));
        }
        fprintf(f, "*Element, type=%s\n", elType.c_str());
        for(int i=0;i<nFaces;++i)
        {
            fprintf(f, "%d, %d, %d, %d\n", i+1, face2vertices(i,0)+1, face2vertices(i,1)+1, face2vertices(i,2)+1);
        }
        fclose(f);
    }
    
    static void writeSet(const std::vector<int> elementIndices, const std::string filename, const bool forceStartAtOne=true)
    {
        const int nElements = (int)elementIndices.size();
        FILE * f = fopen((filename+".inp").c_str(), "w");
        int cnt=0; // max 16 per line
        for(const auto & idx : elementIndices)
        {
            fprintf(f, "%d", (forceStartAtOne ? idx+1 : idx));
            cnt+=1;
            if(cnt != nElements)
            {
                fprintf(f, ", ");
                if(cnt%16==0) fprintf(f, "\n");
            }
        }
        fclose(f);
    }
};

#endif /* WriteINP_hpp */
