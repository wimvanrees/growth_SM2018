//
//  DCSCurrentConfigurationData.hpp
//  Elasticity
//
//  Created by Wim van Rees on 1/18/17.
//  Copyright © 2017 Wim van Rees. All rights reserved.
//

#ifndef DCSCurrentConfigurationData_hpp
#define DCSCurrentConfigurationData_hpp

#include "common.hpp"
#include "ConfigurationData.hpp"
#include "ExtendedTriangleInfo.hpp"
#include "CreateExtendedTriangleInfos.hpp"
#include "TopologyData.hpp"
#include "BoundaryConditionsData.hpp"

struct DCSCurrentConfigurationData
{
    std::vector<ExtendedTriangleInfo> ext_tinfo;
    
    // default constructor
    DCSCurrentConfigurationData()
    {}
    
    void update(const DCSConfigurationData & meshdata, const TopologyData & topo, const BoundaryConditionsData & boundaryConditions)
    {
        const int nFaces = topo.getNumberOfFaces();
        if((int)ext_tinfo.size() != nFaces)
            ext_tinfo.resize(nFaces);
        
        //        CreateExtendedTriangleInfos::create(topo, meshdata, boundaryConditions, ext_tinfo);
        {
            CreateExtendedTriangleInfos_tbb<0> create_infos_0(topo, meshdata, boundaryConditions, ext_tinfo);
            tbb::parallel_for(tbb::blocked_range<int>(0,nFaces),create_infos_0, tbb::auto_partitioner());
        }
        {
            const int nEdges = topo.getNumberOfEdges();
            CreateExtendedTriangleInfos_tbb<1> create_infos_1(topo, meshdata, boundaryConditions, ext_tinfo);
            tbb::parallel_for(tbb::blocked_range<int>(0,nEdges),create_infos_1, tbb::auto_partitioner());
        }
    }
    
    virtual void clear()
    {
        ext_tinfo.clear();
    }
        
    std::vector<int> getEdgeIndices(const TopologyData & topo, const std::vector<int> & face_indices) const
    {
        std::set<int> edge_indices_set; // use set to get the unique list
        const size_t nFaces_subset = face_indices.size();
        for(size_t i=0;i<nFaces_subset;++i)
        {
            const int f_idx = face_indices[i];
            edge_indices_set.insert(topo.face2edges(f_idx,0));
            edge_indices_set.insert(topo.face2edges(f_idx,1));
            edge_indices_set.insert(topo.face2edges(f_idx,2));
        }
        // convert set to vector for easy iteration
        const std::vector<int> edge_indices( edge_indices_set.begin(), edge_indices_set.end() );
        return edge_indices;
    }
    
    void updateSubSet(const DCSConfigurationData & meshdata, const TopologyData & topo, const BoundaryConditionsData & boundaryConditions, const std::vector<int> & face_indices)
    {
        assert((int)ext_tinfo.size() == topo.getNumberOfFaces());
        
        // first get the corresponding edge indices
        const std::vector<int> edge_indices = getEdgeIndices(topo, face_indices);
        
        // then we do our thing
        {
            const int nFaces_subset = (int)face_indices.size();
            CreateExtendedTriangleInfos_tbb_subset<0> create_infos_0(topo, meshdata, boundaryConditions, ext_tinfo, face_indices, edge_indices);
            tbb::parallel_for(tbb::blocked_range<int>(0,nFaces_subset),create_infos_0, tbb::auto_partitioner());
        }
        {
            const int nEdges_subset = (int)edge_indices.size();
            CreateExtendedTriangleInfos_tbb_subset<1> create_infos_1(topo, meshdata, boundaryConditions, ext_tinfo, face_indices, edge_indices);
            tbb::parallel_for(tbb::blocked_range<int>(0,nEdges_subset),create_infos_1, tbb::auto_partitioner());
        }
        
    }
};

#endif /* DCSCurrentConfigurationData_hpp */
