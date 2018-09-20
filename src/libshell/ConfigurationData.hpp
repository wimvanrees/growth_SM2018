//
//  ConfigurationData.hpp
//  Elasticity
//
//  Created by Wim van Rees on 23/04/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef ConfigurationData_hpp
#define ConfigurationData_hpp

#ifdef USELIBLBFGS
#include "lbfgs.h"
#endif

#include "common.hpp"
#include "TopologyData.hpp"
#include "TriangleInfo.hpp"
#include "ExtendedTriangleInfo.hpp"
#include "BoundaryConditionsData.hpp"

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"

// data management class : did I basically create a std::vector ?
struct ConfigurationDataBase
{
    int nVariables;
    Real * dataptr;
    
    // default constructor
    ConfigurationDataBase():
    dataptr(nullptr)
    {}
    
    void init(const int nVariables_in, const Real initVal = 0)
    {
        nVariables = nVariables_in;
        allocate_data();
        for(int i=0;i<nVariables;++i) dataptr[i] = initVal;
    }
    
    void init(const Eigen::Ref<const Eigen::VectorXd> dataptr_in)
    {
        nVariables = dataptr_in.rows();
        allocate_data();
        for(int i=0;i<nVariables;++i) dataptr[i] = dataptr_in(i);
    }
    
    virtual void clear()
    {
        deallocate_data();
    }
    
    // assignment operator
    ConfigurationDataBase & operator= (const ConfigurationDataBase & other)
    {
        if (this != &other) // protect against invalid self-assignment
        {
            nVariables = other.nVariables;
            
            // allocate
            if(dataptr == nullptr || nVariables != other.nVariables)
                allocate_data();
            
            // copy
            std::copy(other.dataptr, other.dataptr + nVariables, dataptr);
        }
        // by convention, always return *this
        return *this;
    }
    
    // copy constructor
    ConfigurationDataBase(const ConfigurationDataBase &other):
    dataptr(nullptr)
    {
        nVariables = other.nVariables;
        allocate_data();
        
        // copy
        std::copy(other.dataptr, other.dataptr + nVariables, dataptr);
    }
    
    friend void swap(ConfigurationDataBase& first, ConfigurationDataBase& second)
    {
        // http://stackoverflow.com/questions/3279543/what-is-the-copy-and-swap-idiom
        using std::swap;
        
        // store the number of variables for each of the two
        const int nVariables_1 = first.nVariables;
        const int nVariables_2 = second.nVariables;
        
        // now we can swap the data
        swap(first.dataptr, second.dataptr);

        // swap the variables
        first.nVariables = nVariables_2;
        second.nVariables = nVariables_1;
    }

    void deallocate_data()
    {
        if(dataptr==nullptr) return;
        
#ifdef USELBFGS
        lbfgs_free(dataptr);
#else
        delete [] dataptr;
#endif
        dataptr = nullptr;
    }
    
    void allocate_data()
    {
        deallocate_data();
#ifdef USELBFGS
        dataptr  = lbfgs_malloc(nVariables);
#else
        dataptr = new Real[nVariables];
#endif
    }
    
    virtual ~ConfigurationDataBase()
    {
        deallocate_data();
    }
    
    Eigen::Ref<Eigen::VectorXd> getData()
    {
        return Eigen::Map<Eigen::VectorXd>(dataptr, nVariables);
    }
    
    Real * getDataPointer() const
    {
        return dataptr;
    }
    
    int getNumberOfDataVariables() const
    {
        return nVariables;
    }
};

// a structure for DCS data (static or dynamic or anything else)
// it basically adds the wrappers and a couple of accessors
struct DCSConfigurationData : ConfigurationDataBase
{
    // the eigen wrappers
    Eigen::Map<Eigen::MatrixXd> vertexdata;
    Eigen::Map<Eigen::VectorXd> edgedata;
    
    // default constructor
    DCSConfigurationData():
    ConfigurationDataBase(),
    vertexdata(nullptr,0,0), edgedata(nullptr,0)
    {}
    
    // assignment operator
    DCSConfigurationData & operator= (const DCSConfigurationData & other)
    {
        if (this != &other) // protect against invalid self-assignment
        {
            const int nVertices = other.vertexdata.rows();
            const int nEdges = other.edgedata.rows();
            
            ConfigurationDataBase::operator=(other);
            mapArrays(nVertices, nEdges);
        }
        // by convention, always return *this
        return *this;
    }
    
    // copy constructor
    DCSConfigurationData(const DCSConfigurationData &other):
    ConfigurationDataBase(other),
    vertexdata(nullptr,0,0), edgedata(nullptr,0)
    {
        const int nVertices = other.vertexdata.rows();;
        const int nEdges = other.edgedata.rows();
        mapArrays(nVertices, nEdges);
    }
    
    friend void swap(DCSConfigurationData& first, DCSConfigurationData& second)
    {
        // http://stackoverflow.com/questions/3279543/what-is-the-copy-and-swap-idiom
        using std::swap;
        
        // store the number of vertexdata/edgedirectors for each of the two
        const int nV_1 = first.vertexdata.rows();
        const int nE_1 = first.edgedata.rows();
        
        const int nV_2 = second.vertexdata.rows();
        const int nE_2 = second.edgedata.rows();
        
        swap(static_cast<ConfigurationDataBase&>(first), static_cast<ConfigurationDataBase&>(second));
        
        first.mapArrays(nV_2, nE_2);
        second.mapArrays(nV_1, nE_1);
    }
    
    // initialize the eigen mappers
    void mapArrays(const int nVertices, const int nEdges)
    {
        assert(nVariables == (3*nVertices + nEdges));
        new (&vertexdata) Eigen::Map<Eigen::MatrixXd>(dataptr, nVertices, 3);
        new (&edgedata) Eigen::Map<Eigen::VectorXd>(dataptr + 3*nVertices, nEdges);
    }
    
    virtual void init(const int nVertices, const int nEdges)
    {
        Eigen::MatrixXd vertexdata(nVertices,3);
        vertexdata.setZero();
        init(vertexdata, nEdges);
    }
    
    virtual void init(const Eigen::MatrixXd & vertexdata_in, const int nEdges)
    {
        Eigen::VectorXd edgedata_in(nEdges);
        edgedata_in.setZero();
        init(vertexdata_in, edgedata_in);
    }
    
    virtual void init(const Eigen::MatrixXd & vertexdata_in, const Eigen::VectorXd & edgedata_in)
    {
        const int nVertices = vertexdata_in.rows();
        const int nEdges = edgedata_in.rows();
        
        // initialize the data array
        ConfigurationDataBase::init(3*nVertices + nEdges);
        
        // map our guys
        mapArrays(nVertices, nEdges);
        
        // assign the variables
        vertexdata = vertexdata_in;
        edgedata = edgedata_in;
    }
    
    void clear() override
    {
        new (&vertexdata) Eigen::Map<Eigen::MatrixXd>(nullptr, 0, 0);
        new (&edgedata) Eigen::Map<Eigen::VectorXd>(nullptr, 0);
        
        ConfigurationDataBase::clear();
    }
    
    TriangleInfo getTriangleInfo(const TopologyData & topo, const int face_idx) const
    {
        const int idx_e0 = topo.face2edges(face_idx,0);
        const int idx_e1 = topo.face2edges(face_idx,1);
        const int idx_e2 = topo.face2edges(face_idx,2);
        
        // vertex indices
        const int idx_v0 = topo.face2vertices(face_idx,0);
        const int idx_v1 = topo.face2vertices(face_idx,1);
        const int idx_v2 = topo.face2vertices(face_idx,2);
        
        // vertex vectors
        const Eigen::Vector3d v0 = vertexdata.row(idx_v0);
        const Eigen::Vector3d v1 = vertexdata.row(idx_v1);
        const Eigen::Vector3d v2 = vertexdata.row(idx_v2);
        
        TriangleInfo info(v0,v1,v2);
        info.face_idx = face_idx;
        info.idx_v0 = idx_v0;
        info.idx_v1 = idx_v1;
        info.idx_v2 = idx_v2;
        info.idx_e0 = idx_e0;
        info.idx_e1 = idx_e1;
        info.idx_e2 = idx_e2;
        
        return info;
    }
    
    int getNumberOfVertices() const
    {
        return vertexdata.rows();
    }
    
    int getNumberOfEdges() const
    {
        return edgedata.rows();
    }
    
    Eigen::Ref<Eigen::MatrixXd> getVertexData()
    {
        return vertexdata;
    }
    
    Eigen::Ref<Eigen::VectorXd> getEdgeData()
    {
        return edgedata;
    }
    
    const Eigen::Ref<const Eigen::MatrixXd> getVertexData() const
    {
        return vertexdata;
    }
    
    const Eigen::Ref<const Eigen::VectorXd> getEdgeData() const
    {
        return edgedata;
    }
};

#endif /* ConfigurationData_hpp */
