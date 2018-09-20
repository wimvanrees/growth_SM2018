//
//  StretchingOperator_Parametric.cpp
//  Elasticity
//
//  Created by Wim van Rees on 1/6/17.
//  Copyright © 2017 Wim van Rees. All rights reserved.
//

#include "StretchingOperator_Parametric.hpp"

#include "EnergyHelper_Parametric.hpp"
#include "EnergyHelper_Parametric_Ortho.hpp"
#include "MergePerFaceQuantities.hpp"
#include "TriangleInfo.hpp"

#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"

template<typename tMesh, typename tMaterialType, MeshLayer layer, bool withGradient>
struct ComputeStretching_Parametric
{
    typedef typename tMesh::tReferenceConfigData tReferenceConfigData;
    typedef typename tMesh::tCurrentConfigData tCurrentConfigData;
    
    const MaterialProperties<tMaterialType>  & material_properties;
    const tCurrentConfigData & currentState;
    const tReferenceConfigData & restState;
    Eigen::Ref<Eigen::MatrixXd> per_face_gradients;
    
    Real energySum;

    void zeroAllSums()
    {
        energySum = 0;
    }
    
    ComputeStretching_Parametric(const MaterialProperties<tMaterialType>  & mat_props_in, const tCurrentConfigData & currentState_in, const tReferenceConfigData & restState_in, Eigen::Ref<Eigen::MatrixXd> per_face_gradients_in):
    material_properties(mat_props_in),
    currentState(currentState_in),
    restState(restState_in),
    per_face_gradients(per_face_gradients_in)
    {
        zeroAllSums();
    }
    
    // split constructor (dont need copy constructor)
    ComputeStretching_Parametric(const ComputeStretching_Parametric & c, tbb::split):
    material_properties(c.material_properties),
    currentState(c.currentState),
    restState(c.restState),
    per_face_gradients(c.per_face_gradients)
    {
        zeroAllSums();
    }
    
    void join(const ComputeStretching_Parametric & j)
    {
        // join the energy
        energySum += j.energySum;
    }
    
    template<bool U = withGradient, MeshLayer L = layer> typename std::enable_if<U && L==single, void>::type
    addGradients(const ExtendedTriangleInfo & info, SaintVenantEnergy<U, tMaterialType, L> & SV_energy)
    {
        SV_energy.compute_gradients();
        
        const int i = info.face_idx;
        
        for(int j=0;j<3;++j)
        {
            per_face_gradients(i, 3*0+j) = SV_energy.gradv0_aa(j);
            per_face_gradients(i, 3*1+j) = SV_energy.gradv1_aa(j);
            per_face_gradients(i, 3*2+j) = SV_energy.gradv2_aa(j);
        }
    }
    
    template<bool U = withGradient, MeshLayer L = layer> typename std::enable_if<U && L!=single, void>::type
    addGradients(const ExtendedTriangleInfo & info, SaintVenantEnergy<U, tMaterialType, L> & SV_energy)
    {
        SV_energy.compute_gradients();
        
        const int i = info.face_idx;
        
        for(int j=0;j<3;++j)
        {
            per_face_gradients(i, 3*0+j) = SV_energy.gradv0_aa(j);
            per_face_gradients(i, 3*1+j) = SV_energy.gradv1_aa(j);
            per_face_gradients(i, 3*2+j) = SV_energy.gradv2_aa(j);
        }
    }
    
    
    template<bool U = withGradient, MeshLayer L = layer> typename std::enable_if<!U, void>::type
    addGradients(const ExtendedTriangleInfo & info, SaintVenantEnergy<U, tMaterialType, L> & SV_energy)
    {
        _unused(info);
        _unused(SV_energy);
        // do nothing
    }
    
    template<MeshLayer L = layer> typename std::enable_if<L==single, Eigen::Matrix2d>::type
    getFirstFundamentalForm(const int i) const
    {
        return this->restState.getFirstFundamentalForm(i);
    }
    
    template<MeshLayer L = layer> typename std::enable_if<L!=single, Eigen::Matrix2d>::type
    getFirstFundamentalForm(const int i) const
    {
        return this->restState.template getFirstFundamentalForm<L>(i);
    }
    
    template<typename tM = tMaterialType> typename std::enable_if<std::is_same<tM, Material_Isotropic>::value, void>::type
    processOneTriangle(const ExtendedTriangleInfo & info, Real & energy)
    {
        const int i = info.face_idx;
        
        const tM matprop = this->material_properties.getFaceMaterial(i);
        const Eigen::Matrix2d aform_bar = getFirstFundamentalForm(i);
        const Eigen::Matrix2d bform_bar = this->restState.getSecondFundamentalForm(i);
        
        SaintVenantEnergy<withGradient, tM, layer> SV_energy(matprop, aform_bar, bform_bar, info);
        
        SV_energy.compute();
        energy += SV_energy.stretching_energy;
        
        addGradients(info, SV_energy);
    }
    
    template<typename tM = tMaterialType> typename std::enable_if<std::is_same<tM, Material_Orthotropic>::value, void>::type
    processOneTriangle(const ExtendedTriangleInfo & info, Real & energy)
    {
        const int i = info.face_idx;
        
        const tM matprop = this->material_properties.getFaceMaterial(i);
        const Eigen::Matrix2d aform_bar = getFirstFundamentalForm(i);
        const Eigen::Matrix2d aform_bar_d1 = matprop.get_aform_bar_d1();
        const Eigen::Matrix2d aform_bar_d2 = aform_bar - aform_bar_d1;
        const Eigen::Matrix2d bform_bar = this->restState.getSecondFundamentalForm(i);
        
        SaintVenantEnergy<withGradient, tM, layer> SV_energy(matprop, aform_bar_d1, aform_bar_d2, bform_bar, info);
        
        SV_energy.compute();
        energy += SV_energy.stretching_energy;
        
        addGradients(info, SV_energy);
    }
    
    
    
    virtual void operator()(const tbb::blocked_range<int> & face_range)
    {
        Real energy = 0.0;
        for (int i=face_range.begin(); i != face_range.end(); ++i)
        {
            const ExtendedTriangleInfo & info = this->currentState.getTriangleInfo(i);
            processOneTriangle(info, energy);
        }
        
        this->energySum += energy;
    }
};


template<typename tMesh, typename tMaterialType, MeshLayer layer>
struct ComputeStretching_Parametric_Subset : ComputeStretching_Parametric<tMesh, tMaterialType, layer, false>
{
    typedef typename tMesh::tReferenceConfigData tReferenceConfigData;
    typedef typename tMesh::tCurrentConfigData tCurrentConfigData;
    
    const std::vector<int> & face_indices;
    
    ComputeStretching_Parametric_Subset(const MaterialProperties<tMaterialType> & mat_props_in, const tCurrentConfigData & currentState_in, const tReferenceConfigData & restState_in, Eigen::Ref<Eigen::MatrixXd> per_face_gradients_in, const std::vector<int> & face_indices_in):
    ComputeStretching_Parametric<tMesh, tMaterialType, layer, false>(mat_props_in, currentState_in, restState_in, per_face_gradients_in),
    face_indices(face_indices_in)
    {}
    
    // split constructor (dont need copy constructor)
    ComputeStretching_Parametric_Subset(const ComputeStretching_Parametric_Subset & c, tbb::split):
    ComputeStretching_Parametric<tMesh, tMaterialType, layer, false>(c),
    face_indices(c.face_indices)
    {
        // for some reason the base class split constructor is not called (I think it just calls the default constructor?)
        // --> need to zero the sums by hand
        
        this->zeroAllSums();
    }
    
    void join(const ComputeStretching_Parametric_Subset & j)
    {
        // join the energy
        this->energySum += j.energySum;
    }
    
    void operator()(const tbb::blocked_range<int> & face_range) override
    {
        Real energy = 0.0;
        
        for (int idx=face_range.begin(); idx != face_range.end(); ++idx)
        {
            const int i = face_indices[idx];
            
            const ExtendedTriangleInfo & info = this->currentState.getTriangleInfo(i);
            this->processOneTriangle(info, energy);
        }
        
        this->energySum += energy;
    }
};

template<typename tMesh, typename tMaterialType, MeshLayer layer>
Real StretchingOperator_Parametric<tMesh, tMaterialType, layer>::computeAll(const tMesh & mesh, Eigen::Ref<Eigen::MatrixXd> gradient_vertices, Eigen::Ref<Eigen::VectorXd> gradient_edges, const bool computeGradient) const
{
    const auto & currentState = mesh.getCurrentConfiguration();
    const auto & restState = mesh.getRestConfiguration();

    const int nFaces = mesh.getNumberOfFaces();
//    const int grainSize = 10;
    
    Eigen::MatrixXd per_face_gradients;
    
    if(not computeGradient)
    {
        this->profiler.push_start("compute energy / gradient");
        
        ComputeStretching_Parametric<tMesh, tMaterialType, layer, false> compute_tbb(material_properties, currentState, restState, per_face_gradients);
        tbb::parallel_reduce(tbb::blocked_range<int>(0,nFaces), compute_tbb, tbb::auto_partitioner());
        
        this->profiler.pop_stop();
        lastEnergy = compute_tbb.energySum;
        
    }
    else
    {
        per_face_gradients.resize(nFaces, 21);
        per_face_gradients.setZero(); // so that entries for bending are not undefined
        
        this->profiler.push_start("compute energy / gradient");
        
        ComputeStretching_Parametric<tMesh, tMaterialType, layer, true> compute_tbb(material_properties, currentState, restState,  per_face_gradients);
        tbb::parallel_reduce(tbb::blocked_range<int>(0,nFaces), compute_tbb, tbb::auto_partitioner());
        
        this->profiler.pop_stop();
        
        lastEnergy = compute_tbb.energySum;
        
        // merge
        const auto & topo = mesh.getTopology();
        const auto & boundaryConditions = mesh.getBoundaryConditions();
        
        const int nVertices = currentState.getNumberOfVertices();
        this->profiler.push_start("merge vertices");
        MergeGradVertices mergevertices(topo, boundaryConditions, per_face_gradients, gradient_vertices);
        tbb::parallel_reduce(tbb::blocked_range<int>(0,nVertices),mergevertices, tbb::auto_partitioner());
        this->profiler.pop_stop();
        
        const int nEdges = topo.getNumberOfEdges();
        this->profiler.push_start("merge edges");
        MergeGradEdges mergeedges(topo, boundaryConditions, per_face_gradients, gradient_edges);
        tbb::parallel_reduce(tbb::blocked_range<int>(0,nEdges),mergeedges, tbb::auto_partitioner());
        this->profiler.pop_stop();
    }
    
    return lastEnergy;
}

template<typename tMesh, typename tMaterialType, MeshLayer layer>
Real StretchingOperator_Parametric<tMesh, tMaterialType, layer>::computeSubsetEnergies(const tMesh & mesh, const std::vector<int> & face_indices) const
{
    const auto & currentState = mesh.getCurrentConfiguration();
    const auto & restState = mesh.getRestConfiguration();
    
    const int nFaces = (int)face_indices.size();
    Eigen::MatrixXd per_face_gradients;
    
    this->profiler.push_start("compute subset energy");
    ComputeStretching_Parametric_Subset<tMesh, tMaterialType, layer> compute_tbb(material_properties, currentState, restState, per_face_gradients, face_indices);
    tbb::parallel_reduce(tbb::blocked_range<int>(0,nFaces), compute_tbb, tbb::auto_partitioner());
    this->profiler.pop_stop();
    
    lastEnergy = compute_tbb.energySum;
    
    return lastEnergy;
}


// explicit instantiations
#include "Mesh.hpp"

// ISOTROPIC MATRIAL
template class StretchingOperator_Parametric<Mesh, Material_Isotropic, single>;
template class ComputeStretching_Parametric<Mesh, Material_Isotropic, single, true>;
template class ComputeStretching_Parametric<Mesh, Material_Isotropic, single, false>;

template class StretchingOperator_Parametric<BilayerMesh, Material_Isotropic, bottom>;
template class ComputeStretching_Parametric<BilayerMesh, Material_Isotropic, bottom, true>;
template class ComputeStretching_Parametric<BilayerMesh, Material_Isotropic, bottom, false>;

template class StretchingOperator_Parametric<BilayerMesh, Material_Isotropic, top>;
template class ComputeStretching_Parametric<BilayerMesh, Material_Isotropic, top, true>;
template class ComputeStretching_Parametric<BilayerMesh, Material_Isotropic, top, false>;

// ORTHOTROPIC MATERIAL
template class StretchingOperator_Parametric<Mesh, Material_Orthotropic, single>;
template class ComputeStretching_Parametric<Mesh, Material_Orthotropic, single, true>;
template class ComputeStretching_Parametric<Mesh, Material_Orthotropic, single, false>;

template class StretchingOperator_Parametric<BilayerMesh, Material_Orthotropic, bottom>;
template class ComputeStretching_Parametric<BilayerMesh, Material_Orthotropic, bottom, true>;
template class ComputeStretching_Parametric<BilayerMesh, Material_Orthotropic, bottom, false>;

template class StretchingOperator_Parametric<BilayerMesh, Material_Orthotropic, top>;
template class ComputeStretching_Parametric<BilayerMesh, Material_Orthotropic, top, true>;
template class ComputeStretching_Parametric<BilayerMesh, Material_Orthotropic, top, false>;

