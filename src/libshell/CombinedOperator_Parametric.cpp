//
//  CombinedOperator_Parametric.cpp
//  Elasticity
//
//  Created by Wim van Rees on 1/6/17.
//  Copyright © 2017 Wim van Rees. All rights reserved.
//

#include "CombinedOperator_Parametric.hpp"

#include "EnergyHelper_Parametric.hpp"
#include "EnergyHelper_Parametric_Ortho.hpp"
#include "MergePerFaceQuantities.hpp"
#include "TriangleInfo.hpp"

#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"

template<typename tMesh, typename tMaterialType, MeshLayer layer, bool withGradient>
struct ComputeCombined_Parametric
{
    typedef typename tMesh::tReferenceConfigData tReferenceConfigData;
    typedef typename tMesh::tCurrentConfigData tCurrentConfigData;
    
    const MaterialProperties<tMaterialType>  & material_properties;
    const tCurrentConfigData & currentState;
    const tReferenceConfigData & restState;
    Eigen::Ref<Eigen::MatrixXd> per_face_gradients;
    
    Real energySum_aa;
    Real energySum_bb;
    Real energySum_ab;

    void zeroAllSums()
    {
        energySum_aa = 0;
        energySum_bb = 0;
        energySum_ab = 0;
    }
    
    ComputeCombined_Parametric(const MaterialProperties<tMaterialType>  & mat_props_in, const tCurrentConfigData & currentState_in, const tReferenceConfigData & restState_in, Eigen::Ref<Eigen::MatrixXd> per_face_gradients_in):
    material_properties(mat_props_in),
    currentState(currentState_in),
    restState(restState_in),
    per_face_gradients(per_face_gradients_in)
    {
        zeroAllSums();
    }
    
    // split constructor (dont need copy constructor)
    ComputeCombined_Parametric(const ComputeCombined_Parametric & c, tbb::split):
    material_properties(c.material_properties),
    currentState(c.currentState),
    restState(c.restState),
    per_face_gradients(c.per_face_gradients)
    {
        zeroAllSums();
    }
    
    void join(const ComputeCombined_Parametric & j)
    {
        // join the energy
        energySum_aa += j.energySum_aa;
        energySum_bb += j.energySum_bb;
        energySum_ab += j.energySum_ab;
    }
    
    template<bool U = withGradient, MeshLayer L = layer> typename std::enable_if<U && L==single, void>::type
    addGradients(const ExtendedTriangleInfo & info, SaintVenantEnergy<U, tMaterialType, L> & SV_energy)
    {
        SV_energy.compute_gradients();
        
        const int i = info.face_idx;
        
        for(int j=0;j<3;++j)
        {
            per_face_gradients(i, 3*0+j) = SV_energy.gradv0_aa(j) + SV_energy.gradv0_bb(j);
            per_face_gradients(i, 3*1+j) = SV_energy.gradv1_aa(j) + SV_energy.gradv1_bb(j);
            per_face_gradients(i, 3*2+j) = SV_energy.gradv2_aa(j) + SV_energy.gradv2_bb(j);
        }
        
        // opposite vertex to each edge (if it exists)
        if(info.other_faces[0] != nullptr)
            for(int j=0;j<3;++j) per_face_gradients(i, 3*3+j) = SV_energy.gradv_other_e0_bb(j);
        
        if(info.other_faces[1] != nullptr)
            for(int j=0;j<3;++j) per_face_gradients(i, 3*4+j) = SV_energy.gradv_other_e1_bb(j);
        
        if(info.other_faces[2] != nullptr)
            for(int j=0;j<3;++j) per_face_gradients(i, 3*5+j) = SV_energy.gradv_other_e2_bb(j);
        
        // theta
        {
            per_face_gradients(i, 3*6+0) = SV_energy.gradphi_e0_bb;
            per_face_gradients(i, 3*6+1) = SV_energy.gradphi_e1_bb;
            per_face_gradients(i, 3*6+2) = SV_energy.gradphi_e2_bb;
        }
    }
    
    template<bool U = withGradient, MeshLayer L = layer> typename std::enable_if<U && L!=single, void>::type
    addGradients(const ExtendedTriangleInfo & info, SaintVenantEnergy<U, tMaterialType, L> & SV_energy)
    {
        SV_energy.compute_gradients();
        
        const int i = info.face_idx;
        
        for(int j=0;j<3;++j)
        {
            per_face_gradients(i, 3*0+j) = SV_energy.gradv0_aa(j) + SV_energy.gradv0_bb(j) + SV_energy.gradv0_ab(j);
            per_face_gradients(i, 3*1+j) = SV_energy.gradv1_aa(j) + SV_energy.gradv1_bb(j) + SV_energy.gradv1_ab(j);
            per_face_gradients(i, 3*2+j) = SV_energy.gradv2_aa(j) + SV_energy.gradv2_bb(j) + SV_energy.gradv2_ab(j);
        }
        
        // opposite vertex to each edge (if it exists)
        if(info.other_faces[0] != nullptr)
            for(int j=0;j<3;++j) per_face_gradients(i, 3*3+j) = SV_energy.gradv_other_e0_bb(j) + SV_energy.gradv_other_e0_ab(j);
        
        if(info.other_faces[1] != nullptr)
            for(int j=0;j<3;++j) per_face_gradients(i, 3*4+j) = SV_energy.gradv_other_e1_bb(j) + SV_energy.gradv_other_e1_ab(j);
        
        if(info.other_faces[2] != nullptr)
            for(int j=0;j<3;++j) per_face_gradients(i, 3*5+j) = SV_energy.gradv_other_e2_bb(j) + SV_energy.gradv_other_e2_ab(j);
        
        // theta
        {
            per_face_gradients(i, 3*6+0) = SV_energy.gradphi_e0_bb + SV_energy.gradphi_e0_ab;
            per_face_gradients(i, 3*6+1) = SV_energy.gradphi_e1_bb + SV_energy.gradphi_e1_ab;
            per_face_gradients(i, 3*6+2) = SV_energy.gradphi_e2_bb + SV_energy.gradphi_e2_ab;
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
    processOneTriangle(const ExtendedTriangleInfo & info, Real & energy_aa, Real & energy_bb, Real & energy_ab)
    {
        const int i = info.face_idx;
        
        const tM matprop = this->material_properties.getFaceMaterial(i);
        const Eigen::Matrix2d aform_bar = getFirstFundamentalForm(i);
        const Eigen::Matrix2d bform_bar = this->restState.getSecondFundamentalForm(i);
        
        SaintVenantEnergy<withGradient, tM, layer> SV_energy(matprop, aform_bar, bform_bar, info);
        
        SV_energy.compute();
        energy_aa += SV_energy.stretching_energy;
        energy_bb += SV_energy.bending_energy;
        if(layer != single) energy_ab += SV_energy.mixed_energy_ab;
        
        addGradients(info, SV_energy);
    }
    
    template<typename tM = tMaterialType> typename std::enable_if<std::is_same<tM, Material_Orthotropic>::value, void>::type
    processOneTriangle(const ExtendedTriangleInfo & info, Real & energy_aa, Real & energy_bb, Real & energy_ab)
    {
        const int i = info.face_idx;
        
        const tM matprop = this->material_properties.getFaceMaterial(i);
        const Eigen::Matrix2d aform_bar = getFirstFundamentalForm(i);
        const Eigen::Matrix2d aform_bar_d1 = matprop.get_aform_bar_d1();
        const Eigen::Matrix2d aform_bar_d2 = aform_bar - aform_bar_d1;
        const Eigen::Matrix2d bform_bar = this->restState.getSecondFundamentalForm(i);
        
        SaintVenantEnergy<withGradient, tM, layer> SV_energy(matprop, aform_bar_d1, aform_bar_d2, bform_bar, info);
        
        SV_energy.compute();
        energy_aa += SV_energy.stretching_energy;
        energy_bb += SV_energy.bending_energy;
        if(layer != single) energy_ab += SV_energy.mixed_energy_ab;
        
        addGradients(info, SV_energy);
    }
    
    
    
    virtual void operator()(const tbb::blocked_range<int> & face_range)
    {
        Real energy_aa = 0.0;
        Real energy_bb = 0.0;
        Real energy_ab = 0.0;
        for (int i=face_range.begin(); i != face_range.end(); ++i)
        {
            const ExtendedTriangleInfo & info = this->currentState.getTriangleInfo(i);
            processOneTriangle(info, energy_aa, energy_bb, energy_ab);
        }
        
        this->energySum_aa += energy_aa;
        this->energySum_bb += energy_bb;
        this->energySum_ab += energy_ab;
    }
};

template<typename tMesh, typename tMaterialType, MeshLayer layer>
struct ComputeCombined_Parametric_Subset : ComputeCombined_Parametric<tMesh, tMaterialType, layer, false>
{
    typedef typename tMesh::tReferenceConfigData tReferenceConfigData;
    typedef typename tMesh::tCurrentConfigData tCurrentConfigData;
    
    const std::vector<int> & face_indices;
    
    ComputeCombined_Parametric_Subset(const MaterialProperties<tMaterialType> & mat_props_in, const tCurrentConfigData & currentState_in, const tReferenceConfigData & restState_in, Eigen::Ref<Eigen::MatrixXd> per_face_gradients_in, const std::vector<int> & face_indices_in):
    ComputeCombined_Parametric<tMesh, tMaterialType, layer, false>(mat_props_in, currentState_in, restState_in, per_face_gradients_in),
    face_indices(face_indices_in)
    {}
    
    // split constructor (dont need copy constructor)
    ComputeCombined_Parametric_Subset(const ComputeCombined_Parametric_Subset & c, tbb::split):
    ComputeCombined_Parametric<tMesh, tMaterialType, layer, false>(c),
    face_indices(c.face_indices)
    {
        // for some reason the base class split constructor is not called (I think it just calls the default constructor?)
        // --> need to zero the sums by hand
        
        this->zeroAllSums();
    }
    
    void join(const ComputeCombined_Parametric_Subset & j)
    {
        // join the energy
        this->energySum_aa += j.energySum_aa;
        this->energySum_bb += j.energySum_bb;
        this->energySum_ab += j.energySum_ab;
    }
    
    void operator()(const tbb::blocked_range<int> & face_range) override
    {
        Real energy_aa = 0.0;
        Real energy_bb = 0.0;
        Real energy_ab = 0.0;
        
        for (int idx=face_range.begin(); idx != face_range.end(); ++idx)
        {
            const int i = face_indices[idx];
            
            const ExtendedTriangleInfo & info = this->currentState.getTriangleInfo(i);
            this->processOneTriangle(info, energy_aa, energy_bb, energy_ab);
        }
        
        this->energySum_aa += energy_aa;
        this->energySum_bb += energy_bb;
        this->energySum_ab += energy_ab;
    }
};

template<typename tMesh, typename tMaterialType, MeshLayer layer>
struct ComputeCombined_Parametric_PerFace : ComputeCombined_Parametric<tMesh, tMaterialType, layer, false>
{
    typedef typename tMesh::tReferenceConfigData tReferenceConfigData;
    typedef typename tMesh::tCurrentConfigData tCurrentConfigData;
    
    Eigen::Ref<Eigen::MatrixXd> per_face_energies;
    
    ComputeCombined_Parametric_PerFace(const MaterialProperties<tMaterialType> & mat_props_in, const tCurrentConfigData & currentState_in, const tReferenceConfigData & restState_in, Eigen::Ref<Eigen::MatrixXd> per_face_gradients_in, Eigen::Ref<Eigen::MatrixXd> per_face_energies_in):
    ComputeCombined_Parametric<tMesh, tMaterialType, layer, false>(mat_props_in, currentState_in, restState_in, per_face_gradients_in),
    per_face_energies(per_face_energies_in)
    {}
    
    // split constructor (dont need copy constructor)
    ComputeCombined_Parametric_PerFace(const ComputeCombined_Parametric_PerFace & c, tbb::split):
    ComputeCombined_Parametric<tMesh, tMaterialType, layer, false>(c),
    per_face_energies(c.per_face_energies)
    {
        // for some reason the base class split constructor is not called (I think it just calls the default constructor?)
        // --> need to zero the sums by hand
        
        this->zeroAllSums();
    }
    
    void join(const ComputeCombined_Parametric_PerFace & j)
    {
        // join the energy
        this->energySum_aa += j.energySum_aa;
        this->energySum_bb += j.energySum_bb;
        this->energySum_ab += j.energySum_ab;
    }
    
    void operator()(const tbb::blocked_range<int> & face_range) override
    {
        Real energy_aa = 0.0;
        Real energy_bb = 0.0;
        Real energy_ab = 0.0;
        
        for (int i=face_range.begin(); i != face_range.end(); ++i)
        {
            Real delta_energy_aa = 0.0;
            Real delta_energy_bb = 0.0;
            Real delta_energy_ab = 0.0;

            const ExtendedTriangleInfo & info = this->currentState.getTriangleInfo(i);
            this->processOneTriangle(info, delta_energy_aa, delta_energy_bb, delta_energy_ab);
            
            // store the energies per face
            per_face_energies(i,0) = delta_energy_aa;
            per_face_energies(i,1) = delta_energy_bb;
            per_face_energies(i,2) = delta_energy_ab;
            
            energy_aa += delta_energy_aa;
            energy_bb += delta_energy_bb;
            energy_ab += delta_energy_ab;
        }
        
        this->energySum_aa += energy_aa;
        this->energySum_bb += energy_bb;
        this->energySum_ab += energy_ab;
    }
};


// only for monolayer isotropic now
template<typename tMesh, typename tMaterialType, MeshLayer layer>
struct ComputeCombined_Parametric_StressPerFace
{
    typedef typename tMesh::tReferenceConfigData tReferenceConfigData;
    typedef typename tMesh::tCurrentConfigData tCurrentConfigData;
 
    const MaterialProperties<tMaterialType>  & material_properties;
    const tCurrentConfigData & currentState;
    const tReferenceConfigData & restState;
    
    Eigen::Ref<Eigen::MatrixXd> per_face_stretch_stress;
    Eigen::Ref<Eigen::MatrixXd> per_face_bend_stress;
    
    ComputeCombined_Parametric_StressPerFace(const MaterialProperties<tMaterialType> & mat_props_in, const tCurrentConfigData & currentState_in, const tReferenceConfigData & restState_in, Eigen::Ref<Eigen::MatrixXd> per_face_stretch_stress_in, Eigen::Ref<Eigen::MatrixXd> per_face_bend_stress_in):
    material_properties(mat_props_in),
    currentState(currentState_in),
    restState(restState_in),
    per_face_stretch_stress(per_face_stretch_stress_in),
    per_face_bend_stress(per_face_bend_stress_in)
    {}
    
    // split constructor (dont need copy constructor)
    ComputeCombined_Parametric_StressPerFace(const ComputeCombined_Parametric_StressPerFace & c, tbb::split):
    material_properties(c.material_properties),
    currentState(c.currentState),
    restState(c.restState),
    per_face_stretch_stress(c.per_face_stretch_stress),
    per_face_bend_stress(c.per_face_bend_stress)
    {
    }
    
    void join(const ComputeCombined_Parametric_StressPerFace & )
    {
    }
    
    template<MeshLayer L = layer, typename tM = tMaterialType> typename std::enable_if<L==single && std::is_same<tM, Material_Isotropic>::value, void>::type
    processOneTriangle(const ExtendedTriangleInfo & info, Eigen::Matrix2d & stress_stretch, Eigen::Matrix2d & stress_bend)
    {
        const int i = info.face_idx;
        
        const tM matprop = this->material_properties.getFaceMaterial(i);
        const Eigen::Matrix2d aform_bar = this->restState.getFirstFundamentalForm(i);
        const Eigen::Matrix2d bform_bar = this->restState.getSecondFundamentalForm(i);
        
        const Eigen::Matrix2d aform_bar_inv = aform_bar.inverse();
        
        const Eigen::Matrix2d aform = info.computeFirstFundamentalForm();
        const Eigen::Matrix2d bform = info.computeSecondFundamentalForm();
        
        const Real alpha = 2*matprop.getStVenantFactor1();
        const Real beta = matprop.getStVenantFactor2();
        
        // here we go
        const Eigen::Matrix2d eps_s = 0.5 * aform_bar_inv * (aform - aform_bar);
        const Eigen::Matrix2d eps_b = - aform_bar_inv * (bform - bform_bar);
        
        stress_stretch = alpha * eps_s.trace() * Eigen::Matrix2d::Identity() + 2 * beta * eps_s;
        stress_bend = alpha * eps_b.trace() * Eigen::Matrix2d::Identity() + 2 * beta * eps_b;
    }
    
    template<MeshLayer L = layer, typename tM = tMaterialType> typename std::enable_if<L!=single || !std::is_same<tM, Material_Isotropic>::value, void>::type
    processOneTriangle(const ExtendedTriangleInfo & , Eigen::Matrix2d & stress_stretch, Eigen::Matrix2d & stress_bend)
    {
        // no implementation : set to zero
        stress_stretch.setZero();
        stress_bend.setZero();
    }
    
    void operator()(const tbb::blocked_range<int> & face_range)
    {
        
        for (int i=face_range.begin(); i != face_range.end(); ++i)
        {
            const ExtendedTriangleInfo & info = this->currentState.getTriangleInfo(i);
            
            Eigen::Matrix2d stress_stretch, stress_bend;
            processOneTriangle(info, stress_stretch, stress_bend);
            
            // store the energies per face
            per_face_stretch_stress(i,0) = stress_stretch(0,0);
            per_face_stretch_stress(i,1) = stress_stretch(0,1); // symmetric
            per_face_stretch_stress(i,2) = stress_stretch(1,1);
            
            per_face_bend_stress(i,0) = stress_bend(0,0);
            per_face_bend_stress(i,1) = stress_bend(0,1); // symmetric
            per_face_bend_stress(i,2) = stress_bend(1,1);
        }
    }
};


// only for monolayer isotropic now
template<typename tMesh, typename tMaterialType, MeshLayer layer>
struct ComputeCombined_Parametric_StrainPerFace
{
    typedef typename tMesh::tReferenceConfigData tReferenceConfigData;
    typedef typename tMesh::tCurrentConfigData tCurrentConfigData;
    
    const MaterialProperties<tMaterialType>  & material_properties;
    const tCurrentConfigData & currentState;
    const tReferenceConfigData & restState;
    
    Eigen::Ref<Eigen::MatrixXd> per_face_stretch_strain;
    Eigen::Ref<Eigen::MatrixXd> per_face_bend_strain;
    
    ComputeCombined_Parametric_StrainPerFace(const MaterialProperties<tMaterialType> & mat_props_in, const tCurrentConfigData & currentState_in, const tReferenceConfigData & restState_in, Eigen::Ref<Eigen::MatrixXd> per_face_stretch_strain_in, Eigen::Ref<Eigen::MatrixXd> per_face_bend_strain_in):
    material_properties(mat_props_in),
    currentState(currentState_in),
    restState(restState_in),
    per_face_stretch_strain(per_face_stretch_strain_in),
    per_face_bend_strain(per_face_bend_strain_in)
    {}
    
    // split constructor (dont need copy constructor)
    ComputeCombined_Parametric_StrainPerFace(const ComputeCombined_Parametric_StrainPerFace & c, tbb::split):
    material_properties(c.material_properties),
    currentState(c.currentState),
    restState(c.restState),
    per_face_stretch_strain(c.per_face_stretch_strain),
    per_face_bend_strain(c.per_face_bend_strain)
    {
    }
    
    void join(const ComputeCombined_Parametric_StrainPerFace & )
    {
    }
    
    template<MeshLayer L = layer, typename tM = tMaterialType> typename std::enable_if<L==single && std::is_same<tM, Material_Isotropic>::value, void>::type
    processOneTriangle(const ExtendedTriangleInfo & info, Eigen::Matrix2d & strain_stretch, Eigen::Matrix2d & strain_bend)
    {
        const int i = info.face_idx;
        
        const Eigen::Matrix2d aform_bar = this->restState.getFirstFundamentalForm(i);
        const Eigen::Matrix2d bform_bar = this->restState.getSecondFundamentalForm(i);
        
        const Eigen::Matrix2d aform_bar_inv = aform_bar.inverse();
        
        const Eigen::Matrix2d aform = info.computeFirstFundamentalForm();
        const Eigen::Matrix2d bform = info.computeSecondFundamentalForm();
        
        // here we go
        strain_stretch = 0.5 * aform_bar_inv * (aform - aform_bar);
        strain_bend = - aform_bar_inv * (bform - bform_bar);
    }
    
    template<MeshLayer L = layer, typename tM = tMaterialType> typename std::enable_if<L!=single || !std::is_same<tM, Material_Isotropic>::value, void>::type
    processOneTriangle(const ExtendedTriangleInfo & , Eigen::Matrix2d & strain_stretch, Eigen::Matrix2d & strain_bend)
    {
        // no implementation : set to zero
        strain_stretch.setZero();
        strain_bend.setZero();
    }
    
    void operator()(const tbb::blocked_range<int> & face_range)
    {
        
        for (int i=face_range.begin(); i != face_range.end(); ++i)
        {
            const ExtendedTriangleInfo & info = this->currentState.getTriangleInfo(i);
            
            Eigen::Matrix2d strain_stretch, strain_bend;
            processOneTriangle(info, strain_stretch, strain_bend);
            
            // store the energies per face
            per_face_stretch_strain(i,0) = strain_stretch(0,0);
            per_face_stretch_strain(i,1) = strain_stretch(0,1);
            per_face_stretch_strain(i,2) = strain_stretch(1,0);
            per_face_stretch_strain(i,3) = strain_stretch(1,1);
            
            per_face_bend_strain(i,0) = strain_bend(0,0);
            per_face_bend_strain(i,1) = strain_bend(0,1);
            per_face_bend_strain(i,2) = strain_bend(1,0);
            per_face_bend_strain(i,3) = strain_bend(1,1);
        }
    }
};


template<typename tMesh, typename tMaterialType, MeshLayer layer>
Real CombinedOperator_Parametric<tMesh, tMaterialType, layer>::computeAll(const tMesh & mesh, Eigen::Ref<Eigen::MatrixXd> gradient_vertices, Eigen::Ref<Eigen::VectorXd> gradient_edges, const bool computeGradient) const
{
    const auto & currentState = mesh.getCurrentConfiguration();
    const auto & restState = mesh.getRestConfiguration();

    const int nFaces = mesh.getNumberOfFaces();
//    const int grainSize = 10;
    
    Eigen::MatrixXd per_face_gradients;
    
    if(not computeGradient)
    {
        this->profiler.push_start("compute energy / gradient");
        
        ComputeCombined_Parametric<tMesh, tMaterialType, layer, false> compute_tbb(material_properties, currentState, restState, per_face_gradients);
        tbb::parallel_reduce(tbb::blocked_range<int>(0,nFaces), compute_tbb, tbb::auto_partitioner());
        
        this->profiler.pop_stop();
        lastEnergy_aa = compute_tbb.energySum_aa;
        lastEnergy_bb = compute_tbb.energySum_bb;
        lastEnergy_ab = compute_tbb.energySum_ab;
        
    }
    else
    {
        per_face_gradients.resize(nFaces, 21);
        
        this->profiler.push_start("compute energy / gradient");
        
        ComputeCombined_Parametric<tMesh, tMaterialType, layer, true> compute_tbb(material_properties, currentState, restState,  per_face_gradients);
        tbb::parallel_reduce(tbb::blocked_range<int>(0,nFaces), compute_tbb, tbb::auto_partitioner());
        
        this->profiler.pop_stop();
        
        lastEnergy_aa = compute_tbb.energySum_aa;
        lastEnergy_bb = compute_tbb.energySum_bb;
        lastEnergy_ab = compute_tbb.energySum_ab;
        
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
    
    return lastEnergy_aa + lastEnergy_bb + lastEnergy_ab;
}

template<typename tMesh, typename tMaterialType, MeshLayer layer>
Real CombinedOperator_Parametric<tMesh, tMaterialType, layer>::computeSubsetEnergies(const tMesh & mesh, const std::vector<int> & face_indices) const
{
    const auto & currentState = mesh.getCurrentConfiguration();
    const auto & restState = mesh.getRestConfiguration();
    
    const int nFaces = (int)face_indices.size();
    Eigen::MatrixXd per_face_gradients;
    
    this->profiler.push_start("compute subset energy");
    ComputeCombined_Parametric_Subset<tMesh, tMaterialType, layer> compute_tbb(material_properties, currentState, restState, per_face_gradients, face_indices);
    tbb::parallel_reduce(tbb::blocked_range<int>(0,nFaces), compute_tbb, tbb::auto_partitioner());
    this->profiler.pop_stop();
    
    lastEnergy_aa = compute_tbb.energySum_aa;
    lastEnergy_bb = compute_tbb.energySum_bb;
    lastEnergy_ab = compute_tbb.energySum_ab;
    
    return lastEnergy_aa + lastEnergy_bb + lastEnergy_ab;
}

template<typename tMesh, typename tMaterialType, MeshLayer layer>
Real CombinedOperator_Parametric<tMesh, tMaterialType, layer>::computePerFaceEnergies(const tMesh & mesh, Eigen::Ref<Eigen::MatrixXd> per_face_energies) const
{
    const auto & currentState = mesh.getCurrentConfiguration();
    const auto & restState = mesh.getRestConfiguration();
    const int nFaces = mesh.getNumberOfFaces();
    
    assert(per_face_energies.rows() == nFaces);
    assert(per_face_energies.cols() == 3);
    
    Eigen::MatrixXd per_face_gradients;
    
    this->profiler.push_start("compute per-face energy");
    ComputeCombined_Parametric_PerFace<tMesh, tMaterialType, layer> compute_tbb(material_properties, currentState, restState, per_face_gradients, per_face_energies);
    tbb::parallel_reduce(tbb::blocked_range<int>(0,nFaces), compute_tbb, tbb::auto_partitioner());
    this->profiler.pop_stop();
    
    lastEnergy_aa = compute_tbb.energySum_aa;
    lastEnergy_bb = compute_tbb.energySum_bb;
    lastEnergy_ab = compute_tbb.energySum_ab;
    
    return lastEnergy_aa + lastEnergy_bb + lastEnergy_ab;
}

template<typename tMesh, typename tMaterialType, MeshLayer layer>
void CombinedOperator_Parametric<tMesh, tMaterialType, layer>::computePerFaceStresses(const tMesh & mesh, Eigen::Ref<Eigen::MatrixXd> per_face_stretch_stress, Eigen::Ref<Eigen::MatrixXd> per_face_bend_stress) const
{
    const auto & currentState = mesh.getCurrentConfiguration();
    const auto & restState = mesh.getRestConfiguration();
    const int nFaces = mesh.getNumberOfFaces();
    
    assert(per_face_stretch_stress.rows() == nFaces);
    assert(per_face_stretch_stress.cols() == 3);
    assert(per_face_bend_stress.rows() == nFaces);
    assert(per_face_bend_stress.cols() == 3);
    
    per_face_stretch_stress.setZero();
    per_face_bend_stress.setZero();
    
    this->profiler.push_start("compute per-face stress");
    ComputeCombined_Parametric_StressPerFace<tMesh, tMaterialType, layer> compute_tbb(material_properties, currentState, restState, per_face_stretch_stress, per_face_bend_stress);
    tbb::parallel_reduce(tbb::blocked_range<int>(0,nFaces), compute_tbb, tbb::auto_partitioner());
    this->profiler.pop_stop();
    
    // done
}

template<typename tMesh, typename tMaterialType, MeshLayer layer>
void CombinedOperator_Parametric<tMesh, tMaterialType, layer>::computePerFaceStrains(const tMesh & mesh, Eigen::Ref<Eigen::MatrixXd> per_face_stretch_strain, Eigen::Ref<Eigen::MatrixXd> per_face_bend_strain) const
{
    const auto & currentState = mesh.getCurrentConfiguration();
    const auto & restState = mesh.getRestConfiguration();
    const int nFaces = mesh.getNumberOfFaces();
    
    assert(per_face_stretch_strain.rows() == nFaces);
    assert(per_face_stretch_strain.cols() == 4);
    assert(per_face_bend_strain.rows() == nFaces);
    assert(per_face_bend_strain.cols() == 4);
    
    per_face_stretch_strain.setZero();
    per_face_bend_strain.setZero();
    
    this->profiler.push_start("compute per-face stress");
    ComputeCombined_Parametric_StrainPerFace<tMesh, tMaterialType, layer> compute_tbb(material_properties, currentState, restState, per_face_stretch_strain, per_face_bend_strain);
    tbb::parallel_reduce(tbb::blocked_range<int>(0,nFaces), compute_tbb, tbb::auto_partitioner());
    this->profiler.pop_stop();
    
    // done
}

// explicit instantiations
#include "Mesh.hpp"

// ISOTROPIC MATRIAL
template class CombinedOperator_Parametric<Mesh, Material_Isotropic, single>;
template class ComputeCombined_Parametric<Mesh, Material_Isotropic, single, true>;
template class ComputeCombined_Parametric<Mesh, Material_Isotropic, single, false>;
template class ComputeCombined_Parametric_PerFace<Mesh, Material_Isotropic, single>;
template class ComputeCombined_Parametric_StressPerFace<Mesh, Material_Isotropic, single>;
template class ComputeCombined_Parametric_StrainPerFace<Mesh, Material_Isotropic, single>;


template class CombinedOperator_Parametric<BilayerMesh, Material_Isotropic, bottom>;
template class ComputeCombined_Parametric<BilayerMesh, Material_Isotropic, bottom, true>;
template class ComputeCombined_Parametric<BilayerMesh, Material_Isotropic, bottom, false>;
template class ComputeCombined_Parametric_PerFace<BilayerMesh, Material_Isotropic, bottom>;
template class ComputeCombined_Parametric_StressPerFace<BilayerMesh, Material_Isotropic, bottom>;
template class ComputeCombined_Parametric_StrainPerFace<BilayerMesh, Material_Isotropic, bottom>;

template class CombinedOperator_Parametric<BilayerMesh, Material_Isotropic, top>;
template class ComputeCombined_Parametric<BilayerMesh, Material_Isotropic, top, true>;
template class ComputeCombined_Parametric<BilayerMesh, Material_Isotropic, top, false>;
template class ComputeCombined_Parametric_PerFace<BilayerMesh, Material_Isotropic, top>;
template class ComputeCombined_Parametric_StressPerFace<BilayerMesh, Material_Isotropic, top>;
template class ComputeCombined_Parametric_StrainPerFace<BilayerMesh, Material_Isotropic, top>;


// ORTHOTROPIC MATERIAL
template class CombinedOperator_Parametric<Mesh, Material_Orthotropic, single>;
template class ComputeCombined_Parametric<Mesh, Material_Orthotropic, single, true>;
template class ComputeCombined_Parametric<Mesh, Material_Orthotropic, single, false>;

template class CombinedOperator_Parametric<BilayerMesh, Material_Orthotropic, bottom>;
template class ComputeCombined_Parametric<BilayerMesh, Material_Orthotropic, bottom, true>;
template class ComputeCombined_Parametric<BilayerMesh, Material_Orthotropic, bottom, false>;

template class CombinedOperator_Parametric<BilayerMesh, Material_Orthotropic, top>;
template class ComputeCombined_Parametric<BilayerMesh, Material_Orthotropic, top, true>;
template class ComputeCombined_Parametric<BilayerMesh, Material_Orthotropic, top, false>;

