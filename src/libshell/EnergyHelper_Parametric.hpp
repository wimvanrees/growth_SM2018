//
//  EnergyHelper_Parametric.hpp
//  Elasticity
//
//  Created by Wim van Rees on 04/06/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef EnergyHelper_Parametric_hpp
#define EnergyHelper_Parametric_hpp

#include "common.hpp"

#include "ExtendedTriangleInfo.hpp"
#include "QuadraticForms.hpp"
#include <type_traits> // for std::conditional

struct EnergyHelper_Parametric
{
    const Real material_prefac_1, material_prefac_2;
    
    const Eigen::Matrix2d aform_bar;
    const Eigen::Matrix2d aform_bar_inv;
    const Real area;
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    
    EnergyHelper_Parametric(const Material_Isotropic & matprop, const Eigen::Matrix2d & aform_bar):
    material_prefac_1(matprop.getStVenantFactor1()),
    material_prefac_2(matprop.getStVenantFactor2()),
    aform_bar(aform_bar),
    aform_bar_inv(aform_bar.inverse()),
    area(0.5*std::sqrt(aform_bar.determinant()))
    {
    }
};

struct StrainData
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    
    Eigen::Matrix2d strain;
};

struct StrainData_WithGradientVerts : StrainData
{
    std::array<Eigen::Matrix2d, 3> grad_strain_v0;
    std::array<Eigen::Matrix2d, 3> grad_strain_v1;
    std::array<Eigen::Matrix2d, 3> grad_strain_v2;
};

struct StrainData_WithGradient : StrainData_WithGradientVerts
{
    std::array<Eigen::Matrix2d, 3> grad_strain_v_other_e0;
    std::array<Eigen::Matrix2d, 3> grad_strain_v_other_e1;
    std::array<Eigen::Matrix2d, 3> grad_strain_v_other_e2;
    
    Eigen::Matrix2d grad_strain_e0, grad_strain_e1, grad_strain_e2;
};

template<bool withGradient>
struct StrainData_Stretching : std::conditional<withGradient, StrainData_WithGradientVerts, StrainData>::type
{
    const EnergyHelper_Parametric & helper;
    const FirstFundamentalForm<withGradient> & firstFF;
    
    StrainData_Stretching(const EnergyHelper_Parametric & helper_in, const FirstFundamentalForm<withGradient> & firstFF_in):
    helper(helper_in),
    firstFF(firstFF_in)
    {
    }
    
    void compute()
    {
        this->strain = helper.aform_bar_inv * firstFF.form - Eigen::Matrix2d::Identity();
    }
    
    template<bool U = withGradient> typename std::enable_if<U, void>::type
    compute_grad()
    {
        for(int d=0;d<3;++d)
        {
            this->grad_strain_v0[d] = helper.aform_bar_inv * (Eigen::Matrix2d() << firstFF.gradform.gradv0_11(d), firstFF.gradform.gradv0_12(d), firstFF.gradform.gradv0_12(d), firstFF.gradform.gradv0_22(d)).finished();
            this->grad_strain_v1[d] = helper.aform_bar_inv * (Eigen::Matrix2d() << firstFF.gradform.gradv1_11(d), firstFF.gradform.gradv1_12(d), firstFF.gradform.gradv1_12(d), firstFF.gradform.gradv1_22(d)).finished();
            this->grad_strain_v2[d] = helper.aform_bar_inv * (Eigen::Matrix2d() << firstFF.gradform.gradv2_11(d), firstFF.gradform.gradv2_12(d), firstFF.gradform.gradv2_12(d), firstFF.gradform.gradv2_22(d)).finished();
        }
    }
};

template<bool withGradient>
struct StrainData_Bending : std::conditional<withGradient, StrainData_WithGradient, StrainData>::type
{
    const EnergyHelper_Parametric & helper;
    const Eigen::Matrix2d bform_bar;
    const SecondFundamentalForm<withGradient> & secondFF;
    
    StrainData_Bending(const EnergyHelper_Parametric & helper_in, const SecondFundamentalForm<withGradient> & secondFF_in):
    helper(helper_in),
    bform_bar(Eigen::Matrix2d::Constant(0)),
    secondFF(secondFF_in)
    {
    }
    
    StrainData_Bending(const EnergyHelper_Parametric & helper_in,  const Eigen::Matrix2d & bform_bar_in, const SecondFundamentalForm<withGradient> & secondFF_in):
    helper(helper_in),
    bform_bar(bform_bar_in),
    secondFF(secondFF_in)
    {
    }
    
    void compute()
    {
        this->strain = helper.aform_bar_inv * (secondFF.form - bform_bar);
    }
    
    template<bool U = withGradient> typename std::enable_if<U, void>::type
    compute_grad()
    {
        for(int d=0;d<3;++d)
        {
            this->grad_strain_v0[d] = helper.aform_bar_inv * (Eigen::Matrix2d() << secondFF.gradform.gradv0_11(d), secondFF.gradform.gradv0_12(d), secondFF.gradform.gradv0_12(d), secondFF.gradform.gradv0_22(d)).finished();
            this->grad_strain_v1[d] = helper.aform_bar_inv * (Eigen::Matrix2d() << secondFF.gradform.gradv1_11(d), secondFF.gradform.gradv1_12(d), secondFF.gradform.gradv1_12(d), secondFF.gradform.gradv1_22(d)).finished();
            this->grad_strain_v2[d] = helper.aform_bar_inv * (Eigen::Matrix2d() << secondFF.gradform.gradv2_11(d), secondFF.gradform.gradv2_12(d), secondFF.gradform.gradv2_12(d), secondFF.gradform.gradv2_22(d)).finished();
            
            this->grad_strain_v_other_e0[d] = helper.aform_bar_inv * (Eigen::Matrix2d() << secondFF.gradform.gradv_other_e0_11(d), secondFF.gradform.gradv_other_e0_12(d), secondFF.gradform.gradv_other_e0_12(d), secondFF.gradform.gradv_other_e0_22(d)).finished();
            this->grad_strain_v_other_e1[d] = helper.aform_bar_inv * (Eigen::Matrix2d() << secondFF.gradform.gradv_other_e1_11(d), secondFF.gradform.gradv_other_e1_12(d), secondFF.gradform.gradv_other_e1_12(d), secondFF.gradform.gradv_other_e1_22(d)).finished();
            this->grad_strain_v_other_e2[d] = helper.aform_bar_inv * (Eigen::Matrix2d() << secondFF.gradform.gradv_other_e2_11(d), secondFF.gradform.gradv_other_e2_12(d), secondFF.gradform.gradv_other_e2_12(d), secondFF.gradform.gradv_other_e2_22(d)).finished();
        }
        
        this->grad_strain_e0 = helper.aform_bar_inv * (Eigen::Matrix2d() << secondFF.gradform.gradphi_e0_11, secondFF.gradform.gradphi_e0_12, secondFF.gradform.gradphi_e0_12, secondFF.gradform.gradphi_e0_22).finished();
        this->grad_strain_e1 = helper.aform_bar_inv * (Eigen::Matrix2d() << secondFF.gradform.gradphi_e1_11, secondFF.gradform.gradphi_e1_12, secondFF.gradform.gradphi_e1_12, secondFF.gradform.gradphi_e1_22).finished();
        this->grad_strain_e2 = helper.aform_bar_inv * (Eigen::Matrix2d() << secondFF.gradform.gradphi_e2_11, secondFF.gradform.gradphi_e2_12, secondFF.gradform.gradphi_e2_12, secondFF.gradform.gradphi_e2_22).finished();
    }
};

struct EnergyNormData
{
    const EnergyHelper_Parametric & helper;
    
    Real trace, trace_sq;
    
    EnergyNormData(const EnergyHelper_Parametric & helper_in):
    helper(helper_in)
    {}
    
    inline Real evaluateMaterialNorm()
    {
        return (helper.material_prefac_1*trace*trace + helper.material_prefac_2*trace_sq)*helper.area;
    }
};

struct EnergyNormData_WithGradient_Verts : EnergyNormData
{
    Eigen::Matrix3d grad_trace_v; // col 0 --> v0, col 1 --> v1, col 2 --> v2
    Eigen::Matrix3d grad_tracesq_v; // same as above
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    
    EnergyNormData_WithGradient_Verts(const EnergyHelper_Parametric & helper_in):
    EnergyNormData(helper_in)
    {}
    
    inline Eigen::Vector3d evaluateGradMaterialNorm(const int idx)
    {
        return Eigen::Vector3d((helper.material_prefac_1 * 2.0 * trace * grad_trace_v.col(idx) + helper.material_prefac_2 * grad_tracesq_v.col(idx))*helper.area);
    }
};

struct EnergyNormData_WithGradient : EnergyNormData_WithGradient_Verts
{
    Eigen::Matrix3d grad_trace_v_other_e; // col 0 --> e0, col 1 --> e1, col 2 --> e2
    Eigen::Matrix3d grad_tracesq_v_other_e; // same as above
    
    Eigen::Vector3d grad_trace_e, grad_tracesq_e; // 0 --> e0, 1 --> e1, 2 --> e2
    
    EnergyNormData_WithGradient(const EnergyHelper_Parametric & helper_in):
    EnergyNormData_WithGradient_Verts(helper_in)
    {}
    
    inline Eigen::Vector3d evaluateGradMaterialNormOther(const int idx)
    {
        return Eigen::Vector3d((helper.material_prefac_1 * 2.0 * trace * grad_trace_v_other_e.col(idx) + helper.material_prefac_2 * grad_tracesq_v_other_e.col(idx))*helper.area);
    }
    
    inline Real evaluateGradMaterialNormPhi(const int idx)
    {
        return (helper.material_prefac_1 * 2.0 * trace * grad_trace_e(idx) + helper.material_prefac_2 * grad_tracesq_e(idx))*helper.area;
    }
};

template<bool withGradient>
struct EnergyNormData_Stretching : std::conditional<withGradient, EnergyNormData_WithGradient_Verts, EnergyNormData>::type
{
    typedef typename std::conditional<withGradient, EnergyNormData_WithGradient_Verts, EnergyNormData>::type EnergyNormData_StretchingBase;
    
    StrainData_Stretching<withGradient> & straindata;
    
    EnergyNormData_Stretching(const EnergyHelper_Parametric & helper_in, StrainData_Stretching<withGradient> & straindata_in):
    EnergyNormData_StretchingBase(helper_in),
    straindata(straindata_in)
    {}
    
    void compute()
    {
        straindata.compute();
        this->trace = straindata.strain.trace();
        this->trace_sq = (straindata.strain * straindata.strain).trace();
    }
    
    template<bool U = withGradient> typename std::enable_if<U, void>::type
    compute_grad()
    {
        straindata.compute_grad();
        for(int d=0;d<3;++d)
        {
            this->grad_trace_v(d,0) = (straindata.grad_strain_v0[d]).trace();
            this->grad_trace_v(d,1) = (straindata.grad_strain_v1[d]).trace();
            this->grad_trace_v(d,2) = (straindata.grad_strain_v2[d]).trace();
            
            this->grad_tracesq_v(d,0) = (straindata.strain * straindata.grad_strain_v0[d] + straindata.grad_strain_v0[d] * straindata.strain).trace();
            this->grad_tracesq_v(d,1) = (straindata.strain * straindata.grad_strain_v1[d] + straindata.grad_strain_v1[d] * straindata.strain).trace();
            this->grad_tracesq_v(d,2) = (straindata.strain * straindata.grad_strain_v2[d] + straindata.grad_strain_v2[d] * straindata.strain).trace();
        }
    }
};

template<bool withGradient>
struct EnergyNormData_Bending : std::conditional<withGradient, EnergyNormData_WithGradient, EnergyNormData>::type
{
    typedef typename std::conditional<withGradient, EnergyNormData_WithGradient, EnergyNormData>::type EnergyNormData_BendingBase;
    
    StrainData_Bending<withGradient> & straindata;
    
    EnergyNormData_Bending(const EnergyHelper_Parametric & helper_in, StrainData_Bending<withGradient> & straindata_in):
    EnergyNormData_BendingBase(helper_in),
    straindata(straindata_in)
    {}
    
    void compute()
    {
        straindata.compute();
        this->trace = straindata.strain.trace();
        this->trace_sq = (straindata.strain * straindata.strain).trace();
    }
    
    template<bool U = withGradient> typename std::enable_if<U, void>::type
    compute_grad()
    {
        straindata.compute_grad();
        for(int d=0;d<3;++d)
        {
            this->grad_trace_v(d,0) = (straindata.grad_strain_v0[d]).trace();
            this->grad_trace_v(d,1) = (straindata.grad_strain_v1[d]).trace();
            this->grad_trace_v(d,2) = (straindata.grad_strain_v2[d]).trace();
            
            this->grad_tracesq_v(d,0) = (straindata.strain * straindata.grad_strain_v0[d] + straindata.grad_strain_v0[d] * straindata.strain).trace();
            this->grad_tracesq_v(d,1) = (straindata.strain * straindata.grad_strain_v1[d] + straindata.grad_strain_v1[d] * straindata.strain).trace();
            this->grad_tracesq_v(d,2) = (straindata.strain * straindata.grad_strain_v2[d] + straindata.grad_strain_v2[d] * straindata.strain).trace();
            
            this->grad_trace_v_other_e(d,0) = (straindata.grad_strain_v_other_e0[d]).trace();
            this->grad_trace_v_other_e(d,1) = (straindata.grad_strain_v_other_e1[d]).trace();
            this->grad_trace_v_other_e(d,2) = (straindata.grad_strain_v_other_e2[d]).trace();
            
            this->grad_tracesq_v_other_e(d,0) = (straindata.strain * straindata.grad_strain_v_other_e0[d] + straindata.grad_strain_v_other_e0[d] * straindata.strain).trace();
            this->grad_tracesq_v_other_e(d,1) = (straindata.strain * straindata.grad_strain_v_other_e1[d] + straindata.grad_strain_v_other_e1[d] * straindata.strain).trace();
            this->grad_tracesq_v_other_e(d,2) = (straindata.strain * straindata.grad_strain_v_other_e2[d] + straindata.grad_strain_v_other_e2[d] * straindata.strain).trace();
        }
        
        this->grad_trace_e(0) = (straindata.grad_strain_e0).trace();
        this->grad_trace_e(1) = (straindata.grad_strain_e1).trace();
        this->grad_trace_e(2) = (straindata.grad_strain_e2).trace();
        
        this->grad_tracesq_e(0) = (straindata.strain * straindata.grad_strain_e0 + straindata.grad_strain_e0 * straindata.strain).trace();
        this->grad_tracesq_e(1) = (straindata.strain * straindata.grad_strain_e1 + straindata.grad_strain_e1 * straindata.strain).trace();
        this->grad_tracesq_e(2) = (straindata.strain * straindata.grad_strain_e2 + straindata.grad_strain_e2 * straindata.strain).trace();
    }
};



struct EnergyNormData_BilayerBase
{
    const EnergyHelper_Parametric & helper;
    
    Real trace_12;
    
    EnergyNormData_BilayerBase(const EnergyHelper_Parametric & helper_in):
    helper(helper_in)
    {}
};

struct EnergyNormData_BilayerBase_WithGradient : EnergyNormData_BilayerBase
{
    Eigen::Matrix3d grad_trace12_v;
    Eigen::Matrix3d grad_trace12_v_other_e;
    Eigen::Vector3d grad_trace12_e;
    
    
    EnergyNormData_BilayerBase_WithGradient(const EnergyHelper_Parametric & helper_in):
    EnergyNormData_BilayerBase(helper_in)
    {}
};

template<bool withGradient>
struct EnergyNormData_Bilayer : std::conditional<withGradient, EnergyNormData_BilayerBase_WithGradient, EnergyNormData_BilayerBase>::type
{
    typedef typename std::conditional<withGradient, EnergyNormData_BilayerBase_WithGradient, EnergyNormData_BilayerBase>::type tEnergyNormData_BilayerBase;
    
    const EnergyNormData_Stretching<withGradient> & energydata_stretching;
    const EnergyNormData_Bending<withGradient> & energydata_bending;
    
    EnergyNormData_Bilayer(const EnergyHelper_Parametric & helper_in, const EnergyNormData_Stretching<withGradient> & energydata_stretching_in, const EnergyNormData_Bending<withGradient> & energydata_bending_in):
    tEnergyNormData_BilayerBase(helper_in),
    energydata_stretching(energydata_stretching_in),
    energydata_bending(energydata_bending_in)
    {}
    
    void compute()
    {
        const StrainData_Stretching<withGradient> & straindata_stretching = energydata_stretching.straindata;
        const StrainData_Bending<withGradient> & straindata_bending = energydata_bending.straindata;
        
        // assume energydata_stretching and energydata_bending have already been computed
        this->trace_12 = (straindata_stretching.strain * straindata_bending.strain).trace();
    }
    
    template<bool U = withGradient> typename std::enable_if<U, void>::type
    compute_grad()
    {
        const StrainData_Stretching<withGradient> & straindata_stretching = energydata_stretching.straindata;
        const StrainData_Bending<withGradient> & straindata_bending = energydata_bending.straindata;
        
        // assume energydata_stretching and energydata_bending have already been computed
        
        for(int d=0;d<3;++d)
        {
            this->grad_trace12_v(d,0) = (straindata_stretching.strain * straindata_bending.grad_strain_v0[d] + straindata_stretching.grad_strain_v0[d] * straindata_bending.strain).trace();
            this->grad_trace12_v(d,1) = (straindata_stretching.strain * straindata_bending.grad_strain_v1[d] + straindata_stretching.grad_strain_v1[d] * straindata_bending.strain).trace();
            this->grad_trace12_v(d,2) = (straindata_stretching.strain * straindata_bending.grad_strain_v2[d] + straindata_stretching.grad_strain_v2[d] * straindata_bending.strain).trace();
            
            // no 'other' contribution from straindata_stretching
            this->grad_trace12_v_other_e(d,0) = (straindata_stretching.strain * straindata_bending.grad_strain_v_other_e0[d]).trace();
            this->grad_trace12_v_other_e(d,1) = (straindata_stretching.strain * straindata_bending.grad_strain_v_other_e1[d]).trace();
            this->grad_trace12_v_other_e(d,2) = (straindata_stretching.strain * straindata_bending.grad_strain_v_other_e2[d]).trace();
        }
        
        // no 'phi' contribution from straindata_stretching
        this->grad_trace12_e(0) = (straindata_stretching.strain * straindata_bending.grad_strain_e0).trace();
        this->grad_trace12_e(1) = (straindata_stretching.strain * straindata_bending.grad_strain_e1).trace();
        this->grad_trace12_e(2) = (straindata_stretching.strain * straindata_bending.grad_strain_e2).trace();
    }
    
    inline Real evaluateMaterialNorm()
    {
        return (this->helper.material_prefac_1*energydata_stretching.trace*energydata_bending.trace + this->helper.material_prefac_2*this->trace_12)*this->helper.area;
    }
    
    template<bool U = withGradient> typename std::enable_if<U, Eigen::Vector3d>::type
    evaluateGradMaterialNorm(const int idx)
    {
        return Eigen::Vector3d((this->helper.material_prefac_1 * (energydata_stretching.trace * energydata_bending.grad_trace_v.col(idx) + energydata_bending.trace * energydata_stretching.grad_trace_v.col(idx)) + this->helper.material_prefac_2 * this->grad_trace12_v.col(idx))*this->helper.area);
    }
    
    template<bool U = withGradient> typename std::enable_if<U, Eigen::Vector3d>::type
    evaluateGradMaterialNormOther(const int idx)
    {
        return Eigen::Vector3d((this->helper.material_prefac_1 * energydata_stretching.trace * energydata_bending.grad_trace_v_other_e.col(idx) + this->helper.material_prefac_2 * this->grad_trace12_v_other_e.col(idx))*this->helper.area);
    }
    
    template<bool U = withGradient> typename std::enable_if<U, Real>::type
    evaluateGradMaterialNormPhi(const int idx)
    {
        return (this->helper.material_prefac_1 * energydata_stretching.trace * energydata_bending.grad_trace_e(idx) + this->helper.material_prefac_2 * this->grad_trace12_e(idx))*this->helper.area;
    }
};


struct SaintVenantEnergyData
{
    Real stretching_energy, bending_energy, mixed_energy_ab;
};


struct SaintVenantEnergyData_WithGradient : SaintVenantEnergyData
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    
    Eigen::Vector3d gradv0_aa, gradv1_aa, gradv2_aa;
    Eigen::Vector3d gradv0_bb, gradv1_bb, gradv2_bb;
    Eigen::Vector3d gradv_other_e0_bb, gradv_other_e1_bb, gradv_other_e2_bb;
    
    Real gradphi_e0_bb, gradphi_e1_bb, gradphi_e2_bb;
};

struct SaintVenantEnergyData_WithGradient_Bilayer : SaintVenantEnergyData_WithGradient
{
    Eigen::Vector3d gradv0_ab, gradv1_ab, gradv2_ab;
    Eigen::Vector3d gradv_other_e0_ab, gradv_other_e1_ab, gradv_other_e2_ab;
    
    Real gradphi_e0_ab, gradphi_e1_ab, gradphi_e2_ab;
};

template<bool withGradient, typename tMaterialType, MeshLayer layer>
struct SaintVenantEnergy : std::conditional<withGradient,
typename std::conditional<layer == single, SaintVenantEnergyData_WithGradient, SaintVenantEnergyData_WithGradient_Bilayer>::type, SaintVenantEnergyData>::type
{
};


template<bool withGradient, MeshLayer layer>
struct SaintVenantEnergy<withGradient, Material_Isotropic, layer> : std::conditional<withGradient,
typename std::conditional<layer == single, SaintVenantEnergyData_WithGradient, SaintVenantEnergyData_WithGradient_Bilayer>::type, SaintVenantEnergyData>::type
{
    const Real h_aa, h_bb, h_ab;
    EnergyHelper_Parametric helper;
    
    FirstFundamentalForm<withGradient> firstFF;
    SecondFundamentalForm<withGradient> secondFF;
    
    StrainData_Stretching<withGradient> straindata_stretching;
    StrainData_Bending<withGradient> straindata_bending;

    EnergyNormData_Stretching<withGradient> energydata_stretching;
    EnergyNormData_Bending<withGradient> energydata_bending;
    EnergyNormData_Bilayer<withGradient> energydata_bilayer;
    
    static Real compute_h_aa(const Real thickness) // static because it has to be there before class is constructed
    {
        return (layer == single ? 1.0 : 0.5)*thickness/4.0;
    }
    static Real compute_h_bb(const Real thickness)
    {
        return (layer == single ? 1.0 : 0.5)*std::pow(thickness,3)/12;
    }
    
    static Real compute_h_ab(const Real thickness)
    {
        return (layer == single ? 0.0 : (layer == bottom ? +1.0 : -1.0))*std::pow(thickness,2)/8.0;
    }
    
    SaintVenantEnergy(const Material_Isotropic & mat_prop_in, const Eigen::Matrix2d & aform_bar_in, const ExtendedTriangleInfo & info_in):
    h_aa(compute_h_aa(mat_prop_in.getThickness())),
    h_bb(compute_h_bb(mat_prop_in.getThickness())),
    h_ab(compute_h_ab(mat_prop_in.getThickness())),
    helper(mat_prop_in, aform_bar_in),
    firstFF(info_in),
    secondFF(info_in),
    straindata_stretching(helper, firstFF),
    straindata_bending(helper, secondFF),
    energydata_stretching(helper, straindata_stretching),
    energydata_bending(helper, straindata_bending),
    energydata_bilayer(helper, energydata_stretching, energydata_bending)
    {}
    
    SaintVenantEnergy(const Material_Isotropic & mat_prop_in, const Eigen::Matrix2d & aform_bar_in, const Eigen::Matrix2d & bform_bar_in, const ExtendedTriangleInfo & info_in):
    h_aa(compute_h_aa(mat_prop_in.getThickness())),
    h_bb(compute_h_bb(mat_prop_in.getThickness())),
    h_ab(compute_h_ab(mat_prop_in.getThickness())),
    helper(mat_prop_in, aform_bar_in),
    firstFF(info_in),
    secondFF(info_in),
    straindata_stretching(helper, firstFF),
    straindata_bending(helper, bform_bar_in, secondFF),
    energydata_stretching(helper, straindata_stretching),
    energydata_bending(helper, straindata_bending),
    energydata_bilayer(helper, energydata_stretching, energydata_bending)
    {}
    
    void compute()
    {
        // compute the strains
        energydata_stretching.compute();
        energydata_bending.compute();
        if(layer != single) energydata_bilayer.compute();
        
        // evaluate the energy (area factor inside evaluateNorm)
        this->stretching_energy = h_aa * energydata_stretching.evaluateMaterialNorm();
        this->bending_energy = h_bb * energydata_bending.evaluateMaterialNorm();
        if(layer != single) this->mixed_energy_ab = h_ab * energydata_bilayer.evaluateMaterialNorm();
    }
    
    // conditional gradient function for single and bilayer
    template<bool U = withGradient> typename std::enable_if<U, void>::type
    compute_gradients_single()
    {
        // compute the strain gradients
        energydata_stretching.compute_grad();
        energydata_bending.compute_grad();
        
        // evaluate
        this->gradv0_aa = h_aa * energydata_stretching.evaluateGradMaterialNorm(0);
        this->gradv1_aa = h_aa * energydata_stretching.evaluateGradMaterialNorm(1);
        this->gradv2_aa = h_aa * energydata_stretching.evaluateGradMaterialNorm(2);
        
        this->gradv0_bb = h_bb * energydata_bending.evaluateGradMaterialNorm(0);
        this->gradv1_bb = h_bb * energydata_bending.evaluateGradMaterialNorm(1);
        this->gradv2_bb = h_bb * energydata_bending.evaluateGradMaterialNorm(2);
        
        this->gradv_other_e0_bb = h_bb * energydata_bending.evaluateGradMaterialNormOther(0);
        this->gradv_other_e1_bb = h_bb * energydata_bending.evaluateGradMaterialNormOther(1);
        this->gradv_other_e2_bb = h_bb * energydata_bending.evaluateGradMaterialNormOther(2);
        
        this->gradphi_e0_bb = h_bb * energydata_bending.evaluateGradMaterialNormPhi(0);
        this->gradphi_e1_bb = h_bb * energydata_bending.evaluateGradMaterialNormPhi(1);
        this->gradphi_e2_bb = h_bb * energydata_bending.evaluateGradMaterialNormPhi(2);
    }
    
    // conditional gradient function for single layer only
    template<bool U = withGradient, MeshLayer L = layer> typename std::enable_if<U && L==single, void>::type
    compute_gradients()
    {
        compute_gradients_single();
    }
    
    // conditional gradient function for bilayer only
    template<bool U = withGradient, MeshLayer L = layer> typename std::enable_if<U && L!=single, void>::type
    compute_gradients()
    {
        compute_gradients_single();
        
        energydata_bilayer.compute_grad();
        
        this->gradv0_ab = h_ab * energydata_bilayer.evaluateGradMaterialNorm(0);
        this->gradv1_ab = h_ab * energydata_bilayer.evaluateGradMaterialNorm(1);
        this->gradv2_ab = h_ab * energydata_bilayer.evaluateGradMaterialNorm(2);
        
        this->gradv_other_e0_ab = h_ab * energydata_bilayer.evaluateGradMaterialNormOther(0);
        this->gradv_other_e1_ab = h_ab * energydata_bilayer.evaluateGradMaterialNormOther(1);
        this->gradv_other_e2_ab = h_ab * energydata_bilayer.evaluateGradMaterialNormOther(2);
        
        this->gradphi_e0_ab = h_ab * energydata_bilayer.evaluateGradMaterialNormPhi(0);
        this->gradphi_e1_ab = h_ab * energydata_bilayer.evaluateGradMaterialNormPhi(1);
        this->gradphi_e2_ab = h_ab * energydata_bilayer.evaluateGradMaterialNormPhi(2);
    }
};



#endif /* EnergyHelper_Parametric_hpp */
