//
//  EnergyHelper_Parametric_Ortho.hpp
//  Elasticity
//
//  Created by Wim van Rees on 1/19/17.
//  Copyright © 2017 Wim van Rees. All rights reserved.
//

#ifndef EnergyHelper_Parametric_Ortho_hpp
#define EnergyHelper_Parametric_Ortho_hpp

#include "common.hpp"
#include "ExtendedTriangleInfo.hpp"
#include "QuadraticForms.hpp"
#include <type_traits> // for std::conditional

struct EnergyHelper_Parametric_Ortho
{
    const Real material_prefac_1_d1;
    const Real material_prefac_1_d2;
    const Real material_prefac_1_d12;
    const Real material_prefac_2_d1;
    const Real material_prefac_2_d2;
    const Real material_prefac_2_d12;
    
    const Eigen::Matrix2d aform_bar;
    const Eigen::Matrix2d aform_bar_inv;
    const std::array<Eigen::Matrix2d, 2> decomp;
    const Real area;
    
    EnergyHelper_Parametric_Ortho(const Material_Orthotropic & matprop, const Eigen::Matrix2d & aform_bar_d1, const Eigen::Matrix2d & aform_bar_d2):
    material_prefac_1_d1(matprop.getStVenantFactor_1_d1()),
    material_prefac_1_d2(matprop.getStVenantFactor_1_d2()),
    material_prefac_1_d12(matprop.getStVenantFactor_1_d12()),
    material_prefac_2_d1(matprop.getStVenantFactor_2_d1()),
    material_prefac_2_d2(matprop.getStVenantFactor_2_d2()),
    material_prefac_2_d12(matprop.getStVenantFactor_2_d12()),
    aform_bar(aform_bar_d1 + aform_bar_d2),
    aform_bar_inv(aform_bar.inverse()),
    decomp{(aform_bar_inv * aform_bar_d1),(aform_bar_inv * aform_bar_d2)},
    area(0.5*std::sqrt(aform_bar.determinant()))
    {
    }
};

struct StrainData_Ortho
{
    std::array<Eigen::Matrix2d, 2> strain;
};

struct StrainData_Ortho_WithGradientVerts : StrainData_Ortho
{
    std::array<std::array<Eigen::Matrix2d, 3>, 2> grad_strain_v0;
    std::array<std::array<Eigen::Matrix2d, 3>, 2> grad_strain_v1;
    std::array<std::array<Eigen::Matrix2d, 3>, 2> grad_strain_v2;
};

struct StrainData_Ortho_WithGradient : StrainData_Ortho_WithGradientVerts
{
    std::array<std::array<Eigen::Matrix2d, 3>, 2> grad_strain_v_other_e0;
    std::array<std::array<Eigen::Matrix2d, 3>, 2> grad_strain_v_other_e1;
    std::array<std::array<Eigen::Matrix2d, 3>, 2> grad_strain_v_other_e2;
    
    std::array<Eigen::Matrix2d, 2> grad_strain_e0, grad_strain_e1, grad_strain_e2;
};

template<bool withGradient>
struct StrainData_Ortho_Stretching : std::conditional<withGradient, StrainData_Ortho_WithGradientVerts, StrainData_Ortho>::type
{
    const EnergyHelper_Parametric_Ortho & helper;
    const FirstFundamentalForm<withGradient> & firstFF;
    
    StrainData_Ortho_Stretching(const EnergyHelper_Parametric_Ortho & helper_in, const FirstFundamentalForm<withGradient> & firstFF_in):
    helper(helper_in),
    firstFF(firstFF_in)
    {
    }
    
    void compute()
    {
        for(int c=0;c<2;++c)
            this->strain[c] = helper.aform_bar_inv * (firstFF.form - helper.aform_bar) * helper.decomp[c];
    }
    
    template<bool U = withGradient> typename std::enable_if<U, void>::type
    compute_grad()
    {
        for(int c=0;c<2;++c)
            for(int d=0;d<3;++d)
            {
                this->grad_strain_v0[c][d] = helper.aform_bar_inv * (Eigen::Matrix2d() << firstFF.gradform.gradv0_11(d), firstFF.gradform.gradv0_12(d), firstFF.gradform.gradv0_12(d), firstFF.gradform.gradv0_22(d)).finished() * helper.decomp[c];
                this->grad_strain_v1[c][d] = helper.aform_bar_inv * (Eigen::Matrix2d() << firstFF.gradform.gradv1_11(d), firstFF.gradform.gradv1_12(d), firstFF.gradform.gradv1_12(d), firstFF.gradform.gradv1_22(d)).finished() * helper.decomp[c];
                this->grad_strain_v2[c][d] = helper.aform_bar_inv * (Eigen::Matrix2d() << firstFF.gradform.gradv2_11(d), firstFF.gradform.gradv2_12(d), firstFF.gradform.gradv2_12(d), firstFF.gradform.gradv2_22(d)).finished() * helper.decomp[c];
            }
    }
};


template<bool withGradient>
struct StrainData_Ortho_Bending : std::conditional<withGradient, StrainData_Ortho_WithGradient, StrainData_Ortho>::type
{
    const EnergyHelper_Parametric_Ortho & helper;
        const Eigen::Matrix2d bform_bar;
    const SecondFundamentalForm<withGradient> & secondFF;
    
    StrainData_Ortho_Bending(const EnergyHelper_Parametric_Ortho & helper_in, const SecondFundamentalForm<withGradient> & secondFF_in):
    helper(helper_in),
    bform_bar(Eigen::Matrix2d::Constant(0)),
    secondFF(secondFF_in)
    {
    }
    
    StrainData_Ortho_Bending(const EnergyHelper_Parametric_Ortho & helper_in, const Eigen::Matrix2d & bform_bar_in, const SecondFundamentalForm<withGradient> & secondFF_in):
    helper(helper_in),
    bform_bar(bform_bar_in),
    secondFF(secondFF_in)
    {
    }
    
    void compute()
    {
        for(int c=0;c<2;++c)
            this->strain[c] = helper.aform_bar_inv * (secondFF.form - bform_bar) * helper.decomp[c];

    }
    
    template<bool U = withGradient> typename std::enable_if<U, void>::type
    compute_grad()
    {
        for(int c=0;c<2;++c)
        {
            for(int d=0;d<3;++d)
            {
                this->grad_strain_v0[c][d] = helper.aform_bar_inv * (Eigen::Matrix2d() << secondFF.gradform.gradv0_11(d), secondFF.gradform.gradv0_12(d), secondFF.gradform.gradv0_12(d), secondFF.gradform.gradv0_22(d)).finished() * helper.decomp[c];
                this->grad_strain_v1[c][d] = helper.aform_bar_inv * (Eigen::Matrix2d() << secondFF.gradform.gradv1_11(d), secondFF.gradform.gradv1_12(d), secondFF.gradform.gradv1_12(d), secondFF.gradform.gradv1_22(d)).finished() * helper.decomp[c];
                this->grad_strain_v2[c][d] = helper.aform_bar_inv * (Eigen::Matrix2d() << secondFF.gradform.gradv2_11(d), secondFF.gradform.gradv2_12(d), secondFF.gradform.gradv2_12(d), secondFF.gradform.gradv2_22(d)).finished() * helper.decomp[c];
                
                this->grad_strain_v_other_e0[c][d] = helper.aform_bar_inv * (Eigen::Matrix2d() << secondFF.gradform.gradv_other_e0_11(d), secondFF.gradform.gradv_other_e0_12(d), secondFF.gradform.gradv_other_e0_12(d), secondFF.gradform.gradv_other_e0_22(d)).finished() * helper.decomp[c];
                this->grad_strain_v_other_e1[c][d] = helper.aform_bar_inv * (Eigen::Matrix2d() << secondFF.gradform.gradv_other_e1_11(d), secondFF.gradform.gradv_other_e1_12(d), secondFF.gradform.gradv_other_e1_12(d), secondFF.gradform.gradv_other_e1_22(d)).finished() * helper.decomp[c];
                this->grad_strain_v_other_e2[c][d] = helper.aform_bar_inv * (Eigen::Matrix2d() << secondFF.gradform.gradv_other_e2_11(d), secondFF.gradform.gradv_other_e2_12(d), secondFF.gradform.gradv_other_e2_12(d), secondFF.gradform.gradv_other_e2_22(d)).finished() * helper.decomp[c];
            }
            
            this->grad_strain_e0[c] = helper.aform_bar_inv * (Eigen::Matrix2d() << secondFF.gradform.gradphi_e0_11, secondFF.gradform.gradphi_e0_12, secondFF.gradform.gradphi_e0_12, secondFF.gradform.gradphi_e0_22).finished() * helper.decomp[c];
            this->grad_strain_e1[c] = helper.aform_bar_inv * (Eigen::Matrix2d() << secondFF.gradform.gradphi_e1_11, secondFF.gradform.gradphi_e1_12, secondFF.gradform.gradphi_e1_12, secondFF.gradform.gradphi_e1_22).finished() * helper.decomp[c];
            this->grad_strain_e2[c] = helper.aform_bar_inv * (Eigen::Matrix2d() << secondFF.gradform.gradphi_e2_11, secondFF.gradform.gradphi_e2_12, secondFF.gradform.gradphi_e2_12, secondFF.gradform.gradphi_e2_22).finished() * helper.decomp[c];
        }
    }
};

struct EnergyNormData_Ortho
{
    const EnergyHelper_Parametric_Ortho & helper;
    
    std::array<Real, 2> trace;
    std::array<Real, 2> trace_sq;
    Real trace_12;
    
    EnergyNormData_Ortho(const EnergyHelper_Parametric_Ortho & helper_in):
    helper(helper_in)
    {}
    
    inline Real evaluateMaterialNorm()
    {
        const Real fac1_d1  = helper.material_prefac_1_d1 * trace[0] * trace[0];
        const Real fac1_d2  = helper.material_prefac_1_d2 * trace[1] * trace[1];
        const Real fac1_d12 = helper.material_prefac_1_d12 * 2.0 * trace[0] * trace[1]; // (trace0 * trace1 + trace1 * trace0)

        const Real fac2_d1 = helper.material_prefac_2_d1 * trace_sq[0];
        const Real fac2_d2 = helper.material_prefac_2_d2 * trace_sq[1];
        const Real fac2_d12 = helper.material_prefac_2_d12 * 2.0 * trace_12; // (trace (0 * 1) + trace (1 * 0))
        
        return (fac1_d1 + fac1_d2 + fac1_d12 + fac2_d1 + fac2_d2 + fac2_d12)*helper.area;
    }
};

struct EnergyNormData_Ortho_WithGradient_Verts : EnergyNormData_Ortho
{
    std::array<Eigen::Matrix3d, 2> grad_trace_v; // col 0 --> v0, col 1 --> v1, col 2 --> v2
    std::array<Eigen::Matrix3d, 2> grad_tracesq_v; // col 0 --> v0, col 1 --> v1, col 2 --> v2
    Eigen::Matrix3d grad_trace12_v; // same as above
    
    EnergyNormData_Ortho_WithGradient_Verts(const EnergyHelper_Parametric_Ortho & helper_in):
    EnergyNormData_Ortho(helper_in)
    {}
    
    inline Eigen::Vector3d evaluateGradMaterialNorm(const int idx)
    {
        const Eigen::Vector3d grad_fac1_d1 = helper.material_prefac_1_d1 * 2.0 * trace[0] * grad_trace_v[0].col(idx);
        const Eigen::Vector3d grad_fac1_d2 = helper.material_prefac_1_d2 * 2.0 * trace[1] * grad_trace_v[1].col(idx);
        const Eigen::Vector3d grad_fac1_d12 = helper.material_prefac_1_d12 * 2.0 * (trace[0] * grad_trace_v[1].col(idx) + trace[1] * grad_trace_v[0].col(idx));
        
        const Eigen::Vector3d grad_fac2_d1 = helper.material_prefac_2_d1 * grad_tracesq_v[0].col(idx);
        const Eigen::Vector3d grad_fac2_d2 = helper.material_prefac_2_d2 * grad_tracesq_v[1].col(idx);
        const Eigen::Vector3d grad_fac2_d12 = helper.material_prefac_2_d12 * 2.0 * grad_trace12_v.col(idx);
        
        return Eigen::Vector3d( (grad_fac1_d1 + grad_fac1_d2 + grad_fac1_d12 + grad_fac2_d1 + grad_fac2_d2 + grad_fac2_d12)*helper.area );
    }
};

struct EnergyNormData_Ortho_WithGradient : EnergyNormData_Ortho_WithGradient_Verts
{
    std::array<Eigen::Matrix3d, 2> grad_trace_v_other_e; // col 0 --> e0, col 1 --> e1, col 2 --> e2
    std::array<Eigen::Matrix3d, 2> grad_tracesq_v_other_e; // col 0 --> e0, col 1 --> e1, col 2 --> e2
    Eigen::Matrix3d grad_trace12_v_other_e; // same as above
    
    std::array<Eigen::Vector3d, 2> grad_trace_e;
    std::array<Eigen::Vector3d, 2> grad_tracesq_e;
    Eigen::Vector3d grad_trace12_e; // 0 --> e0, 1 --> e1, 2 --> e2
    
    EnergyNormData_Ortho_WithGradient(const EnergyHelper_Parametric_Ortho & helper_in):
    EnergyNormData_Ortho_WithGradient_Verts(helper_in)
    {}
    
    inline Eigen::Vector3d evaluateGradMaterialNormOther(const int idx)
    {
        const Eigen::Vector3d grad_fac1_d1 = helper.material_prefac_1_d1 * 2.0 * trace[0] * grad_trace_v_other_e[0].col(idx);
        const Eigen::Vector3d grad_fac1_d2 = helper.material_prefac_1_d2 * 2.0 * trace[1] * grad_trace_v_other_e[1].col(idx);
        const Eigen::Vector3d grad_fac1_d12 = helper.material_prefac_1_d12 * 2.0 * (trace[0] * grad_trace_v_other_e[1].col(idx) + trace[1] * grad_trace_v_other_e[0].col(idx));
        
        const Eigen::Vector3d grad_fac2_d1 = helper.material_prefac_2_d1 * grad_tracesq_v_other_e[0].col(idx);
        const Eigen::Vector3d grad_fac2_d2 = helper.material_prefac_2_d2 * grad_tracesq_v_other_e[1].col(idx);
        const Eigen::Vector3d grad_fac2_d12 = helper.material_prefac_2_d12 * 2.0 * grad_trace12_v_other_e.col(idx);
        
        return Eigen::Vector3d( (grad_fac1_d1 + grad_fac1_d2 + grad_fac1_d12 + grad_fac2_d1 + grad_fac2_d2 + grad_fac2_d12)*helper.area );
    }
    
    inline Real evaluateGradMaterialNormPhi(const int idx)
    {
        const Real grad_fac1_d1 = helper.material_prefac_1_d1 * 2.0 * trace[0] * grad_trace_e[0](idx);
        const Real grad_fac1_d2 = helper.material_prefac_1_d2 * 2.0 * trace[1] * grad_trace_e[1](idx);
        const Real grad_fac1_d12 = helper.material_prefac_1_d12 * 2.0 * (trace[0] * grad_trace_e[1](idx) + trace[1] * grad_trace_e[0](idx));
        
        const Real grad_fac2_d1 = helper.material_prefac_2_d1 * grad_tracesq_e[0](idx);
        const Real grad_fac2_d2 = helper.material_prefac_2_d2 * grad_tracesq_e[1](idx);
        const Real grad_fac2_d12 = helper.material_prefac_2_d12 * 2.0 * grad_trace12_e(idx);
        
        return (grad_fac1_d1 + grad_fac1_d2 + grad_fac1_d12 + grad_fac2_d1 + grad_fac2_d2 + grad_fac2_d12)*helper.area;
    }
};

template<bool withGradient>
struct EnergyNormData_Ortho_Stretching : std::conditional<withGradient, EnergyNormData_Ortho_WithGradient_Verts, EnergyNormData_Ortho>::type
{
    typedef typename std::conditional<withGradient, EnergyNormData_Ortho_WithGradient_Verts, EnergyNormData_Ortho>::type EnergyNormData_Ortho_StretchingBase;
    
    StrainData_Ortho_Stretching<withGradient> & straindata;
    
    EnergyNormData_Ortho_Stretching(const EnergyHelper_Parametric_Ortho & helper_in, StrainData_Ortho_Stretching<withGradient> & straindata_in):
    EnergyNormData_Ortho_StretchingBase(helper_in),
    straindata(straindata_in)
    {}
    
    void compute()
    {
        straindata.compute();
        for(int c=0;c<2;++c)
        {
            this->trace[c] = straindata.strain[c].trace();
            this->trace_sq[c] = (straindata.strain[c] * straindata.strain[c]).trace();
        }

        this->trace_12 = (straindata.strain[0] * straindata.strain[1]).trace();
    }
    
    template<bool U = withGradient> typename std::enable_if<U, void>::type
    compute_grad()
    {
        straindata.compute_grad();
        for(int c=0;c<2;++c)
        {
            for(int d=0;d<3;++d)
            {
                this->grad_trace_v[c](d,0) = (straindata.grad_strain_v0[c][d]).trace();
                this->grad_trace_v[c](d,1) = (straindata.grad_strain_v1[c][d]).trace();
                this->grad_trace_v[c](d,2) = (straindata.grad_strain_v2[c][d]).trace();
                
                this->grad_tracesq_v[c](d,0) = (straindata.strain[c] * straindata.grad_strain_v0[c][d] + straindata.grad_strain_v0[c][d] * straindata.strain[c]).trace();
                this->grad_tracesq_v[c](d,1) = (straindata.strain[c] * straindata.grad_strain_v1[c][d] + straindata.grad_strain_v1[c][d] * straindata.strain[c]).trace();
                this->grad_tracesq_v[c](d,2) = (straindata.strain[c] * straindata.grad_strain_v2[c][d] + straindata.grad_strain_v2[c][d] * straindata.strain[c]).trace();
            }
        }
        for(int d=0;d<3;++d)
        {
            this->grad_trace12_v(d,0) = (straindata.strain[0] * straindata.grad_strain_v0[1][d] + straindata.grad_strain_v0[0][d] * straindata.strain[1]).trace();
            this->grad_trace12_v(d,1) = (straindata.strain[0] * straindata.grad_strain_v1[1][d] + straindata.grad_strain_v1[0][d] * straindata.strain[1]).trace();
            this->grad_trace12_v(d,2) = (straindata.strain[0] * straindata.grad_strain_v2[1][d] + straindata.grad_strain_v2[0][d] * straindata.strain[1]).trace();
        }
    }
};

template<bool withGradient>
struct EnergyNormData_Ortho_Bending : std::conditional<withGradient, EnergyNormData_Ortho_WithGradient, EnergyNormData_Ortho>::type
{
    typedef typename std::conditional<withGradient, EnergyNormData_Ortho_WithGradient, EnergyNormData_Ortho>::type EnergyNormData_Ortho_BendingBase;
    
    StrainData_Ortho_Bending<withGradient> & straindata;
    
    EnergyNormData_Ortho_Bending(const EnergyHelper_Parametric_Ortho & helper_in, StrainData_Ortho_Bending<withGradient> & straindata_in):
    EnergyNormData_Ortho_BendingBase(helper_in),
    straindata(straindata_in)
    {}
    
    void compute()
    {
        straindata.compute();
        for(int c=0;c<2;++c)
        {
            this->trace[c] = straindata.strain[c].trace();
            this->trace_sq[c] = (straindata.strain[c] * straindata.strain[c]).trace();
        }
        
        this->trace_12 = (straindata.strain[0] * straindata.strain[1]).trace();
    }
    
    template<bool U = withGradient> typename std::enable_if<U, void>::type
    compute_grad()
    {
        straindata.compute_grad();
        for(int c=0;c<2;++c)
        {
            for(int d=0;d<3;++d)
            {
                this->grad_trace_v[c](d,0) = (straindata.grad_strain_v0[c][d]).trace();
                this->grad_trace_v[c](d,1) = (straindata.grad_strain_v1[c][d]).trace();
                this->grad_trace_v[c](d,2) = (straindata.grad_strain_v2[c][d]).trace();
                
                this->grad_tracesq_v[c](d,0) = (straindata.strain[c] * straindata.grad_strain_v0[c][d] + straindata.grad_strain_v0[c][d] * straindata.strain[c]).trace();
                this->grad_tracesq_v[c](d,1) = (straindata.strain[c] * straindata.grad_strain_v1[c][d] + straindata.grad_strain_v1[c][d] * straindata.strain[c]).trace();
                this->grad_tracesq_v[c](d,2) = (straindata.strain[c] * straindata.grad_strain_v2[c][d] + straindata.grad_strain_v2[c][d] * straindata.strain[c]).trace();
                
                this->grad_trace_v_other_e[c](d,0) = (straindata.grad_strain_v_other_e0[c][d]).trace();
                this->grad_trace_v_other_e[c](d,1) = (straindata.grad_strain_v_other_e1[c][d]).trace();
                this->grad_trace_v_other_e[c](d,2) = (straindata.grad_strain_v_other_e2[c][d]).trace();
                
                this->grad_tracesq_v_other_e[c](d,0) = (straindata.strain[c] * straindata.grad_strain_v_other_e0[c][d] + straindata.grad_strain_v_other_e0[c][d] * straindata.strain[c]).trace();
                this->grad_tracesq_v_other_e[c](d,1) = (straindata.strain[c] * straindata.grad_strain_v_other_e1[c][d] + straindata.grad_strain_v_other_e1[c][d] * straindata.strain[c]).trace();
                this->grad_tracesq_v_other_e[c](d,2) = (straindata.strain[c] * straindata.grad_strain_v_other_e2[c][d] + straindata.grad_strain_v_other_e2[c][d] * straindata.strain[c]).trace();
                
            }
            
            this->grad_trace_e[c](0) = (straindata.grad_strain_e0[c]).trace();
            this->grad_trace_e[c](1) = (straindata.grad_strain_e1[c]).trace();
            this->grad_trace_e[c](2) = (straindata.grad_strain_e2[c]).trace();
            
            this->grad_tracesq_e[c](0) = (straindata.strain[c] * straindata.grad_strain_e0[c] + straindata.grad_strain_e0[c] * straindata.strain[c]).trace();
            this->grad_tracesq_e[c](1) = (straindata.strain[c] * straindata.grad_strain_e1[c] + straindata.grad_strain_e1[c] * straindata.strain[c]).trace();
            this->grad_tracesq_e[c](2) = (straindata.strain[c] * straindata.grad_strain_e2[c] + straindata.grad_strain_e2[c] * straindata.strain[c]).trace();
        }
        
        
        for(int d=0;d<3;++d)
        {
            this->grad_trace12_v(d,0) = (straindata.strain[0] * straindata.grad_strain_v0[1][d] + straindata.grad_strain_v0[0][d] * straindata.strain[1]).trace();
            this->grad_trace12_v(d,1) = (straindata.strain[0] * straindata.grad_strain_v1[1][d] + straindata.grad_strain_v1[0][d] * straindata.strain[1]).trace();
            this->grad_trace12_v(d,2) = (straindata.strain[0] * straindata.grad_strain_v2[1][d] + straindata.grad_strain_v2[0][d] * straindata.strain[1]).trace();
            
            this->grad_trace12_v_other_e(d,0) = (straindata.strain[0] * straindata.grad_strain_v_other_e0[1][d] + straindata.grad_strain_v_other_e0[0][d] * straindata.strain[1]).trace();
            this->grad_trace12_v_other_e(d,1) = (straindata.strain[0] * straindata.grad_strain_v_other_e1[1][d] + straindata.grad_strain_v_other_e1[0][d] * straindata.strain[1]).trace();
            this->grad_trace12_v_other_e(d,2) = (straindata.strain[0] * straindata.grad_strain_v_other_e2[1][d] + straindata.grad_strain_v_other_e2[0][d] * straindata.strain[1]).trace();
        }

        this->grad_trace12_e(0) = (straindata.strain[0] * straindata.grad_strain_e0[1] + straindata.grad_strain_e0[0] * straindata.strain[1]).trace();
        this->grad_trace12_e(1) = (straindata.strain[0] * straindata.grad_strain_e1[1] + straindata.grad_strain_e1[0] * straindata.strain[1]).trace();
        this->grad_trace12_e(2) = (straindata.strain[0] * straindata.grad_strain_e2[1] + straindata.grad_strain_e2[0] * straindata.strain[1]).trace();
    }
};

// WIM STOPPED HERE : NEED TO ADD grad_tracesq BELOW AND THEN FIX THE EVALUATENORM METHODS -- THEN BILAYER!!

struct EnergyNormData_Ortho_BilayerBase
{
    const EnergyHelper_Parametric_Ortho & helper;
    
    std::array<std::array<Real, 2>, 2> trace_12; // {{d1d1, d1d2},{d2d1, d2d2}}
    
    EnergyNormData_Ortho_BilayerBase(const EnergyHelper_Parametric_Ortho & helper_in):
    helper(helper_in)
    {}
};

struct EnergyNormData_Ortho_BilayerBase_WithGradient : EnergyNormData_Ortho_BilayerBase
{
    std::array<std::array<Eigen::Matrix3d, 2>, 2> grad_trace12_v;
    std::array<std::array<Eigen::Matrix3d, 2>, 2> grad_trace12_v_other_e;
    std::array<std::array<Eigen::Vector3d, 2>, 2> grad_trace12_e;
    
    
    EnergyNormData_Ortho_BilayerBase_WithGradient(const EnergyHelper_Parametric_Ortho & helper_in):
    EnergyNormData_Ortho_BilayerBase(helper_in)
    {}
};


template<bool withGradient>
struct EnergyNormData_Ortho_Bilayer : std::conditional<withGradient, EnergyNormData_Ortho_BilayerBase_WithGradient, EnergyNormData_Ortho_BilayerBase>::type
{
    typedef typename std::conditional<withGradient, EnergyNormData_Ortho_BilayerBase_WithGradient, EnergyNormData_Ortho_BilayerBase>::type tEnergyNormData_Ortho_BilayerBase;
    
    const EnergyNormData_Ortho_Stretching<withGradient> & energydata_stretching;
    const EnergyNormData_Ortho_Bending<withGradient> & energydata_bending;
    
    EnergyNormData_Ortho_Bilayer(const EnergyHelper_Parametric_Ortho & helper_in, const EnergyNormData_Ortho_Stretching<withGradient> & energydata_stretching_in, const EnergyNormData_Ortho_Bending<withGradient> & energydata_bending_in):
    tEnergyNormData_Ortho_BilayerBase(helper_in),
    energydata_stretching(energydata_stretching_in),
    energydata_bending(energydata_bending_in)
    {}
    
    void compute()
    {
        const StrainData_Ortho_Stretching<withGradient> & straindata_stretching = energydata_stretching.straindata;
        const StrainData_Ortho_Bending<withGradient> & straindata_bending = energydata_bending.straindata;
        
        // assume energydata_stretching and energydata_bending have already been computed
        for(int c1=0;c1<2;++c1)
            for(int c2=0;c2<2;++c2)
                this->trace_12[c1][c2] = (straindata_stretching.strain[c1] * straindata_bending.strain[c2]).trace();
    }
    
    template<bool U = withGradient> typename std::enable_if<U, void>::type
    compute_grad()
    {
        const StrainData_Ortho_Stretching<withGradient> & straindata_stretching = energydata_stretching.straindata;
        const StrainData_Ortho_Bending<withGradient> & straindata_bending = energydata_bending.straindata;
        
        // assume energydata_stretching and energydata_bending have already been computed
        for(int c1=0;c1<2;++c1)
            for(int c2=0;c2<2;++c2)
            {
                for(int d=0;d<3;++d)
                {
                    this->grad_trace12_v[c1][c2](d,0) = (straindata_stretching.strain[c1] * straindata_bending.grad_strain_v0[c2][d] + straindata_stretching.grad_strain_v0[c1][d] * straindata_bending.strain[c2]).trace();
                    this->grad_trace12_v[c1][c2](d,1) = (straindata_stretching.strain[c1] * straindata_bending.grad_strain_v1[c2][d] + straindata_stretching.grad_strain_v1[c1][d] * straindata_bending.strain[c2]).trace();
                    this->grad_trace12_v[c1][c2](d,2) = (straindata_stretching.strain[c1] * straindata_bending.grad_strain_v2[c2][d] + straindata_stretching.grad_strain_v2[c1][d] * straindata_bending.strain[c2]).trace();
                    
                    // no 'other' contribution from straindata_stretching
                    this->grad_trace12_v_other_e[c1][c2](d,0) = (straindata_stretching.strain[c1] * straindata_bending.grad_strain_v_other_e0[c2][d]).trace();
                    this->grad_trace12_v_other_e[c1][c2](d,1) = (straindata_stretching.strain[c1] * straindata_bending.grad_strain_v_other_e1[c2][d]).trace();
                    this->grad_trace12_v_other_e[c1][c2](d,2) = (straindata_stretching.strain[c1] * straindata_bending.grad_strain_v_other_e2[c2][d]).trace();
                }
        
                // no 'phi' contribution from straindata_stretching
                this->grad_trace12_e[c1][c2](0) = (straindata_stretching.strain[c1] * straindata_bending.grad_strain_e0[c2]).trace();
                this->grad_trace12_e[c1][c2](1) = (straindata_stretching.strain[c1] * straindata_bending.grad_strain_e1[c2]).trace();
                this->grad_trace12_e[c1][c2](2) = (straindata_stretching.strain[c1] * straindata_bending.grad_strain_e2[c2]).trace();
            }
    }
    
    inline Real evaluateMaterialNorm()
    {
        const Real fac1_d1  = this->helper.material_prefac_1_d1 * energydata_stretching.trace[0] * energydata_bending.trace[0];
        const Real fac1_d2  = this->helper.material_prefac_1_d2 * energydata_stretching.trace[1] * energydata_bending.trace[1];
        const Real fac1_d12 = this->helper.material_prefac_1_d12 * (energydata_stretching.trace[0] * energydata_bending.trace[1] + energydata_stretching.trace[1] * energydata_bending.trace[0]); // (trace0 * trace1 + trace1 * trace0)
        
        const Real fac2_d1 = this->helper.material_prefac_2_d1 * this->trace_12[0][0];
        const Real fac2_d2 = this->helper.material_prefac_2_d2 * this->trace_12[1][1];
        const Real fac2_d12 = this->helper.material_prefac_2_d12 * (this->trace_12[0][1] + this->trace_12[1][0]);
        
        return (fac1_d1 + fac1_d2 + fac1_d12 + fac2_d1 + fac2_d2 + fac2_d12)*this->helper.area;
    }
    
    template<bool U = withGradient> typename std::enable_if<U, Eigen::Vector3d>::type
    evaluateGradMaterialNorm(const int idx)
    {
        const Eigen::Vector3d grad_fac1_d1 = this->helper.material_prefac_1_d1 * (energydata_stretching.trace[0] * energydata_bending.grad_trace_v[0].col(idx) + energydata_stretching.grad_trace_v[0].col(idx) * energydata_bending.trace[0]);
        const Eigen::Vector3d grad_fac1_d2 = this->helper.material_prefac_1_d2 * (energydata_stretching.trace[1] * energydata_bending.grad_trace_v[1].col(idx) + energydata_stretching.grad_trace_v[1].col(idx) * energydata_bending.trace[1]);
        
        const Eigen::Vector3d grad_fac1_d12 = this->helper.material_prefac_1_d12 * (energydata_stretching.trace[0] * energydata_bending.grad_trace_v[1].col(idx) + energydata_stretching.grad_trace_v[0].col(idx) * energydata_bending.trace[1] + energydata_stretching.trace[1] * energydata_bending.grad_trace_v[0].col(idx) + energydata_stretching.grad_trace_v[1].col(idx) * energydata_bending.trace[0]);
        
        const Eigen::Vector3d grad_fac2_d1 = this->helper.material_prefac_2_d1 * this->grad_trace12_v[0][0].col(idx);
        const Eigen::Vector3d grad_fac2_d2 = this->helper.material_prefac_2_d2 * this->grad_trace12_v[1][1].col(idx);
        const Eigen::Vector3d grad_fac2_d12 = this->helper.material_prefac_2_d12 * (this->grad_trace12_v[0][1].col(idx) + this->grad_trace12_v[1][0].col(idx));
        
        return Eigen::Vector3d( (grad_fac1_d1 + grad_fac1_d2 + grad_fac1_d12 + grad_fac2_d1 + grad_fac2_d2 + grad_fac2_d12)*this->helper.area );
    }
    
    template<bool U = withGradient> typename std::enable_if<U, Eigen::Vector3d>::type
    evaluateGradMaterialNormOther(const int idx)
    {
        const Eigen::Vector3d grad_fac1_d1 = this->helper.material_prefac_1_d1 *
        (energydata_stretching.trace[0] * energydata_bending.grad_trace_v_other_e[0].col(idx));
        const Eigen::Vector3d grad_fac1_d2 = this->helper.material_prefac_1_d2 *
        (energydata_stretching.trace[1] * energydata_bending.grad_trace_v_other_e[1].col(idx));
        
        const Eigen::Vector3d grad_fac1_d12 = this->helper.material_prefac_1_d12 *
        (energydata_stretching.trace[0] * energydata_bending.grad_trace_v_other_e[1].col(idx) +
         energydata_stretching.trace[1] * energydata_bending.grad_trace_v_other_e[0].col(idx));
        
        const Eigen::Vector3d grad_fac2_d1 = this->helper.material_prefac_2_d1 * this->grad_trace12_v_other_e[0][0].col(idx);
        const Eigen::Vector3d grad_fac2_d2 = this->helper.material_prefac_2_d2 * this->grad_trace12_v_other_e[1][1].col(idx);
        const Eigen::Vector3d grad_fac2_d12 = this->helper.material_prefac_2_d12 * (this->grad_trace12_v_other_e[0][1].col(idx) + this->grad_trace12_v_other_e[1][0].col(idx));
        
        return Eigen::Vector3d( (grad_fac1_d1 + grad_fac1_d2 + grad_fac1_d12 + grad_fac2_d1 + grad_fac2_d2 + grad_fac2_d12)*this->helper.area );
    }
    
    template<bool U = withGradient> typename std::enable_if<U, Real>::type
    evaluateGradMaterialNormPhi(const int idx)
    {
        const Real grad_fac1_d1 = this->helper.material_prefac_1_d1 *
        (energydata_stretching.trace[0] * energydata_bending.grad_trace_e[0](idx));
        const Real grad_fac1_d2 = this->helper.material_prefac_1_d2 *
        (energydata_stretching.trace[1] * energydata_bending.grad_trace_e[1](idx));
        
        const Real grad_fac1_d12 = this->helper.material_prefac_1_d12 *
        (energydata_stretching.trace[0] * energydata_bending.grad_trace_e[1](idx) +
         energydata_stretching.trace[1] * energydata_bending.grad_trace_e[0](idx));
        
        const Real grad_fac2_d1 = this->helper.material_prefac_2_d1 * this->grad_trace12_e[0][0](idx);
        const Real grad_fac2_d2 = this->helper.material_prefac_2_d2 * this->grad_trace12_e[1][1](idx);
        const Real grad_fac2_d12 = this->helper.material_prefac_2_d12 * (this->grad_trace12_e[0][1](idx) + this->grad_trace12_e[1][0](idx));
        
        return (grad_fac1_d1 + grad_fac1_d2 + grad_fac1_d12 + grad_fac2_d1 + grad_fac2_d2 + grad_fac2_d12)*this->helper.area;
    }
};

#include "EnergyHelper_Parametric.hpp"

template<bool withGradient, MeshLayer layer>
struct SaintVenantEnergy<withGradient, Material_Orthotropic, layer> : std::conditional<withGradient,
typename std::conditional<layer == single, SaintVenantEnergyData_WithGradient, SaintVenantEnergyData_WithGradient_Bilayer>::type, SaintVenantEnergyData>::type
{
    const Real h_aa, h_bb, h_ab;
    EnergyHelper_Parametric_Ortho helper;
    
    FirstFundamentalForm<withGradient> firstFF;
    SecondFundamentalForm<withGradient> secondFF;
    
    StrainData_Ortho_Stretching<withGradient> straindata_stretching;
    StrainData_Ortho_Bending<withGradient> straindata_bending;
    
    EnergyNormData_Ortho_Stretching<withGradient> energydata_stretching;
    EnergyNormData_Ortho_Bending<withGradient> energydata_bending;
    EnergyNormData_Ortho_Bilayer<withGradient> energydata_bilayer;
    
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
    
    SaintVenantEnergy(const Material_Orthotropic & mat_prop_in, const Eigen::Matrix2d & aform_bar_d1_in, const Eigen::Matrix2d & aform_bar_d2_in, const ExtendedTriangleInfo & info_in):
    h_aa(compute_h_aa(mat_prop_in.getThickness())),
    h_bb(compute_h_bb(mat_prop_in.getThickness())),
    h_ab(compute_h_ab(mat_prop_in.getThickness())),
    helper(mat_prop_in, aform_bar_d1_in, aform_bar_d2_in),
    firstFF(info_in),
    secondFF(info_in),
    straindata_stretching(helper, firstFF),
    straindata_bending(helper, secondFF),
    energydata_stretching(helper, straindata_stretching),
    energydata_bending(helper, straindata_bending),
    energydata_bilayer(helper, energydata_stretching, energydata_bending)
    {}
    
    SaintVenantEnergy(const Material_Orthotropic & mat_prop_in, const Eigen::Matrix2d & aform_bar_d1_in, const Eigen::Matrix2d & aform_bar_d2_in, const Eigen::Matrix2d & bform_bar_in, const ExtendedTriangleInfo & info_in):
    h_aa(compute_h_aa(mat_prop_in.getThickness())),
    h_bb(compute_h_bb(mat_prop_in.getThickness())),
    h_ab(compute_h_ab(mat_prop_in.getThickness())),
    helper(mat_prop_in, aform_bar_d1_in, aform_bar_d2_in),
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

#endif /* EnergyHelper_Parametric_Ortho_hpp */
