//
//  MaterialProperties.hpp
//  Elasticity
//
//  Created by Wim van Rees on 1/6/17.
//  Copyright © 2017 Wim van Rees. All rights reserved.
//

#ifndef MaterialProperties_h
#define MaterialProperties_h

#include "common.hpp"

struct Material_Isotropic
{
    Real Young; // we need to change Young in the inverse-design growth process where we use Young as design parameter    
    const Real Poisson;
    Real thickness; // we need to change thickness in Derivative_FD_Inverse::computeGradient_thickness
    const Real density;
    
    Material_Isotropic(const Real Y, const Real P, const Real t, const Real rho = -1):
    Young(Y),
    Poisson(P),
    thickness(t),
    density(rho)
    {}
    
    Real getThickness() const
    {
        return thickness;
    }
    
    Real getDensity() const
    {
        return density;
    }

    Real getYoung() const
    {
        return Young;
    }

    Real getStVenantFactor1() const
    {
        return 0.5 * Young*Poisson/(1.0-Poisson*Poisson); // alpha/2
    }
    
    Real getStVenantFactor2() const
    {
        return  0.5 * Young/(1.0+Poisson); // beta
    }
};

struct Material_Orthotropic
{
    Real Young_d1, Young_d2, G_12; // we need to be able to change Young since there are no methods for array input of Young in Ortho_Array
    Real nu_d1, nu_d2;
    Real thickness;  // we need to be able to change thickness since there are no methods for array input of thickness in Ortho_Array
    Real density;
    Eigen::Matrix2d aform_bar_d1; // aform_bar_d2 = aform_bar - aform_bar_d1
    
    Material_Orthotropic(const Real Y1, const Real Y2, const Real G, const Real nu_d1, const Real t, const Real rho = -1):
    Young_d1(Y1),
    Young_d2(Y2),
    G_12(G),
    nu_d1(nu_d1),
    nu_d2(nu_d1*Young_d2/Young_d1),
    thickness(t),
    density(rho)
    {}
    
    Eigen::Matrix2d get_aform_bar_d1() const
    {
        return aform_bar_d1;
    }
    
    Real getThickness() const
    {
        return thickness;
    }
    
    Real getDensity() const
    {
        return density;
    }
    
    Real getStVenantFactor_1_d1() const
    {
        return 0.5 * nu_d1*Young_d1 / (1.0 - nu_d1 * nu_d2);
    }
    
    Real getStVenantFactor_1_d2() const
    {
        return 0.5 * nu_d1*Young_d2 / (1.0 - nu_d1 * nu_d2);
         // == 0.5 * Y1 * nu_2 / (1.0 - nu_d1 * nu_d2)
    }
    
    Real getStVenantFactor_1_d12() const
    {
        return 0.25*(Young_d1*nu_d2 + Young_d2*nu_d1) / (1.0 - nu_d1 * nu_d2); // extra 1/2 because we split the 2ab term into ab + ba
    }
    
    Real getStVenantFactor_2_d1() const
    {
        return 0.5 * (1.0 - nu_d1) * Young_d1 / (1.0 - nu_d1 * nu_d2);
    }
    
    Real getStVenantFactor_2_d2() const
    {
        return 0.5 * (1.0 - nu_d1) * Young_d2 / (1.0 - nu_d1 * nu_d2);
    }
    
    Real getStVenantFactor_2_d12() const
    {
        return G_12; // for isotropic : this would be E/(2(1+nu)), same as factor_2 above
    }
};



template<typename tMaterialType>
struct MaterialProperties
{
    virtual tMaterialType getFaceMaterial(const int face_idx) const = 0;
};


struct MaterialProperties_Iso_Constant : MaterialProperties<Material_Isotropic>
{
    typedef Material_Isotropic tMaterialType;
    Material_Isotropic material; // only one for the entire mesh
    
    MaterialProperties_Iso_Constant(const Real Y, const Real P, const Real t, const Real rho = -1):
    material(Y, P, t, rho)
    {}

    Material_Isotropic getFaceMaterial(const int ) const
    {
        return material;
    }
};

struct MaterialProperties_Iso_Array : MaterialProperties<Material_Isotropic>
{
    typedef Material_Isotropic tMaterialType;
    std::vector<Material_Isotropic> materials;
    
    
    MaterialProperties_Iso_Array(const Real E_in, const Real P_in, const Eigen::VectorXd & t_in)
    {
        const int nFaces = t_in.rows();
        materials.reserve(nFaces);
        for(int i=0;i<nFaces;++i)
            materials.emplace_back(E_in, P_in, t_in(i));
    }
    
    MaterialProperties_Iso_Array(const Eigen::VectorXd & E_in, const Real P_in, const Eigen::VectorXd & t_in)
    {
        const int nFaces = t_in.rows();
        assert(E_in.rows() == nFaces);
        
        materials.reserve(nFaces);
        for(int i=0;i<nFaces;++i)
            materials.emplace_back(E_in(i), P_in, t_in(i));
    }
    
    MaterialProperties_Iso_Array(const Eigen::VectorXd & E_in, const Real P_in, const Real t_in)
    {
        const int nFaces = E_in.rows();
        
        materials.reserve(nFaces);
        for(int i=0;i<nFaces;++i)
            materials.emplace_back(E_in(i), P_in, t_in);
    }
    
    MaterialProperties_Iso_Array(const Eigen::VectorXd & E_in, const Real P_in, const Real t_in, const Eigen::VectorXd & rho_in)
    {
        const int nFaces = E_in.rows();
        
        materials.reserve(nFaces);
        for(int i=0;i<nFaces;++i)
            materials.emplace_back(E_in(i), P_in, t_in, rho_in(i));
    }
    
    Material_Isotropic getFaceMaterial(const int face_idx) const
    {
        return materials[face_idx];
    }
};

struct MaterialProperties_Ortho_Constant : MaterialProperties<Material_Orthotropic>
{
    typedef Material_Orthotropic tMaterialType;
    Material_Orthotropic material; // only one for the entire mesh
    
    MaterialProperties_Ortho_Constant(const Real Y1_in, const Real Y2_in, const Real G_in, const Real nu_d1_in, const Real t_in, const Real rho_in = -1):
    material(Y1_in, Y2_in, G_in, nu_d1_in, t_in, rho_in)
    {
    }
    
    Material_Orthotropic getFaceMaterial(const int ) const
    {
        return material;
    }
};

struct MaterialProperties_Ortho_Array : MaterialProperties<Material_Orthotropic>
{
    typedef Material_Orthotropic tMaterialType;
    std::vector<Material_Orthotropic> materials;
    
    MaterialProperties_Ortho_Array(const int nFaces, const Real Y1_in, const Real Y2_in, const Real G_in, const Real nu_d1_in, const Real t_in, const Real rho_in = -1)
    {
        materials.reserve(nFaces);
        for(int i=0;i<nFaces;++i)
            materials.emplace_back(Y1_in, Y2_in, G_in, nu_d1_in, t_in, rho_in);
    }
    
    Material_Orthotropic getFaceMaterial(const int face_idx) const
    {
        return materials[face_idx];
    }
};

#endif /* MaterialProperties_h */
