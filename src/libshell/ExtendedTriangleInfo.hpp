//
//  ExtendedTriangleInfo.hpp
//  Elasticity
//
//  Created by Wim van Rees on 1/2/17.
//  Copyright © 2017 Wim van Rees. All rights reserved.
//

#ifndef ExtendedTriangleInfo_hpp
#define ExtendedTriangleInfo_hpp

#include "common.hpp"

struct QuadraticFormGradientData_Verts
{
    Eigen::Vector3d gradv0_11, gradv0_12, gradv0_22;
    Eigen::Vector3d gradv1_11, gradv1_12, gradv1_22;
    Eigen::Vector3d gradv2_11, gradv2_12, gradv2_22;
};

struct QuadraticFormGradientData : QuadraticFormGradientData_Verts
{
    int idx_v_other_e0;
    int idx_v_other_e1;
    int idx_v_other_e2;
    
    bool isInterior_e0, isInterior_e1, isInterior_e2;
    
    Real gradphi_e0_11, gradphi_e0_12, gradphi_e0_22;
    Real gradphi_e1_11, gradphi_e1_12, gradphi_e1_22;
    Real gradphi_e2_11, gradphi_e2_12, gradphi_e2_22;
    
    Eigen::Vector3d gradv_other_e0_11, gradv_other_e0_12, gradv_other_e0_22;
    Eigen::Vector3d gradv_other_e1_11, gradv_other_e1_12, gradv_other_e1_22;
    Eigen::Vector3d gradv_other_e2_11, gradv_other_e2_12, gradv_other_e2_22;
};

struct ExtendedTriangleInfo
{
    int face_idx;
    
    int idx_v0, idx_v1, idx_v2; /*!< vertex indices */
    int idx_e0, idx_e1, idx_e2; /*!< edge indices */
    
    /*! vertex locations */
    Eigen::Vector3d v0;
    Eigen::Vector3d v1;
    Eigen::Vector3d v2;
    
    /*! edge vectors */
    Eigen::Vector3d e0;
    Eigen::Vector3d e1;
    Eigen::Vector3d e2;
    
    Real double_face_area;
    Eigen::Vector3d edgelength, sign, theta, phi, alpha, height, cosGamma;
    Eigen::Vector3d face_normal;
    Eigen::Matrix3d edge_normals;
    
    
    Eigen::Vector3b v0_fixed, v1_fixed, v2_fixed;
    bool e0_clamped, e1_clamped, e2_clamped;
    Eigen::Vector3d clamped_normal_e0, clamped_normal_e1, clamped_normal_e2;
    
    std::array<ExtendedTriangleInfo *, 3> other_faces;
    Eigen::Vector3i other_face_edge_idx;
    
    template<int vidx>
    Eigen::Matrix3d get_gradv_normal() const
    {
        const int vidx_plus = (vidx + 1)%3;
        return 1.0/height(vidx_plus) * edge_normals.col(vidx_plus) * face_normal.transpose();
    }
    
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    
    ExtendedTriangleInfo():
    other_faces({nullptr, nullptr, nullptr})
    {}
    
    void setTheta(const int eidx, const Real theta_val)
    {
        theta(eidx) = theta_val;
        alpha(eidx) = 0.5*theta(eidx) + sign(eidx)*phi(eidx);
        
    }
    
    template<int vidx, int hidx>
    Eigen::Vector3d get_gradv_height() const
    {
        // can we not just use double_face_area / length(i) and take derivative of that??
        if((vidx+1)%3 == hidx) // compile-time conditional
            return -edge_normals.col(hidx);
        
        const int gammaidx = ((hidx+1)%3) == (vidx+1)%3 ? (hidx+2)%3 : (hidx+1)%3;
        
        return height(hidx)/height((vidx+1)%3) * cosGamma(gammaidx) * edge_normals.col(hidx);
    }
    
/*
    void get_gradv_height0(Eigen::Vector3d & gradv0, Eigen::Vector3d & gradv1, Eigen::Vector3d & gradv2) const
    {
        gradv0 =  height(0)/height(1) * cosGamma(2) * edge_normals.col(0);
        gradv1 =  height(0)/height(2) * cosGamma(1) * edge_normals.col(0);
        gradv2 = -edge_normals.col(0);
    }
    
    void get_gradv_height1(Eigen::Vector3d & gradv0, Eigen::Vector3d & gradv1, Eigen::Vector3d & gradv2) const
    {
        gradv0 = -edge_normals.col(1);
        gradv1 =  height(1)/height(2) * cosGamma(0) * edge_normals.col(1);
        gradv2 =  height(1)/height(0) * cosGamma(2) * edge_normals.col(1);
    }
    
    void get_gradv_height2(Eigen::Vector3d & gradv0, Eigen::Vector3d & gradv1, Eigen::Vector3d & gradv2) const
    {
        gradv0 =  height(2)/height(1) * cosGamma(0) * edge_normals.col(2);
        gradv1 = -edge_normals.col(2);
        gradv2 =  height(2)/height(0) * cosGamma(1) * edge_normals.col(2);
    }
*/
    
    
    template<int eidx>
    bool get_gradv_theta(Eigen::Vector3d & gradv0, Eigen::Vector3d & gradv1, Eigen::Vector3d & gradv2, Eigen::Vector3d & other_gradv2, int & other_vidx) const
    {
        if(other_faces[eidx] == nullptr)
        {
            gradv0.setZero();
            gradv1.setZero();
            gradv2.setZero();
            other_gradv2.setZero();
            other_vidx = -1;
            return false;
        }
        
        const ExtendedTriangleInfo & other = *other_faces[eidx];
        
        // these are two triangles that share an edge denoted by eidx_me (0,1,2) for me and eidx_other (0,1,2) for the other
        
        // renumber them so that for each, edge0 is the shared edge
        const int idx_0 = (eidx + 0)%3;
        const int idx_1 = (eidx + 1)%3;
        const int idx_2 = (eidx + 2)%3;
        
        // number the other edge indices also consistently
        const int other_eidx = other_face_edge_idx(eidx);
        const int other_idx_0 = (other_eidx + 0)%3;
        const int other_idx_1 = (other_eidx + 1)%3;
        const int other_idx_2 = (other_eidx + 2)%3;
        
        /*
         // after renumbering we can draw the following picture:
         // the primes denote 'other'
         // gi = cos(gamma_i)
         // shared edge is e0 = v1-v0 and also e0' = v1' - v0'
         // the heights are the opposites to their edges
         
                 v2'
         
                 /\
                /  \
          e2'  / g0'\ e1'
              /      \
             /        \
         v0'/ g1'  g2' \ v1'
            ------------
         v1 \ g2   g1  / v0
             \        /
              \      /
           e1  \ g0 /  e2
                \  /
                 \/
         
                 v2
         */
        
        // compute the gradients as new variables but this should be optimized away by the compiler (hopefully - else just have to write messier code)
        const Eigen::Vector3d gradv0rel = cosGamma(idx_2)/height(idx_1) * face_normal + other.cosGamma(other_idx_1)/other.height(other_idx_2) * other.face_normal;
        const Eigen::Vector3d gradv1rel = cosGamma(idx_1)/height(idx_2) * face_normal + other.cosGamma(other_idx_2)/other.height(other_idx_1) * other.face_normal;
        const Eigen::Vector3d gradv2rel = -1.0/height(idx_0) * face_normal;
        const Eigen::Vector3d other_gradv2rel = -1.0/other.height(other_idx_0) * other.face_normal;
        
        // assign the gradients to the correct vertices (reverse the renumbering)
        // if eidx == 0 --> v0rel -> v0, v1rel -> v1, v2rel -> v2
        // if eidx == 1 --> v0rel -> v1, v1rel -> v2, v2rel -> v0
        // if eidx == 2 --> v0rel -> v2, v1rel -> v0, v2rel -> v1
        
        gradv0 = (eidx == 0 ? gradv0rel : (eidx == 1 ? gradv2rel : gradv1rel));
        gradv1 = (eidx == 0 ? gradv1rel : (eidx == 1 ? gradv0rel : gradv2rel));
        gradv2 = (eidx == 0 ? gradv2rel : (eidx == 1 ? gradv1rel : gradv0rel));
        other_gradv2 = other_gradv2rel;
        
        // assign the vertex index of v2rel of the other triangle
        other_vidx = (other_eidx == 0 ? other.idx_v2 : (other_eidx == 1 ? other.idx_v0 : other.idx_v1));
        
        return true;
    }
    
    Eigen::Vector3d computeFaceCenter() const
    {
        return (v0 + v1 + v2)/3.0;
    }
    
    template<int vidx>
    Eigen::Vector3d computeGradFaceCenter() const
    {
        return 1.0/3.0*Eigen::Vector3d::Constant(1);
    }
    
    Eigen::Vector3d computeFaceNormalFromDirectors() const
    {
        const Eigen::Vector3d midedge_e0 = std::cos(alpha(0))*face_normal + std::sin(alpha(0))*edge_normals.col(0);
        const Eigen::Vector3d midedge_e1 = std::cos(alpha(1))*face_normal + std::sin(alpha(1))*edge_normals.col(1);
        const Eigen::Vector3d midedge_e2 = std::cos(alpha(2))*face_normal + std::sin(alpha(2))*edge_normals.col(2);
        
        return ((midedge_e0 + midedge_e1 + midedge_e2)/3.0).normalized();
        
    }
    
    Eigen::Matrix2d computeFirstFundamentalForm() const
    {
        const Real e1_dot_e1 = edgelength(1)*edgelength(1);
        const Real e1_dot_e2 = e1.dot(e2);
        const Real e2_dot_e2 = edgelength(2)*edgelength(2);
        return (Eigen::Matrix2d() << e1_dot_e1, e1_dot_e2, e1_dot_e2, e2_dot_e2).finished();
    }
    
    Eigen::Matrix2d computeSecondFundamentalForm() const
    {
        const Real s_alpha_e0 = std::sin(alpha(0));
        const Real s_alpha_e1 = std::sin(alpha(1));
        const Real s_alpha_e2 = std::sin(alpha(2));
        
        // exception for clamped edges
        const Real n2_dot_e1 = (e2_clamped ? clamped_normal_e2.dot(e1) : +height(2)*s_alpha_e2);
        const Real n0_dot_e1 = (e0_clamped ? clamped_normal_e0.dot(e1) : -height(0)*s_alpha_e0);
        //        const Real n1_dot_e1 = 0.0;//(e1_clamped ? clamped_normal_e1.dot(e1) :  0.0); // clamped_normal_e1 should be normal to e1
        const Real n0_dot_e2 = (e0_clamped ? clamped_normal_e0.dot(e2) : -n0_dot_e1);
        const Real n1_dot_e2 = (e1_clamped ? clamped_normal_e1.dot(e2) : -height(1)*s_alpha_e1);
        
        return (Eigen::Matrix2d() << 2.0*(n0_dot_e1 - n2_dot_e1) , -2.0*n0_dot_e1, -2.0*n0_dot_e1, 2.0*(n1_dot_e2 - n0_dot_e2)).finished();
    }
    
    QuadraticFormGradientData_Verts computeGradFirstFundamentalForm() const
    {
        QuadraticFormGradientData_Verts gradFirstFF;
        
        // gradient with respect to vertices
        
        // derivatives of e1.dot(e1)
        gradFirstFF.gradv0_11.setZero();
        gradFirstFF.gradv1_11 = -2.0*e1;
        gradFirstFF.gradv2_11 = +2.0*e1;
        
        // derivatives of e1.dot(e2)
        gradFirstFF.gradv0_12 =  e1;
        gradFirstFF.gradv1_12 = -e2;
        gradFirstFF.gradv2_12 =  e2 - e1;
        
        // derivatives of e2.dot(e2)
        gradFirstFF.gradv0_22 = +2.0*e2;
        gradFirstFF.gradv1_22.setZero();
        gradFirstFF.gradv2_22 = -2.0*e2;
        
        return gradFirstFF;
    }
    
    QuadraticFormGradientData computeGradSecondFundamentalForm() const
    {
        const Real s_alpha_e0 = std::sin(alpha(0));
        const Real s_alpha_e1 = std::sin(alpha(1));
        const Real s_alpha_e2 = std::sin(alpha(2));
        
        QuadraticFormGradientData gradSecondFF;
        
        // gradient with respect to edge directors
        gradSecondFF.gradphi_e0_11 = -2*height(0)*std::cos(alpha(0))*sign(0);
        gradSecondFF.gradphi_e0_12 = +2*height(0)*std::cos(alpha(0))*sign(0);
        gradSecondFF.gradphi_e0_22 = -2*height(0)*std::cos(alpha(0))*sign(0);
        
        gradSecondFF.gradphi_e1_11 =  0.0;
        gradSecondFF.gradphi_e1_12 =  0.0;
        gradSecondFF.gradphi_e1_22 = -2*height(1)*std::cos(alpha(1))*sign(1);
        
        gradSecondFF.gradphi_e2_11 = -2*height(2)*std::cos(alpha(2))*sign(2);
        gradSecondFF.gradphi_e2_12 =  0.0;
        gradSecondFF.gradphi_e2_22 =  0.0;
     
        // gradient with respect to vertices
        // gradient of the triangle heights
        const Eigen::Vector3d gradv0_h0 = get_gradv_height<0,0>();
        const Eigen::Vector3d gradv1_h0 = get_gradv_height<1,0>();
        const Eigen::Vector3d gradv2_h0 = get_gradv_height<2,0>();
        
        const Eigen::Vector3d gradv0_h1 = get_gradv_height<0,1>();
        const Eigen::Vector3d gradv1_h1 = get_gradv_height<1,1>();
        const Eigen::Vector3d gradv2_h1 = get_gradv_height<2,1>();
        
        const Eigen::Vector3d gradv0_h2 = get_gradv_height<0,2>();
        const Eigen::Vector3d gradv1_h2 = get_gradv_height<1,2>();
        const Eigen::Vector3d gradv2_h2 = get_gradv_height<2,2>();
        
        // gradient of dihedral angles
        Eigen::Vector3d gradv0_theta_e0, gradv1_theta_e0, gradv2_theta_e0;
        Eigen::Vector3d gradv0_theta_e1, gradv1_theta_e1, gradv2_theta_e1;
        Eigen::Vector3d gradv0_theta_e2, gradv1_theta_e2, gradv2_theta_e2;
        Eigen::Vector3d gradv_other_theta_e0, gradv_other_theta_e1, gradv_other_theta_e2;
        
        gradSecondFF.isInterior_e0 = get_gradv_theta<0>(gradv0_theta_e0, gradv1_theta_e0, gradv2_theta_e0, gradv_other_theta_e0, gradSecondFF.idx_v_other_e0);
        gradSecondFF.isInterior_e1 = get_gradv_theta<1>(gradv0_theta_e1, gradv1_theta_e1, gradv2_theta_e1, gradv_other_theta_e1, gradSecondFF.idx_v_other_e1);
        gradSecondFF.isInterior_e2 = get_gradv_theta<2>(gradv0_theta_e2, gradv1_theta_e2, gradv2_theta_e2, gradv_other_theta_e2, gradSecondFF.idx_v_other_e2);
        
        // gradient of the bform entries
        const Real h0_cos_alpha0 = height(0)*std::cos(alpha(0));
        const Real h1_cos_alpha1 = height(1)*std::cos(alpha(1));
        const Real h2_cos_alpha2 = height(2)*std::cos(alpha(2));
        
        const Eigen::Vector3d gradv0_n2_dot_e1 = (e2_clamped ? Eigen::Vector3d::Zero() : Eigen::Vector3d(gradv0_h2*s_alpha_e2 + 0.5*h2_cos_alpha2*gradv0_theta_e2));
        const Eigen::Vector3d gradv1_n2_dot_e1 = (e2_clamped ? Eigen::Vector3d(-clamped_normal_e2) : Eigen::Vector3d(gradv1_h2*s_alpha_e2 + 0.5*h2_cos_alpha2*gradv1_theta_e2));
        const Eigen::Vector3d gradv2_n2_dot_e1 = (e2_clamped ?  clamped_normal_e2 : Eigen::Vector3d(gradv2_h2*s_alpha_e2 + 0.5*h2_cos_alpha2*gradv2_theta_e2));
        const Eigen::Vector3d gradv_other_e2_n2_dot_e1 = (e2_clamped ? Eigen::Vector3d::Zero() : Eigen::Vector3d(0.5*h2_cos_alpha2*gradv_other_theta_e2));
        
        const Eigen::Vector3d gradv0_n0_dot_e1 = (e0_clamped ? Eigen::Vector3d::Zero() : Eigen::Vector3d(-gradv0_h0*s_alpha_e0 - 0.5*h0_cos_alpha0*gradv0_theta_e0));
        const Eigen::Vector3d gradv1_n0_dot_e1 = (e0_clamped ? Eigen::Vector3d(-clamped_normal_e0) : Eigen::Vector3d(-gradv1_h0*s_alpha_e0 - 0.5*h0_cos_alpha0*gradv1_theta_e0));
        const Eigen::Vector3d gradv2_n0_dot_e1 = (e0_clamped ?  clamped_normal_e0 : Eigen::Vector3d(-gradv2_h0*s_alpha_e0 - 0.5*h0_cos_alpha0*gradv2_theta_e0));
        const Eigen::Vector3d gradv_other_e0_n0_dot_e1 = (e0_clamped ?  Eigen::Vector3d::Zero() : Eigen::Vector3d(-0.5*h0_cos_alpha0*gradv_other_theta_e0));
        
        const Eigen::Vector3d gradv0_n0_dot_e2 = (e0_clamped ?  clamped_normal_e0 : Eigen::Vector3d(-gradv0_n0_dot_e1));
        const Eigen::Vector3d gradv1_n0_dot_e2 = (e0_clamped ? Eigen::Vector3d::Zero() :  Eigen::Vector3d(-gradv1_n0_dot_e1));
        const Eigen::Vector3d gradv2_n0_dot_e2 = (e0_clamped ?  Eigen::Vector3d(-clamped_normal_e0) : Eigen::Vector3d(-gradv2_n0_dot_e1));
        const Eigen::Vector3d gradv_other_e0_n0_dot_e2 = (e0_clamped ?  Eigen::Vector3d::Zero() : Eigen::Vector3d(-gradv_other_e0_n0_dot_e1));
        
        const Eigen::Vector3d gradv0_n1_dot_e2 = (e1_clamped ?  clamped_normal_e1 : Eigen::Vector3d(-gradv0_h1*s_alpha_e1 - 0.5*h1_cos_alpha1*gradv0_theta_e1));
        const Eigen::Vector3d gradv1_n1_dot_e2 = (e1_clamped ? Eigen::Vector3d::Zero() :  Eigen::Vector3d(-gradv1_h1*s_alpha_e1 - 0.5*h1_cos_alpha1*gradv1_theta_e1));
        const Eigen::Vector3d gradv2_n1_dot_e2 = (e1_clamped ?  Eigen::Vector3d(-clamped_normal_e1) : Eigen::Vector3d(-gradv2_h1*s_alpha_e1 - 0.5*h1_cos_alpha1*gradv2_theta_e1));
        const Eigen::Vector3d gradv_other_e1_n1_dot_e2 = (e1_clamped ?  Eigen::Vector3d::Zero() : Eigen::Vector3d(- 0.5*h1_cos_alpha1*gradv_other_theta_e1));
        
        
        gradSecondFF.gradv0_11 = -2*(gradv0_n2_dot_e1 - gradv0_n0_dot_e1);
        gradSecondFF.gradv0_12 = -2*gradv0_n0_dot_e1;
        gradSecondFF.gradv0_22 = -2*(gradv0_n0_dot_e2 - gradv0_n1_dot_e2);
        
        gradSecondFF.gradv1_11 = -2*(gradv1_n2_dot_e1 - gradv1_n0_dot_e1);
        gradSecondFF.gradv1_12 = -2*gradv1_n0_dot_e1;
        gradSecondFF.gradv1_22 = -2*(gradv1_n0_dot_e2 - gradv1_n1_dot_e2);
        
        gradSecondFF.gradv2_11 = -2*(gradv2_n2_dot_e1 - gradv2_n0_dot_e1);
        gradSecondFF.gradv2_12 = -2*gradv2_n0_dot_e1;
        gradSecondFF.gradv2_22 = -2*(gradv2_n0_dot_e2 - gradv2_n1_dot_e2);
        
        gradSecondFF.gradv_other_e0_11 = +2*gradv_other_e0_n0_dot_e1;
        gradSecondFF.gradv_other_e1_11 = Eigen::Vector3d::Zero();
        gradSecondFF.gradv_other_e2_11 = -2*gradv_other_e2_n2_dot_e1;
        
        gradSecondFF.gradv_other_e0_12 = -2*gradv_other_e0_n0_dot_e1;
        gradSecondFF.gradv_other_e1_12 = Eigen::Vector3d::Zero();
        gradSecondFF.gradv_other_e2_12 = Eigen::Vector3d::Zero();
        
        gradSecondFF.gradv_other_e0_22 = -2*gradv_other_e0_n0_dot_e2;
        gradSecondFF.gradv_other_e1_22 = +2*gradv_other_e1_n1_dot_e2;
        gradSecondFF.gradv_other_e2_22 = Eigen::Vector3d::Zero();

        return gradSecondFF;
    }
};


#endif /* ExtendedTriangleInfo_h */
