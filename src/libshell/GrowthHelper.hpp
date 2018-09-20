//
//  GrowthHelper.hpp
//  Elasticity
//
//  Created by Wim van Rees on 8/17/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef GrowthHelper_hpp
#define GrowthHelper_hpp

#include "common.hpp"
#include "TriangleInfo.hpp"

#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>


class DecomposedGrowthState
{
protected:
    Real s1, s2;
    Eigen::Vector3d v1, v2;
    
    void decompose_metric(const Eigen::MatrixXd & rxy, const Eigen::Matrix2d & metric)
    {
        // decompute metric
        const Eigen::Matrix2d a_i = rxy.transpose() * rxy;
        
        // compute delta_metric
        const Eigen::Matrix2d delta_a = a_i.inverse() * metric;
        
        // eigensolver
        Eigen::EigenSolver<Eigen::Matrix2d> es(delta_a);
        const Eigen::Vector2d evals = es.eigenvalues().real();
        const Eigen::Matrix2d evecs = es.eigenvectors().real();
        
        // evecs is the V-matrix (generalized eigen-vectors) but not appropriately normalized yet -- normalize
        // const Eigen::Matrix2d shouldBeI = evecs.transpose() * a_i * evecs;
        // const Eigen::Matrix2d V = evecs * shouldBeI.inverse().sqrt();
        
        const Eigen::MatrixXd v12_evec = rxy * evecs;
        Eigen::Vector3d v1_evec = v12_evec.col(0).normalized();
        Eigen::Vector3d v2_evec = v12_evec.col(1).normalized();
        
        // note: for everything else here we assume that v1 and v2 are orthogonal
        // however if s1==s2 --> the values of v1 and v2 are arbitrary and so from whatever we do above, they might turn out to not be orthogonal
        // so, in this case we orthogonalize them by hand (they are already normalized)
        // but since we can still have them equal, instead we set them equal to rxy_1 and rxy_2 first and orthogonalize second
        if(std::abs(evals(0) - evals(1)) < 1e-12)
        {
            v1_evec = rxy.col(0);
            v2_evec = rxy.col(1);
            
            v1_evec.normalize();
            v2_evec = v2_evec - v1_evec.dot(v2_evec)*v1_evec;
            v2_evec.normalize();
        }
        
        // finally we assign
        s1 = std::sqrt(evals(0));
        s2 = std::sqrt(evals(1));
        v1 = v1_evec;
        v2 = v2_evec;
        
        if(s2 > s1) // make sure that s1 is the biggest (does not matter really, just convention -- will still give the same matrix)
        {
            std::swap(s1,s2);
            v1.swap(v2);
        }
#ifndef NDEBUG
        if(not checkVectorValidity())
        {
            std::cout << "decomposition is not valid : the eigenvectors are not unit length and/or not orthogonal" << std::endl;
            print();
            printf("computed eigenvalues are %10.10e \t %10.10e\n", evals(0), evals(1));
            printf("computed evec1 is %10.10e \t %10.10e\n", evecs(0,0), evecs(1,0));
            printf("computed evec2 is %10.10e \t %10.10e\n", evecs(0,1), evecs(1,1));
            printf("input rxy_1 is %10.10e \t %10.10e \t %10.10e\n", rxy(0,0), rxy(1,0), rxy(2,0));
            printf("input rxy_2 is %10.10e \t %10.10e \t %10.10e\n", rxy(0,1), rxy(1,1), rxy(2,1));
            printf("input metric is %10.10e \t %10.10e \t %10.10e\n", metric(0,0), metric(0,1), metric(1,1));
            
            assert(checkVectorValidity());
        }
#endif
    }
    
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    
    DecomposedGrowthState(const Real s1_in, const Real s2_in, const Eigen::Vector3d & v1_in, const Eigen::Vector3d & v2_in):
    s1(s1_in),
    s2(s2_in),
    v1(v1_in),
    v2(v2_in)
    {}
    
    DecomposedGrowthState(const Eigen::MatrixXd & rxy, const Eigen::Matrix2d & metric)
    {
        decompose_metric(rxy, metric);
    }
    
    void changeMetric(const Eigen::MatrixXd & rxy, const Eigen::Matrix2d & metric)
    {
        decompose_metric(rxy, metric);
    }
    
    bool checkVectorValidity(const Real tol = 1e-9) const
    {
        const bool unit_v1 = std::abs(v1.dot(v1) - 1) < tol;
        const bool unit_v2 = std::abs(v2.dot(v2) - 1) < tol;
        const bool ortho = std::abs(v1.dot(v2)) < tol;
        return (unit_v1 and unit_v2 and ortho);
    }
    
    void print() const
    {
        printf("s1/s2 : %10.10e \t %10.10e \t\t v1 : %10.10e, %10.10e, %10.10e \t v2 : %10.10e, %10.10e, %10.10e\n",s1,s2,v1(0),v1(1),v1(2),v2(0),v2(1),v2(2));
    }
    
    Eigen::Matrix2d computeMetric(const Eigen::MatrixXd & rxy) const
    {
        const Eigen::Matrix2d a_i = rxy.transpose() * rxy;
        const Eigen::Matrix2d Lsq = (Eigen::Matrix2d() << s1*s1,0,0,s2*s2).finished();
        Eigen::MatrixXd v12(3,2);
        v12 << v1,v2;
        const Eigen::Matrix2d V = (v12.transpose() * rxy).inverse();
        const Eigen::Matrix2d a_f = a_i * V * Lsq * V.transpose() * a_i;
        return a_f;
    }
    
    Real get_s1() const {return s1;}
    Real get_s2() const {return s2;}
    Eigen::Vector3d get_v1() const {return v1;}
    Eigen::Vector3d get_v2() const {return v2;}
};

class GrowthState
{
protected:
    const Eigen::MatrixXd rxy_base;
    const Eigen::Matrix2d a_final;
    const DecomposedGrowthState decomposed_final;
    Eigen::Matrix2d a_init;
    DecomposedGrowthState decomposed_init;
    
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    
    GrowthState(const Eigen::MatrixXd & rxy_base_in, const Eigen::Matrix2d & a_final_in):
    rxy_base(rxy_base_in),
    a_final(a_final_in),
    decomposed_final(rxy_base, a_final),
    a_init(rxy_base.transpose() * rxy_base),
    decomposed_init(rxy_base, a_init)
    {
    }
    
    GrowthState(const Eigen::MatrixXd & rxy_base_in, const Eigen::Matrix2d & a_init_in, const Eigen::Matrix2d & a_final_in):
    rxy_base(rxy_base_in),
    a_final(a_final_in),
    decomposed_final(rxy_base, a_final),
    a_init(a_init_in),
    decomposed_init(rxy_base, a_init)
    {
    }
    
    
    void changeInitialState(const Eigen::Matrix2d & a_init_in)
    {
        a_init = a_init_in;
        decomposed_init.changeMetric(rxy_base, a_init);
    }
    
    Eigen::Matrix2d interpolate_from_iso(const Real t) const
    {
        // we start from isotropic : just interpolate the growth factors
        const Real s1_interp = (1.0 - t) * decomposed_init.get_s1() + t*decomposed_final.get_s1();
        const Real s2_interp = (1.0 - t) * decomposed_init.get_s2() + t*decomposed_final.get_s2();
        
        const DecomposedGrowthState decomposed_intermediate(s1_interp, s2_interp, decomposed_final.get_v1(), decomposed_final.get_v2());
        return decomposed_intermediate.computeMetric(rxy_base);
    }
    
    Eigen::Matrix2d interpolate_from_iso_logeucl(const Real t) const
    {
        // we start from isotropic : just interpolate the growth factors
        const Real s1_interp = std::exp( (1.0 - t) * std::log(decomposed_init.get_s1()) + t*std::log(decomposed_final.get_s1()) );
        const Real s2_interp = std::exp( (1.0 - t) * std::log(decomposed_init.get_s2()) + t*std::log(decomposed_final.get_s2()) );
        
        const DecomposedGrowthState decomposed_intermediate(s1_interp, s2_interp, decomposed_final.get_v1(), decomposed_final.get_v2());
        return decomposed_intermediate.computeMetric(rxy_base);
    }
    
    Eigen::Matrix2d interpolate(const Real t) const
    {
        if(std::abs(decomposed_init.get_s1() - decomposed_init.get_s2()) < 1e-12) return interpolate_from_iso(t);
        
        // t between 0 and 1
        return (1.0 - t)*a_init + t*a_final;
    }
    
    Eigen::Matrix2d interpolateLogEucl(const Real t) const
    {
        if(std::abs(decomposed_init.get_s1() - decomposed_init.get_s2()) < 1e-12) return interpolate_from_iso_logeucl(t);
        
        // t between 0 and 1
        return ((1.0 - t)*a_init.log() + t * a_final.log()).exp();
    }
    
    Eigen::Matrix2d getInitGrowthMetric() const
    {
        return a_init;
    }
    
    Eigen::Matrix2d getFinalGrowthMetric() const
    {
        return a_final;
    }
    
    const DecomposedGrowthState & getDecomposedInitState() const
    {
        return decomposed_init;
    }
    
    const DecomposedGrowthState & getDecomposedFinalState() const
    {
        return decomposed_final;
    }
};



template<typename tMesh>
struct GrowthHelper
{
    static void anglesToVectors(const Eigen::Ref<const Eigen::VectorXd> angles, Eigen::Ref<Eigen::MatrixXd> vectors, const Real phase = 0)
    {
        const int nFaces = angles.rows();
        for(int i=0;i<nFaces;++i)
        {
            const Real phi = angles(i);
            const Eigen::Vector3d dir = (Eigen::Vector3d() <<  std::cos(phi + phase), std::sin(phi + phase), 0).finished();
            for(int d=0;d<3;++d)
                vectors(i,d) = dir(d);
        }
    }
    
    static void anglesToVectors(const Eigen::Ref<const Eigen::VectorXd> angles, const Eigen::Ref<const Eigen::VectorXd> lengths, Eigen::Ref<Eigen::MatrixXd> vectors, const Real phase = 0)
    {
        const int nFaces = angles.rows();
        for(int i=0;i<nFaces;++i)
        {
            const Real phi = angles(i);
            const Eigen::Vector3d dir = (Eigen::Vector3d() <<  std::cos(phi + phase), std::sin(phi + phase), 0).finished();
            for(int d=0;d<3;++d)
                vectors(i,d) = lengths(i)*dir(d);
        }
    }
    
    static void vectorsToAngles(const Eigen::Ref<const Eigen::MatrixXd> dir1, const Eigen::Ref<const Eigen::MatrixXd> dir2, Eigen::Ref<Eigen::VectorXd> angles, Eigen::Ref<Eigen::VectorXd> growthRates1, Eigen::Ref<Eigen::VectorXd> growthRates2)
    {
        const int nFaces = dir1.rows();
        for(int i=0;i<nFaces;++i)
        {
            const Eigen::Vector3d & dir1_i = dir1.row(i);
            const Eigen::Vector3d & dir2_i = dir2.row(i);
            
            const Real mag_dir1_i = dir1_i.norm();
            const Real mag_dir2_i = dir2_i.norm();

            const Real angle_i = std::atan2(dir1_i(1), dir1_i(0));
            
            // fill the arrays
            angles(i) = angle_i;
            growthRates1(i) = mag_dir1_i;
            growthRates2(i) = mag_dir2_i;

        }
    }
    static void computeAbarsIsoGrowth(const tMesh & mesh, const Eigen::Ref<const Eigen::VectorXd> growthRates, tVecMat2d & aforms)
    {
        const int nFaces = mesh.getNumberOfFaces();
        const auto & topo = mesh.getTopology();
        const auto & reststate = mesh.getRestConfiguration(); // does not matter which one at this point
        
        aforms.resize(nFaces);
        
        for(int i=0;i<nFaces;++i)
        {
            const TriangleInfo info = reststate.getTriangleInfoLite(topo, i);
            const Real growth = std::pow(1 + growthRates(i),2);
            
            const Eigen::Vector3d dir_p = (Eigen::Vector3d() << 1, 0, 0).finished();
            const Eigen::Vector3d dir_o = (Eigen::Vector3d() << 0, 1, 0).finished();
            
            const Real e1_p = (info.e1).dot(dir_p);
            const Real e1_o = (info.e1).dot(dir_o);
            
            const Real e2_p = (info.e2).dot(dir_p);
            const Real e2_o = (info.e2).dot(dir_o);
            
            const Real a11 = (e1_p*e1_p + e1_o*e1_o)*growth;
            const Real a12 = (e1_p*e2_p + e1_o*e2_o)*growth;
            const Real a22 = (e2_p*e2_p + e2_o*e2_o)*growth;
            
            aforms[i](0,0) = a11;
            aforms[i](0,1) = aforms[i](1,0) = a12;
            aforms[i](1,1) = a22;
        }
    }
    
    static void computeAbarsIsoGrowth_Gradient(const tMesh & mesh, const Eigen::Ref<const Eigen::VectorXd> growthRates, const Eigen::Ref<const Eigen::MatrixXd> gradEng_abar, Eigen::Ref<Eigen::VectorXd> gradEng_rate)
    {
        // compute gradient abar wrt growth
        const int nFaces = mesh.getNumberOfFaces();
        const auto & topo = mesh.getTopology();
        const auto & reststate = mesh.getRestConfiguration(); // does not matter which one at this point
        
        for(int i=0;i<nFaces;++i)
        {
            const TriangleInfo info = reststate.getTriangleInfoLite(topo, i);
            const Real dgrowth = 2.0*(1 + growthRates(i));
            
            const Eigen::Vector3d dir_p = (Eigen::Vector3d() << 1, 0, 0).finished();
            const Eigen::Vector3d dir_o = (Eigen::Vector3d() << 0, 1, 0).finished();
            
            const Real e1_p = (info.e1).dot(dir_p);
            const Real e1_o = (info.e1).dot(dir_o);
            
            const Real e2_p = (info.e2).dot(dir_p);
            const Real e2_o = (info.e2).dot(dir_o);
            
            const Real da11 = (e1_p*e1_p + e1_o*e1_o)*dgrowth;
            const Real da12 = (e1_p*e2_p + e1_o*e2_o)*dgrowth;
            const Real da22 = (e2_p*e2_p + e2_o*e2_o)*dgrowth;
            
            gradEng_rate(i) = gradEng_abar(i,0)*da11 + gradEng_abar(i,1)*da12 + gradEng_abar(i,2)*da22;
        }
    }
    
    static void computeAbarsOrthoGrowth(const tMesh & mesh, const Eigen::Ref<const Eigen::VectorXd> growthAngles, const Eigen::Ref<const Eigen::VectorXd> growthRates_p, const Eigen::Ref<const Eigen::VectorXd> growthRates_o, tVecMat2d & aforms)
    {
        const int nFaces = mesh.getNumberOfFaces();
        const auto & topo = mesh.getTopology();
        const auto & reststate = mesh.getRestConfiguration(); // does not matter which one at this point
        
        aforms.resize(nFaces);
        
        for(int i=0;i<nFaces;++i)
        {
            const TriangleInfo info = reststate.getTriangleInfoLite(topo, i);
            
            const Real growth_p = std::pow(1 + growthRates_p(i),2);
            const Real growth_o = std::pow(1 + growthRates_o(i),2);
            const Real phi = growthAngles(i);
            
            const Eigen::Vector3d dir_p = (Eigen::Vector3d() <<  std::cos(phi), std::sin(phi), 0).finished();
            const Eigen::Vector3d dir_o = (Eigen::Vector3d() << -std::sin(phi), std::cos(phi), 0).finished();

//          new edges are the following
//          const Eigen::Vector3d e0_new = (1.0 + growthRate_p)*(e0.dot(dir_p))*dir_p + (1.0 + growthRate_o)*(e0.dot(dir_o))*dir_o;
//          const Eigen::Vector3d e1_new = (1.0 + growthRate_p)*(e1.dot(dir_p))*dir_p + (1.0 + growthRate_o)*(e1.dot(dir_o))*dir_o;
//          const Eigen::Vector3d e2_new = (1.0 + growthRate_p)*(e2.dot(dir_p))*dir_p + (1.0 + growthRate_o)*(e2.dot(dir_o))*dir_o;
//          a11 = e1_new.dot(e1_new), a12 = e1_new.dot(e2_new), a22 = e2_new.dot(e2_new)
//          also : e0_new + e1_new + e2_new = 0
            
            const Real e1_p = (info.e1).dot(dir_p);
            const Real e1_o = (info.e1).dot(dir_o);
            
            const Real e2_p = (info.e2).dot(dir_p);
            const Real e2_o = (info.e2).dot(dir_o);
            
            const Real a11 = e1_p*e1_p*growth_p + e1_o*e1_o*growth_o;
            const Real a12 = e1_p*e2_p*growth_p + e1_o*e2_o*growth_o;
            const Real a22 = e2_p*e2_p*growth_p + e2_o*e2_o*growth_o;
            
            aforms[i](0,0) = a11;
            aforms[i](0,1) = aforms[i](1,0) = a12;
            aforms[i](1,1) = a22;
        }
    }
    
    static void computeAbarsOrthoGrowth_Gradient(const tMesh & mesh, const Eigen::Ref<const Eigen::VectorXd> growthAngles, const Eigen::Ref<const Eigen::VectorXd> growthRates_p, const Eigen::Ref<const Eigen::VectorXd> growthRates_o, const Eigen::Ref<const Eigen::MatrixXd> gradEng_abar, Eigen::Ref<Eigen::VectorXd> gradEng_angle, Eigen::Ref<Eigen::VectorXd> gradEng_rate_p, Eigen::Ref<Eigen::VectorXd> gradEng_rate_o)
    {
        // compute gradient abar wrt growth
        const int nFaces = mesh.getNumberOfFaces();
        const auto & topo = mesh.getTopology();
        const auto & reststate = mesh.getRestConfiguration(); // does not matter which one at this point
        
        for(int i=0;i<nFaces;++i)
        {
            const TriangleInfo info = reststate.getTriangleInfoLite(topo, i);
            
            const Real growth_p = std::pow(1 + growthRates_p(i),2);
            const Real growth_o = std::pow(1 + growthRates_o(i),2);
            const Real phi = growthAngles(i);

            const Eigen::Vector3d dir_p = (Eigen::Vector3d() <<  std::cos(phi), std::sin(phi), 0).finished();
            const Eigen::Vector3d dir_o = (Eigen::Vector3d() << -std::sin(phi), std::cos(phi), 0).finished();

            // derivatives
            const Real dgrowth_p = 2.0*(1 + growthRates_p(i));
            const Real dgrowth_o = 2.0*(1 + growthRates_o(i));
            
            const Eigen::Vector3d ddir_p = (Eigen::Vector3d() << -std::sin(phi),  std::cos(phi), 0).finished();
            const Eigen::Vector3d ddir_o = (Eigen::Vector3d() << -std::cos(phi), -std::sin(phi), 0).finished();
            
            
            const Real e1_p = (info.e1).dot(dir_p);
            const Real e1_o = (info.e1).dot(dir_o);
            
            const Real e2_p = (info.e2).dot(dir_p);
            const Real e2_o = (info.e2).dot(dir_o);
            
            const Real e1_dp = (info.e1).dot(ddir_p);
            const Real e1_do = (info.e1).dot(ddir_o);
            
            const Real e2_dp = (info.e2).dot(ddir_p);
            const Real e2_do = (info.e2).dot(ddir_o);

            const Real da11_dphi = 2.0*e1_p*e1_dp*growth_p + 2.0*e1_o*e1_do*growth_o;
            const Real da12_dphi = (e1_p*e2_dp + e1_dp*e2_p)*growth_p + (e1_o*e2_do + e1_do*e2_o)*growth_o;
            const Real da22_dphi = 2.0*e2_p*e2_dp*growth_p + 2.0*e2_o*e2_do*growth_o;
            
            const Real da11_dgrowth_p = e1_p*e1_p*dgrowth_p;
            const Real da12_dgrowth_p = e1_p*e2_p*dgrowth_p;
            const Real da22_dgrowth_p = e2_p*e2_p*dgrowth_p;
            
            const Real da11_dgrowth_o = e1_o*e1_o*dgrowth_o;
            const Real da12_dgrowth_o = e1_o*e2_o*dgrowth_o;
            const Real da22_dgrowth_o = e2_o*e2_o*dgrowth_o;
            
            gradEng_angle(i)  = gradEng_abar(i,0)*da11_dphi      + gradEng_abar(i,1)*da12_dphi      + gradEng_abar(i,2)*da22_dphi;
            gradEng_rate_p(i) = gradEng_abar(i,0)*da11_dgrowth_p + gradEng_abar(i,1)*da12_dgrowth_p + gradEng_abar(i,2)*da22_dgrowth_p;
            gradEng_rate_o(i) = gradEng_abar(i,0)*da11_dgrowth_o + gradEng_abar(i,1)*da12_dgrowth_o + gradEng_abar(i,2)*da22_dgrowth_o;
        }
    }
    
    
    static void computeAbarsOrthoGrowth_GradientVertices(const tMesh & mesh, const Eigen::Ref<const Eigen::VectorXd> growthAngles, const Eigen::Ref<const Eigen::VectorXd> growthRates_p, const Eigen::Ref<const Eigen::VectorXd> growthRates_o, const Eigen::Ref<const Eigen::MatrixXd> gradEng_abar, Eigen::Ref<Eigen::MatrixXd> gradEng_verts)
    {
        // compute gradient abar wrt growth
        const int nFaces = mesh.getNumberOfFaces();
        const auto & topo = mesh.getTopology();
        const auto & reststate = mesh.getRestConfiguration(); // does not matter which one at this point
        
        for(int i=0;i<nFaces;++i)
        {
            const TriangleInfo info = reststate.getTriangleInfoLite(topo, i);
            
            const Real growth_p = std::pow(1 + growthRates_p(i),2);
            const Real growth_o = std::pow(1 + growthRates_o(i),2);
            const Real phi = growthAngles(i);
            
            const Eigen::Vector3d dir_p = (Eigen::Vector3d() <<  std::cos(phi), std::sin(phi), 0).finished();
            const Eigen::Vector3d dir_o = (Eigen::Vector3d() << -std::sin(phi), std::cos(phi), 0).finished();
            
            // e1 = v2 - v1
            // e2 = v0 - v2
            const Real e1_p = (info.e1).dot(dir_p);
            const Real e1_o = (info.e1).dot(dir_o);
            
            const Real e2_p = (info.e2).dot(dir_p);
            const Real e2_o = (info.e2).dot(dir_o);
            
//            const Real a11 = e1_p*e1_p*growth_p + e1_o*e1_o*growth_o;
//            const Real a12 = e1_p*e2_p*growth_p + e1_o*e2_o*growth_o;
//            const Real a22 = e2_p*e2_p*growth_p + e2_o*e2_o*growth_o;
            
            const Eigen::Vector3d gradv0_a12 = e1_p*growth_p*dir_p + e1_o*growth_o*dir_o;
            const Eigen::Vector3d gradv0_a22 = 2.0*(e2_p*growth_p*dir_p + e2_o*growth_o*dir_o);
            
            const Eigen::Vector3d gradv1_a11 = -2.0*(e1_p*growth_p*dir_p + e1_o*growth_o*dir_o);
            const Eigen::Vector3d gradv1_a12 = -e2_p*growth_p*dir_p - e2_o*growth_o*dir_o;
            
            const Eigen::Vector3d gradv2_a11 = 2.0*(e1_p*growth_p*dir_p + e1_o*growth_o*dir_o);
            const Eigen::Vector3d gradv2_a12 = growth_p*(e2_p - e1_p)*dir_p + growth_o*(e2_o - e1_o)*dir_o;
            const Eigen::Vector3d gradv2_a22 = -2.0*(e2_p*growth_p*dir_p + e2_o*growth_o*dir_o);
            
            for(int d=0;d<2;++d)
            {
                gradEng_verts(info.idx_v0,d) += gradEng_abar(i,1)*gradv0_a12(d) + gradEng_abar(i,2)*gradv0_a22(d);
                gradEng_verts(info.idx_v1,d) += gradEng_abar(i,0)*gradv1_a11(d) + gradEng_abar(i,1)*gradv1_a12(d);
                gradEng_verts(info.idx_v2,d) += gradEng_abar(i,0)*gradv2_a11(d) + gradEng_abar(i,1)*gradv2_a12(d) + gradEng_abar(i,2)*gradv2_a22(d);
            }
            
        }
    }

    
    
    static void computeAbarsOrthoGrowthShell(const tMesh & mesh, const Eigen::Ref<const Eigen::MatrixXd> growthdirs_1, const Eigen::Ref<const Eigen::VectorXd> growthRates_1, const Eigen::Ref<const Eigen::VectorXd> growthRates_2, tVecMat2d & aforms)
    {
        const int nFaces = mesh.getNumberOfFaces();
        const auto & topo = mesh.getTopology();
        const auto & reststate = mesh.getRestConfiguration(); // does not matter which one at this point
        
        aforms.resize(nFaces);
        
        for(int i=0;i<nFaces;++i)
        {
            const TriangleInfo info = reststate.getTriangleInfoLite(topo, i);
            
            const Real growth_1 = std::pow(1 + growthRates_1(i),2);
            const Real growth_2 = std::pow(1 + growthRates_2(i),2);
            const Eigen::Vector3d dir_1 = growthdirs_1.row(i).normalized(); // should be normalized already but do it again just to be sure
            
            // get the face normal to compute dir_2
            const Eigen::Vector3d facenormal = (info.e2).cross(info.e0).normalized();
            const Eigen::Vector3d dir_2 = (facenormal.cross(dir_1)).normalized(); // so that dir_1 cross dir_2 = facenormal
            
            const Real e1_1 = (info.e1).dot(dir_1);
            const Real e1_2 = (info.e1).dot(dir_2);
            
            const Real e2_1 = (info.e2).dot(dir_1);
            const Real e2_2 = (info.e2).dot(dir_2);
            
            const Real a11 = e1_1*e1_1*growth_1 + e1_2*e1_2*growth_2;
            const Real a12 = e1_1*e2_1*growth_1 + e1_2*e2_2*growth_2;
            const Real a22 = e2_1*e2_1*growth_1 + e2_2*e2_2*growth_2;
            
            aforms[i](0,0) = a11;
            aforms[i](0,1) = aforms[i](1,0) = a12;
            aforms[i](1,1) = a22;
        }
    }
    
    
    
    
    
    
    
    
    static void computeAbarsIsoGrowthPerEdge(const tMesh & mesh, const Eigen::Ref<const Eigen::VectorXd> growthRates, tVecMat2d & aforms)
    {
        // now the growth rates are prescribed per-edge
        const int nFaces = mesh.getNumberOfFaces();
        const auto & topo = mesh.getTopology();
        const auto & reststate = mesh.getRestConfiguration(); // does not matter which one at this point
        
        aforms.resize(nFaces);
        
        for(int i=0;i<nFaces;++i)
        {
            const TriangleInfo info = reststate.getTriangleInfoLite(topo, i);
            const Real growth_e0 = std::pow(1 + growthRates(info.idx_e0), 2);
            const Real growth_e1 = std::pow(1 + growthRates(info.idx_e1), 2);
            const Real growth_e2 = std::pow(1 + growthRates(info.idx_e2), 2);
            
            const Real e0_magsq = (info.e0).dot(info.e0);
            const Real e1_magsq = (info.e1).dot(info.e1);
            const Real e2_magsq = (info.e2).dot(info.e2);
            
            const Real a11 = growth_e1 * e1_magsq;
            const Real a12 = -0.5*( growth_e1 * e1_magsq + growth_e2 * e2_magsq - growth_e0 * e0_magsq ); // checked : this is correct
            const Real a22 = growth_e2 * e2_magsq;
            
            aforms[i](0,0) = a11;
            aforms[i](0,1) = aforms[i](1,0) = a12;
            aforms[i](1,1) = a22;
        }
    }
    
    static void computeAbarsIsoGrowthPerEdge_Gradient(const tMesh & mesh, const Eigen::Ref<const Eigen::VectorXd> growthRates, const Eigen::Ref<const Eigen::MatrixXd> gradEng_abar, Eigen::Ref<Eigen::VectorXd> gradEng_rate)
    {
        // compute gradient abar wrt growth
        const int nFaces = mesh.getNumberOfFaces();
        const auto & topo = mesh.getTopology();
        const auto & reststate = mesh.getRestConfiguration(); // does not matter which one at this point
        
        for(int i=0;i<nFaces;++i)
        {
            const TriangleInfo info = reststate.getTriangleInfoLite(topo, i);
            
            const Real dgrowth_e0 = 2 * (1 + growthRates(info.idx_e0));
            const Real dgrowth_e1 = 2 * (1 + growthRates(info.idx_e1));
            const Real dgrowth_e2 = 2 * (1 + growthRates(info.idx_e2));
            
            const Real e0_magsq = (info.e0).dot(info.e0);
            const Real e1_magsq = (info.e1).dot(info.e1);
            const Real e2_magsq = (info.e2).dot(info.e2);
            
            const Real da11_e1 = dgrowth_e1 * e1_magsq;

            const Real da12_e0 =  0.5 * dgrowth_e0 * e0_magsq;
            const Real da12_e1 = -0.5 * dgrowth_e1 * e1_magsq;
            const Real da12_e2 = -0.5 * dgrowth_e2 * e2_magsq;

            const Real da22_e2 = dgrowth_e2 * e2_magsq;

            gradEng_rate(info.idx_e0) += gradEng_abar(i,1)*da12_e0;
            gradEng_rate(info.idx_e1) += gradEng_abar(i,0)*da11_e1 + gradEng_abar(i,1)*da12_e1;
            gradEng_rate(info.idx_e2) += gradEng_abar(i,1)*da12_e2 + gradEng_abar(i,2)*da22_e2;
        }
    }
    
    
    
    
    
    static void computeAbarsIsoGrowthPerVertex(const tMesh & mesh, const Eigen::Ref<const Eigen::VectorXd> growthRates, tVecMat2d & aforms)
    {
        // now the growth rates are prescribed per-edge
        const int nFaces = mesh.getNumberOfFaces();
        const auto & topo = mesh.getTopology();
        const auto & reststate = mesh.getRestConfiguration(); // does not matter which one at this point
        
        aforms.resize(nFaces);
        
        for(int i=0;i<nFaces;++i)
        {
            const TriangleInfo info = reststate.getTriangleInfoLite(topo, i);
            const Real growth_v0 = 1 + growthRates(info.idx_v0);
            const Real growth_v1 = 1 + growthRates(info.idx_v1);
            const Real growth_v2 = 1 + growthRates(info.idx_v2);
            
            const Real growth_e1 = 0.5*(growth_v2 + growth_v1);
            const Real growth_e2 = 0.5*(growth_v0 + growth_v2);
            
            const Eigen::Vector3d dir_p = (Eigen::Vector3d() << 1, 0, 0).finished();
            const Eigen::Vector3d dir_o = (Eigen::Vector3d() << 0, 1, 0).finished();
            
            const Real e1_p = (info.e1).dot(dir_p);
            const Real e1_o = (info.e1).dot(dir_o);
            
            const Real e2_p = (info.e2).dot(dir_p);
            const Real e2_o = (info.e2).dot(dir_o);
            
            const Real e1_p_growth = e1_p * growth_e1;
            const Real e1_o_growth = e1_o * growth_e1;
            const Real e2_p_growth = e2_p * growth_e2;
            const Real e2_o_growth = e2_o * growth_e2;
            
            const Real a11 = e1_p_growth*e1_p_growth + e1_o_growth*e1_o_growth;
            const Real a12 = e1_p_growth*e2_p_growth + e1_o_growth*e2_o_growth;
            const Real a22 = e2_p_growth*e2_p_growth + e2_o_growth*e2_o_growth;
            
            aforms[i](0,0) = a11;
            aforms[i](0,1) = aforms[i](1,0) = a12;
            aforms[i](1,1) = a22;
        }
    }
    
    static void computeAbarsIsoGrowthPerVertex_Gradient(const tMesh & mesh, const Eigen::Ref<const Eigen::VectorXd> growthRates, const Eigen::Ref<const Eigen::MatrixXd> gradEng_abar, Eigen::Ref<Eigen::VectorXd> gradEng_rate)
    {
        // compute gradient abar wrt growth
        const int nFaces = mesh.getNumberOfFaces();
        const auto & topo = mesh.getTopology();
        const auto & reststate = mesh.getRestConfiguration(); // does not matter which one at this point
        
        for(int i=0;i<nFaces;++i)
        {
            const TriangleInfo info = reststate.getTriangleInfoLite(topo, i);
            const Real growth_v0 = 1 + growthRates(info.idx_v0);
            const Real growth_v1 = 1 + growthRates(info.idx_v1);
            const Real growth_v2 = 1 + growthRates(info.idx_v2);
            
            const Real growth_e1 = 0.5*(growth_v2 + growth_v1);
            const Real growth_e2 = 0.5*(growth_v0 + growth_v2);
            
            const Eigen::Vector3d dir_p = (Eigen::Vector3d() << 1, 0, 0).finished();
            const Eigen::Vector3d dir_o = (Eigen::Vector3d() << 0, 1, 0).finished();
            
            const Real e1_p = (info.e1).dot(dir_p);
            const Real e1_o = (info.e1).dot(dir_o);
            
            const Real e2_p = (info.e2).dot(dir_p);
            const Real e2_o = (info.e2).dot(dir_o);
            
            const Real e1_p_growth = e1_p * growth_e1;
            const Real e1_o_growth = e1_o * growth_e1;
            const Real e2_p_growth = e2_p * growth_e2;
            const Real e2_o_growth = e2_o * growth_e2;
            
            const Real da11_e1 = 2.0*(e1_p_growth * e1_p + e1_o_growth * e1_o);
            const Real da12_e1 = e1_p * e2_p_growth + e1_o * e2_o_growth;
            const Real da12_e2 = e1_p_growth * e2_p + e1_o_growth * e2_o;
            const Real da22_e2 = 2.0*(e2_p_growth * e2_p + e2_o_growth * e2_o);
            
            const Real da12_v0 = 0.5*da12_e2;
            const Real da22_v0 = 0.5*da22_e2;

            const Real da11_v1 = 0.5*da11_e1;
            const Real da12_v1 = 0.5*da12_e1;
            
            const Real da11_v2 = 0.5*da11_e1;
            const Real da12_v2 = 0.5*(da12_e1 + da12_e2);
            const Real da22_v2 = 0.5*da22_e2;
            
            gradEng_rate(info.idx_v0) += gradEng_abar(i,1)*da12_v0 + gradEng_abar(i,2)*da22_v0;
            gradEng_rate(info.idx_v1) += gradEng_abar(i,0)*da11_v1 + gradEng_abar(i,1)*da12_v1;
            gradEng_rate(info.idx_v2) += gradEng_abar(i,0)*da11_v2 + gradEng_abar(i,1)*da12_v2 + gradEng_abar(i,2)*da22_v2;
        }
    }
};

#endif /* GrowthHelper_hpp */
