//
//  Test_EnergyOperatorDCS_Gradient.cpp
//  Elasticity
//
//  Created by Wim van Rees on 8/20/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#include "gtest/gtest.h"

#include "common.hpp"
#include "Timer.hpp"
#include "Geometry.hpp"
#include "Mesh.hpp"
#include "MaterialProperties.hpp"

#include "ComputeBasicGradients.hpp"
#include "Derivative_FD.hpp"

#include "CombinedOperator_Parametric.hpp"

template <typename tM, typename tMat, template<typename, typename, MeshLayer> class tS, MeshLayer tL>
struct TypeDefinitions
{
    typedef tM tMesh;
    typedef tMat tMaterialProperties;
    typedef typename tMat::tMaterialType tMaterialType;
    typedef tS<tM, tMaterialType, tL> tEnergyOperatorDCS;
    static const MeshLayer layer = tL;
};

template <typename tCombinedTypes>
class Test_EnergyOperatorDCS_Gradient : public ::testing::Test
{
public:
    typedef typename tCombinedTypes::tMesh tMesh;
    typedef typename tCombinedTypes::tMaterialProperties tMaterialProperties;
    typedef typename tCombinedTypes::tEnergyOperatorDCS tEnergyOperatorDCS;
    
protected:
    tMesh mesh;
    
    // create a mesh
    virtual void initGrid()
    {
        // use a shell for non-zero bbar
        const Real edgeLength = 0.25;
        SphericalShell_Lambert geometry(1.0, 1.0, 0.0, edgeLength, true);
        geometry.setQuiet();
        mesh.init(geometry);
    }
    
    void applyRandomDeformation(const Real deltaV, const Real deltaE)
    {
        std::mt19937 gen;
        gen.seed(42);
        std::uniform_real_distribution<Real> distV(-deltaV, deltaV);
        std::uniform_real_distribution<Real> distE(-deltaE, deltaE);
        auto perturb_v = [&](Eigen::Vector3d in)
        {
            Eigen::Vector3d retval;
            retval << in(0)+distV(gen), in(1)+distV(gen), 0;//in(2)+dist(gen);
            return retval;
        };
        auto perturb_e = [&](Real in)
        {
            return in + distE(gen);
        };
        
        mesh.changeVertices(perturb_v);
        mesh.changeEdgeDirectors(perturb_e);
    }
    
    std::pair<bool, bool> testGradient(const tEnergyOperatorDCS & meshoperator)
    {
        const int nVertices = mesh.getNumberOfVertices();
        const int nEdges = mesh.getNumberOfEdges();
        
        Eigen::MatrixXd gradient_vertices_exact(nVertices,3);
        Eigen::VectorXd gradient_edges_exact(nEdges);
        gradient_vertices_exact.setZero();
        gradient_edges_exact.setZero();
        
        Timer<std::chrono::milliseconds> timer;
        timer.start();
        
        mesh.updateDeformedConfiguration();
        meshoperator.compute(mesh, gradient_vertices_exact, gradient_edges_exact);
        timer.stop();
       
        // compute the approximate gradient
        Eigen::MatrixXd gradient_vertices_approx(nVertices,3);
        Eigen::VectorXd gradient_edges_approx(nEdges);
        
        const Real eps = 1e-4;
        const int order=3; // 4th order
        
        Derivative_FD<tMesh, tEnergyOperatorDCS> derivative_FD(order, mesh, meshoperator, eps);
        
        // compute gradient using finite differences
        timer.reset();
        timer.start();

        derivative_FD.computeGradient(gradient_vertices_approx, gradient_edges_approx);
        
        timer.stop();
        
//        std::cout << "Approx gradient (order=" << order+1 << ") took " << timer.totTime << " ms" << std::endl;

        const Real L1Error_V =(gradient_vertices_exact - gradient_vertices_approx).lpNorm<1>()/gradient_vertices_approx.rows();
        const Real maxError_V = (gradient_vertices_exact - gradient_vertices_approx).lpNorm<Eigen::Infinity>();
//        std::cout << "=== avg/max gradient error vertices = " << L1Error_V << " " << maxError_V << std::endl;
        
        const Real L1Error_E =(gradient_edges_exact - gradient_edges_approx).lpNorm<1>()/gradient_edges_approx.rows();
        const Real maxError_E = (gradient_edges_exact - gradient_edges_approx).lpNorm<Eigen::Infinity>();
//        std::cout << "=== avg/max gradient error edges = " << L1Error_E << " " << maxError_E << std::endl;
  
        _unused(L1Error_V);
        _unused(L1Error_E);
        
        const bool verticesOK = maxError_V < 1e-10;
        const bool edgesOK = maxError_E < 1e-10;

        return std::make_pair(verticesOK, edgesOK);
//        if(maxError_V > 1e-6)
//        {
//            std::cout << "V EXACT\n" << gradient_vertices_exact << std::endl;
//            std::cout << "V APPROX\n" << gradient_vertices_approx << std::endl;
//        }
//        
//        if(maxError_E > 1e-6)
//        {
//            std::cout << "E EXACT\n" << gradient_edges_exact << std::endl;
//            std::cout << "E APPROX\n" << gradient_edges_approx << std::endl;
//        }
        
    }
    
    virtual tMaterialProperties getMaterialProperties() = 0;
    
public:
    Test_EnergyOperatorDCS_Gradient()
    {
    }
    
    bool testGradient()
    {
        initGrid();
        applyRandomDeformation(0.005, 0.005);
        
        tMaterialProperties matprop = getMaterialProperties();
        tEnergyOperatorDCS engOp(matprop);
        std::pair<bool, bool> errors = testGradient(engOp);
        
        return (errors.first==true and errors.second==true);
    }
};

template <typename tCombinedTypes>
class Test_EnergyOperatorDCS_Gradient_Iso : public Test_EnergyOperatorDCS_Gradient<tCombinedTypes>
{
public:
    typedef typename tCombinedTypes::tMaterialProperties tMaterialProperties;
protected:
    tMaterialProperties getMaterialProperties() override
    {
        const Real E = 1.0;
        const Real nu = 0.3;
        const Real h = 0.01;
        tMaterialProperties matprop(E, nu, h);
        return matprop;
    }
};

template <typename tCombinedTypes>
class Test_EnergyOperatorDCS_Gradient_Ortho : public Test_EnergyOperatorDCS_Gradient<tCombinedTypes>
{
public:
    typedef typename tCombinedTypes::tMaterialProperties tMaterialProperties;
protected:
    tMaterialProperties getMaterialProperties() override
    {
        const Real E1 = 1.0;
        const Real E2 = 2.0;
        const Real nu_1 = 0.3;
        const Real G = 1.5;
        const Real h = 0.01;
        tMaterialProperties matprop(E1, E2, G, nu_1, h);
        
        // set aform_bar_d1
        {
            const int nFaces = this->mesh.getNumberOfFaces();
            const auto & restState = this->mesh.getRestConfiguration();
            for(int i=0;i<nFaces;++i)
            {
                const TriangleInfo info = restState.getTriangleInfoLite(this->mesh.getTopology(), i);
                // compute face normal
                const Eigen::Vector3d restface_normal = ((info.e2).cross(info.e0)).normalized();
                
                // compute orthogonal directions for this face
                const Eigen::Vector3d facectr = info.computeFaceCenter();
                const Eigen::Vector3d orthodir_d1 = (Eigen::Vector3d() << -facectr(1), facectr(0), 0).finished(); // azimuthal
                const Eigen::Vector3d a1 = (orthodir_d1 - (orthodir_d1.dot(restface_normal))*restface_normal).normalized(); // a1
//                const Eigen::Vector3d a2 = (restface_normal.cross(a1)).normalized(); // a2 --> no need to do this 
                
                // create the aform_bar
                // e1.dot(e1) = (e1.dot(a1))^2 + (e1.dot(a2))^2
                const Real a11_d1 = std::pow((info.e1).dot(a1), 2);
                const Real a12_d1 = (info.e1).dot(a1) * (info.e2).dot(a1);
                const Real a22_d1 = std::pow((info.e2).dot(a1), 2);
                
                matprop.material.aform_bar_d1 << a11_d1, a12_d1, a12_d1, a22_d1;
            }
        }
        return matprop;
    }
};



// The list of types we want to test (DCS).
typedef ::testing::Types<
TypeDefinitions<Mesh, MaterialProperties_Iso_Constant, CombinedOperator_Parametric, single>,
TypeDefinitions<BilayerMesh, MaterialProperties_Iso_Constant, CombinedOperator_Parametric, top>,
TypeDefinitions<BilayerMesh, MaterialProperties_Iso_Constant, CombinedOperator_Parametric, bottom>
> Test_EnergyOperatorDCS_GradientImplementations_Iso;

typedef ::testing::Types<
TypeDefinitions<Mesh, MaterialProperties_Ortho_Constant, CombinedOperator_Parametric, single>,
TypeDefinitions<BilayerMesh, MaterialProperties_Ortho_Constant, CombinedOperator_Parametric, top>,
TypeDefinitions<BilayerMesh, MaterialProperties_Ortho_Constant, CombinedOperator_Parametric, bottom>
> Test_EnergyOperatorDCS_GradientImplementations_Ortho;

// define the testcase (DCS)
TYPED_TEST_CASE(Test_EnergyOperatorDCS_Gradient_Iso, Test_EnergyOperatorDCS_GradientImplementations_Iso);
TYPED_TEST_CASE(Test_EnergyOperatorDCS_Gradient_Ortho, Test_EnergyOperatorDCS_GradientImplementations_Ortho);

TYPED_TEST(Test_EnergyOperatorDCS_Gradient_Iso, IsGradientCorrectDCS_Iso) {
    EXPECT_TRUE(this->testGradient());
}

TYPED_TEST(Test_EnergyOperatorDCS_Gradient_Ortho, IsGradientCorrectDCS_Ortho) {
    EXPECT_TRUE(this->testGradient());
}
