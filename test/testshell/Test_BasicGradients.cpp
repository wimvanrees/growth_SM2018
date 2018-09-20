//
//  Test_BasicGradients.cpp
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

#include "ComputeBasicGradients.hpp"
#include "Derivative_FD.hpp"


template <typename tM, int tIdx, template<typename, int> class tS>
struct TypeDefinitions
{
    typedef tM tMesh;
    typedef tS<tM, tIdx> tTestGradient;
};

template <typename tCombinedTypes>
class Test_BasicGradients : public ::testing::Test
{
public:
    typedef typename tCombinedTypes::tMesh tMesh;
    typedef typename tCombinedTypes::tTestGradient tTestGradient;
    
protected:
    tMesh mesh;
    
    // create a mesh
    virtual void initGrid()
    {
        // use a shell for non-zero bbar
        const Real edgeLength = 0.25;
        SphericalShell_Lambert geometry(1.0, 1.0, 0.0, edgeLength, false); // no boundary conditions here
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
            retval << in(0)+distV(gen), in(1)+distV(gen), in(2)+distV(gen);
            return retval;
        };
        auto perturb_e = [&](Real in)
        {
            return in + distE(gen);
        };
        
        mesh.changeVertices(perturb_v);
        mesh.changeEdgeDirectors(perturb_e);
    }
    
    std::pair<bool, bool> testGradient()
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
        tTestGradient gradOp;
        const Real eng = gradOp.compute(mesh, gradient_vertices_exact, gradient_edges_exact);
        timer.stop();
        std::cout << "[----------] === energy = " << eng << std::endl;
        
        // compute the approximate gradient
        Eigen::MatrixXd gradient_vertices_approx(nVertices,3);
        Eigen::VectorXd gradient_edges_approx(nEdges);

        const Real eps = 1e-4;
        const int order=3; // 4th order
        
        Derivative_FD<tMesh, tTestGradient> derivative_FD(order, mesh, gradOp, eps);
        
        // compute gradient using finite differences
        timer.reset();
        timer.start();
        
        derivative_FD.computeGradient(gradient_vertices_approx, gradient_edges_approx);
        
        timer.stop();
        
//        std::cout << "Approx gradient (order=" << order+1 << ") took " << timer.totTime << " ms" << std::endl;
        
        const Real L1Error_V =(gradient_vertices_exact - gradient_vertices_approx).lpNorm<1>()/gradient_vertices_approx.rows();
        const Real maxError_V = (gradient_vertices_exact - gradient_vertices_approx).lpNorm<Eigen::Infinity>();
        std::cout << "[----------] === avg/max gradient error vertices = " << L1Error_V << " " << maxError_V << std::endl;
        
        const Real L1Error_E =(gradient_edges_exact - gradient_edges_approx).lpNorm<1>()/gradient_edges_approx.rows();
        const Real maxError_E = (gradient_edges_exact - gradient_edges_approx).lpNorm<Eigen::Infinity>();
        std::cout << "[----------] === avg/max gradient error edges = " << L1Error_E << " " << maxError_E << std::endl;
  
        _unused(L1Error_V);
        _unused(L1Error_E);
        
        const bool verticesOK = maxError_V < 1e-9;
        const bool edgesOK = maxError_E < 1e-9;

        if(maxError_V > 1e-6)
        {
            std::cout << "V EXACT\n" << gradient_vertices_exact << std::endl;
            std::cout << "V APPROX\n" << gradient_vertices_approx << std::endl;
        }
//
//        if(maxError_E > 1e-6)
//        {
//            std::cout << "E EXACT\n" << gradient_edges_exact << std::endl;
//            std::cout << "E APPROX\n" << gradient_edges_approx << std::endl;
//        }

        return std::make_pair(verticesOK, edgesOK);	        
    }
    
public:
    Test_BasicGradients()
    {
    }
    
    bool test()
    {
        initGrid();
        applyRandomDeformation(0.005, 0.005);
        
        std::pair<bool, bool> errors = testGradient();
        
        return (errors.first==true and errors.second==true);
    }
};


// The list of types we want to test
typedef ::testing::Types<
TypeDefinitions<Mesh, 0, TestBasicGradients::Test_GradHeight>,
TypeDefinitions<Mesh, 1, TestBasicGradients::Test_GradHeight>,
TypeDefinitions<Mesh, 2, TestBasicGradients::Test_GradHeight>,
TypeDefinitions<Mesh, 0, TestBasicGradients::Test_GradTheta>,
TypeDefinitions<Mesh, 1, TestBasicGradients::Test_GradTheta>,
TypeDefinitions<Mesh, 2, TestBasicGradients::Test_GradTheta>
> Test_BasicGradientsImplementations;


// define the testcase (DCS)
TYPED_TEST_CASE(Test_BasicGradients, Test_BasicGradientsImplementations);

TYPED_TEST(Test_BasicGradients, IsBasicGradientCorrect) {
    EXPECT_TRUE(this->test());
}
