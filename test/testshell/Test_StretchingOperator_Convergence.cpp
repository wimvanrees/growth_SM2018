//
//  Test_StretchingOperator_Convergence.cpp
//  Elasticity
//
//  Created by Wim van Rees on 2/24/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#include "gtest/gtest.h"

#include "Test_Helpers_Convergence.hpp"
#include "common.hpp"
#include "Geometry.hpp"
#include "Mesh.hpp"
#include "MaterialProperties.hpp"
#include "CombinedOperator_Parametric.hpp"

// write a templatized test case. we have two templates: mesh and operator type, but google test only provides interface for one : use a structure that is itself templated that we can use
template <typename tM, template<typename, typename, MeshLayer> class tS>
struct TypeDefinitions
{
    typedef tM tMesh;
    typedef tS<tM, Material_Isotropic, single> tStretchingOperator;
};

template <typename tCombinedTypes>
class Test_StretchingOperatorDCS_Convergence : public ::testing::Test
{
public:
    typedef typename tCombinedTypes::tMesh tMesh;
    typedef typename tCombinedTypes::tStretchingOperator tStretchingOperator;
    
protected:
    tMesh mesh;
    const Real halfSideLength;

    // create a method
    virtual void initGrid(const Real maxRelTriangleArea)
    {
        // create a plate
        RectangularPlate geometry(halfSideLength, halfSideLength, maxRelTriangleArea, {true, true}, {true, true});
        geometry.setQuiet();
        mesh.init(geometry);
    }
    
    void applyDeformation(const Real w0, const Real c)
    {
        const Real a = halfSideLength;
        auto add_displacement = [a, w0, c](Eigen::Vector3d in)
        {
            const Real x = in(0);
            const Real y = in(1);
            const Real w = w0 * std::cos(M_PI*x/(2*a))*std::cos(M_PI*y/(2*a));
            const Real u = c * std::sin(M_PI*x/a)*std::cos(M_PI*y/(2*a));
            const Real v = c * std::sin(M_PI*y/a)*std::cos(M_PI*x/(2*a));
            Eigen::Vector3d retval;
            retval << x + u, y + v, in(2) - w;
            return retval;
        };
        mesh.changeVertices(add_displacement);
        mesh.updateDeformedConfiguration();
    }
    
    Real getExactEnergy(const Real Young, const Real nu, const Real thickness, const Real w0, const Real c)
    {
        // some abbreviations
        const Real a = halfSideLength;
        const Real E = Young;
        const Real h = thickness;
        
        // timoshenko : for nu = 0.25 (p 420 equation h - large deflections of plates)
        //    const Real exactEng = E*h/(7.5) * ( 5./64*std::pow(M_PI*w0,4)/(a*a) - 17./6 * std::pow(M_PI*w0,2)*c/a + c*c*(35*M_PI*M_PI/4.0 + 80./9) );
        
        // instead general expression for any nu
        const Real exactEng = E*h*( (-45*std::pow(M_PI,4)*std::pow(w0,4) - 384*a*c*std::pow(M_PI,2)*std::pow(w0,2)*(-5 + 3*nu) + 64*std::pow(a,2)*std::pow(c,2)*(9*std::pow(M_PI,2)*(-9 + nu) - 64*(1 + nu)))/(4608.*std::pow(a,2)*(-1 + std::pow(nu,2))) ) ;
        return exactEng;
        
    }
    
    std::pair<int, Real> computeEnergy(const Real maxRelTriangleArea, Real & exactEnergy)
    {
        // init
        initGrid(maxRelTriangleArea);
        const int nFaces = mesh.getNumberOfFaces();
        
        const Real E = 1e6;
        const Real nu = 0.3;
        const Real h = 0.001;
        
        // deform
//        const Real w0 = 0.1; // free parameter
//        const Real c = 0.05;//

//        const Real q = 1e-6;// loading : equilibrium configuration gives rise to this energy:
//        const Real exactw0_prefac = (8*std::pow(2,2.0/3.0))/(std::pow(M_PI,2)*std::pow((-45*std::pow(M_PI,2)*(-9 + nu) - 64*(20 + nu*(-35 + 9*nu)))/((-1 + std::pow(nu,2))*(9*std::pow(M_PI,2)*(-9 + nu) - 64*(1 + nu))),1.0/3.0));
//        const Real w0 = exactw0_prefac * halfSideLength * std::pow(halfSideLength*q/(h*E), 1.0/3.0);
        const Real w0 = 5*h;
        const Real c = (3*std::pow(M_PI,2)*std::pow(w0,2)*(-5 + 3*nu))/(halfSideLength*(9*std::pow(M_PI,2)*(-9 + nu) - 64*(1 + nu)));
        
        applyDeformation(w0, c);
        
        // compute energy
        MaterialProperties_Iso_Constant matprop(E, nu, h);
        tStretchingOperator stretching(matprop);
        stretching.compute(mesh);
        // exact the stretching energy from the energy array
        Real computedEnergy = -1;
        std::vector<std::pair<std::string, Real>> computedEnergies;
        stretching.addEnergy(computedEnergies);
        for(size_t i=0;i<computedEnergies.size();++i)
            if((computedEnergies[i].first).find("stretching")!=std::string::npos)
            {
                computedEnergy = computedEnergies[i].second;
                break;
            }
        
        exactEnergy = getExactEnergy(E, nu, h, w0, c);
        
        return std::make_pair(nFaces, computedEnergy);
    }
    
public:
    Test_StretchingOperatorDCS_Convergence():
    halfSideLength(1.0)
    {
    }
    
    bool testConvergence()
    {
        const std::vector<Real> relAreas = {1.0/16, 1.0/32, 1.0/64, 1.0/128, 1.0/256, 1.0/512, 1.0/1024, 1.0/2048, 1.0/4096};//, 1.0/8192, 1.0/16384, 1.0/32768, 1.0/65536, 1.0/131072};
        
        // compute all energies
        const size_t nRuns = relAreas.size();
        std::vector<std::pair<int, Real>> results(nRuns);
        Real exactEnergy;
        for(size_t i=0;i<nRuns;++i)
            results[i] = computeEnergy(relAreas[i], exactEnergy);
        
        // compute the slopes
        const std::vector<Real> slopes = Test_Helpers::Convergence::computeConvergence(results, exactEnergy, true);
        
        // test the slopes
        const int targetSlope = 1.0; // linear in the area
        const Real acceptableError = 0.2;// accept anything between 0.9 and 1.1
        const bool acceptedConvergence = Test_Helpers::Convergence::testConvergence(slopes, targetSlope, acceptableError);
//        if(not acceptedConvergence)
//            computeConvergence(results, exactEnergy, true);
        
        return acceptedConvergence;
    }
};

// The list of types we want to test.
typedef ::testing::Types<
TypeDefinitions<Mesh, CombinedOperator_Parametric>
> Test_StretchingOperatorDCS_ConvergenceImplementations;

// define the testcase
TYPED_TEST_CASE(Test_StretchingOperatorDCS_Convergence, Test_StretchingOperatorDCS_ConvergenceImplementations);

TYPED_TEST(Test_StretchingOperatorDCS_Convergence, IsConverged) {
    // Inside the test body, you can refer to the type parameter by
    // TypeParam, and refer to the fixture class by TestFixture.  We
    // don't need them in this example.
        EXPECT_TRUE(this->testConvergence());
}


