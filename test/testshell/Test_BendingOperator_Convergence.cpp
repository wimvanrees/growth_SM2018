//
//  Test_BendingOperator_Convergence.cpp
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
template <typename tM, template<typename, typename, MeshLayer> class tB>
struct TypeDefinitions
{
    typedef tM tMesh;
    typedef tB<tM, Material_Isotropic, single> tBendingOperator;
};

template <typename tCombinedTypes>
class Test_BendingOperatorDCS_Convergence : public ::testing::Test
{
public:
    typedef typename tCombinedTypes::tMesh tMesh;
    typedef typename tCombinedTypes::tBendingOperator tBendingOperator;
    
protected:
    tMesh mesh;
    const Real radius;
    
    // create a method
    virtual void initGrid(const int nPointsAlongBoundary)
    {
        // create a plate
        CircularPlate geometry(radius, nPointsAlongBoundary, true);
        geometry.setQuiet();
        mesh.init(geometry);
    }
    
    void applyDeformation(const Real w0)
    {
        const Real a = radius;
        
        // vertex positions
        auto displacement = [&](Eigen::Vector3d vertex) -> Eigen::Vector3d
        {
            const Real x = vertex(0);
            const Real y = vertex(1);
            const Real radiussq = x*x + y*y;
            
            const Real w = w0*std::pow(1.0 - radiussq/(a*a), 2);
            
            Eigen::Vector3d retval;
            retval << vertex(0), vertex(1), vertex(2)-w;
            return retval;
        };

        mesh.changeVertices(displacement);
        
        // edge directors
        auto circular_normal = [&](const Eigen::Vector3d pos) -> Eigen::Vector3d
        {
            const Real x = pos(0);
            const Real y = pos(1);
            const Real dwdx = -2.0*w0*(1.0 - (x*x+y*y)/(a*a)) * 2*x/(a*a) ;
            const Real dwdy = -2.0*w0*(1.0 - (x*x+y*y)/(a*a)) * 2*y/(a*a) ;
            
            Eigen::Vector3d n;
            n << dwdx , dwdy, +1.0;
            return n.normalized();
        };
        
#if 1==0
        mesh.setEdgeDirectors(circular_normal);
#else
        {
            const int nEdges = mesh.getNumberOfEdges();
            Eigen::MatrixXd edgeNormals(nEdges,3);
            
            const auto edge2vertices = mesh.getTopology().getEdge2Vertices();
            const auto restvertices = mesh.getRestConfiguration().getVertices();
            
            for(int i=0;i<nEdges;++i)
            {
                const int idx_v0 = edge2vertices(i,0);
                const int idx_v1 = edge2vertices(i,1);
                
                const Eigen::Vector3d v0 = restvertices.row(idx_v0);
                const Eigen::Vector3d v1 = restvertices.row(idx_v1);
                
                const Eigen::Vector3d midedge_pos = 0.5*(v0 + v1);
                const Eigen::Vector3d midedge_normal = circular_normal(midedge_pos);
                for(int d=0;d<3;++d)
                    edgeNormals(i,d) = midedge_normal(d);
            }
            
            mesh.getCurrentConfiguration().setEdgeDirectors(mesh.getTopology(), edgeNormals);
        }
#endif

        mesh.updateDeformedConfiguration();
    }
    
    Real getExactEnergy(const Real Young, const Real nu, const Real thickness, const Real w0)
    {
        // some abbreviations
        const Real a = radius;
        const Real E = Young;
        const Real h = thickness;
        
        const Real D = E*h*h*h/(12*(1.0-nu*nu));
        const Real exactEng = 32*M_PI/3.0 * w0*w0/(a*a) * D; // from timoshenko page 400
        
        return exactEng;
        
    }
    
    std::pair<int, Real> computeEnergy(const int nPointsAlongBoundary, Real & exactEnergy)
    {
        // init
        initGrid(nPointsAlongBoundary);
        const int nFaces = mesh.getNumberOfFaces();
        
        const Real E = 1e6;
        const Real nu = 0.3;//0.05;
        const Real h = 0.001*radius;
        
        // deform
        const Real w0 = 5*h; // free parameter
        applyDeformation(w0);
        
        // compute energy        
        MaterialProperties_Iso_Constant matprop(E, nu, h);
        tBendingOperator bending(matprop);
        bending.compute(mesh);
        // extract the bending energy from the energy array
        Real computedEnergy = -1;
        std::vector<std::pair<std::string, Real>> computedEnergies;
        bending.addEnergy(computedEnergies);
        for(size_t i=0;i<computedEnergies.size();++i)
            if((computedEnergies[i].first).find("bending")!=std::string::npos)
            {
                computedEnergy = computedEnergies[i].second;
                break;
            }
        
        exactEnergy = getExactEnergy(E, nu, h, w0);
        
        return std::make_pair(nFaces, computedEnergy);
    }
    
public:
    Test_BendingOperatorDCS_Convergence():
    radius(1.0)
    {
    }
    
    bool testConvergence()
    {
        // add sqrt(2) increments in between to get nfaces increasing by 2 instead of 4 between successive guys
        // number of faces ~ N^2/(sqrt(3)*pi)
        std::vector<int> nPoints = {16, 23, 32, 45, 64, 91, 128, 181, 256, 362, 512};//, 724, 1024, 1448, 2048};
        
        // compute all energies
        const size_t nRuns = nPoints.size();
        std::vector<std::pair<int, Real>> results(nRuns);
        Real exactEnergy;
        for(size_t i=0;i<nRuns;++i)
            results[i] = computeEnergy(nPoints[i], exactEnergy);
        
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
> Test_BendingOperatorDCS_ConvergenceImplementations;

// define the testcase
TYPED_TEST_CASE(Test_BendingOperatorDCS_Convergence, Test_BendingOperatorDCS_ConvergenceImplementations);

TYPED_TEST(Test_BendingOperatorDCS_Convergence, IsConverged) {
    // Inside the test body, you can refer to the type parameter by
    // TypeParam, and refer to the fixture class by TestFixture.  We
    // don't need them in this example.
    EXPECT_TRUE(this->testConvergence());
}

