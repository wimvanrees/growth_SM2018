#include "gtest/gtest.h"

#include "common.hpp"
#include "Geometry.hpp"
#include "Mesh.hpp"
#include "MaterialProperties.hpp"
#include "CombinedOperator_Parametric.hpp"

TEST(test_energy_staticeq, energy_is_zero_in_equilibrium)
{
    // set up a simple geometry
    const Real Lx = 1.0;
    const Real Ly = 1.0;
    const Real relArea = 0.005;
    RectangularPlate geometry(Lx, Ly, relArea, {false,false}, {false,false});
    geometry.setQuiet();
    
    // create a mesh
    Mesh mesh;
    mesh.init(geometry, false);
    
    // set up an energy operator
    const Real Young = 1.0;
    const Real poisson = 0.3;
    const Real thickness = 0.01*Lx;
    MaterialProperties_Iso_Constant matprop(Young, poisson, thickness);
    CombinedOperator_Parametric<Mesh, Material_Isotropic, single> engOp(matprop);
    
    // compute energy
    const Real energy = engOp.compute(mesh);

    // make sure the energy is close to
    const int nVertices = mesh.getNumberOfVertices();
    const Real tol =  Young*thickness*nVertices*std::numeric_limits<Real>::epsilon();
    
    std::cout << "[----------] energy = \t " << energy << "\t" << tol << std::endl;
    ASSERT_NEAR(energy, 0.0, tol);
}
