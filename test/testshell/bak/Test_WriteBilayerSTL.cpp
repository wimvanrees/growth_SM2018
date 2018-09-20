//
//  Test_WriteBilayerSTL.cpp
//  Elasticity
//
//  Created by Wim van Rees on 7/29/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#include "Test_WriteBilayerSTL.hpp"
#include "Geometry.hpp"
#include "WriteSTLVolumetric.hpp"
#include "WriteSTLVolumetricCurved.hpp"

void Test_WriteBilayerSTL::init()
{
    // initialize geometry
    const Real Lx = 0.5;
    const Real Ly = 5.0;
//    const Real edgeLength = 0.1;
//    RectangularPlate_RightAngle geometry(Lx, Ly, edgeLength, false, false);
    const Real relArea = 0.01;
    RectangularPlate geometry(Lx, Ly, relArea, {false,false}, {false,false});
    
    mesh.init(geometry, false);
}

void Test_WriteBilayerSTL::run()
{
//    test_twotriangles();
//    test_annulus();
    test_writeSingleLayerCurved();
}

void Test_WriteBilayerSTL::test_twotriangles()
{
    TwoTriangles geometry;
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi face2vertices;
    Eigen::MatrixXb vertices_bc;
    geometry.get(vertices,face2vertices,vertices_bc);
    const int nFaces = face2vertices.rows();
    Eigen::VectorXi faces_STL(nFaces);
    for(int i=0;i<nFaces;++i)
        faces_STL(i) = i;
    
    WriteSTLVolumetric::write(face2vertices, vertices, faces_STL, 0.01, "twotriangles");
}

void Test_WriteBilayerSTL::test_annulus()
{
    const Real innerR = 0.1;
    const Real outerR = 0.5;
    const Real edgeLength = 0.025;
    AnnulusPlate geometry(outerR, innerR, edgeLength);
    
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi face2vertices;
    Eigen::MatrixXb vertices_bc;
    geometry.get(vertices,face2vertices,vertices_bc);
    const int nFaces = face2vertices.rows();
    Eigen::VectorXi faces_STL(nFaces/2);
    for(int i=0;i<nFaces/2;++i)
        faces_STL(i) = i;
    WriteSTLVolumetric::write(face2vertices, vertices, faces_STL, 0.01,  "annulus_1");
    for(int i=0;i<nFaces/2;++i)
        faces_STL(i) = i+nFaces/2;
    WriteSTLVolumetric::write(face2vertices, vertices, faces_STL, 0.01,  "annulus_2");
    
}


void Test_WriteBilayerSTL::test_writeSingleLayerCurved()
{
    // create a cylinder
    const Real radius = 1.0;
    const Real length = 4.0;
    const Real edgeLength = 0.1;
    CylinderTube geometry(radius, length, edgeLength);
    mesh.init(geometry);
    
    // create random thickness function
    std::mt19937 gen;
    std::uniform_real_distribution<Real> distr(0.1, 0.5);
    Eigen::VectorXd thickness(mesh.getNumberOfFaces());
    for(int i=0;i<mesh.getNumberOfFaces();++i)
        thickness(i) = distr(gen);
    
    WriteSTLVolumetricCurved::write(mesh.getTopology(), mesh.getCurrentConfiguration(), thickness, "cylinder_varthick");
}