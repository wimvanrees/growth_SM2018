//
//  ComputeHausdorffDistance.cpp
//  Elasticity
//
//  Created by Wim van Rees on 12/28/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#include "ComputeHausdorffDistance.hpp"

#include <igl/hausdorff.h>
#include <Eigen/Geometry>
#include "WriteVTK.hpp"


Real ComputeHausdorffDistance::compute(const Eigen::Ref<const Eigen::MatrixXd> vertices_A, const Eigen::Ref<const Eigen::MatrixXd> vertices_B, const Eigen::Ref<const Eigen::MatrixXi> faces, const Real rescale) const
{
    // compute the least-squares transformation between two point sets
    //https://eigen.tuxfamily.org/dox/group__Geometry__Module.html#gab3f5a82a24490b936f8694cf8fef8e60
    
    const Eigen::MatrixXd trafo = Eigen::umeyama(vertices_A.transpose(), vertices_B.transpose(), rescale);
    const Eigen::Matrix3d rotmat = trafo.block<3,3>(0,0);
    const Eigen::Vector3d transv = trafo.block<3,1>(0,3);
    // now we have the rotation and scaling matrices to transform vertices_A into vertices_B
    
    // apply the transform
    const int nVertices = vertices_A.rows();
    Eigen::MatrixXd vertices_A_trafo(nVertices,3);
    for(int i=0;i<nVertices;++i)
    {
        const Eigen::Vector3d vertA = vertices_A.row(i);
        const Eigen::Vector3d vertA_trafo = rotmat * vertA + transv;
        vertices_A_trafo.row(i) = vertA_trafo;
    }
    
    if(bDump)
    {
        WriteVTK writerA(vertices_A, faces);
        writerA.write("haussdorf_A");
        
        WriteVTK writerB(vertices_B, faces);
        writerB.write("haussdorf_B");
        
        WriteVTK writerA_trafo(vertices_A_trafo, faces);
        writerA_trafo.write("haussdorf_A_trafo");
    }
    
    // compute the distance
    Real retval;
//    igl::hausdorff<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd, Eigen::MatrixXi, Real>(vertices_A_trafo, faces, vertices_B, faces, retval);
    
    // IGL does not handle Eigen::Ref well -- need to define a reference explicitly
    const Eigen::MatrixXd & ref_vertices_B = vertices_B;
    const Eigen::MatrixXi & ref_faces = faces;
    igl::hausdorff(vertices_A_trafo, ref_faces, ref_vertices_B, ref_faces, retval);
    
    return retval;
}
