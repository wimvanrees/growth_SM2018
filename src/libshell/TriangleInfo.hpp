//
//  TriangleInfo.hpp
//  Elasticity
//
//  Created by Wim van Rees on 04/04/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef TriangleInfo_hpp
#define TriangleInfo_hpp

#include "common.hpp"

/*! \struct TriangleInfo
 * \brief Store indices, vertex and edge vectors for a single triangle
 *
 * The ordering of the vertices (v0,v1,v2) corresponds to the face2vertices(face_idx,0:2) of the mesh
 * The ordering of the edges is e0 = v1-v0 ; e1 = v2-v1 ; e2 = v0-v2  ; (see constructor)
 */
struct TriangleInfo
{
    int face_idx;
    
    int idx_v0, idx_v1, idx_v2; /*!< vertex indices */
    int idx_e0, idx_e1, idx_e2; /*!< edge indices */
    
    /*! vertex locations */
    const Eigen::Vector3d v0;
    const Eigen::Vector3d v1;
    const Eigen::Vector3d v2;
    
    /*! edge vectors */
    const Eigen::Vector3d e0;
    const Eigen::Vector3d e1;
    const Eigen::Vector3d e2;
    
    TriangleInfo(const Eigen::Vector3d v0, const Eigen::Vector3d v1, const Eigen::Vector3d v2):
    v0(v0),v1(v1),v2(v2),e0(v1-v0),e1(v2-v1),e2(v0-v2)
    {}
    
    Eigen::Vector3d computeFaceCenter() const
    {
        return (v0 + v1 + v2)/3.0;
    }
};

#endif /* TriangleInfo_hpp */
