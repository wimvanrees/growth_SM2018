//
//  Triangulate_TRIANGLE.hpp
//  Elasticity
//
//  Created by Wim van Rees on 3/31/17.
//  Copyright © 2017 Wim van Rees. All rights reserved.
//

#ifndef Triangulate_TRIANGLE_hpp
#define Triangulate_TRIANGLE_hpp

#include "common.hpp"

namespace TRIANGLE
{

    // Triangulate the interior of a polygon using the triangle library.
    //
    // Inputs:
    //   V #V by 2 list of 2D vertex positions
    //   E #E by 2 list of vertex ids forming unoriented edges of the boundary of the polygon
    //   H #H by 2 coordinates of points contained inside holes of the polygon
    //   A #A by 3 : first two are coordinates of points contained inside attribute regions of the polygon. third is area constraint (we do not use the attribute property)
    //   flags  string of options pass to triangle (see triangle documentation)
    // Outputs:
    //   V2  #V2 by 2  coordinates of the vertives of the generated triangulation
    //   F2  #F2 by 3  list of indices forming the faces of the generated triangulation
    //
    // TODO: expose the option to prevent Steiner points on the boundary
    //

    // this is modified from the IGL method
    inline void triangulate(
                                const Eigen::MatrixXd& V,
                                const Eigen::MatrixXi& E,
                                const Eigen::MatrixXd& H,
                                const Eigen::MatrixXd& A,
                                const std::string flags,
                                Eigen::MatrixXd& V2,
                                Eigen::MatrixXi& F2)
    {
        using namespace std;
        using namespace Eigen;
        
        // Prepare the flags
        string full_flags = flags + "pzB";
        
        // Prepare the input struct
        triangulateio in;
        
        assert(V.cols() == 2);
        
        in.numberofpoints = V.rows();
        in.pointlist = (double*)calloc(in.numberofpoints*2,sizeof(double));
        for (unsigned i=0;i<V.rows();++i)
            for (unsigned j=0;j<2;++j)
                in.pointlist[i*2+j] = V(i,j);
        
        in.numberofpointattributes = 0;
        in.pointmarkerlist = (int*)calloc(in.numberofpoints,sizeof(int));
        for (unsigned i=0;i<V.rows();++i)
            in.pointmarkerlist[i] = 1;
        
        in.trianglelist = NULL;
        in.numberoftriangles = 0;
        in.numberofcorners = 0;
        in.numberoftriangleattributes = 0;
        in.triangleattributelist = NULL;
        
        in.numberofsegments = E.rows();
        in.segmentlist = (int*)calloc(in.numberofsegments*2,sizeof(int));
        for (unsigned i=0;i<E.rows();++i)
            for (unsigned j=0;j<2;++j)
                in.segmentlist[i*2+j] = E(i,j);
        in.segmentmarkerlist = (int*)calloc(in.numberofsegments,sizeof(int));
        for (unsigned i=0;i<E.rows();++i)
            in.segmentmarkerlist[i] = 1;
        
        in.numberofholes = H.rows();
        in.holelist = (double*)calloc(in.numberofholes*2,sizeof(double));
        for (unsigned i=0;i<H.rows();++i)
            for (unsigned j=0;j<2;++j)
                in.holelist[i*2+j] = H(i,j);

        in.numberofregions = A.rows();
        in.regionlist = (double*)calloc(in.numberofregions*4,sizeof(double));
        for (unsigned i=0;i<A.rows();++i)
        {
            // positions
            for (unsigned j=0;j<2;++j)
                in.regionlist[i*4+j] = A(i,j);
            // attributes
            in.regionlist[i*4+2] = 0;
            // area
            in.regionlist[i*4+3] = A(i,2);
        }
        
        // Prepare the output struct
        triangulateio out;
        
        out.pointlist = NULL;
        out.trianglelist = NULL;
        out.segmentlist = NULL;
        
        // Call triangle
        ::triangulate(const_cast<char*>(full_flags.c_str()), &in, &out, 0);
        
        // Return the mesh
        V2.resize(out.numberofpoints,2);
        for (unsigned i=0;i<V2.rows();++i)
            for (unsigned j=0;j<2;++j)
                V2(i,j) = out.pointlist[i*2+j];
        
        F2.resize(out.numberoftriangles,3);
        for (unsigned i=0;i<F2.rows();++i)
            for (unsigned j=0;j<3;++j)
                F2(i,j) = out.trianglelist[i*3+j];
        
        // Cleanup in
        free(in.pointlist);
        free(in.pointmarkerlist);
        free(in.segmentlist);
        free(in.segmentmarkerlist);
        free(in.holelist);
        free(in.regionlist);
        
        // Cleanup out
        free(out.pointlist);
        free(out.trianglelist);
        free(out.segmentlist);        
    }
}

#endif /* Triangulate_TRIANGLE_hpp */
