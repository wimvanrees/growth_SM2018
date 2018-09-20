//
//  common.hpp
//  Elasticity
//
//  Created by Wim van Rees on 2/15/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef common_hpp
#define common_hpp

#include <iostream>
#include <cassert>
#include <random>
#include <limits>
#include <utility>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <map>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <array>
#include <set>
#include <memory>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/StdVector>

//#include "ArgumentParser.h"

typedef double Real;
typedef unsigned long long tUint;
typedef std::vector<Eigen::Matrix2d, Eigen::aligned_allocator<Eigen::Matrix2d>> tVecMat2d;
typedef std::vector<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d>> tVecMat3d;
enum MeshLayer {single, bottom, top};

//#define USETANTHETA

// little helper to avoid unused warnings when something is called in release mode and variable only checked in debug
#define _unused(x) ((void)(x))

namespace rnd
{
    /*! global random number generator */    
    static std::mt19937 gen;
}


namespace helpers
{
    
    /*! convert an integer to a string with a specified number of zeros */
    inline std::string ToString(int value,int digitsCount)
    {
        std::ostringstream os;
        os<<std::setfill('0')<<std::setw(digitsCount)<<value;
        return os.str();
    }
    
    /*! check if two numbers are close (within some accuracy) */
    inline bool isclose(const Real val1, const Real val2, const bool printOnly=false, const Real rtol=1e-5, const Real atol=1e-8)
    {
        const Real diff = std::abs(val1-val2) ;
        bool retval = true;
        if(diff > (atol + rtol * std::abs(val2)))
            retval = false;
        if(diff > (atol + rtol * std::abs(val1)))
            retval = false;

        if(printOnly && not retval)
        {
            std::cout << "PROBLEM : " << val1 << "\t" << val2 << "\t" << diff << std::endl;
            retval = true;
        }
        
        return retval;
    }
 
    
    inline bool isPointInProjectedTriangle(const Real x, const Real y, const Eigen::Vector3d & v0, const Eigen::Vector3d & v1, const Eigen::Vector3d & v2)
    {
        // dont use the third dimension
        const Real xp[3] = {v0(0),v1(0),v2(0)};
        const Real yp[3] = {v0(1),v1(1),v2(1)};
        
        const int npoly = 3; // triangle 
        
        int i,j,c=0;
        for (i = 0, j = npoly-1; i < npoly; j = i++) {
            if ((((yp[i] <= y) && (y < yp[j])) ||
                 ((yp[j] <= y) && (y < yp[i]))) &&
                (x < (xp[j] - xp[i]) * (y - yp[i]) / (yp[j] - yp[i]) + xp[i]))
                c = !c;
        }
        return c;
    }
    
    inline Eigen::Vector3d getBaryCentricWeights2D(const Eigen::Vector3d & v0, const Eigen::Vector3d & v1, const Eigen::Vector3d & v2, const Real x, const Real y)
    {
        assert(isPointInProjectedTriangle(x, y, v0, v1, v2));
        
        Eigen::Vector3d me;
        me << x, y, 0;
        // compute my point in barycentric coordinates of undeformed configuration
        // Compute barycentric coordinates (u, v, w) for
        // point p with respect to triangle (a, b, c)
        // http://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
        const Eigen::Vector3d tmp0 = v1 - v0;
        const Eigen::Vector3d tmp1 = v2 - v0;
        const Eigen::Vector3d tmp2 = me - v0;
        const Real d00 = tmp0.dot(tmp0);
        const Real d01 = tmp0.dot(tmp1);
        const Real d11 = tmp1.dot(tmp1);
        const Real d20 = tmp2.dot(tmp0);
        const Real d21 = tmp2.dot(tmp1);
        const Real denom = d00 * d11 - d01 * d01;
        const Real lambda1 = (d11 * d20 - d01 * d21) / denom;
        const Real lambda2 = (d00 * d21 - d01 * d20) / denom;
        const Real lambda0 = 1.0 - lambda1 - lambda2;
        // if me == v0 --> lambda0 = 1
        // if me == v1 --> lambda1 = 1
        // if me == v2 --> lambda2 = 1
        
        Eigen::Vector3d retval;
        retval << lambda0, lambda1, lambda2;
        return retval;
    }
    
    template<int component>
    inline Real getDisplacementOfPointInRestTriangle(const Eigen::Vector3d & rv0, const Eigen::Vector3d & rv1, const Eigen::Vector3d & rv2, const Eigen::Vector3d & v0, const Eigen::Vector3d & v1, const Eigen::Vector3d & v2, const Real x, const Real y)
    {

        const Eigen::Vector3d lambdas = getBaryCentricWeights2D(rv0, rv1, rv2, x, y);
        
        const Real displ_0 = (v0 - rv0)(component);
        const Real displ_1 = (v1 - rv1)(component);
        const Real displ_2 = (v2 - rv2)(component);
        
        // get my displacement
        return -(lambdas(0)*displ_0 + lambdas(1)*displ_1 + lambdas(2)*displ_2);
    }
    
    
    inline Real interpolateOverTriangle(const Eigen::Vector3d & v0, const Eigen::Vector3d & v1, const Eigen::Vector3d & v2, const Real x, const Real y)
    {
        // linear interpolation over the triangle
        const Real A = v0(1)*(v1(2)-v2(2)) + v1(1)*(v2(2)-v0(2)) + v2(1)*(v0(2)-v1(2));
        const Real B = v0(2)*(v1(0)-v2(0)) + v1(2)*(v2(0)-v0(0)) + v2(2)*(v0(0)-v1(0));
        const Real C = v0(0)*(v1(1)-v2(1)) + v1(0)*(v2(1)-v0(1)) + v2(0)*(v0(1)-v1(1));
        const Real D = v0(0)*(v1(1)*v2(2)-v1(2)*v2(1)) + v1(0)*(v2(1)*v0(2)-v2(2)*v0(1)) + v2(0)*(v0(1)*v1(2)-v0(2)*v1(1));
        
        const Real ptZ = (D-A*x-B*y)/C;

        return ptZ;
    }
    
    inline Real computeDistanceToLineSegment(const Eigen::Vector2d & l1, const Eigen::Vector2d & l2, const Eigen::Vector2d & pt)
    {
        //http://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
        const Real Lsq = (l1-l2).squaredNorm();
        if(Lsq < std::numeric_limits<Real>::epsilon())
            return (l1-pt).norm();
        
        const Real t =  std::max(0.0, std::min(1.0, (pt-l1).dot(l2-l1) / Lsq));
        const Eigen::Vector2d projection = l1 + t * (l2 - l1);
        return (projection - pt).norm();
    }
    
    inline int getIntersectionPointBetweenVectors(const Eigen::Vector2d & pt1, const Eigen::Vector2d & n1, const Eigen::Vector2d & pt2, const Eigen::Vector2d & n2, Eigen::Vector2d & result)
    {
        // return 0 if intersection is found
        // return 1 if an intersection is found, but beyond the extent of pt1+n1 and pt2+n2
        // return 2 if no intersection is found (parallel vectors)
        
        const Eigen::Vector2d n1_hat = n1.normalized();
        const Eigen::Vector2d n2_hat = n2.normalized();
        
        const Real pt1_vec1 = pt1.dot(n1_hat);
        const Real pt1_vec2 = pt1.dot(n2_hat);
        const Real pt2_vec1 = pt2.dot(n1_hat);
        const Real pt2_vec2 = pt2.dot(n2_hat);
        const Real vec1_vec2 = n1_hat.dot(n2_hat);
        
        // deal with orthogonal
        if(std::abs(vec1_vec2 - 1) < std::numeric_limits<Real>::epsilon())
        {
            if((pt1-pt2).norm()  < std::numeric_limits<Real>::epsilon())
            {
                result = pt1;
                return 0;
            }
            else
            {
                return 2;
            }
        }
        
        const Real beta1 = (pt2_vec1 + (pt1_vec2 - pt2_vec2) * vec1_vec2 - pt1_vec1) / (1.0 - std::pow(vec1_vec2,2));
        const Real beta2 = (pt1_vec2 + (pt2_vec1 - pt1_vec1) * vec1_vec2 - pt2_vec2) / (1.0 - std::pow(vec1_vec2,2));
        
        //const Eigen::Vector2d retval1 = pt1 + beta1*n1_hat;
        const Eigen::Vector2d retval2 = pt2 + beta2*n2_hat;
        
        // set the intersection point
        result = retval2;
        
        const bool inRange_1 = (beta1 >= 0 && beta1 <= n1.norm());
        const bool inRange_2 = (beta2 >= 0 && beta2 <= n2.norm());
        return (inRange_1 && inRange_2) ? 0 : 1;
    }
    
    inline std::string removeExtension(const std::string input, const std::string extension)
    {
        // remove extension from input. Assume extension is given as (eg) .vtp (including the period)
        // if not found, return input unaltered
        const std::size_t found_ext = input.rfind(extension);
        const std::string output = (found_ext != std::string::npos ? input.substr(0, found_ext) : input);
        return output;
    }
    
    inline void catastrophe(std::string s, const char* file, const int line, const bool isFatal=true)    // write ``error: s and exit program
    {
        std::cerr << "Something went wrong, error: " << s << std::endl;
        std::cerr << "Called from file : " << file << " at line number " << line << std::endl;
        if(isFatal)
        {
            std::cout << "Exiting..." << std::endl;
            std::exit(1);
        }
    }
    
    template<class Matrix>
    inline void write_matrix_binary(const std::string & filename, const Matrix& matrix)
    {
        // from https://stackoverflow.com/a/25389481
        std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
        typename Matrix::Index rows=matrix.rows(), cols=matrix.cols();
        out.write((char*) (&rows), sizeof(typename Matrix::Index));
        out.write((char*) (&cols), sizeof(typename Matrix::Index));
        out.write((char*) matrix.data(), rows*cols*sizeof(typename Matrix::Scalar) );
        out.close();
    }
    
    template<class Matrix>
    inline void read_matrix_binary(const std::string & filename, Matrix& matrix)
    {
        // from https://stackoverflow.com/a/25389481
        std::ifstream in(filename, std::ios::in | std::ios::binary);
        typename Matrix::Index rows=0, cols=0;
        in.read((char*) (&rows),sizeof(typename Matrix::Index));
        in.read((char*) (&cols),sizeof(typename Matrix::Index));
        matrix.resize(rows, cols);
        in.read( (char *) matrix.data() , rows*cols*sizeof(typename Matrix::Scalar) );
        in.close();
    }
    
    template<class Matrix>
    inline void write_vecmat_binary(const std::string & filename, const std::vector<Matrix, Eigen::aligned_allocator<Matrix>> & vecmat)
    {
        // from https://stackoverflow.com/a/25389481
        std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
        const size_t vecsize = vecmat.size();
        if(vecsize==0) return;
        typename Matrix::Index rows=vecmat[0].rows(), cols=vecmat[0].cols();
        out.write((char*) &vecsize, sizeof(vecsize));
        out.write((char*) (&rows), sizeof(typename Matrix::Index));
        out.write((char*) (&cols), sizeof(typename Matrix::Index));
        out.write((char*)&vecmat[0], vecsize * rows*cols*sizeof(typename Matrix::Scalar));
        out.close();
    }
    
    template<class Matrix>
    inline void read_vecmat_binary(const std::string & filename, std::vector<Matrix, Eigen::aligned_allocator<Matrix>> & vecmat)
    {
        std::ifstream in(filename, std::ios::in | std::ios::binary);
        size_t vecsize = 0;
        typename Matrix::Index rows=0, cols=0;
        
        in.read((char*) (&vecsize),sizeof(vecsize));
        assert(vecsize > 0);
        
        in.read((char*) (&rows),sizeof(typename Matrix::Index));
        in.read((char*) (&cols),sizeof(typename Matrix::Index));
        
        assert(rows == Matrix::RowsAtCompileTime);
        assert(cols == Matrix::ColsAtCompileTime);
        
        vecmat.resize(vecsize);
        in.read( (char *) &vecmat[0] , vecsize * rows*cols*sizeof(typename Matrix::Scalar) );
        in.close();
    }
}

namespace Eigen
{
    /*! typedef vector of booleans */
    typedef Matrix<bool, 3, 1> Vector3b;
    
    /*! typedef vector of booleans */
    typedef Matrix<bool, Dynamic, 1> VectorXb;
    typedef Matrix<bool, Dynamic, Dynamic> MatrixXb;
}

#endif /* common_hpp */
