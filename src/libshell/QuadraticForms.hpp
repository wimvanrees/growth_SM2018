//
//  QuadraticForms.hpp
//  Elasticity
//
//  Created by Wim van Rees on 1/5/17.
//  Copyright © 2017 Wim van Rees. All rights reserved.
//

#ifndef QuadraticForms_hpp
#define QuadraticForms_hpp

#include "common.hpp"
#include "ExtendedTriangleInfo.hpp"

struct QuadraticForm
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    
    Eigen::Matrix2d form;
};

struct QuadraticForm_WithGradient_Verts : QuadraticForm
{
    QuadraticFormGradientData_Verts gradform;    
};

struct QuadraticForm_WithGradient : QuadraticForm
{
    QuadraticFormGradientData gradform;    
};

template<bool withGradient>
struct FirstFundamentalForm : std::conditional<withGradient, QuadraticForm_WithGradient_Verts, QuadraticForm>::type
{
    FirstFundamentalForm(const ExtendedTriangleInfo & info)
    {
        set(info);
    }
    
    template<bool U = withGradient> typename std::enable_if<!U, void>::type
    set(const ExtendedTriangleInfo & info)
    {
        this->form = info.computeFirstFundamentalForm();
    }
    
    template<bool U = withGradient> typename std::enable_if<U, void>::type
    set(const ExtendedTriangleInfo & info)
    {
        this->form = info.computeFirstFundamentalForm();
        this->gradform = info.computeGradFirstFundamentalForm();
    }
};

template<bool withGradient>
struct SecondFundamentalForm : std::conditional<withGradient, QuadraticForm_WithGradient, QuadraticForm>::type
{
    SecondFundamentalForm(const ExtendedTriangleInfo & info)
    {
        set(info);
    }
    
    template<bool U = withGradient> typename std::enable_if<!U, void>::type
    set(const ExtendedTriangleInfo & info)
    {
        this->form = info.computeSecondFundamentalForm();
    }
    
    template<bool U = withGradient> typename std::enable_if<U, void>::type
    set(const ExtendedTriangleInfo & info)
    {
        this->form = info.computeSecondFundamentalForm();
        this->gradform = info.computeGradSecondFundamentalForm();
    }
};






#endif /* QuadraticForms_hpp */
