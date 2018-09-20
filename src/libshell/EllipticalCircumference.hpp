//
//  EllipticalCircumference.hpp
//  Elasticity
//
//  Created by Wim van Rees on 9/17/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef EllipticalCircumference_hpp
#define EllipticalCircumference_hpp

#ifdef USEGSL
#include <gsl/gsl_sf_ellint.h>
#endif

struct EllipticalCircumference
{
    const Real radiusX;
    const Real radiusY;
    const Real ecc;
    const Real circum_rel;
    
    EllipticalCircumference(const Real rX, const Real rY):
    radiusX(rX),
    radiusY(rY),
    ecc(std::sqrt(1.0 - std::pow(radiusY/radiusX, 2))),
    circum_rel(getArcLengthRel(2.0*M_PI))
    {
        assert(rX >= 0);
        assert(rY >= 0);
        assert(rX >= rY);
#ifndef USEGSL
        std::cout << "PROBLEM : To use EllipticalCircumference : need to compile with GSL library\n";
#endif
    }
    
    
    Real eval_ellint_E(const Real phi) const
    {
#ifdef USEGSL
        return gsl_sf_ellint_E(phi,ecc, 0);
#else
        return phi; // closest assumption for circle
#endif
    }
    
    Real getCircumference() const
    {
        return radiusX * circum_rel;
    }
    
    Real getArcLengthRel(const Real theta) const
    {
        return eval_ellint_E(theta);
    }
    
    Real getArcLengthRelDerivative(const Real theta) const
    {
        return std::sqrt(1.0 - std::pow(ecc*std::sin(theta),2));
    }
    
    Real getThetaFromArclength(const Real arclength) const
    {
        const Real arclength_rel = arclength / radiusX;
        
        // solve Newton's equations
        const Real tol = 1e-6;
        const int nMaxIters = 100;
        
        Real theta = arclength_rel; // for circle this would be correct
        
        Real err = 2*tol;
        int iter = 0;
        
        while(err > tol)
        {
            const Real fx = getArcLengthRel(theta) - arclength_rel; // want to have s(theta_sol) = arclength_rel
            const Real dfx = getArcLengthRelDerivative(theta);
            
            if(dfx < std::numeric_limits<Real>::epsilon()) break;
            
            const Real diff = fx/dfx;
            theta -= diff;
            
            err = std::abs(diff);
            
            iter += 1;
            if(iter > nMaxIters) break;
        }
        
        return theta;
    }
    
    std::pair<Real, Real> getPointFromArcLength(const Real arclength) const
    {
        // find out theta
        const Real theta = getThetaFromArclength(arclength);
        
        // get the point
        const Real ptX = radiusX*std::cos(theta);
        const Real ptY = radiusY*std::sin(theta);
        
        return std::make_pair(ptX, ptY);
    }
};

#endif /* EllipticalCircumference_hpp */
