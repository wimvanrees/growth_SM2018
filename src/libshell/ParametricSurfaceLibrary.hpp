//
//  ParametricSurfaceLibrary.hpp
//  Elasticity
//
//  Created by Wim van Rees on 7/7/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef ParametricSurfaceLibrary_hpp
#define ParametricSurfaceLibrary_hpp

#include "common.hpp"
#include "JacobiEllipticFunctions.hpp"

struct ParametricSurface
{
    virtual Eigen::Vector3d operator()(const Real uu, const Real vv) const = 0;
    virtual Eigen::Vector3d get_xu(const Real uu, const Real vv) const = 0;
    virtual Eigen::Vector3d get_xv(const Real uu, const Real vv) const = 0;
    virtual Real getMeanCurvature(const Real uu, const Real vv) const = 0;
    virtual Real getGaussCurvature(const Real uu, const Real vv) const = 0;
    
    virtual std::pair<Real,Real> getExtent_U() const = 0;
    virtual std::pair<Real,Real> getExtent_V() const = 0;
    
    std::tuple<Real, Real, Real> get_FFF(const Real uu, const Real vv, const Eigen::Vector3d & e1, const Eigen::Vector3d & e2) const
    {
        const Eigen::Vector3d xu = this->get_xu(uu, vv);
        const Eigen::Vector3d xv = this->get_xv(uu, vv);
        
        // [e1x e1y] [xu.xu xu.xv; xu.xv xv.xv] [e1x e1y]
        // [e1x e1y] [e1x * xu.xu + e1y*xu.xv ; e1x * xu.xv + e1y*xv.xv]
        // e1x*e1x * xu.xu + 2 * e1x*e1y * xu.xv + e1h
        const Real xu_dot_xu = xu.dot(xu);
        const Real xu_dot_xv = xu.dot(xv);
        const Real xv_dot_xv = xv.dot(xv);
        
        const Real I_E = xu_dot_xu*e1(0)*e1(0) + 2*xu_dot_xv*e1(0)*e1(1) + xv_dot_xv*e1(1)*e1(1);
        const Real I_F = xu_dot_xu*e1(0)*e2(0) + xu_dot_xv*(e1(0)*e2(1) + e1(1)*e2(0)) + xv_dot_xv*e1(1)*e2(1);
        const Real I_G = xu_dot_xu*e2(0)*e2(0) + 2*xu_dot_xv*e2(0)*e2(1) + xv_dot_xv*e2(1)*e2(1);
        
        return std::make_tuple(I_E, I_F, I_G);
    }
    
    virtual Eigen::Vector3d getNormal(const Real uu, const Real vv) const
    {
        const Eigen::Vector3d xu = get_xu(uu,vv);
        const Eigen::Vector3d xv = get_xv(uu,vv);
        const Eigen::Vector3d xu_cross_xv = xu.cross(xv);
        return xu_cross_xv.normalized();
    }
};

struct ParametricCylinder : ParametricSurface
{
    const Real radius, length, theta;
    
    ParametricCylinder(const Real R, const Real L, const Real th = 2.0*M_PI):
    radius(R),
    length(L),
    theta(th)
    {}
    
    Eigen::Vector3d operator()(const Real uu, const Real vv) const override
    {
        Eigen::Vector3d retval;
        retval << vv, radius*std::cos(uu), radius*std::sin(uu);
        return retval;
    }
    
    Eigen::Vector3d get_xu(const Real uu, const Real ) const override
    {
        Eigen::Vector3d retval;
        retval << 0, -radius*std::sin(uu), radius*std::cos(uu);
        return retval;
    }
    
    Eigen::Vector3d get_xv(const Real , const Real ) const override
    {
        Eigen::Vector3d retval;
        retval << 1,0,0;
        return retval;
    }
    
    Real getMeanCurvature(const Real , const Real ) const override
    {
        return 0.5/radius;
    }
    
    Real getGaussCurvature(const Real, const Real ) const override
    {
        return 0.0;
    }
    
    std::pair<Real,Real> getExtent_U() const override
    {
        return std::make_pair(0, theta);
    }
    
    std::pair<Real,Real> getExtent_V() const override
    {
        return std::make_pair(0, length);
    }
};


struct ParametricIsoScaledSurface : ParametricSurface
{
    const Real scaleFac;
    
    ParametricIsoScaledSurface(const Real scaleFac):
    scaleFac(scaleFac)
    {}
    
    Eigen::Vector3d operator()(const Real uu, const Real vv) const override
    {
        Eigen::Vector3d retval;
        retval << scaleFac*uu, scaleFac*vv, 0;
        return retval;
    }
    
    Eigen::Vector3d get_xu(const Real , const Real ) const override
    {
        Eigen::Vector3d retval;
        retval << scaleFac, 0, 0;
        return retval;
    }
    
    Eigen::Vector3d get_xv(const Real , const Real ) const override
    {
        Eigen::Vector3d retval;
        retval << 0, scaleFac, 0;
        return retval;
    }
    
    Real getMeanCurvature(const Real , const Real ) const override
    {
        return 0.0;
    }
    
    Real getGaussCurvature(const Real, const Real ) const override
    {
        return 0.0;
    }
    
    std::pair<Real,Real> getExtent_U() const override
    {
        return std::make_pair(-1, 1);
    }
    
    std::pair<Real,Real> getExtent_V() const override
    {
        return std::make_pair(-1, 1);
    }
};


struct ParametricCone : ParametricSurface
{
    const Real baseRadius, height;
    
    ParametricCone(const Real R, const Real H):
    baseRadius(R),
    height(H)
    {}
    
    Eigen::Vector3d operator()(const Real uu, const Real vv) const override
    {
        Eigen::Vector3d retval;
        const Real rad = std::sqrt(uu*uu + vv*vv);
        const Real myH = height*(1.0 - rad/baseRadius);
        retval << uu, vv, myH;
        return retval;
    }
    
    Eigen::Vector3d get_xu(const Real uu, const Real vv) const override
    {
        Eigen::Vector3d retval;
        const Real rad = std::sqrt(uu*uu + vv*vv);
        retval << 1, 0, -height*uu/(baseRadius*rad);
        return retval;
    }
    
    Eigen::Vector3d get_xv(const Real uu, const Real vv) const override
    {
        Eigen::Vector3d retval;
        const Real rad = std::sqrt(uu*uu + vv*vv);
        retval << 0, 1, -height*vv/(baseRadius*rad);
        return retval;
    }
    
    Real getMeanCurvature(const Real , const Real ) const override
    {
        return -1; // not yet implemented
    }
    
    Real getGaussCurvature(const Real, const Real ) const override
    {
        return 0.0;
    }
    
    std::pair<Real,Real> getExtent_U() const override
    {
        return std::make_pair(-baseRadius, baseRadius);
    }
    
    std::pair<Real,Real> getExtent_V() const override
    {
        return std::make_pair(-baseRadius, baseRadius);
    }
};

struct ParametricOmega : ParametricSurface
{
    const Real height;
    const Real innerRadius, outerRadius_left, outerRadius_rght;
    const Real innerAngle, outerAngle_left, outerAngle_rght;
    
    Eigen::Vector2d origin_left, origin_rght;
    
    ParametricOmega(const Real height, const Real innerRadius, const Real outerRadius, const Real innerAngle, const Real outerAngle):
    height(height),
    innerRadius(innerRadius),
    outerRadius_left(outerRadius),
    outerRadius_rght(outerRadius),
    innerAngle(innerAngle),
    outerAngle_left(outerAngle),
    outerAngle_rght(outerAngle)
    {
        // middle origin is at (0,0)

        // compute the location of the connection point on the left
        const Eigen::Vector2d pt_left = (Eigen::Vector2d() << innerRadius*std::sin(-0.5*innerAngle), innerRadius*std::cos(-0.5*innerAngle)).finished();
        const Eigen::Vector2d pt_rght = (Eigen::Vector2d() << innerRadius*std::sin(+0.5*innerAngle), innerRadius*std::cos(+0.5*innerAngle)).finished();
        
        // compute the direction vectors between origin and these points
        const Eigen::Vector2d dir_left = pt_left.normalized();
        const Eigen::Vector2d dir_rght = pt_rght.normalized();
        
        // compute the origins of the outer circles : they fall along the same direction vector, just further out
        origin_left = pt_left + dir_left*outerRadius_left;
        origin_rght = pt_rght + dir_rght*outerRadius_rght;
    }
    
    Eigen::Vector3d getNormal(const Real uu, const Real vv) const override
    {
        // my own normal function so that it is consistently outwards
        const Eigen::Vector3d surface_pt = this->operator()(uu,vv);
        Real nX,nY;
        
        if(uu < outerAngle_left)
        {
            // left branch
            nX = origin_left(0) - surface_pt(0);
            nY = origin_left(1) - surface_pt(2);
        }
        else if(uu < outerAngle_left + innerAngle)
        {
            // middle branch
            nX = surface_pt(0);
            nY = surface_pt(2);
        }
        else
        {
            // right branch
            nX = origin_rght(0) - surface_pt(0);
            nY = origin_rght(1) - surface_pt(2);
        }
        Eigen::Vector3d retval;
        retval << nX,0,nY;
        return retval.normalized();
    }
    
    Eigen::Vector3d operator()(const Real uu, const Real vv) const override
    {
        Real myX, myY;
        if(uu < outerAngle_left)
        {
            // left branch
            const Real uu_rel = outerAngle_left - 0.5*innerAngle - uu; // between outerPhi - 0.5*innerPhi and -0.5*innerPhi
            myX = origin_left(0) - outerRadius_left*std::sin(uu_rel);
            myY = origin_left(1) - outerRadius_left*std::cos(uu_rel);
        }
        else if(uu < outerAngle_left + innerAngle)
        {
            // middle branch
            const Real uu_rel = uu - outerAngle_left - 0.5*innerAngle; // between -0.5*innerAngle and +0.5*innerAngle
            myX = innerRadius*std::sin(uu_rel);
            myY = innerRadius*std::cos(uu_rel);
        }
        else
        {
            // right branch
            const Real uu_rel = outerAngle_left + 1.5*innerAngle - uu; // between +0.5*innerAngle and +0.5*innerAngle - outerAngle_rght
            myX = origin_rght(0) - outerRadius_rght*std::sin(uu_rel);
            myY = origin_rght(1) - outerRadius_rght*std::cos(uu_rel);
        }
        
        Eigen::Vector3d retval;
        retval << myX,vv,myY;
        return retval;
    }
    
    Eigen::Vector3d get_xu(const Real uu, const Real ) const override
    {
        Real dux, duy;
        if(uu < outerAngle_left)
        {
            // left branch
            const Real uu_rel = outerAngle_left - 0.5*innerAngle - uu; // between outerPhi - 0.5*innerPhi and -0.5*innerPhi
            dux = -outerRadius_left*std::cos(uu_rel) * -1;
            duy = +outerRadius_left*std::sin(uu_rel) * -1;
        }
        else if(uu < outerAngle_left + outerAngle_rght)
        {
            // middle branch
            const Real uu_rel = uu - outerAngle_left - 0.5*innerAngle; // between -0.5*innerAngle and +0.5*innerAngle
            dux =  innerRadius*std::cos(uu_rel);
            duy = -innerRadius*std::sin(uu_rel);
        }
        else
        {
            // right branch
            const Real uu_rel = outerAngle_left + 1.5*innerAngle - uu; // between +0.5*innerAngle to +0.5*innerAngle - outerAngle_rght
            dux = -outerRadius_rght*std::cos(uu_rel) * -1;
            duy = +outerRadius_rght*std::sin(uu_rel) * -1;
        }
        
        Eigen::Vector3d retval;
        retval << dux,0.0,duy;
        return retval;
    }
    
    Eigen::Vector3d get_xv(const Real , const Real ) const override
    {
        Eigen::Vector3d retval;
        retval << 0,1,0;
        return retval;
    }
    
    Real getMeanCurvature(const Real , const Real ) const override
    {
        std::cout << "mean curvature not implemented for parametric omega\n";
        return -1;
    }
    
    Real getGaussCurvature(const Real, const Real ) const override
    {
        std::cout << "gauss curvature not implemented for parametric omega\n";
        return -1;
    }
    
    std::pair<Real,Real> getExtent_U() const override
    {
        return std::make_pair(0, outerAngle_left + innerAngle + outerAngle_rght);
    }
    
    std::pair<Real,Real> getExtent_V() const override
    {
        return std::make_pair(0, height);
    }
};

struct ParametricHelicoid : ParametricSurface
{
    const Real cval;
    const Real height;
    
    ParametricHelicoid(const Real cval_in, const Real height_in):
    cval(cval_in),
    height(height_in)
    {}
    
    Eigen::Vector3d operator()(const Real uu, const Real vv) const override
    {
        Eigen::Vector3d retval;
        retval << uu*std::cos(vv), uu*std::sin(vv), cval*vv;
        return retval;
    }
    
    Eigen::Vector3d get_xu(const Real , const Real vv) const override
    {
        Eigen::Vector3d retval;
        retval << std::cos(vv), std::sin(vv), 0.0;
        return retval;
    }
    
    Eigen::Vector3d get_xv(const Real uu, const Real vv) const override
    {
        Eigen::Vector3d retval;
        retval << -uu*std::sin(vv), uu*std::cos(vv), cval;
        return retval;
    }
    
    Real getMeanCurvature(const Real , const Real ) const override
    {
        return 0.0;
    }
    
    Real getGaussCurvature(const Real uu, const Real ) const override
    {
        return  -cval*cval/std::pow(cval*cval + uu*uu,2);
    }
    
    std::pair<Real,Real> getExtent_U() const override
    {
        return std::make_pair(-0.5*height/cval, 0.5*height/cval);
    }
    
    std::pair<Real,Real> getExtent_V() const override
    {
        return std::make_pair(0, 2*M_PI);
    }
    
};

struct ParametricCatenoid : ParametricSurface
{
    const Real cval;
    const Real height;
    
    ParametricCatenoid(const Real cval_in, const Real height_in):
    cval(cval_in), // inner radius = c
    height(height_in) // height determines outer radius: Rmax = c*cosh(h/(2*c))
    {}
    

    Eigen::Vector3d operator()(const Real uu, const Real vv) const override
    {
        Eigen::Vector3d retval;
        const Real tx = cval * std::cosh(vv/cval)*std::cos(uu);
        const Real ty = vv;
        const Real tz = cval * std::cosh(vv/cval)*std::sin(uu);
        retval << tx, ty, tz;
        return retval;
    }
    
    Eigen::Vector3d get_xu(const Real uu, const Real vv) const override
    {
        const Eigen::Vector3d xu = (Eigen::Vector3d() <<
                                    -cval*std::cosh(vv/cval)*std::sin(uu),
                                    0.0,
                                    cval*std::cosh(vv/cval)*std::cos(uu)
                                    ).finished();
        return xu;
    }
    
    Eigen::Vector3d get_xv(const Real uu, const Real vv) const override
    {
        const Eigen::Vector3d xv = (Eigen::Vector3d() <<
                                    std::cos(uu)*std::sinh(vv/cval),
                                    1.0,
                                    std::sin(uu)*std::sinh(vv/cval)
                                    ).finished();
        return xv;
    }
    
    Real getMeanCurvature(const Real , const Real ) const override
    {
        return 0.0;
    }
    
    Real getGaussCurvature(const Real , const Real vv) const override
    {
        return  -1.0/(cval*cval)*1.0/std::pow(std::cosh(vv/cval), 4.0);
    }
    
    std::pair<Real,Real> getExtent_U() const override
    {
        return std::make_pair(0, 2.0*M_PI);
    }
    
    std::pair<Real,Real> getExtent_V() const override
    {
        return std::make_pair(-0.5*height, 0.5*height);
    }
};

struct ParametricSphericalShell : ParametricSurface
{
    const Real radius;
    const Real halfOpeningAngle;
    const bool scaleBaseToOne;
    
    ParametricSphericalShell(const Real radius_in, const Real angle_in, const bool scaleBaseToOne_in = false):
    radius(radius_in),
    halfOpeningAngle(angle_in),
    scaleBaseToOne(scaleBaseToOne_in)
    {}
    
    Eigen::Vector3d operator()(const Real uu, const Real vv) const override
    {
        // stereographic projection
        const Real urel = uu/radius;
        const Real vrel = vv/radius;
        
        const Real scaleFac = scaleBaseToOne ? (1+radius*radius)/(2*radius) : 1.0; // so that xmax/xmin are within +1/-1
        
        const Real denum = 1 + urel*urel + vrel*vrel;
        const Real x = scaleFac * 2*urel/denum;
        const Real y = scaleFac * 2*vrel/denum;
        const Real z = scaleFac * ( -1 + urel*urel + vrel*vrel)/denum;
        
//        //const Real diskRadius = radius*std::sin(halfOpeningAngle); // radius of base disk
//        
//        const Real plane_height = radius * std::cos(halfOpeningAngle); // distance from plane to origin
//        
//        const Real inplane_radius = std::sqrt(uu*uu + vv*vv);
//        const Real distance_to_center = std::sqrt( std::pow(inplane_radius,2) + std::pow(plane_height,2) );
//        
//        const Real deltaR = radius  - distance_to_center;
//        
//        const Real x_over_z = uu / plane_height;
//        const Real y_over_z = vv / plane_height;
//        const Real deltaX = deltaR * x_over_z / (1.0 + std::pow(x_over_z,2)); // sin(atan(x/(Rcosalpha)))
//        const Real deltaY = deltaR * y_over_z / (1.0 + std::pow(y_over_z,2)); // sin(atan(y/(Rcosalpha)))
//        const Real deltaZ = std::sqrt(deltaR*deltaR - deltaX*deltaX - deltaY*deltaY);
//        
        Eigen::Vector3d retval;
//        retval << uu + deltaX, vv + deltaY, deltaZ - plane_height;
        retval << x,y,z;
        return retval;
    }
    
    Eigen::Vector3d get_xu(const Real uu, const Real vv) const override
    {
        const Real urel = uu/radius;
        const Real vrel = vv/radius;
        
        const Real scaleFac = scaleBaseToOne ? (1+radius*radius)/(2*radius) : 1.0; // so that xmax/xmin are within +1/-1
        
        const Real denum = std::pow(1.0 + urel*urel + vrel*vrel,2);
        
        const Real xu_x = scaleFac * 2*(1 - urel*urel + vrel*vrel) / denum;
        const Real xu_y = scaleFac * -4*urel*vrel / denum;
        const Real xu_z = scaleFac * 4*urel/ denum;
        
        const Eigen::Vector3d xu = (Eigen::Vector3d() << xu_x, xu_y, xu_z).finished();
        return xu;
    }
    
    Eigen::Vector3d get_xv(const Real uu, const Real vv) const override
    {
        const Real urel = uu/radius;
        const Real vrel = vv/radius;
        
        const Real scaleFac = scaleBaseToOne ? (1+radius*radius)/(2*radius) : 1.0; // so that xmax/xmin are within +1/-1
        
        const Real denum = std::pow(1.0 + urel*urel + vrel*vrel,3);
        
        const Real xv_x = scaleFac * -4*urel*vrel / denum;
        const Real xv_y = scaleFac * 2*(1 + urel*urel - vrel*vrel) / denum;
        const Real xv_z = scaleFac * 4*vrel/ denum;
        
        const Eigen::Vector3d xv = (Eigen::Vector3d() << xv_x, xv_y, xv_z).finished();
        return xv;
    }
    
    Real getMeanCurvature(const Real , const Real ) const override
    {
        return 1.0/radius;
    }
    
    Real getGaussCurvature(const Real , const Real ) const override
    {
        return  1.0/(radius*radius);
    }
    
    std::pair<Real,Real> getExtent_U() const override
    {
        return std::make_pair(-radius, radius);
    }
    
    std::pair<Real,Real> getExtent_V() const override
    {
        return std::make_pair(-radius, radius);
    }
};





struct ParametricPeirceHemisphere : ParametricSurface
{
    const Real radius;
    const Real Kk;
    JacobiEllipticFunctions jacobi;
    
    ParametricPeirceHemisphere(const Real radius_in):
    radius(radius_in),
    Kk(1.8540746773013719184)
    {}
    
    
    // some helper methods
    
    
    Eigen::Vector2d hemi_trafo(const Eigen::Vector2d uv) const
    {
        // extract the center diamond from a square with dimensions [0,2K] x [-K, K]
        // the new extents are K*[ 1 - sqrt(2)/2, 1 + sqrt(2)/2] x K*[-sqrt(2)/2, sqrt(2)/2]
        Eigen::Vector2d retval;
        retval <<
        0.5*std::sqrt(2.0)*( (uv(0) - Kk) + uv(1) ) + Kk,
        0.5*std::sqrt(2.0)*(-(uv(0) - Kk) + uv(1) );
        return retval;
    }
    
    Eigen::Vector2d inverse_peirce(const Eigen::Vector2d uv) const
    {
        // use the jacobi elliptic function to transform the coordinates of a square to polar coordinates of a disk (still defined on square)
        std::complex<Real> snc, cnc, dnc;
        jacobi.getJacobiComplex(uv(0), uv(1), 0.5, snc, cnc, dnc);
        
        // invert
        const Real rad = 2.0*std::atan(std::abs(cnc))/(0.5*M_PI); // divide by pi/2 to make the r coordinate go from 0 to 1
        const Real theta = std::arg(cnc);
        Eigen::Vector2d retval;
        retval << rad, theta;
        return retval;
    }
    
    Eigen::Vector2d polar2cart(const Eigen::Vector2d uv) const
    {
        // convert polar coordinates to cartesian coordinates
        Eigen::Vector2d retval;
        retval << uv(0)*std::cos(uv(1)), uv(0)*std::sin(uv(1));
        return retval;
    }
    
    Eigen::Vector3d forward_stereographic(const Eigen::Vector2d uv) const
    {
        // convert cartesian coordinates to a 3D surface using the stereographic projection
        const Real urel = uv(0);
        const Real vrel = uv(1);
        
        const Real denum = 1.0 + urel*urel + vrel*vrel;
        const Real x = 2*urel/denum;
        const Real y = 2*vrel/denum;
        const Real z = (urel*urel + vrel*vrel - 1)/denum;
        
        Eigen::Vector3d retval;
        retval << x,y,z;
        return retval;
    }
    
    Eigen::Vector3d operator()(const Real uu, const Real vv) const override
    {
        Eigen::Vector2d uv;
        uv << uu, vv;
        
        const Eigen::Vector2d uv_trafo = hemi_trafo(uv);
        const Eigen::Vector2d uv_disk_polar = inverse_peirce(uv_trafo);
        const Eigen::Vector2d uv_disk_cart = polar2cart(uv_disk_polar);
        const Eigen::Vector3d retval = forward_stereographic(uv_disk_cart);
        return retval;
    }
    
    // override the getNormals method directly so no need for the xu and xv methods
    Eigen::Vector3d getNormal(const Real uu, const Real vv) const override
    {
        Eigen::Vector2d uv;
        uv << uu, vv;
        const Eigen::Vector2d uv_trafo = hemi_trafo(uv);
        const Eigen::Vector2d uv_disk_polar = inverse_peirce(uv_trafo);
        const Eigen::Vector2d uv_disk_cart = polar2cart(uv_disk_polar);
        const Eigen::Vector3d XYZ = forward_stereographic(uv_disk_cart);
        // we are in the negative hemisphere
        return XYZ.normalized();
    }
    
    Eigen::Vector3d get_xu(const Real uu, const Real vv) const override
    {
//        helpers::catastrophe("should not need the derivative of the peirce hemisphere",__FILE__, __LINE__, false);
//        return Eigen::Vector3d();                
        
        // finite difference
        const Real eps = 1e-8;
        const Eigen::Vector3d val_p = this->operator()(uu + eps, vv);
        const Eigen::Vector3d val_m = this->operator()(uu - eps, vv);
        const Eigen::Vector3d retval = 0.5/eps * (val_p - val_m);
        return retval;
    }
    Eigen::Vector3d get_xv(const Real uu, const Real vv) const override
    {
//        helpers::catastrophe("should not need the derivative of the peirce hemisphere",__FILE__, __LINE__, false);
//        return Eigen::Vector3d();
        
        // finite difference
        const Real eps = 1e-8;
        const Eigen::Vector3d val_p = this->operator()(uu, vv + eps);
        const Eigen::Vector3d val_m = this->operator()(uu, vv - eps);
        const Eigen::Vector3d retval = 0.5/eps * (val_p - val_m);
        return retval;
    }
    
    Real getMeanCurvature(const Real , const Real ) const override
    {
        return 1.0/radius;
    }
    
    Real getGaussCurvature(const Real , const Real ) const override
    {
        return  1.0/(radius*radius);
    }
    
    std::pair<Real,Real> getExtent_U() const override
    {
        return std::make_pair(Kk - 0.5*std::sqrt(2.)*Kk, Kk + 0.5*std::sqrt(2.)*Kk);
    }
    
    std::pair<Real,Real> getExtent_V() const override
    {
        return std::make_pair(-0.5*std::sqrt(2.)*Kk, 0.5*std::sqrt(2.)*Kk);
    }
};


struct ParametricPeirceHemisphere_v2 : ParametricSurface
{
    const Real radius;
    
    ParametricPeirceHemisphere_v2(const Real radius_in):
    radius(radius_in)
    {}
    
    Eigen::Vector3d operator()(const Real uu, const Real vv) const override
    {
        const Real nu = 0.25*M_PI;
        
        // interpret uu, vv as longitude / latitude
        // http://geographiclib.sourceforge.net/1.41/jacobi.html
        
        const Real xx = radius * std::sin(uu) * std::sqrt( 1.0 - std::pow(std::sin(nu) * std::sin(vv),2) );
        const Real yy = radius * std::sin(vv) * std::sqrt( 1.0 - std::pow(std::cos(nu) * std::sin(uu),2) );
        const Real zz = radius * std::cos(uu) * std::cos(vv);
        const Eigen::Vector3d retval = (Eigen::Vector3d() << xx,yy,zz).finished();
        
        return retval;
    }
    
    // override the getNormals method directly so no need for the xu and xv methods
    Eigen::Vector3d getNormal(const Real uu, const Real vv) const override
    {
        const Eigen::Vector3d pos = this->operator()(uu,vv);
        return pos.normalized();
    }
    
    Eigen::Vector3d get_xu(const Real , const Real ) const override
    {
        helpers::catastrophe("should not need the derivative of the peirce hemisphere",__FILE__, __LINE__, false);
        return Eigen::Vector3d();
    }
    Eigen::Vector3d get_xv(const Real , const Real ) const override
    {
        helpers::catastrophe("should not need the derivative of the peirce hemisphere",__FILE__, __LINE__, false);
        return Eigen::Vector3d();
    }
    
    Real getMeanCurvature(const Real , const Real ) const override
    {
        return 1.0/radius;
    }
    
    Real getGaussCurvature(const Real , const Real ) const override
    {
        return  1.0/(radius*radius);
    }
    
    std::pair<Real,Real> getExtent_U() const override
    {
        return std::make_pair(-0.5*M_PI, 0.5*M_PI);
    }
    
    std::pair<Real,Real> getExtent_V() const override
    {
        return std::make_pair(-0.5*M_PI, 0.5*M_PI);
    }
};


struct ParametricEnneper : ParametricSurface
{
    const int nFac;
    const Real rExtent;
    
    ParametricEnneper(const int nFac, const Real rExtent=1):
    nFac(nFac),
    rExtent(rExtent)
    {}
    
    Eigen::Vector3d operator()(const Real uu, const Real vv) const override
    {
        // shortcuts
        const Real r = uu;
        const Real t = vv; // theta
        const int n = nFac;
        
        const Real r_pow_2n = std::pow(r,2*n);
        
        const Real x = r*std::cos(t) + (r_pow_2n*std::cos(t - 2*n*t))/(r - 2*n*r);
        const Real y = r*std::sin(t) + (r_pow_2n*std::sin(t - 2*n*t))/(r - 2*n*r);
        const Real z = 2*std::pow(r,n)*std::cos(n*t)/n;
     
        Eigen::Vector3d retval;
        retval << x,y,z;
        return retval;
    }
    
    Eigen::Vector3d get_xu(const Real uu, const Real vv) const override
    {
        // shortcuts
        const Real r = uu;
        const Real t = vv; // theta
        const int n = nFac;
        
        const Real dx = std::cos(t) - std::pow(r, 2*n-2)*std::cos(t - 2*n*t);
        const Real dy = std::sin(t) - std::pow(r, 2*n-2)*std::sin(t - 2*n*t);
        const Real dz = 2*std::pow(r, n-1)*std::cos(n*t);
        
        const Eigen::Vector3d xu = (Eigen::Vector3d() << dx, dy, dz ).finished();
        return xu;
    }
    
    Eigen::Vector3d get_xv(const Real uu, const Real vv) const override
    {
        // shortcuts
        const Real r = uu;
        const Real t = vv; // theta
        const int n = nFac;
        
        const Real dx = -r*std::sin(t) - std::pow(r, 2*n-1)*std::sin(t - 2*n*t);
        const Real dy =  r*std::cos(t) + std::pow(r, 2*n-1)*std::cos(t - 2*n*t);
        const Real dz = -2*std::pow(r, n)*std::sin(n*t);
        
        const Eigen::Vector3d xv = (Eigen::Vector3d() << dx, dy, dz ).finished();
        return xv;
    }
    
    Real getMeanCurvature(const Real , const Real ) const override
    {
        return 0.0; // minimal surface
    }
    
    Real getGaussCurvature(const Real , const Real ) const override
    {
        std::cout << "gaussian curvature for generalized enneper surfaces not yet implemented\n";
        return 0;
    }
    
    std::pair<Real,Real> getExtent_U() const override
    {
        return std::make_pair(0, rExtent);
    }
    
    std::pair<Real,Real> getExtent_V() const override
    {
        return std::make_pair(-M_PI, M_PI);
    }
};


struct ParametricSphere : ParametricSurface
{
    const Real radius;
    
    ParametricSphere(const Real radius_in):
    radius(radius_in)
    {}
    
    Eigen::Vector3d operator()(const Real uu, const Real vv) const override
    {
        Eigen::Vector3d retval;
        const Real tx = radius*std::cos(uu)*std::sin(vv);
        const Real ty = radius*std::sin(uu)*std::sin(vv);
        const Real tz = radius*std::cos(vv);
        retval << tx, ty, tz;
        return retval;
    }
    
    Eigen::Vector3d get_xu(const Real uu, const Real vv) const override
    {
        const Eigen::Vector3d xu = (Eigen::Vector3d() <<
                                    -radius*std::sin(uu)*std::sin(vv),
                                     radius*std::cos(uu)*std::sin(vv),
                                    0.0
                                    ).finished();
        return xu;
    }
    
    Eigen::Vector3d get_xv(const Real uu, const Real vv) const override
    {
        const Eigen::Vector3d xv = (Eigen::Vector3d() <<
                                    radius*std::cos(uu)*std::cos(vv),
                                    radius*std::sin(uu)*std::cos(vv),
                                    -radius*std::sin(vv)
                                    ).finished();
        return xv;
    }
    
    Real getMeanCurvature(const Real , const Real ) const override
    {
        return 1.0/radius;
    }
    
    Real getGaussCurvature(const Real , const Real ) const override
    {
        return  1.0/(radius*radius);
    }
    
    std::pair<Real,Real> getExtent_U() const override
    {
        return std::make_pair(0.0, 2.0*M_PI);
    }
    
    std::pair<Real,Real> getExtent_V() const override
    {
        return std::make_pair(0.0, M_PI);
    }
};



struct ParametricPseudoSphere : ParametricSurface
{
    const Real minU;
    const Real maxU;
    ParametricPseudoSphere(const Real minU_in, const Real maxU_in):
    minU(minU_in), // >= 0 (but if it's 0, the area at u=0 will go to zero also --> maybe make it 0.1 or sth)
    maxU(maxU_in) // should be bigger than zero, but it will rapidly result in crazy surface... something between 1 and 3 is reasonable
    {}
    
    Eigen::Vector3d operator()(const Real uu, const Real vv) const override
    {
        Eigen::Vector3d retval;
        const Real sech_u = 1.0/std::cosh(uu);
        
        const Real tx = sech_u * std::cos(vv);
        const Real ty = sech_u * std::sin(vv);
        const Real tz = uu - std::tanh(uu);
        retval << tx, ty, tz;
        return retval;
    }
    
    Eigen::Vector3d get_xu(const Real uu, const Real vv) const override
    {
        const Real sech_u = 1.0/std::cosh(uu);
        const Real tanh_u = std::tanh(uu);
        
        const Real tx_u = -std::cos(vv)*sech_u*tanh_u;
        const Real ty_u = -std::sin(vv)*sech_u*tanh_u;
        const Real tz_u = tanh_u*tanh_u;
        
        Eigen::Vector3d retval;
        retval << tx_u,ty_u,tz_u;
        return retval;
    }
    
    Eigen::Vector3d get_xv(const Real uu, const Real vv) const override
    {
        const Real sech_u = 1.0/std::cosh(uu);
        const Real tx_v = -std::sin(vv)*sech_u;
        const Real ty_v =  std::cos(vv)*sech_u;
        const Real tz_v = 0.0;
        
        Eigen::Vector3d retval;
        retval << tx_v, ty_v, tz_v;
        return retval;
    }
    
    Real getMeanCurvature(const Real uu, const Real ) const override
    {
        const Real csch_u = 1.0/std::sinh(uu);
        return 0.5*(std::sinh(uu) - csch_u);
    }
    
    Real getGaussCurvature(const Real , const Real ) const override
    {
        return -1;
    }
    
    std::pair<Real,Real> getExtent_U() const override
    {
        return std::make_pair(minU, maxU);
    }
    
    std::pair<Real,Real> getExtent_V() const override
    {
        return std::make_pair(0.0, 2.0*M_PI);
    }
};


struct ParametricSinusoidalSurface : ParametricSurface
{
    const Real amplitude;
    
    ParametricSinusoidalSurface(const Real amplitude):
    amplitude(amplitude)
    {}
    
    Eigen::Vector3d operator()(const Real uu, const Real vv) const override
    {
        const Real xx = std::sin(2.0*M_PI*uu);
        const Real yy = std::cos(4.0*M_PI*vv);
        const Real xx_yy = xx*yy;
        const Real expfac = std::exp(-(uu*uu + vv*vv)/0.5);
        
        Eigen::Vector3d retval;
        retval << uu, vv, amplitude*expfac*xx_yy;
        
        return retval;
    }
    
    Eigen::Vector3d get_xu(const Real uu, const Real vv) const override
    {
        const Real expFac = std::exp(-(uu*uu + vv*vv)/0.5);
        const Real dzdu = amplitude*2.0*expFac*std::cos(4.0*M_PI*vv)*(M_PI*std::cos(2.0*M_PI*uu) - 2.0*uu*std::sin(2.0*M_PI*uu));
        
        const Eigen::Vector3d xu = (Eigen::Vector3d() <<
                                    1.0,
                                    0.0,
                                    dzdu
                                    ).finished();
        return xu;
    }
    
    Eigen::Vector3d get_xv(const Real uu, const Real vv) const override
    {
        const Real expFac = std::exp(-(uu*uu + vv*vv)/0.5);
        const Real dzdv = -amplitude*4.0*expFac*std::sin(2.0*M_PI*uu)*(M_PI*std::sin(4.0*M_PI*vv) + vv*std::cos(4.0*M_PI*vv));
        
        const Eigen::Vector3d xv = (Eigen::Vector3d() <<
                                    0.0,
                                    1.0,
                                    dzdv
                                    ).finished();
        return xv;
    }
    
    Real getMeanCurvature(const Real , const Real ) const override
    {
        std::cout << "mean curvature not implemented for parametric sinusoidal surface\n";
        return -1;
    }
    
    Real getGaussCurvature(const Real , const Real ) const override
    {
        std::cout << "gauss curvature not implemented for parametric sinusoidal surface\n";
        return -1;
    }
    
    std::pair<Real,Real> getExtent_U() const override
    {
        return std::make_pair(-1.0, 1.0);
    }
    
    std::pair<Real,Real> getExtent_V() const override
    {
        return std::make_pair(-1.0, 1.0);
    }
};

struct ParametricFunnel : ParametricSurface
{
    const Real uBound_lower, uBound_upper, afac;
    
    ParametricFunnel(const Real lowerRadius, const Real upperRadius, const Real height):
    uBound_lower(lowerRadius),
    uBound_upper(upperRadius),
    afac(height / (std::log(upperRadius) - std::log(lowerRadius)))
    {}
    
    Eigen::Vector3d operator()(const Real uu, const Real vv) const override
    {
        const Real tx = uu*std::cos(vv);
        const Real ty = uu*std::sin(vv);;
        const Real tz = afac*std::log(uu);
        
        Eigen::Vector3d retval;
        retval << tx, ty, tz;
        return retval;
    }
    
    Eigen::Vector3d get_xu(const Real uu, const Real vv) const override
    {
        const Real xu_x = std::cos(vv);
        const Real xu_y = std::sin(vv);
        const Real xu_z = afac/uu;
        
        Eigen::Vector3d retval;
        retval << xu_x, xu_y, xu_z;
        return retval;
    }
    
    Eigen::Vector3d get_xv(const Real uu, const Real vv) const override
    {
        const Real xv_x = -uu*std::sin(vv);
        const Real xv_y =  uu*std::cos(vv);
        const Real xv_z =  0.0;
        
        Eigen::Vector3d retval;
        retval << xv_x, xv_y, xv_z;
        return retval;
    }
    
    Real getMeanCurvature(const Real uu, const Real ) const override
    {
        return std::pow(afac,3)/(2.0 * uu * std::pow(afac*afac + uu*uu, 1.5));
    }
    
    Real getGaussCurvature(const Real uu, const Real ) const override
    {
        return -std::pow(afac,2)/std::pow(afac*afac + uu*uu, 2);
    }
    
    std::pair<Real,Real> getExtent_U() const override
    {
        return std::make_pair(uBound_lower, uBound_upper);
    }
    
    std::pair<Real,Real> getExtent_V() const override
    {
        return std::make_pair(0, 2.0*M_PI);
    }
};


struct ParametricCallaLily : ParametricSurface
{
    Eigen::Vector3d operator()(const Real uu, const Real vv) const override
    {
       
        // ready to go
        const Real expFac = std::exp(4.0*uu/(8.0 + uu));
        
        
        const Real ty = -0.025*(20.0*expFac + M_PI + 2*vv)*std::cos(vv); // flip y by minus one so the normals for undeformed and deformed point in same direction
        const Real tz = 0.025*(20.0*expFac + M_PI + 2*vv)*std::sin(vv);
        const Real tx = 32./3 * (2 + uu)/(8 + uu);
        
        Eigen::Vector3d retval;
        retval << tx, ty, tz;
        return retval;
    }
    
    Eigen::Vector3d get_xu(const Real uu, const Real vv) const override
    {
        const Real expFac = std::exp(4.0*uu/(8.0 + uu));
        
        const Real xu_y = -(16*expFac*std::cos(vv))/std::pow(8 + uu,2);
        const Real xu_z = (16*expFac*std::sin(vv))/std::pow(8 + uu,2);
        const Real xu_x = 64/std::pow(8 + uu,2);
        
        Eigen::Vector3d retval;
        retval << xu_x, xu_y, xu_z;
        return retval;
    }
    
    Eigen::Vector3d get_xv(const Real uu, const Real vv) const override
    {
        
        const Real expFac = std::exp(4.0*uu/(8.0 + uu));
        
        const Real xv_y =  -0.025*(2.0*std::cos(vv) - (20.0*expFac + M_PI + 2*vv)*std::sin(vv));
        const Real xv_z =  0.025*(2.0*std::sin(vv) + (20.0*expFac + M_PI + 2*vv)*std::cos(vv));
        const Real xv_x =  0.0;
        
        Eigen::Vector3d retval;
        retval << xv_x, xv_y, xv_z;
        return retval;
    }
    
    Real getMeanCurvature(const Real uu, const Real ) const override
    {
        // WIM : not sure if this is correct still after taking into account the actual shape
        return -0.5 * 128*std::exp(-4*uu/(uu+8))/std::pow(std::exp(8*uu/(uu+8)) + 16, 1.5);
    }
    
    Real getGaussCurvature(const Real uu, const Real ) const override
    {
        // WIM : not sure if this is correct still after taking into account the actual shape
        return -64.0/std::pow(std::exp(8*uu/(uu+8)) + 16, 2);
    }
    
    std::pair<Real,Real> getExtent_U() const override
    {
        std::cout << "should not be here (getExtent_U) - custom mapping needed" << std::endl;
        return std::make_pair(-2.0, 184/25.); // actual mapping function is more complex : do not need this
    }
    
    std::pair<Real,Real> getExtent_V() const override
    {
        std::cout << "should not be here (getExtent_V) - custom mapping needed" << std::endl;
        return std::make_pair(0, 1);
    }
};




struct ParametricHyperbolicParaboloid : ParametricSurface
{
    const Real afac,bfac;
    
    ParametricHyperbolicParaboloid(const Real afac_in, const Real bfac_in):
    afac(afac_in),
    bfac(bfac_in)
    {}
    
    Eigen::Vector3d operator()(const Real uu, const Real vv) const override
    {
        Eigen::Vector3d retval;
        retval << uu, vv, std::pow(uu/afac,2) - std::pow(vv/bfac,2);
        return retval;
    }
    
    Eigen::Vector3d get_xu(const Real uu, const Real ) const override
    {
        const Eigen::Vector3d xu = (Eigen::Vector3d() <<
                                    1.0,
                                    0.0,
                                    2.0*uu/(afac*afac)
                                    ).finished();
        return xu;
    }
    
    Eigen::Vector3d get_xv(const Real , const Real vv) const override
    {
        const Eigen::Vector3d xv = (Eigen::Vector3d() <<
                                    0.0,
                                    1.0,
                                    -2.0*vv/(bfac*bfac)
                                    ).finished();
        return xv;
    }
    
    Real getMeanCurvature(const Real uu, const Real vv) const override
    {
        const Real ufac = 4.0*uu*uu/std::pow(afac,4);
        const Real vfac = 4.0*vv*vv/std::pow(bfac,4);

        const Real asq = std::pow(afac,2);
        const Real bsq = std::pow(bfac,2);
        
        const Real denum = asq*bsq*std::pow(1.0 + ufac + vfac, 3.0/2.0);
        const Real num = -asq + bsq - 4.0*uu*uu/asq + 4.0*vv*vv/bsq;
        
        return num/denum;
    }
    
    Real getGaussCurvature(const Real uu, const Real vv) const override
    {
        const Real ufac = 4.0*uu*uu/std::pow(afac,4);
        const Real vfac = 4.0*vv*vv/std::pow(bfac,4);
        
        const Real denum = std::pow(afac * bfac * (1.0 + ufac + vfac), 2);
        return -4.0/denum;
    }
    
    std::pair<Real,Real> getExtent_U() const override
    {
        return std::make_pair(-1.0, 1.0);
    }
    
    std::pair<Real,Real> getExtent_V() const override
    {
        return std::make_pair(-1.0, 1.0);
    }
};
#endif /* ParametricSurfaceLibrary_hpp */
