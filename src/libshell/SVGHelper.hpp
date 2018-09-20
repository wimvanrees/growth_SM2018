//
//  SVGHelper.h
//  Elasticity
//
//  Created by Wim van Rees on 8/8/17.
//  Copyright © 2017 Wim van Rees. All rights reserved.
//

#ifndef SVGHelper_hpp
#define SVGHelper_hpp

#include "common.hpp"

namespace SVGHelper
{
    class Segment
    {
    public:
        virtual Eigen::Vector2d eval(const Real t) const = 0;
        
        Real getLength(const int Nt = 1000) const
        {
            Real retval = 0.0;
            const Real dt = 1.0/(Nt-1.);
            Eigen::Vector2d lastPt = this->eval(0);
            
            for(int i=1;i<Nt;++i)
            {
                const Real t = i*dt;
                const Eigen::Vector2d pt = this->eval(t);
                retval += (pt - lastPt).norm();
                lastPt = pt;
            }
            return retval;
        }
        
        virtual ~Segment() = default;
    };

    class Line : public Segment
    {
        const Eigen::Vector2d p0, p1;
        
    public:
        Line(const Eigen::Vector2d & p0_in, const Eigen::Vector2d & p1_in):
        p0(p0_in),
        p1(p1_in)
        {}
        
        Eigen::Vector2d eval(const Real t) const override
        {
            assert(t >= 0);
            assert(t <= 1);
            
            return (1.0 - t)*p0 + t*p1;
        }
    };


    class QuadraticBezier : public Segment
    {
        const Eigen::Vector2d P0, P1, P2;
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        QuadraticBezier(const Eigen::Vector2d & P0_in, const Eigen::Vector2d & P1_in, const Eigen::Vector2d & P2_in):
        P0(P0_in),
        P1(P1_in),
        P2(P2_in)
        {}
        
        Eigen::Vector2d eval(const Real t) const override
        {
            assert(t >= 0);
            assert(t <= 1);
            
            const Real one_minus_t = 1.0 - t;
            return std::pow(one_minus_t, 2)*P0 + 2*one_minus_t*t*P1 + std::pow(t,2)*P2;
        }
    };


    class CubicBezier : public Segment
    {
        const Eigen::Vector2d P0, P1, P2, P3;
        
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        CubicBezier(const Eigen::Vector2d & P0_in, const Eigen::Vector2d & P1_in, const Eigen::Vector2d & P2_in, const Eigen::Vector2d & P3_in):
        P0(P0_in),
        P1(P1_in),
        P2(P2_in),
        P3(P3_in)
        {}
        
        Eigen::Vector2d eval(const Real t) const override
        {
            assert(t >= 0);
            assert(t <= 1);
            
            const Real one_minus_t = 1.0 - t;
            
            return std::pow(one_minus_t, 3)*P0 + 3*std::pow(one_minus_t,2)*t*P1 + 3*one_minus_t*std::pow(t,2)*P2 + std::pow(t,3)*P3;
        }
    };

    class EllipticalArc : public Segment
    {
        // http://ericeastwood.com/blog/25/curves-and-arcs-quadratic-cubic-elliptical-svg-implementations
        // https://www.w3.org/TR/SVG/implnote.html#ArcConversionEndpointToCenter
        
        Eigen::Vector2d P0, P1;
        Eigen::Matrix2d rotMat;
        Eigen::Vector2d rxy;
        Eigen::Vector2d cxy;
        Real theta1, dtheta;
        
        Real computeAngle(const Eigen::Vector2d & u, const Eigen::Vector2d & v) const
        {
            const Real sign = (u(0)*v(1) - u(1)*v(0)) <  0 ? -1 : +1;
            const Real angle = std::acos(u.dot(v)/(u.norm() * v.norm()));
            return sign*angle;
        }
        
        void initializeFromEndpointParametrization(const Real xAxisRotation, const Eigen::Vector2d & rxy_raw, const int large_arc, const int sweep)
        {
            // Following "Conversion from endpoint to center parameterization"
            // http://www.w3.org/TR/SVG/implnote.html#ArcConversionEndpointToCenter
            
            // compute the rotation matrix
            const Real crot = std::cos(xAxisRotation);
            const Real srot = std::sin(xAxisRotation);
            this->rotMat = (Eigen::Matrix2d() <<  crot, srot, -srot, crot).finished();
            
            // compute the correct ellipse radii
            const Eigen::Vector2d rxy_abs = rxy_raw.cwiseAbs();
            const Eigen::Vector2d x1p = 0.5 * rotMat * (P0 - P1);
            
            const Real Lambda = std::pow(x1p(0)/rxy_abs(0),2) + std::pow(x1p(1)/rxy_abs(1),2);
            this->rxy = rxy_abs;
            if(Lambda > 1) this->rxy *= std::sqrt(Lambda);
            
            // compute center
            const Real rx_sq = rxy(0)*rxy(0);
            const Real ry_sq = rxy(1)*rxy(1);
            const Real x1p_sq = x1p(0)*x1p(0);
            const Real y1p_sq = x1p(1)*x1p(1);
            
            const Real fac1_num = rx_sq*ry_sq - rx_sq*y1p_sq - ry_sq*x1p_sq;
            const Real fac2_num = rx_sq*y1p_sq + ry_sq*x1p_sq;
            const Real signfac = (large_arc == sweep ? -1 : +1);
            const Real prefac = signfac*std::sqrt(fac1_num/fac2_num);
            const Eigen::Vector2d cxyp = prefac * (Eigen::Vector2d() << rxy(0)*x1p(1)/rxy(1), -rxy(1)*x1p(0)/rxy(0)).finished();
            this->cxy = rotMat.transpose() * cxyp + 0.5*(P0 + P1);
            
            // compute the angles
            
            const Eigen::Vector2d theta_vec_1 = (Eigen::Vector2d() << 1,0).finished();
            const Eigen::Vector2d theta_vec_2 = (Eigen::Vector2d() << (x1p(0) - cxyp(0))/rxy(0), (x1p(1) - cxyp(1))/rxy(1)).finished();
            const Eigen::Vector2d theta_vec_3 = (Eigen::Vector2d() << (-x1p(0) - cxyp(0))/rxy(0), (-x1p(1) - cxyp(1))/rxy(1)).finished();
            
            this->theta1 = computeAngle(theta_vec_1, theta_vec_2);
            const Real dtheta_raw = computeAngle(theta_vec_2, theta_vec_3);
            this->dtheta = (dtheta_raw > 0 && sweep == 0) ? dtheta_raw - 2*M_PI : ((dtheta_raw < 0 && sweep == 1) ? dtheta_raw + 2*M_PI : dtheta_raw);
        }
        
        
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        EllipticalArc(const Eigen::Vector2d & P0_in, const Eigen::Vector2d & P1_in, const Real xAxisRotation_in, const Eigen::Vector2d & rxy_in, const int large_arc_in, const int sweep_in):
        P0(P0_in),
        P1(P1_in)
        {
            initializeFromEndpointParametrization(xAxisRotation_in, rxy_in, large_arc_in, sweep_in);
        }
        
        Eigen::Vector2d eval(const Real t) const override
        {
            assert(t >= 0);
            assert(t <= 1);
            
            // If the endpoints are identical, then this is equivalent to omitting the elliptical arc segment entirely.
            if( P0.isApprox(P1) ) return P0;
            
            // If rx = 0 or ry = 0 then this arc is treated as a straight line segment joining the endpoints.
            if( rxy.isApprox(Eigen::Vector2d::Zero()) ) return (1.0 - t)*P0 + t*P1;
            
            // else compute the value
            const Real angle = theta1 + t*dtheta;
            const Eigen::Vector2d ell_xy = (Eigen::Vector2d() << rxy(0)*std::cos(angle), rxy(1)*std::sin(angle)).finished();
            
            const Eigen::Vector2d pt = rotMat.transpose() * ell_xy + cxy;
            return pt;
        }
        
    };

    class Test_EllipticalArc
    {
        
    public:
        void test()
        {
            Eigen::Vector2d start,end,rxy;
            start<<0,0;
            end<<1,0;
            rxy<<0.5, 0.75;
            
            const Real rotation = 0.25*M_PI;
            const int large_arc = 1;
            const int sweep = 0;
            
            EllipticalArc arc(start, end, rotation, rxy, large_arc, sweep);
            
            const Real dt = 0.1;
            Real t = 0.0;
            while(t <= 1.0)
            {
                const Eigen::Vector2d pt = arc.eval(t);
                printf("%10.10e : \t %10.10e, %10.10e\n", t, pt(0), pt(1));
                t+=dt;
            }
        }
    };


    class Path
    {
        std::vector<std::pair<std::shared_ptr<Segment>, Real>> segments;
        
    public:
        
        void addCurve(const std::shared_ptr<Segment> & segment, const Real endPath)
        {
            segments.push_back(std::make_pair(segment, endPath));
        }
        
        Eigen::Vector2d eval(const Real T) const
        {
            assert(T>=0);
            assert(T<=1);
            
            int idx = 0;
            Real firstT = 0.0;
            for(const auto & seg : segments)
            {
                if(T <= seg.second) break; // need <= to also keep track of 1
                firstT = seg.second;
                idx++;
            }
            assert(idx < (int)segments.size());
            
            const Real lastT = segments[idx].second;
            
            assert(T >= firstT);
            assert(T <= lastT);
            
            const Real t = (T - firstT)/(lastT - firstT);
            return segments[idx].first->eval(t);
        }
        
        
        Real getLength(const int Nt = 1000) const
        {
            Real retval = 0.0;
            for(const auto & seg : segments)
            {
                retval += seg.first->getLength(Nt);
            }
            return retval;
        }
        
        int getNumberOfSegments() const
        {
            return (int)segments.size();
        }
    };

    struct SVGParser
    {
        static std::vector<Path> parseFile(const std::string fname)
        {
            std::vector<Path> paths;
            
            std::ifstream svgdata(fname);
            if(svgdata.fail())
            {
                helpers::catastrophe("FILE "+fname+" DOES NOT EXIST OR OPEN CORRECTLY", __FILE__, __LINE__);
            }
                                     
            int nPaths;
            svgdata >> nPaths;
            printf("Found %d paths\n\n", nPaths);
            
            for(int i=0;i<nPaths;++i)
            {
                Path path;
                
                int nSegs;
                svgdata >> nSegs;
                printf("found %d segments for path %d\n", nSegs, i);
                
                Real last_tEnd = -1;
                _unused(last_tEnd); // only used in debug mode
                for(int s=0;s<nSegs;++s)
                {
                    Real tStart, tEnd;
                    svgdata >> tStart >> tEnd;
                    
                    if(s==0) assert(std::abs(tStart) < 1e-12);
                    if(s==nSegs-1) assert(std::abs(tEnd - 1) < 1e-12);
                    if(s > 0) assert(std::abs(tStart - last_tEnd) < 1e-12);
                    
                    last_tEnd = tEnd;
                    
                    const bool zeroLength = (std::abs(tStart - tEnd) < 1e-12);

                    int segmentType;
                    svgdata >> segmentType;

                    if(zeroLength) printf("Found a segment with zero length : path %d, seg %d, tStart %10.10e, tEnd %10.10e, type %d\n", i, s, tStart, tEnd, segmentType);
                                            
                    //printf("\tsegment %d is of type %d, and runs from %10.10e to %10.10e %s\n", s, segmentType,tStart,tEnd, (zeroLength ? "-- zero length detected!" : ""));
                    
                    if(segmentType==1) // linear
                    {
                        Eigen::Vector2d P0, P1;
                        svgdata >> P0(0) >> P0(1) >> P1(0) >> P1(1);
                        Line line(P0, P1);
                        if(not zeroLength) path.addCurve(std::make_shared<Line>(line), tEnd);
                    }
                    else if(segmentType==2) // quadratic
                    {
                        Eigen::Vector2d P0, P1, P2;
                        svgdata >> P0(0) >> P0(1) >> P1(0) >> P1(1) >> P2(0) >> P2(1);
                        QuadraticBezier quad(P0, P1, P2);
                        if(not zeroLength) path.addCurve(std::make_shared<QuadraticBezier>(quad), tEnd);
                    }
                    else if(segmentType==3) // cubic
                    {
                        Eigen::Vector2d P0, P1, P2, P3;
                        svgdata >> P0(0) >> P0(1) >> P1(0) >> P1(1) >> P2(0) >> P2(1) >> P3(0) >> P3(1);
                        CubicBezier cube(P0, P1, P2, P3);
                        if(not zeroLength) path.addCurve(std::make_shared<CubicBezier>(cube), tEnd);
                    }
                    else if(segmentType==4) //arc
                    {
                        Eigen::Vector2d P0, P1;
                        Eigen::Vector2d rxy;
                        Real rotation;
                        int large_arc, sweep;
                        svgdata >> P0(0) >> P0(1) >> P1(0) >> P1(1) >> rxy(0) >> rxy(1) >> rotation >> large_arc >> sweep;
                        
                        EllipticalArc arc(P0, P1, rotation, rxy, large_arc, sweep);
                        if(not zeroLength) path.addCurve(std::make_shared<EllipticalArc>(arc), tEnd);
                    }
                }
                
                paths.push_back(path);
                printf("\t(path %d has length %10.10e)\n", i, path.getLength());
            }
            
            return paths;
        }
    };



    class DiscretePath
    {
        Eigen::MatrixXd points;
        std::array<Real, 4> extents;
        
        void fillPoints(const Path & path, const Real dpath)
        {
            const Real length = path.getLength();
            const int Nt = std::max(10, (int)(length/dpath));
            const Real dt = 1.0/(Nt - 1.0);

            points.resize(Nt, 2);
            //#pragma omp parallel for if(Nt > 1000) schedule(guided)
            for(int i=0;i<Nt;++i)
            {
                points.row(i) = path.eval(i*dt);
            }
        }
        
        void computeExtents()
        {
            extents = {1e9, -1e9, 1e9, -1e9};
            
            const int nPoints = getNumberOfPoints();
            
            for(int k=0;k<nPoints;++k)
            {
                for(int d=0;d<2;++d)
                {
                    extents[2*d+0] = std::min(extents[2*d+0], points(k,d));
                    extents[2*d+1] = std::max(extents[2*d+1], points(k,d));
                }
            }
        }
        
    public:
        DiscretePath(const Path & path, const Real dpath)
        {
            fillPoints(path, dpath);
            computeExtents();
        }
        
        void clipPoints(const std::array<Real, 2> xBounds, const std::array<Real, 2> yBounds)
        {
            // copy the old points
            const Eigen::MatrixXd points_all = points;
         
            const int Nt = getNumberOfPoints();
            Eigen::VectorXb isPointInside(Nt);

            // count the number of points after clipping
            int cnt = 0;
            for(int i=0;i<Nt;++i)
            {
                bool isInside = true;
                isInside = (isInside and (points(i,0) >= xBounds[0]));
                isInside = (isInside and (points(i,0) <= xBounds[1]));
                isInside = (isInside and (points(i,1) >= yBounds[0]));
                isInside = (isInside and (points(i,1) <= yBounds[1]));
                
                isPointInside(i) = isInside;
                
                if(isInside) cnt++;
            }

            points.resize(cnt, 2);
            points.setZero();
            
            if(cnt == 0) return; // this path needs to be dropped
            
            cnt = 0;
            for(int i=0;i<Nt;++i)
            {
                if(isPointInside(i))
                {
                    for(int d=0;d<2;++d) points(cnt, d) = points_all(i,d);
                    cnt++;
                }
            }
            
            // recompute extents
            computeExtents();
        }
        
        Eigen::Vector2d eval(const int t_idx) const
        {
            return points.row(t_idx);
        }
        
        const Eigen::Ref<const Eigen::MatrixXd> getPoints() const
        {
            return points;
        }
        
        Real getLength() const
        {
            Real retval = 0.0;
            const int nPoints = getNumberOfPoints();
            
            for(int i=0;i<nPoints-1;++i)
            {
                retval += (points.row(i) - points.row(i+1)).norm();
            }
            return retval;
        }
        
        void translatePoints(const Real deltaX, const Real deltaY)
        {
            const int nPoints = getNumberOfPoints();
            for(int i=0;i<nPoints;++i)
            {
                points(i,0) += deltaX;
                points(i,1) += deltaY;
            }
            
            // recompute extents
            computeExtents();
        }
        
        void rescalePoints(const Real scaleFac)
        {
            const int nPoints = getNumberOfPoints();
            for(int i=0;i<nPoints;++i)
                for(int d=0;d<2;++d)
                    points(i,d) *= scaleFac;
            
            // recompute extents
            computeExtents();
        }
        
        int getNumberOfPoints() const
        {
            return points.rows();
        }
        
        const std::array<Real, 4> & getExtents() const
        {
            return extents;
        }
    };


    class DiscreteSVG
    {
        const Real dpath;
        std::vector<DiscretePath> paths;
        
    public:
        DiscreteSVG(const std::string filename, const Real dpath_in):
        dpath(dpath_in)
        {
            const std::vector<Path> cont_paths = SVGParser::parseFile(filename);
            const int nPaths = (int)cont_paths.size();
            
            paths.reserve(nPaths);
            
            for(int i=0;i<nPaths;++i)
            {
                paths.emplace_back(cont_paths[i], dpath);
                //printf("path %d has continuous length %10.10e, and discrete length %10.10e, and %d points. Cont path first/last = %10.10e, %10.10e \t %10.10e, %10.10e \t Discr path first/last = %10.10e, %10.10e \t %10.10e, %10.10e\n", i, cont_paths[i].getLength(), paths[i].getLength(), paths[i].getNumberOfPoints(), cont_paths[i].eval(0)(0), cont_paths[i].eval(0)(1), cont_paths[i].eval(1)(0), cont_paths[i].eval(1)(1), paths[i].eval(0)(0), paths[i].eval(0)(1), paths[i].eval(paths[i].getNumberOfPoints()-1)(0), paths[i].eval(paths[i].getNumberOfPoints()-1)(1));
            }
            std::cout << "done with creating discrete paths" << std::endl;
        }
        
        std::array<Real, 4> getGlobalExtents() const
        {
            std::array<Real, 4> global_extents = {1e9, -1e9, 1e9, -1e9};
            for(const auto & path : paths)
            {
                const std::array<Real, 4> & path_extents = path.getExtents();
                for(int d=0;d<2;++d)
                {
                    global_extents[2*d+0] = std::min(global_extents[2*d+0], path_extents[2*d+0]);
                    global_extents[2*d+1] = std::max(global_extents[2*d+1], path_extents[2*d+1]);
                }
            }
            return global_extents;
        }
        
        
        void translateRescalePaths(const Real originX, const Real originY, const Real targetSizeX, const Real targetSizeY)
        {
            // first we translate
            {
                std::array<Real, 4> global_extents = getGlobalExtents();
                const Real sizeX = global_extents[1] - global_extents[0];
                const Real sizeY = global_extents[3] - global_extents[2];
                
                const Real deltaX = -global_extents[0] - 0.5*sizeX + originX;
                const Real deltaY = -global_extents[2] - 0.5*sizeY + originY;
                
                for(auto & path : paths)
                    path.translatePoints(deltaX, deltaY);
            }

            // now want to compute one rescaling factor for the entire grid - assume origin is correct
            {
                if(targetSizeX < 0 and targetSizeY < 0) return;
                if(targetSizeX > 0 and targetSizeY > 0)
                {
                    std::cout << "can not rescale paths when both scaling factors are positive - only one global positive scaling factor is allowed " << std::endl;
                    return;
                }
                
                std::array<Real, 4> global_extents = getGlobalExtents();
                const Real scaleFac = (targetSizeY < 0 ? (targetSizeX / (global_extents[1] - global_extents[0])) : (targetSizeY / (global_extents[3] - global_extents[2])));
                
                for(auto & path : paths)
                    path.rescalePoints(scaleFac);
            }
        }
        
        
        void clipAllPoints(const std::array<Real, 2> xBounds, const std::array<Real, 2> yBounds)
        {
            for(auto & path : paths)
                path.clipPoints(xBounds, yBounds);
            
            // remove empty paths from vector
            const int nPaths_0 = (int)paths.size();
            paths.erase(std::remove_if(paths.begin(), paths.end(),[](const DiscretePath& path) { return path.getNumberOfPoints()==0;}), paths.end());
            const int nPaths_1 = (int)paths.size();
            if(nPaths_1 != nPaths_0)
                printf("Removed %d elements from path vector: resized it from %d to %d elements\n", nPaths_0 - nPaths_1, nPaths_0, nPaths_1);
        }
        
        std::pair<Real, Eigen::Vector2d> getInfoForPoint_Filament(const Eigen::Vector2d & xc, const Real sigma, const Real dropSize = 0.5, const Real dropFac = 1.0) const
        {
            // dropsize is a prefactor for sigma (newsigma = dropsize * sigma)
            // dropfac is a prefactor for the weight
            Eigen::Vector2d retdir;
            retdir.setZero();
            int retdir_ref_path_idx = -1;
            
            Real retval = 0.0;
            
            const Real maxSigma = (dropSize > 1 ? dropSize : 1.0)*sigma;
            
            const int nPaths = (int)paths.size();
            for(int i=0;i<nPaths;++i)
            {
                const DiscretePath & path = paths[i];
                
                const int Nt = path.getNumberOfPoints();
                const Eigen::Ref<const Eigen::MatrixXd> points = path.getPoints();
                const std::array<Real, 4> & extents = path.getExtents();
                
                const bool x_out_of_range = (xc(0) < extents[0] - 1.05*maxSigma) or (xc(0) > extents[1] + 1.05*maxSigma);
                const bool y_out_of_range = (xc(1) < extents[2] - 1.05*maxSigma) or (xc(1) > extents[3] + 1.05*maxSigma);
                
                if(x_out_of_range or y_out_of_range) continue;
                
                const Real dt = 1.0/(Nt-1.);
                for(int tt = 0; tt<Nt; ++tt)
                {
                    const Eigen::Vector2d xp = points.row(tt);
                    
                    const Eigen::Vector2d tang_dir_raw =
                    (tt==0 ? ((Eigen::Vector2d)points.row(tt+1) - xp)/dt :
                     (tt==(Nt-1) ? (xp - (Eigen::Vector2d)points.row(tt-1))/dt :
                      ((Eigen::Vector2d)points.row(tt+1) - (Eigen::Vector2d)points.row(tt-1))/(2.0*dt)));
                    const Real tang_dir_length = tang_dir_raw.norm();
                    
                    const bool endPoint = ( (tt == 0) or (tt == Nt-1) ) ? true : false;
                    const Real prefac = (endPoint ? dropFac : 1.0);
                    const Real mySigma = (endPoint ? dropSize : 1.0)*sigma;
                    
                    Real weight = 0;
                    if(tang_dir_length > std::numeric_limits<Real>::epsilon())
                    {
                        const Eigen::Vector2d tang_dir = tang_dir_raw / tang_dir_length;
                        const Eigen::Vector2d norm_dir = (Eigen::Vector2d() << tang_dir(1), -tang_dir(0)).finished();
                        const Real reldist_t = std::abs((xc - xp).dot(tang_dir)) / mySigma;
                        const Real reldist_n = std::abs((xc - xp).dot(norm_dir)) / mySigma;
                        const Real weight_t = prefac*std::max(0.0, 1.0 - reldist_t);
                        const Real weight_n = prefac*std::max(0.0, 1.0 - reldist_n);
                        weight = weight_t * weight_n;
                    }
                    else
                    {
                        const Real reldist_x = std::abs(xc(0) - xp(0))/mySigma;
                        const Real reldist_y = std::abs(xc(1) - xp(1))/mySigma;
                        const Real weight_x = prefac*std::max(0.0, 1.0 - reldist_x);
                        const Real weight_y = prefac*std::max(0.0, 1.0 - reldist_y);
                        weight = weight_x * weight_y;
                    }
                    
                    Real dpath_me; // actual length of my path (because we are likely not arclength parametrized)
                    if(tt==0)
                    {
                        dpath_me = 0.5*((Eigen::Vector2d)points.row(1) - (Eigen::Vector2d)points.row(0)).norm();
                    }
                    else if(tt==(Nt-1))
                    {
                        dpath_me = 0.5*((Eigen::Vector2d)points.row(Nt-1) - (Eigen::Vector2d)points.row(Nt-2)).norm();
                    }
                    else
                    {
                        dpath_me = 0.5*((Eigen::Vector2d)points.row(tt+1) - (Eigen::Vector2d)points.row(tt)).norm() + 0.5*((Eigen::Vector2d)points.row(tt) - (Eigen::Vector2d)points.row(tt-1)).norm();
                    }
                    
                    retval += weight * dpath_me / mySigma;
                    
                    // NOW THE DIRECTION
                    

                    // since our print path is agnostic to the orientation of tangent dir we should make a consistent choice
                    // choose the version that is closest to the current guy
                    
                    Eigen::Vector2d tang_dir;
                    tang_dir.setZero();
                    
                    if(tang_dir_length > std::numeric_limits<Real>::epsilon())
                    {
                        tang_dir = tang_dir_raw / tang_dir_length;
                        
                        if( (retdir_ref_path_idx < 0) or (retdir_ref_path_idx == i) )
                        {
                            // if we haven't set a path yet, take this path as reference (and ignore the sign complications for all other points on this path)
                            retdir_ref_path_idx = i;
                        }
                        else
                        {
                            // need to make a consistent choice between interpolating + and - tang dir (if retdir already > 0)
                            const Real theta_plus = std::atan2(tang_dir(1), tang_dir(0));
                            const Real theta_mins = std::atan2(-tang_dir(1), -tang_dir(0));
                            const Real theta_dir = std::atan2(retdir(1), retdir(0));
                            
                            Real angle_diff_plus = fmod(std::abs(theta_plus - theta_dir) + 2*M_PI, 2*M_PI);
                            angle_diff_plus = std::min(angle_diff_plus, std::abs(angle_diff_plus - 2*M_PI));
                            Real angle_diff_mins = fmod(std::abs(theta_mins - theta_dir) + 2*M_PI, 2*M_PI);
                            angle_diff_mins = std::min(angle_diff_mins, std::abs(angle_diff_mins - 2*M_PI));
                            
                            const bool flipSign = (angle_diff_mins < angle_diff_plus);
                            
                            tang_dir *= (flipSign ? -1 : +1);
                        }
                        
                        retdir += weight * tang_dir * tang_dir_length / mySigma * dt;
                    }
                }
            }
            
            return std::make_pair(retval, retdir.normalized());
        }
        
        
        
        Real getRhoForPoint(const Eigen::Vector2d & xc, const Real sigma) const
        {
            Real retval = 0.0;
            
            const int nPaths = (int)paths.size();
            for(int i=0;i<nPaths;++i)
            {
                const DiscretePath & path = paths[i];
                
                const int Nt = path.getNumberOfPoints();
                const Eigen::Ref<const Eigen::MatrixXd> points = path.getPoints();
                const std::array<Real, 4> & extents = path.getExtents();
                
                const bool x_out_of_range = (xc(0) < extents[0] - 1.05*sigma) or (xc(0) > extents[1] + 1.05*sigma);
                const bool y_out_of_range = (xc(1) < extents[2] - 1.05*sigma) or (xc(1) > extents[3] + 1.05*sigma);
                
                if(x_out_of_range or y_out_of_range) continue;
                
                for(int tt = 0; tt<Nt; ++tt)
                {
                    const Eigen::Vector2d xp = points.row(tt);
                    const Real reldist = (xc - xp).norm() / sigma;
                    retval += std::max(0.0, 1 - reldist);
                }
            }
            
            return retval;
        }
    };
}

#endif /* SVGHelper_h */
