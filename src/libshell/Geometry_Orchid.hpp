//
//  Geometry_Orchid.hpp
//  Elasticity
//
//  Created by Wim van Rees on 8/23/17.
//  Copyright © 2017 Wim van Rees. All rights reserved.
//

#ifndef Geometry_Orchid_hpp
#define Geometry_Orchid_hpp

#include "Geometry.hpp"
#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <functional>

#include "EllipticalCircumference.hpp"


class Geometry_Flower_Base : public PlateGeometry
{
protected:
    typedef std::array<Real, 2> vec2;
    
    //**** SOME HELPER FUNCTIONS ****//
    
    inline std::pair<int, Real> getNumberOfPoints(const Real length, const Real delta_target, const bool include_end = false) const
    {
        const int nPoints = (int)std::ceil(length / delta_target);
        const Real delta = (include_end ? length / (nPoints - 1.) : length / nPoints);
        
        return std::make_pair(nPoints, delta);
    }
    
    inline Real computeLength(const vec2 & pt_1, const vec2 & pt_2) const
    {
        return std::sqrt(std::pow(pt_1[0] - pt_2[0],2) + std::pow(pt_1[1] - pt_2[1],2));
    }
    
    inline void translate(std::vector<vec2> & points, const vec2 deltapos) const
    {
        for(auto & point : points)
        {
            point[0] += deltapos[0];
            point[1] += deltapos[1];
        }
    }
    
    inline void rotate(std::vector<vec2> & points, const Real theta) const
    {
        const int nPoints = (int)points.size();
        const Real ctheta = std::cos(theta);
        const Real stheta = std::sin(theta);
        for(int i=0;i<nPoints;++i)
        {
            const Real xOld = points[i][0];
            const Real yOld = points[i][1];
            points[i][0] = ctheta*xOld + stheta*yOld;
            points[i][1] =-stheta*xOld + ctheta*yOld;
        }
    }
    
    //**** SOME HELPER FUNCTIONS ****//
    
    std::vector<vec2> connectWithStraightLine(const vec2 & start_pt_raw, const vec2 & end_pt, const Real delta, const bool include_start = true, const bool include_end = true) const
    {
        const Real length_raw = computeLength(start_pt_raw, end_pt);
        
        // get the direction of the line
        const vec2 line_dir = {(end_pt[0] - start_pt_raw[0])/length_raw, (end_pt[1] - start_pt_raw[1])/length_raw};
        
        // adjust the starting point
        const vec2 start_pt = {start_pt_raw[0] + (include_start ? 0.0 : delta)*line_dir[0], start_pt_raw[1] + (include_start ? 0.0 : delta)*line_dir[1]};
        
        // now we do as normal
        const Real length = computeLength(start_pt, end_pt);
        const auto line_info = getNumberOfPoints(length, delta, include_end);
        std::vector<vec2> retval(line_info.first);
        
        for(int i=0;i<line_info.first;++i)
        {
            retval[i][0] = start_pt[0] + i*line_info.second * line_dir[0];
            retval[i][1] = start_pt[1] + i*line_info.second * line_dir[1];
        }
        return retval;
    }
    
    std::vector<vec2> connectWithCCWArc(const Real radius, const Real start_theta_raw, const Real end_theta_raw, const Real delta, const bool include_start = true, const bool include_end = true) const
    {
        const Real delta_angle = delta/radius;
        
        const Real start_theta = start_theta_raw + (include_start ? 0.0 : delta_angle);
        const Real end_theta = (end_theta_raw < start_theta ? end_theta_raw + 2*M_PI : end_theta_raw);
        
        const Real delta_theta = end_theta - start_theta;
        const auto arc_info = getNumberOfPoints(delta_theta*radius, delta, include_end);
        std::vector<vec2> retval(arc_info.first);
        
        for(int i=0;i<arc_info.first;++i)
        {
            const Real theta = start_theta + i * arc_info.second/radius;
            retval[i][0] = radius * std::cos(theta);
            retval[i][1] = radius * std::sin(theta);
        }
        
        return retval;
    }
    
    std::vector<vec2> connectWithCCWArc(const Real radius, const vec2 & start_pt, const vec2 & end_pt, const Real delta, const bool include_start = true, const bool include_end = true) const
    {
        const Real theta_start = std::atan2(start_pt[1], start_pt[0]);
        const Real theta_end = std::atan2(end_pt[1], end_pt[0]);
        return connectWithCCWArc(radius, theta_start, theta_end, delta, include_start, include_end);
    }
    
    std::vector<vec2> connectWithCCWEllipse(const Real radiusX, const Real radiusY, const Real start_theta, const Real end_theta_raw, const Real delta) const // always contains starting and ending point
    {
        const Real end_theta = (end_theta_raw < start_theta ? end_theta_raw + 2*M_PI : end_theta_raw);
        
        const Real max_rad = std::max(radiusX, radiusY);
        const Real min_rad = std::min(radiusX, radiusY);
        const bool flip = radiusY > radiusX;
        
        EllipticalCircumference ellipse(max_rad, min_rad);
        
        {
            const Real arclength_start = max_rad * ellipse.getArcLengthRel(start_theta);
            const auto point_start = ellipse.getPointFromArcLength(arclength_start);
            const Real arclength_end = max_rad * ellipse.getArcLengthRel(end_theta);
            const auto point_end = ellipse.getPointFromArcLength(arclength_end);
            printf("start : %10.10e pi,  %10.10e, (%10.10e, %10.10e) \t end : %10.10e pi ,  %10.10e, (%10.10e, %10.10e)\n", start_theta/M_PI, arclength_start, point_start.first, point_start.second, end_theta/M_PI, arclength_end, point_end.first, point_end.second);
            
            const Real start_theta_check = ellipse.getThetaFromArclength(arclength_start);
            const Real end_theta_check = ellipse.getThetaFromArclength(arclength_end);
            printf("check theta start : %10.10e, %10.10e \t\t check theta end : %10.10e, %10.10e\n", start_theta, start_theta_check, end_theta, end_theta_check);
        }
        // compute the arclength at the start
        const Real arclength_start = max_rad * ellipse.getArcLengthRel(start_theta);
        const Real arclength_end = max_rad * ellipse.getArcLengthRel(end_theta);
        const Real arclength_diff = arclength_end - arclength_start;
        
        const auto ellarc_info = getNumberOfPoints(arclength_diff, delta, true);
        std::vector<vec2> retval(ellarc_info.first);
        
        for(int i=0;i<ellarc_info.first;++i)
        {
            const Real arclength = arclength_start + i * ellarc_info.second;
            const auto point = ellipse.getPointFromArcLength(arclength);
            // rotate it by pi/2 in CCW direction if its flipped
            retval[i][0] = flip ? -point.second : point.first;
            retval[i][1] = flip ? point.first : point.second;
        }
        
        printf("OUT :first and last point : %10.10e, %10.10e \t %10.10e, %10.10e\n", retval[0][0], retval[0][1], retval[ellarc_info.first-1][0], retval[ellarc_info.first-1][1]);
        return retval;
    }
    
    std::vector<vec2> connectWithCCWEllipse(const Real radiusX, const Real radiusY, const vec2 & start_pt, const vec2 & end_pt, const Real delta) const
    {
        printf("IN : first and last point : %10.10e, %10.10e \t %10.10e, %10.10e\n", start_pt[0], start_pt[1], end_pt[0], end_pt[1]);
        
        const bool flip = radiusY > radiusX; // if the major/minor axes are flipped: we should rotate the points
        
        const Real theta_start = flip ? std::atan2(-start_pt[0]/radiusX, start_pt[1]/radiusY) : std::atan2(start_pt[1]/radiusY, start_pt[0]/radiusX);
        const Real theta_end =   flip ? std::atan2(-end_pt[0]/radiusX, end_pt[1]/radiusY) : std::atan2(end_pt[1]/radiusY, end_pt[0]/radiusX);
        
        return connectWithCCWEllipse(radiusX, radiusY, theta_start, theta_end, delta);
    }
    
    std::vector<vec2> foldingPetal(const Real petal_height, const Real petal_radius, const Real offset_y, const Real delta, const Real center_r, const Real theta = 0) const
    {
        // want to find the starting point of the arc : solve for the circle-circle intersection
        // solution to FullSimplify[Solve[{x^2 + y^2 == R^2, (x - a)^2 + (y - b)^2 == R2^2}, {x, y}], Assumptions -> {R > 0, a > 0, b > 0, R2 > 0}]
        // where R = center_r, R2 = petal_radius, b = 0.5*petal_height, a = -sqrt(R2^2 - b^2)
        const Real petal_ctr_y =  0.5*petal_height + offset_y;
        const Real petal_ctr_x = -std::sqrt(std::pow(petal_radius,2) - std::pow(0.5*petal_height,2));
        
        // starting point
        const Real x_start = (std::pow(petal_ctr_x,3) + petal_ctr_x*(std::pow(petal_ctr_y,2) + std::pow(center_r,2) - std::pow(petal_radius,2)) + petal_ctr_y*std::sqrt(-std::pow(std::pow(petal_ctr_x,2) + std::pow(petal_ctr_y,2) - std::pow(center_r,2),2) + 2*(std::pow(petal_ctr_x,2) + std::pow(petal_ctr_y,2) + std::pow(center_r,2))*std::pow(petal_radius,2) - std::pow(petal_radius,4)))/(2.*(std::pow(petal_ctr_x,2) + std::pow(petal_ctr_y,2)));
        const Real y_start = (std::pow(petal_ctr_y,3) + petal_ctr_y*(std::pow(petal_ctr_x,2) + std::pow(center_r,2) - std::pow(petal_radius,2)) - petal_ctr_x*std::sqrt(-std::pow(std::pow(petal_ctr_x,2) + std::pow(petal_ctr_y,2) - std::pow(center_r,2),2) + 2*(std::pow(petal_ctr_x,2) + std::pow(petal_ctr_y,2) + std::pow(center_r,2))*std::pow(petal_radius,2) - std::pow(petal_radius,4)))/(2.*(std::pow(petal_ctr_x,2) + std::pow(petal_ctr_y,2)));
        
        // create starting and ending point from perspective of circle
        const vec2 pt_start = {x_start - petal_ctr_x, y_start - petal_ctr_y};
        const vec2 pt_top_1 = {0.0 - petal_ctr_x, petal_height + offset_y - petal_ctr_y};
        
        // other side : petal_ctr_x is flippe and x_start is flipped also
        const vec2 pt_top_2 = {0.0 + petal_ctr_x, petal_height + offset_y - petal_ctr_y};
        const vec2 pt_end = {-x_start + petal_ctr_x, y_start - petal_ctr_y};
        
        std::vector<vec2> part1 = connectWithCCWArc(petal_radius, pt_start, pt_top_1, delta, true, false); // does not include tip
        std::vector<vec2> part2 = connectWithCCWArc(petal_radius, pt_top_2, pt_end, delta, true, true); // includes tip
        
        // shift the first part
        for(auto & pt : part1)
        {
            pt[0] += petal_ctr_x;
            pt[1] += petal_ctr_y;
        }
        
        // shift the second part
        for(auto & pt : part2)
        {
            pt[0] -= petal_ctr_x;
            pt[1] += petal_ctr_y;
        }
        
        // combine
        std::vector<vec2> retval;
        retval.reserve(part1.size() + part2.size());
        retval.insert(retval.end(), part1.begin(), part1.end());
        retval.insert(retval.end(), part2.begin(), part2.end());
        
        // rotate
        if(theta != 0)
            rotate(retval, theta);
        
        return retval;
    }
    
    std::vector<vec2> shortPetal(const Real short_l, const Real short_w, const Real delta, const Real center_r, const Real theta = 0) const
    {
        // want to find the starting point of the ellipse : solve for the intersection
        // solution to FullSimplify[Solve[{x^2 + y^2 == R^2, (x/a)^2 + ((y - b)/b)^2 == 1}, {x, y}],Assumptions -> {R > 0, a > 0, b > 0}]
        // where R = center_r, a = short_w/2, and b = short_l/2
        
        const Real x_start = (short_w*std::sqrt(-2*std::pow(center_r,2) + (std::pow(short_l,2)*(-std::pow(short_w,2) + std::sqrt(std::pow(short_w,4) + 4*std::pow(center_r,2)*(std::pow(short_l,2) - std::pow(short_w,2)))))/(std::pow(short_l,2) - std::pow(short_w,2))))/(std::sqrt(2)*std::sqrt((short_l - short_w)*(short_l + short_w)));
        const Real y_start =  (short_l*(-std::pow(short_w,2) + std::sqrt(std::pow(short_w,4) + 4*std::pow(center_r,2)*(std::pow(short_l,2) - std::pow(short_w,2)))))/(2.*(std::pow(short_l,2) - std::pow(short_w,2)));
        printf("xstart, ystart = %10.10e , %10.10e\n", x_start, y_start);
        
        // create an ellipse from the start to the end with the origin at (0,0)
        const vec2 pt_start = {x_start, y_start - 0.5*short_l};
        const vec2 pt_end = {-x_start, y_start - 0.5*short_l};
        std::vector<vec2> retval = connectWithCCWEllipse(0.5*short_w, 0.5*short_l, pt_start, pt_end, delta);
        
        // shift the ellipse
        for(auto & pt : retval) pt[1] += 0.5*short_l;
        
        // rotate the ellipse
        if(theta != 0)
            rotate(retval, theta);
        
        return retval;
    }
    
    std::vector<vec2> longPetal(const Real long_l, const Real long_w, const Real delta, const Real center_r, const Real theta = 0) const
    {
        // compute the start/end points for each segment
        const vec2 startRight = {0.5*long_w, std::sqrt(center_r * center_r - 0.25*long_w*long_w)};
        const vec2 endLeft = {-0.5*long_w, std::sqrt(center_r * center_r - 0.25*long_w*long_w)};
        
        const vec2 endRight = {0.5*long_w, center_r + long_l};
        const vec2 startLeft = {-0.5*long_w, center_r + long_l};
        const Real length_right = computeLength(startRight, endRight);
        const Real length_left = computeLength(startLeft, endLeft);
        
        // compute number of points
        const auto info_right = getNumberOfPoints(length_right, delta);
        const auto info_left = getNumberOfPoints(length_left, delta);
        const auto info_curved = getNumberOfPoints(0.5*M_PI*long_w, delta);
        
        const int nPoints = info_right.first + info_left.first + info_curved.first;
        std::vector<vec2> retval(nPoints);
        int idx = 0;
        for(int i=0;i<info_right.first;++i, ++idx)
        {
            retval[idx][0] = 0.5*long_w;
            retval[idx][1] = startRight[1] + i*info_right.second;
        }
        for(int i=0;i<info_curved.first;++i, ++idx)
        {
            const Real theta = i * info_curved.second/(0.5*long_w);
            retval[idx][0] = 0.5*long_w*std::cos(theta);
            retval[idx][1] = center_r + long_l + 0.5*long_w*std::sin(theta);
        }
        for(int i=0;i<info_left.first;++i, ++idx)
        {
            retval[idx][0] = -0.5*long_w;
            retval[idx][1] = startLeft[1] - (i+1)*info_left.second;
        }
        
        if(theta != 0)
            rotate(retval, theta);
        
        return retval;
    }
    
    std::vector<vec2> bottomPetal(const Real bot_l, const Real bot_w, const Real delta, const Real center_r, const Real theta=0) const
    {
        const Real bot_r = 0.5*bot_w; // rb in mathematica
        
        const Real sqrtfac = std::sqrt(std::pow(center_r + bot_r,2) - std::pow(std::sqrt(2)*bot_r + 0.5*bot_w, 2));
        
        // get all the circular arc segments (they contain the first and last points already)
        const vec2 center_1 = {-0.5*bot_w, -sqrtfac};
        std::vector<vec2> segment_1 = connectWithCCWArc(bot_r, 0.75*M_PI, M_PI, delta, true, true);
        translate(segment_1, center_1);
        
        const vec2 center_3 = {-0.5*bot_w, -0.5*bot_l - sqrtfac};
        std::vector<vec2> segment_3 = connectWithCCWArc(bot_r, M_PI, 1.5*M_PI, delta, true, true);
        translate(segment_3, center_3);
        
        const vec2 center_5 = {0.0, -bot_l - bot_r - sqrtfac};
        std::vector<vec2> segment_5 = connectWithCCWArc(0.5*bot_w, M_PI, 2*M_PI, delta, true, true);
        translate(segment_5, center_5);
        
        const vec2 center_7 = {0.5*bot_w, -0.5*bot_l - sqrtfac};
        std::vector<vec2> segment_7 = connectWithCCWArc(bot_r, 1.5*M_PI, 2*M_PI, delta, true, true);
        translate(segment_7, center_7);
        
        const vec2 center_9 = {0.5*bot_w, -sqrtfac};
        std::vector<vec2> segment_9 = connectWithCCWArc(bot_r, 0.0, 0.25*M_PI, delta, true, true); // also contains end point
        translate(segment_9, center_9);
        
        // get the straight line segments
        const std::vector<vec2> segment_2 = connectWithStraightLine(segment_1[segment_1.size()-1], segment_3[0], delta, false, false);
        const std::vector<vec2> segment_4 = connectWithStraightLine(segment_3[segment_3.size()-1], segment_5[0], delta, false, false);
        const std::vector<vec2> segment_6 = connectWithStraightLine(segment_5[segment_5.size()-1], segment_7[0], delta, false, false);
        const std::vector<vec2> segment_8 = connectWithStraightLine(segment_7[segment_7.size()-1], segment_9[0], delta, false, false);
        
        std::vector<vec2> retval;
        retval.reserve(segment_1.size() + segment_2.size() + segment_3.size() + segment_4.size() + segment_5.size() + segment_6.size() + segment_7.size() + segment_8.size() + segment_9.size());
        retval.insert(retval.end(), segment_1.begin(), segment_1.end());
        retval.insert(retval.end(), segment_2.begin(), segment_2.end());
        retval.insert(retval.end(), segment_3.begin(), segment_3.end());
        retval.insert(retval.end(), segment_4.begin(), segment_4.end());
        retval.insert(retval.end(), segment_5.begin(), segment_5.end());
        retval.insert(retval.end(), segment_6.begin(), segment_6.end());
        retval.insert(retval.end(), segment_7.begin(), segment_7.end());
        retval.insert(retval.end(), segment_8.begin(), segment_8.end());
        retval.insert(retval.end(), segment_9.begin(), segment_9.end());
        
        if(theta != 0)
            rotate(retval, theta);
        
        return retval;
    }
    
    std::vector<vec2> path_from_petals(const Real center_r, const Real delta, const std::vector<std::reference_wrapper<const std::vector<vec2>> > petals) const
    {
        const int nPetals = (int)petals.size();
        
        std::vector<vec2> retval;
        
        for(int i=0;i<nPetals;++i)
        {
            const int ip = (i+1)%nPetals;
            // add the petal
            const std::vector<vec2> & petal_points = std::remove_reference<const std::vector<vec2>&>::type(petals[i]);
            const std::vector<vec2> & petal_points_next = std::remove_reference<const std::vector<vec2>&>::type(petals[ip]);
            retval.insert( retval.end(), petal_points.begin(), petal_points.end() );
            // add the arc to the next petal
            const std::vector<vec2> arc_points = connectWithCCWArc(center_r, petal_points[petal_points.size() - 1], petal_points_next[0], delta, false, false); // assume petals contain starting and end points
            retval.insert( retval.end(), arc_points.begin(), arc_points.end() );
        }
        
        return retval;
    }
    
public:
    Geometry_Flower_Base()
    {}
};

class Geometry_FoldingFlower : public Geometry_Flower_Base
{
protected:
    const Real petal_height;
    const Real center_r;
    const Real offset_y;
    const Real edgeLength;
    
    Real petal_radius;
    
    void getPlateGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        const Real petal_dtheta = 2*M_PI/5;
        
        // create the flower : CCW
        const std::vector<vec2> petal_1 = foldingPetal(petal_height, petal_radius, offset_y, edgeLength, center_r, 0.0);
        const std::vector<vec2> petal_2 = foldingPetal(petal_height, petal_radius, offset_y, edgeLength, center_r, -petal_dtheta);
        const std::vector<vec2> petal_3 = foldingPetal(petal_height, petal_radius, offset_y, edgeLength, center_r, -2*petal_dtheta);
        const std::vector<vec2> petal_4 = foldingPetal(petal_height, petal_radius, offset_y, edgeLength, center_r, -3*petal_dtheta);
        const std::vector<vec2> petal_5 = foldingPetal(petal_height, petal_radius, offset_y, edgeLength, center_r, -4*petal_dtheta);

        std::vector<std::reference_wrapper<const std::vector<vec2>> > petals;
        petals.push_back(std::ref(petal_1));
        petals.push_back(std::ref(petal_2));
        petals.push_back(std::ref(petal_3));
        petals.push_back(std::ref(petal_4));
        petals.push_back(std::ref(petal_5));
        
        const std::vector<vec2> all_points = path_from_petals(center_r, edgeLength, petals);
        const int nTotalPoints = (int)all_points.size();
        
        // convert to eigen arrays
        Eigen::MatrixXd polygon_vertices(nTotalPoints, 2);
        Eigen::MatrixXi polygon_edges(nTotalPoints, 2);
        
        for(int i=0;i<nTotalPoints;++i)
        {
            polygon_vertices(i,0) = all_points[i][0];
            polygon_vertices(i,1) = all_points[i][1];
            
            polygon_edges(i,0) = i;
            polygon_edges(i,1) = (i+1)%nTotalPoints;
        }
        
        const std::string flags = "q20a"+std::to_string(std::sqrt(3.0)/4.0*edgeLength*edgeLength)+(verbose ? "" : "Q");;
        Triangulate triangulate(polygon_vertices, polygon_edges, flags);
        triangulate.get(vertices, face2vertices, vertices_bc);
    }
    
public:
    Geometry_FoldingFlower(const Real petal_height_in, const Real center_r_in, const Real offset_y_in, const Real edgeLength_in):
    Geometry_Flower_Base(),
    petal_height(petal_height_in),
    center_r(center_r_in),
    offset_y(offset_y_in),
    edgeLength(edgeLength_in)
    {
        petal_radius = petal_height / (2 * std::cos(0.5*M_PI - M_PI/5.0)); //nPetals = 5
    }
    
};
class Geometry_Orchid : public Geometry_Flower_Base
{
protected:
    // long petal parameters
    const Real long_l; // total length
    const Real long_w; // total width
    
    // short petal parameters
    const Real short_l; // total length
    const Real short_w; // total width
    
    // bottom petal parameters
    const Real bot_l;
    // bot_w defined below
    
    // center
    const Real center_r;
    
    // discretization length
    const Real edgeLength;
    
    Real theta, bot_w;
    
    void getPlateGeometry(Eigen::MatrixXd & vertices, Eigen::MatrixXi & face2vertices, Eigen::MatrixXb & vertices_bc) const override
    {
        // create the orchid
        const std::vector<vec2> long_petal_right = longPetal(long_l, long_w, edgeLength, center_r, theta);
        const std::vector<vec2> short_petal_top = shortPetal(short_l, short_w, edgeLength, center_r, 0.0);
        const std::vector<vec2> long_petal_left = longPetal(long_l, long_w, edgeLength, center_r, -theta);
        const std::vector<vec2> short_petal_left = shortPetal(short_l, short_w, edgeLength, center_r, -2*theta);
        const std::vector<vec2> bottom_petal = bottomPetal(bot_l, bot_w, edgeLength, center_r, 0.0);
        const std::vector<vec2> short_petal_right = shortPetal(short_l, short_w, edgeLength, center_r, +2*theta);
        
        std::vector<std::reference_wrapper<const std::vector<vec2>> > petals;
        petals.push_back(std::ref(long_petal_right));
        petals.push_back(std::ref(short_petal_top));
        petals.push_back(std::ref(long_petal_left));
        petals.push_back(std::ref(short_petal_left));
        petals.push_back(std::ref(bottom_petal));
        petals.push_back(std::ref(short_petal_right));
        
        const std::vector<vec2> all_points = path_from_petals(center_r, edgeLength, petals);
        const int nTotalPoints = (int)all_points.size();
        
        // convert to eigen arrays
        Eigen::MatrixXd polygon_vertices(nTotalPoints, 2);
        Eigen::MatrixXi polygon_edges(nTotalPoints, 2);

        for(int i=0;i<nTotalPoints;++i)
        {
            polygon_vertices(i,0) = all_points[i][0];
            polygon_vertices(i,1) = all_points[i][1];
            
            polygon_edges(i,0) = i;
            polygon_edges(i,1) = (i+1)%nTotalPoints;
        }
        
        const std::string flags = "q20a"+std::to_string(std::sqrt(3.0)/4.0*edgeLength*edgeLength)+(verbose ? "" : "Q");;
        Triangulate triangulate(polygon_vertices, polygon_edges, flags);
        triangulate.get(vertices, face2vertices, vertices_bc);
    }
    
public:
    Geometry_Orchid(const Real edgeLength_in):
    Geometry_Flower_Base(),
    long_l(30.0),
    long_w(2.5),
    short_l(20.0),
    short_w(6.0),
    bot_l(6.0),
    center_r(5.0),
    edgeLength(edgeLength_in)
    {
        // derived parameters
        theta = std::asin(long_w/(2.*center_r)) + std::asin((short_w* std::sqrt(1 - std::pow(std::pow(short_l,2) - std::sqrt(std::pow(short_w,4) + std::pow(center_r,2)*(std::pow(short_l,2) - std::pow(short_w,2))),2)/std::pow(std::pow(short_l,2) - std::pow(short_w,2),2)))/center_r);
        bot_w = (1 + 1/std::sqrt(2))*center_r*std::sin(2*std::asin(long_w/(2.*center_r)) + 3*std::asin((short_w*std::sqrt(1 - std::pow(std::pow(short_l,2) - std::sqrt(std::pow(short_w,4) + (std::pow(short_l,2) - std::pow(short_w,2))*std::pow(center_r,2)),2)/ std::pow(std::pow(short_l,2) - std::pow(short_w,2),2)))/center_r));
    }
};


#endif /* Geometry_Orchid_hpp */
