//
//  Sim_MonoLayer_Growth.cpp
//  Elasticity
//
//  Created by Wim van Rees on 8/17/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#include "Sim_MonoLayer_Growth.hpp"
#include "Geometry.hpp"
#include "GrowthHelper.hpp"
#include "CombinedOperator_Parametric.hpp"
#include "MaterialProperties.hpp"

#include <igl/edge_topology.h>

void Sim_MonoLayer_Growth::run()
{
    const std::string runCase = parser.parse<std::string>("-case", "");
    if(runCase == "cone")
        runCone();
    else if(runCase == "edgegrowth")
        runEdgeGrowth();
    else if(runCase == "basic_disk")
        run_basic_disk();
    else
    {
        std::cout << "No valid runCase defined. Options are \n";
        std::cout << "\t -case cone\n";
        std::cout << "\t -case edgegrowth\n";
        std::cout << "\t -case basic_disk\n";
    }
}


void Sim_MonoLayer_Growth::runCone()
{
    tag = "coneGrowth";
    
    const Real radius = 1.0;
    const Real length = 4;
    const Real edgeLength = parser.parse<Real>("-res", 0.05);
    CylinderTube geometry(radius, length, edgeLength);
    mesh.init(geometry);
    
    // set the growth profile
    const int nFaces = mesh.getNumberOfFaces();
    std::vector<GrowthState> growth;
    growth.reserve(nFaces);
    
    {
        const Real cfac = parser.parse<Real>("-cfac", 1.0);
        const Real pfac = parser.parse<Real>("-pfac", 4.0);
        
        const auto restState = mesh.getRestConfiguration();
        for(int i=0;i<nFaces;++i)
        {
            const auto rinfo = restState.getTriangleInfoLite(mesh.getTopology(), i);
            const auto restfacecenter = rinfo.computeFaceCenter();
            
            // get the face position
            const Real fx = restfacecenter(0);
            const Real fy = restfacecenter(1);
            const Real fz = restfacecenter(2);
            
            // compute the growth factors
            const Real radius = std::sqrt(fx*fx + fy*fy);
            const Real s1 = 1.0 + cfac * std::pow(fz / length, pfac); // azimuthal
            const Real s2 = 1.0;
            
            // compute the growth directions in R3
            const Eigen::Vector3d v1 = (Eigen::Vector3d() << -fy/radius, fx/radius, 0.0).finished(); // azimuthal dir
//            const Eigen::Vector3d v2 = (Eigen::Vector3d() <<  0.0, 0.0, 1.0).finished(); // other tangent dir
            const Eigen::Vector3d facenormal = (rinfo.e2).cross(rinfo.e0).normalized();
            const Eigen::Vector3d v2 = (facenormal.cross(v1)).normalized(); // so that dir_1 cross dir_2 = facenormal

            // compute the grown metric through DecomposedGrowthState (wrt current state)
            DecomposedGrowthState growth_state(s1, s2, v1, v2);
            
            Eigen::MatrixXd rxy_base(3,2);
            rxy_base << rinfo.e1, rinfo.e2;
            
            const Eigen::Matrix2d a_final = growth_state.computeMetric(rxy_base);
            //growth.emplace_back(rxy_base, a_final);
            GrowthState growthState(rxy_base, a_final);
            growth.push_back(growthState);

        }
    }
    
    const bool interpolate_logeucl = parser.parse<bool>("-interp_logeucl", false);
    
    assignGrowthToMetric(growth, 0.0, interpolate_logeucl);
    
    dumpOrthoNew(growth, tag+"_init");
    
    // define the material
    const Real E = 1.0;
    const Real nu = 0.3;
    const Real h = 1e-2;
    MaterialProperties_Iso_Constant matprop(E, nu, h);
    CombinedOperator_Parametric<tMesh> engOp(matprop);
    
    dumpOrthoNew(growth, tag+"_final_"+helpers::ToString(0,2));
    
    std::vector<Real> swellingRates = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    const size_t nSwellingRates = swellingRates.size();
    for(size_t i=0;i<nSwellingRates;++i)
    {
        // set growth
        assignGrowthToMetric(growth, swellingRates[i], interpolate_logeucl);
        
        // add noise
        addNoiseToVertices(0.02*edgeLength);
        addNoiseToEdgeDirectors(0.01*M_PI);
        
        // solve minimum
        Real eps = 1e-2;
        minimizeEnergy(engOp, eps);
        
        // dump
        dumpOrthoNew(growth, tag+"_final_"+helpers::ToString((int)i+1,2)); // i+1 since we already dumped 0
    }
}

void Sim_MonoLayer_Growth::assignGrowthToMetric(const std::vector<GrowthState> & growth, const Real t, const bool interp_logeucl)
{
    const int nFaces = mesh.getNumberOfFaces();
    tVecMat2d & firstFF = mesh.getRestConfiguration().getFirstFundamentalForms();
    
    for(int i=0;i<nFaces;++i)
        firstFF[i] = (interp_logeucl ? growth[i].interpolateLogEucl(t) : growth[i].interpolate(t));
}


void Sim_MonoLayer_Growth::runEdgeGrowth()
{
    tag = "test_growthedge";
    
    const Real Lx = 4.0;
    const Real Ly = 1.0;
    
    const Real maxTriangleAreaRel = 0.0001;//0.0005;
    
    RectangularPlate geometry(0.5*Lx, 0.5*Ly, maxTriangleAreaRel, {false, false}, {true, false});
    mesh.init(geometry, true);
    
    // now we create the abar
    const int nFaces = mesh.getNumberOfFaces();
    Eigen::VectorXd growthRates_1(nFaces);
    Eigen::VectorXd growthRates_2(nFaces);
    Eigen::VectorXd growthAngles(nFaces);
    {
        const auto restState = mesh.getRestConfiguration();
        const Real small_l_fac = 40*1e-3;
        
        for(int i=0;i<nFaces;++i)
        {
            const auto restfacecenter = restState.getTriangleInfoLite(mesh.getTopology(), i).computeFaceCenter();
            
            growthAngles(i) = 0.0; // so that dir_1 = dir_p = (1,0,0) and dir_2 = dir_o = (0,1,0)
            const Real L_min_z = 0.5*Ly - restfacecenter(1); // z=-L/2 --> L, z=L/2 --> 0
            growthRates_1(i) = 1.0 / (1.0 +  L_min_z/small_l_fac) - 1.0 / (1.0 +  Ly/small_l_fac); // between 0 (clamped edge) and 1 - 1/(1 + Ly/small_l_fac)
            growthRates_2(i) = 0.0;
        }
    }
    GrowthHelper<tMesh>::computeAbarsOrthoGrowth(mesh, growthAngles, growthRates_1, growthRates_2, mesh.getRestConfiguration().getFirstFundamentalForms());
    
    dumpWithGrowthRates(growthRates_1, tag+"_init");
    
    // put some random perturbations
    if(false)
    {
        std::mt19937 gen;
        gen.seed(42);
        std::uniform_real_distribution<Real> dist(-0.005, 0.005);
        auto perturb_v = [&](Eigen::Vector3d in)
        {
            Eigen::Vector3d retval;
            retval << in(0)+dist(gen), in(1)+dist(gen), in(2)+dist(gen);
            return retval;
        };
        auto perturb_e = [&](Real in)
        {
            return in + dist(gen);
        };
        
        mesh.changeVertices(perturb_v);
        mesh.changeEdgeDirectors(perturb_e);
    }
    else if (true)
    {
        // deform into parabolic surface with origin on edge
        const Real maxHeight = 0.1*Ly;
        
        auto vertices = mesh.getCurrentConfiguration().getVertices();
        const auto & restvertices = mesh.getRestConfiguration().getVertices();
        const int nVertices = mesh.getNumberOfVertices();
        for(int i=0;i<nVertices;++i)
        {
            const Real relY = (0.5 + restvertices(i,1)/Ly); // between 0 and 1
            const Real deformFac = maxHeight/(0.25*Lx*Lx) * relY;

            const Real xSq = std::pow(restvertices(i,0),2);
            const Real height = deformFac*xSq;
            vertices(i,2) = height;
        }
        
    }
    
    dumpWithGrowthRates(growthRates_1, tag+"_deform");
    
    // define the material
    const Real E = 1e8;
    const Real nu = 0.3;
    const Real h = 1e-4;
    MaterialProperties_Iso_Constant matprop(E, nu, h);
    CombinedOperator_Parametric<tMesh> engOp(matprop);
    
    std::vector<Real> swellingRates = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    const size_t nSwellingRates = swellingRates.size();
    for(size_t i=0;i<nSwellingRates;++i)
    {
        // set growth
        tVecMat2d & aforms = mesh.getRestConfiguration().getFirstFundamentalForms();
        GrowthHelper<tMesh>::computeAbarsOrthoGrowth(mesh, growthAngles, swellingRates[i]*growthRates_1, swellingRates[i]*growthRates_2, aforms);
        
        Real eps = 1e-2;
        minimizeEnergy(engOp, eps);
        
        dumpWithGrowthRates(growthRates_1, tag+"_"+helpers::ToString((int)i, 2)+"_final");
    }
}


void Sim_MonoLayer_Growth::run_basic_disk()
{
    tag = "test_basic_disk";
    
    const Real radius = parser.parse<Real>("-R", 1.0);
    const int res = parser.parse<int>("-res", 128);

    const Real innerR = parser.parse<Real>("-innerR", 0.0);
    if(std::abs(innerR) > std::numeric_limits<Real>::epsilon())
    {
        AnnulusPlate geometry(radius, innerR, 2.0*M_PI/res * radius);
        mesh.init(geometry);
    }
    else
    {
        CircularPlate geometry(radius, res, false);
        mesh.init(geometry);
    }
    
    // now we create the abar
    const int nFaces = mesh.getNumberOfFaces();
    
    std::vector<GrowthState> growth;
    growth.reserve(nFaces);
    
    const std::string growthCase = parser.parse<std::string>("-growthcase", "ortho");

    if(growthCase == "ortho")
    {
        const Real growthAngle = parser.parse<Real>("-growthangle", 0.0)*M_PI; // 0 is radial, 1/2 is azimuthal
        
        const Real s1 = parser.parse<Real>("-s1", 1.0) + 1;
        const Real s2 = parser.parse<Real>("-s2", 0.0) + 1;
        
        const auto restState = mesh.getRestConfiguration();
        for(int i=0;i<nFaces;++i)
        {
            const auto rinfo = restState.getTriangleInfoLite(mesh.getTopology(), i);
            const auto restfacecenter = rinfo.computeFaceCenter();
            // get the azimuthal direction
            const Real fx = restfacecenter(0);
            const Real fy = restfacecenter(1);
            
            const Real rad = std::sqrt(fx*fx + fy*fy);
            if(rad < std::numeric_limits<Real>::epsilon())
            {
                helpers::catastrophe("Error : can't have orthotropic growth at the origin -- choose an annulus instead (i.e. set innerR > 0)\n", __FILE__, __LINE__);
            }
            else
            {
                const Eigen::Vector3d azi = (Eigen::Vector3d() << -fy/rad, fx/rad, 0.0).finished(); // azimuthal dir
                const Eigen::Vector3d facenormal = (Eigen::Vector3d() << 0.0, 0.0, 1.0).finished();
                const Eigen::Vector3d radial = (facenormal.cross(azi)).normalized();

                const Eigen::Vector3d v1 =  std::cos(growthAngle)*radial + std::sin(growthAngle)*azi;
                const Eigen::Vector3d v2 = -std::sin(growthAngle)*radial + std::cos(growthAngle)*azi;
                DecomposedGrowthState growth_state(s1, s2, v1, v2);
                
                Eigen::MatrixXd rxy_base(3,2);
                rxy_base << rinfo.e1, rinfo.e2;
                
                const Eigen::Matrix2d a_final = growth_state.computeMetric(rxy_base);
		//                growth.emplace_back(rxy_base, a_final);
	    GrowthState growthState(rxy_base, a_final);
            growth.push_back(growthState);
		
            }
        }
    }
    else if(growthCase == "spherical")
    {
        // sphere mapped to disk according to stereographic projection
        const auto restState = mesh.getRestConfiguration();
        for(int i=0;i<nFaces;++i)
        {
            const auto rinfo = restState.getTriangleInfoLite(mesh.getTopology(), i);
            const auto restfacecenter = rinfo.computeFaceCenter();
            // get the radius (make it a hemisphere -- normalize by disk radius)
            const Real fx = restfacecenter(0)/radius;
            const Real fy = restfacecenter(1)/radius;
            
            const Real rsq = (fx*fx + fy*fy);
            
            const Real s_iso = 2.0/(1.0 + rsq); // 1 at edges, 2 at center
            const Eigen::Vector3d v1 = (Eigen::Vector3d() << 1,0,0).finished(); // does not matter as long as its tangential
            const Eigen::Vector3d v2 = (Eigen::Vector3d() << 0,1,0).finished(); // does not matter as long as its tangential
            
            DecomposedGrowthState growth_state(s_iso, s_iso, v1, v2);
            
            Eigen::MatrixXd rxy_base(3,2);
            rxy_base << rinfo.e1, rinfo.e2;
            
            const Eigen::Matrix2d a_final = growth_state.computeMetric(rxy_base);
	    GrowthState growthState(rxy_base, a_final);
            growth.push_back(growthState);
        }
    }
    else if(growthCase == "validation")
    {
        // sphere mapped to disk according to sin(r) projection
        const auto restState = mesh.getRestConfiguration();
        for(int i=0;i<nFaces;++i)
        {
            const auto rinfo = restState.getTriangleInfoLite(mesh.getTopology(), i);
            const auto restfacecenter = rinfo.computeFaceCenter();
            // get the radius (make it a hemisphere -- normalize by disk radius)
            const Real fx = restfacecenter(0)/radius;
            const Real fy = restfacecenter(1)/radius;
            
            const Real rad = std::sqrt(fx*fx + fy*fy);
            
            if(rad < std::numeric_limits<Real>::epsilon())
            {
                helpers::catastrophe("Error : can't have orthotropic growth at the origin -- choose an annulus instead (i.e. set innerR > 0)\n", __FILE__, __LINE__);
            }
            else
            {
                const Real radsq = rad*rad;
                const Real s_azimuthal = (radsq < std::numeric_limits<Real>::epsilon() ? 1.0 : std::sin(rad)/rad);
                const Real s_radial = 1.0; // no growth
                
                const Eigen::Vector3d v1 = (Eigen::Vector3d() << -fy/rad, fx/rad, 0.0).finished(); // azimuthal dir
                const Eigen::Vector3d facenormal = (Eigen::Vector3d() << 0.0, 0.0, 1.0).finished();
                const Eigen::Vector3d v2 = (facenormal.cross(v1)).normalized(); // radial dir
                
                DecomposedGrowthState growth_state(s_azimuthal, s_radial, v1, v2);
                
                Eigen::MatrixXd rxy_base(3,2);
                rxy_base << rinfo.e1, rinfo.e2;
                
                const Eigen::Matrix2d a_final = growth_state.computeMetric(rxy_base);
                //growth.emplace_back(rxy_base, a_final);
		GrowthState growthState(rxy_base, a_final);
		growth.push_back(growthState);
		
            }
        }
    }
    else
    {
        std::cout << "Pick a valid growth case: should be one of:" << std::endl;
        std::cout << "\t - ortho " << std::endl;
        std::cout << "\t - spherical " << std::endl;
        std::cout << "\t - validation " << std::endl;
        return;
    }
    
    const bool interpolate_logeucl = parser.parse<bool>("-interp_logeucl", false);
    assignGrowthToMetric(growth, 0.0, interpolate_logeucl);
    dumpOrthoNew(growth, tag+"_init");
    
    // define the material
    const Real E = parser.parse<Real>("-E", 1.0);
    const Real nu = (growthCase == "validation" ? 0.5 : 0.3);
    const Real h = parser.parse<Real>("-h", 0.01);
    MaterialProperties_Iso_Constant matprop(E, nu, h);
    CombinedOperator_Parametric<tMesh> engOp(matprop);
    
    dumpOrthoNew(growth, tag+"_final_"+helpers::ToString(0,2));
    
    std::vector<Real> swellingRates = {0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    const size_t nSwellingRates = swellingRates.size();
    for(size_t i=0;i<nSwellingRates;++i)
    {
        // set growth
        assignGrowthToMetric(growth, swellingRates[i], interpolate_logeucl);
        
        // add noise
        addNoiseToVertices_c<2>(0.01*h);
        if(growthCase != "validation") addNoiseToEdgeDirectors(0.01*M_PI);
        
        // solve minimum
        Real eps = 1e-2;
        minimizeEnergy(engOp, eps);
        
        // dump
        dumpOrthoNew(growth, tag+"_final_"+helpers::ToString((int)i+1,2)); // i+1 since we already dumped 0
    }
    
    if(growthCase == "validation")
    {
        // store the final solution stretching and bending energies
        const Real eng_stretch = engOp.getLastStretchingEnergy();
        const Real eng_bend = engOp.getLastBendingEnergy();

        const Real Eref =  4.0*(1 - nu*nu) / M_PI;
        FILE * f = fopen((tag+"_final_energies.dat").c_str(), "w");
        fprintf(f, "# thickness \t stretch eng \t bend eng \n");
        fprintf(f, "%10.10e \t %10.10e \t %10.10e\n", h, eng_stretch * Eref / E, eng_bend * Eref / E);
        fclose(f);
    }
}


void Sim_MonoLayer_Growth::dumpWithGrowthRates(const Eigen::Ref<const Eigen::VectorXd> growthRates, const std::string filename, const bool addCurvature)
{
    WriteVTK writer(mesh.getCurrentConfiguration().getVertices(), mesh.getTopology().getFace2Vertices());
    writer.addScalarFieldToFaces(growthRates, "growthrates");
    
    // add curvature
    if(addCurvature)
    {
        const int nFaces = mesh.getNumberOfFaces();
        Eigen::VectorXd gauss(nFaces);
        Eigen::VectorXd mean(nFaces);
        ComputeCurvatures<tMesh> computeCurvatures;
        computeCurvatures.compute(mesh, gauss, mean);
        writer.addScalarFieldToFaces(gauss, "gauss");
        writer.addScalarFieldToFaces(mean, "mean");
    }
    
    writer.write(filename);
}



void Sim_MonoLayer_Growth::dumpOrthoNew(const std::vector<GrowthState> & growth, const std::string filename, const bool restConfig)
{
    const auto cvertices = restConfig ? mesh.getRestConfiguration().getVertices() : mesh.getCurrentConfiguration().getVertices();
    const auto cface2vertices = mesh.getTopology().getFace2Vertices();
    WriteVTK writer(cvertices, cface2vertices);
    
    const int nFaces = mesh.getNumberOfFaces();
    Eigen::VectorXd growthFacs_1(nFaces), growthFacs_2(nFaces);
    Eigen::MatrixXd growthDirs_1(nFaces, 3), growthDirs_2(nFaces, 3);
    for(int i=0;i<nFaces;++i)
    {
        const DecomposedGrowthState & decomposed = restConfig ? growth[i].getDecomposedInitState() : growth[i].getDecomposedFinalState();
        
        growthFacs_1(i) = decomposed.get_s1();
        growthFacs_2(i) = decomposed.get_s2();
        
        growthDirs_1.row(i) = decomposed.get_v1();
        growthDirs_2.row(i) = decomposed.get_v2();
    }
    
    writer.addScalarFieldToFaces(growthFacs_1, "rate1");
    writer.addScalarFieldToFaces(growthFacs_2, "rate2");
    writer.addVectorFieldToFaces(growthDirs_1, "dir1");
    writer.addVectorFieldToFaces(growthDirs_2, "dir2");
    
    {
        Eigen::VectorXd gauss(nFaces);
        Eigen::VectorXd mean(nFaces);
        ComputeCurvatures<tMesh> computeCurvatures;
        computeCurvatures.compute(mesh, gauss, mean);
        writer.addScalarFieldToFaces(gauss, "gauss");
        writer.addScalarFieldToFaces(mean, "mean");
    }
    
    writer.write(filename);
}

void Sim_MonoLayer_Growth::init()
{
}

