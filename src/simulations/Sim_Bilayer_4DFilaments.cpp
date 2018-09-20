//
//  Sim_Bilayer_4DFilaments.cpp
//  Elasticity
//
//  Created by Wim van Rees on 8/10/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#include "Sim_Bilayer_4DFilaments.hpp"
#include "Geometry.hpp"
#include "MaterialProperties.hpp"
#include "CombinedOperator_Parametric.hpp"
#include "EnergyOperatorList.hpp"

#include "ReadVTK.hpp"
#include "WriteSTLVolumetricCurved.hpp"

#include "SVGHelper.hpp"
#include "Geometry_Orchid.hpp"

#include "Parametrizer.hpp"
#include "HLBFGS_Wrapper_Parametrized.hpp"

template<typename tMesh>
class Parametrizer_FixedZPlane : public Parametrizer<tMesh>
{
public:
    using Parametrizer<tMesh>::mesh; // to avoid having to type this-> all the time
    using Parametrizer<tMesh>::data; // to avoid having to type this-> all the time
    
protected:
    const Real penalizationFac;
    const Real halfBufferWidth;
    const Real plane_z_value;
    const bool flip;
    
    Real computeEnergyPenalization(const Real z)
    {
        if(z > plane_z_value + halfBufferWidth)
            return flip ? penalizationFac : 0;
        if(z < plane_z_value - halfBufferWidth)
            return flip ? 0 : penalizationFac;
        
        const Real zRel = (z - (plane_z_value - halfBufferWidth))/(2*halfBufferWidth);// between 0 (z= (plane_z_value - halfBufferWidth) and 1 (z = (plane_z_value + halfBufferWidth))
        const Real moll = std::cos(M_PI*(zRel-1)); // between -1 (zRel=0) and 1 (zRel=1)
        const Real retval = 0.5*(1 - moll); // between 0 (zRel=1) and 1 (zRel=0)
        return penalizationFac*(flip ? 1 - retval : retval);
    }
    
    Real computeEnergyPenalizationDerivative(const Real z)
    {
        if(z > (plane_z_value + halfBufferWidth) or z < (plane_z_value - halfBufferWidth)) return 0;
        
        const Real zRel = (z - (plane_z_value - halfBufferWidth))/(2*halfBufferWidth);// between 0 (z= (plane_z_value - halfBufferWidth) and 1 (z = (plane_z_value + halfBufferWidth))
        const Real d_zRel = 1.0/(2*halfBufferWidth);
        const Real d_moll = -M_PI*std::sin(M_PI*(zRel-1))*d_zRel;
        const Real d_retval = -0.5*d_moll;
        return penalizationFac*(flip ? -1 : +1)*d_retval;
    }
    
    public:
    
    Parametrizer_FixedZPlane(tMesh & mesh_in, const Real penal_fac_in = 10.0, const Real halfwidth_in = 0.1, const Real z_val_in = 0.0, const bool flip_in = false):
    Parametrizer<tMesh>(mesh_in),
    penalizationFac(penal_fac_in),
    halfBufferWidth(halfwidth_in),
    plane_z_value(z_val_in),
    flip(flip_in)
    {
        assert(halfBufferWidth > 0);
    }
    
    void initSolution(const Eigen::Ref<const Eigen::VectorXd> , const bool ) override
    {
        // do nothing : we are use the mesh data
    }
    
    Real * getDataPointer() const override
    {
        return mesh.getDataPointer();
    }
    
    int getNumberOfVariables() const override
    {
        return 3*mesh.getNumberOfVertices() + mesh.getNumberOfEdges(); // DCS energy minimization
    }
    
    void updateSolution() override
    {
        mesh.updateDeformedConfiguration();
    }
    
    Real computeEnergyContribution() override
    {
        // here we have a function that penalizes the vertex to go to the negative z-plane
        const int nVertices = mesh.getNumberOfVertices();
        const auto vertices = mesh.getCurrentConfiguration().getVertices();
        Real retval = 0.0;
        for(int i=0;i<nVertices;++i)
            retval += computeEnergyPenalization(vertices(i,2));
        
        return retval;
    }
    
    
    void updateGradient(const int nVars, const Eigen::Ref<const Eigen::VectorXd> energyGradient, Real * const grad_ptr) override
    {
        assert(nVars == getNumberOfVariables());
        
        // set grad_ptr equal to the actual gradient
        for(int i=0;i<nVars;++i)
            grad_ptr[i] = energyGradient(i);
        
        // wrapper around the solution of the z-coordinates of the vertices
        const int nVertices = mesh.getNumberOfVertices();
        Eigen::Map<Eigen::VectorXd> gradVertices_z(grad_ptr + 2*nVertices, nVertices);
        const auto vertices = mesh.getCurrentConfiguration().getVertices();
        
        // add the contribution from the mollifier
        for(int i=0;i<nVertices;++i)
            gradVertices_z(i) += computeEnergyPenalizationDerivative(vertices(i,2));
    }
};


void Sim_Bilayer_4DFilaments::init()
{
    // do nothing until we run something
}

void Sim_Bilayer_4DFilaments::run()
{
    const std::string runCase = parser.parse<std::string>("-case", "");
    
    if(runCase == "vtu2stl")
        vtu_to_stl();
    else if(runCase == "catenoid_printpath")
        run_catenoid_printpath();
    else if(runCase == "helicoid_printpath")
        run_helicoid_printpath();
    else if(runCase == "sombrero_printpath")
        run_sombrero_printpath();
    else if(runCase == "logspiral_printpath")
        run_logspiral_printpath();
    else if(runCase == "folding_flower_printpath")
        run_folding_flower_printpath();
    else if(runCase == "orchid_printpath")
        run_orchid_printpath();
    else
    {
        std::cout << "No valid test case defined. Options are \n";
        std::cout << "\t -case vtu2stl\n";
        
        std::cout << "\t -case catenoid_printpath\n";
        std::cout << "\t -case helicoid_printpath\n";
        std::cout << "\t -case sombrero_printpath\n";
        std::cout << "\t -case logspiral_printpath\n";
        std::cout << "\t -case folding_flower_printpath\n";
        std::cout << "\t -case orchid_printpath\n";
        
    }
}

void Sim_Bilayer_4DFilaments::get_density_from_svg(const DensityInputs & densityInputs, Eigen::VectorXd & density_bot, Eigen::VectorXd & density_top, Eigen::MatrixXd & dir_bot, Eigen::MatrixXd & dir_top)
{
    SVGHelper::DiscreteSVG svg_bot(densityInputs.filename_bot, densityInputs.dpath);
    SVGHelper::DiscreteSVG svg_top(densityInputs.filename_top, densityInputs.dpath);

    {
        const std::array<Real, 4> global_extents_bot = svg_bot.getGlobalExtents();
        const std::array<Real, 4> global_extents_top = svg_top.getGlobalExtents();
        printf(" ---- raw input file ---- \n");
        printf("bot : x ranges from %10.10e to %10.10e \t y ranges from %10.10e to %10.10e\n", global_extents_bot[0], global_extents_bot[1], global_extents_bot[2], global_extents_bot[3]);
        printf("top : x ranges from %10.10e to %10.10e \t y ranges from %10.10e to %10.10e\n", global_extents_top[0], global_extents_top[1], global_extents_top[2], global_extents_top[3]);
    }

    
    // rescale
    svg_bot.translateRescalePaths(densityInputs.origin_X, densityInputs.origin_Y, densityInputs.target_size_X, densityInputs.target_size_Y);
    svg_top.translateRescalePaths(densityInputs.origin_X, densityInputs.origin_Y, densityInputs.target_size_X, densityInputs.target_size_Y);
    
    {
        const std::array<Real, 4> global_extents_bot = svg_bot.getGlobalExtents();
        const std::array<Real, 4> global_extents_top = svg_top.getGlobalExtents();
        printf(" ---- after translating/rescaling ---- \n");
        printf("bot : x ranges from %10.10e to %10.10e \t y ranges from %10.10e to %10.10e\n", global_extents_bot[0], global_extents_bot[1], global_extents_bot[2], global_extents_bot[3]);
        printf("top : x ranges from %10.10e to %10.10e \t y ranges from %10.10e to %10.10e\n", global_extents_top[0], global_extents_top[1], global_extents_top[2], global_extents_top[3]);
    }

    
    
    // clip
    svg_bot.clipAllPoints({densityInputs.min_extent_X, densityInputs.max_extent_X},{densityInputs.min_extent_Y, densityInputs.max_extent_Y});
    svg_top.clipAllPoints({densityInputs.min_extent_X, densityInputs.max_extent_X},{densityInputs.min_extent_Y, densityInputs.max_extent_Y});
    
    {
        printf(" ---- after clipping ---- \n");
        const std::array<Real, 4> global_extents_bot = svg_bot.getGlobalExtents();
        const std::array<Real, 4> global_extents_top = svg_top.getGlobalExtents();
        
        printf("bot : x ranges from %10.10e to %10.10e \t y ranges from %10.10e to %10.10e\n", global_extents_bot[0], global_extents_bot[1], global_extents_bot[2], global_extents_bot[3]);
        printf("top : x ranges from %10.10e to %10.10e \t y ranges from %10.10e to %10.10e\n", global_extents_top[0], global_extents_top[1], global_extents_top[2], global_extents_top[3]);
    }

    // fill the arrays
    const int nFaces = mesh.getNumberOfFaces();
    density_bot.resize(nFaces);
    density_top.resize(nFaces);
    dir_bot.resize(nFaces, 3);
    dir_top.resize(nFaces, 3);
    

    for(int i=0;i<nFaces;++i)
    {
        const TriangleInfo info = mesh.getRestConfiguration().getTriangleInfoLite(mesh.getTopology(), i);
        const Eigen::Vector3d restfacecenter = info.computeFaceCenter();
        const Eigen::Vector2d pos = (Eigen::Vector2d() << restfacecenter(0), restfacecenter(1)).finished();
                
        const auto info_bot = svg_bot.getInfoForPoint_Filament(pos, densityInputs.sigma, densityInputs.dropSize, densityInputs.dropFac);
        const auto info_top = svg_top.getInfoForPoint_Filament(pos, densityInputs.sigma, densityInputs.dropSize, densityInputs.dropFac);

        // rescale density using filament_thickness / sigma (and cap at 1)
        density_bot(i) = std::min(1.0, info_bot.first * densityInputs.filament_thickness / densityInputs.sigma);
        density_top(i) = std::min(1.0, info_top.first * densityInputs.filament_thickness / densityInputs.sigma);
        for(int d=0;d<3;++d)
        {
            dir_bot(i,d) = (d==2 ? 0.0 : info_bot.second(d));
            dir_top(i,d) = (d==2 ? 0.0 : info_top.second(d));
        }
    }
    
}

void Sim_Bilayer_4DFilaments::run_catenoid_printpath()
{
    tag = "catenoid_printpath";
    
    const Real Lx = 4;
    const Real Ly = 20;
    const Real relArea = parser.parse<Real>("-res", 0.0005);
    
    // initialize geometry
    RectangularPlate geometry(Lx, Ly, relArea, {false,false}, {false,false});
    mesh.init(geometry, false);
    
    // get densities
    const Real filament_thickness = parser.parse<Real>("-filthickness", 0.250);
    const Real max_filament_spacing = 0.75;

    const Real sigma_prefac = parser.parse<Real>("-sigmafac", 1.25);
    const Real sigma = parser.parse<Real>("-sigma", sigma_prefac*max_filament_spacing); // use sigma = 2 * distance
    const Real dpath = parser.parse<Real>("-dpath", 0.5);
    
    const std::string filename_bot = parser.parse<std::string>("-fname_bot", "/Users/wvanrees/Documents/project_2/plates/4D_printing/path2rho/catenoid_vec_1.dat");
    const std::string filename_top = parser.parse<std::string>("-fname_top", "/Users/wvanrees/Documents/project_2/plates/4D_printing/path2rho/catenoid_vec_2.dat");
    
    DensityInputs densityInputs(filename_bot, filename_top, sigma, filament_thickness, dpath);
    
    // set up the specifics for this geometry that are non-default
    densityInputs.target_size_X = parser.parse<Real>("-rescaleX", 8.);
    densityInputs.min_extent_X = parser.parse<Real>("-minX", -4.);
    densityInputs.max_extent_X = parser.parse<Real>("-maxX", +4.);
    densityInputs.min_extent_Y = parser.parse<Real>("-minY", -20.);
    densityInputs.max_extent_Y = parser.parse<Real>("-maxY", +20.);
    densityInputs.dropSize = parser.parse<Real>("-dropsize", 1.0);
    densityInputs.dropFac = parser.parse<Real>("-dropfac", 1.0);
    
    Eigen::VectorXd density_bot, density_top;
    Eigen::MatrixXd dir_bot, dir_top;
    get_density_from_svg(densityInputs, density_bot, density_top, dir_bot, dir_top);
    setup_general_printpath(filament_thickness, density_bot, density_top, dir_bot, dir_top);
}


void Sim_Bilayer_4DFilaments::run_helicoid_printpath()
{
    tag = "helicoid_printpath";
    
    const Real Lx = 4;
    const Real Ly = 20;
    const Real relArea = parser.parse<Real>("-res", 0.0005);
    
    // initialize geometry
    RectangularPlate geometry(Lx, Ly, relArea, {false,false}, {false,false});
    mesh.init(geometry, false);
    
    // get densities
    const Real filament_thickness = parser.parse<Real>("-filthickness", 0.4);
    const Real max_filament_spacing = 0.75;

    const Real sigma_prefac = parser.parse<Real>("-sigmafac", 1.25);
    const Real sigma = parser.parse<Real>("-sigma", sigma_prefac*max_filament_spacing);// use sigma = 2 * distance
    const Real dpath = parser.parse<Real>("-dpath", 0.5);
    
    const std::string filename_bot = parser.parse<std::string>("-fname_bot", "/Users/wvanrees/Documents/project_2/plates/4D_printing/path2rho/helicoid_vec_1.dat");
    const std::string filename_top = parser.parse<std::string>("-fname_top", "/Users/wvanrees/Documents/project_2/plates/4D_printing/path2rho/helicoid_vec_2.dat");
    
    DensityInputs densityInputs(filename_bot, filename_top, sigma, filament_thickness, dpath);
    
    // set up the specifics for this geometry that are non-default
    densityInputs.target_size_X = parser.parse<Real>("-rescaleX", 8.);
    densityInputs.min_extent_X = parser.parse<Real>("-minX", -4.);
    densityInputs.max_extent_X = parser.parse<Real>("-maxX", +4.);
    densityInputs.min_extent_Y = parser.parse<Real>("-minY", -20.);
    densityInputs.max_extent_Y = parser.parse<Real>("-maxY", +20.);
    densityInputs.dropSize = parser.parse<Real>("-dropsize", 1.0);
    densityInputs.dropFac = parser.parse<Real>("-dropfac", 1.0);
    
    Eigen::VectorXd density_bot, density_top;
    Eigen::MatrixXd dir_bot, dir_top;
    get_density_from_svg(densityInputs, density_bot, density_top, dir_bot, dir_top);
    setup_general_printpath(filament_thickness, density_bot, density_top, dir_bot, dir_top);
}

void Sim_Bilayer_4DFilaments::run_sombrero_printpath()
{
    tag = "sombrero_printpath";
    
    const Real rad = 7.5;
    const int nPts = parser.parse<int>("-res", 256);
    
    // initialize geometry
    CircularPlate geometry(rad, nPts, false);
    mesh.init(geometry, false);
    
    // get densities
    const Real filament_thickness = 0.250;
    const Real max_filament_spacing = (2*M_PI*rad)/(2.*20); // 2*20 spokes around perimeter
    
    const Real sigma_prefac = parser.parse<Real>("-sigmafac", 1.25);
    const Real sigma = parser.parse<Real>("-sigma", sigma_prefac*max_filament_spacing);// use sigma = 2 * distance
    const Real dpath = parser.parse<Real>("-dpath", 0.5);
    
    const std::string filename_bot = parser.parse<std::string>("-fname_bot", "/Users/wvanrees/Documents/project_2/plates/4D_printing/path2rho/sombrero_vec_1.dat");
    const std::string filename_top = parser.parse<std::string>("-fname_top", "/Users/wvanrees/Documents/project_2/plates/4D_printing/path2rho/sombrero_vec_2.dat");
    
    DensityInputs densityInputs(filename_bot, filename_top, sigma, filament_thickness, dpath);
    
    // set up the specifics for this geometry that are non-default
    densityInputs.target_size_Y = parser.parse<Real>("-rescaleY", 15.);
    densityInputs.dropSize = parser.parse<Real>("-dropsize", 1.0);
    densityInputs.dropFac = parser.parse<Real>("-dropfac", 1.0);
    
    Eigen::VectorXd density_bot, density_top;
    Eigen::MatrixXd dir_bot, dir_top;
    get_density_from_svg(densityInputs, density_bot, density_top, dir_bot, dir_top);
    
    {
        const int nFaces = mesh.getNumberOfFaces();
        const Real tmpFac1 = 2.0;
        const Real tmpFac2 = 5.0;
        for(int i=0;i<nFaces;++i)
        {
            const TriangleInfo rinfo = mesh.getRestConfiguration().getTriangleInfoLite(mesh.getTopology(), i);
            Eigen::MatrixXd rxy_base(3,2);
            rxy_base << rinfo.e1, rinfo.e2;
            
            const Eigen::Vector3d v1_top = dir_top.row(i);
            const Eigen::Vector3d v2_top = (Eigen::Vector3d() << v1_top(1), -v1_top(0), 0.0).finished();
            
            DecomposedGrowthState growthState(tmpFac1, tmpFac2, v1_top, v2_top);
            const Eigen::Matrix2d metric = growthState.computeMetric(rxy_base);
            const Eigen::Matrix2d metric_2 = (metric.log()).exp();
            DecomposedGrowthState growthState_2(rxy_base, metric_2);
            
            const Real s1_2 = growthState_2.get_s1();
            const Real s2_2 = growthState_2.get_s2();
            const Eigen::Vector3d v1_2 = growthState_2.get_v1();
            const Eigen::Vector3d v2_2 = growthState_2.get_v2();
            
            if(std::abs(s1_2 - tmpFac1) < std::abs(s2_2 - tmpFac1))
            {
                dir_top.row(i) = v1_2;
            }
            else
            {
                dir_top.row(i) = v2_2;
            }
        }            
    }
    setup_general_printpath(filament_thickness, density_bot, density_top, dir_bot, dir_top);
}


void Sim_Bilayer_4DFilaments::run_logspiral_printpath()
{
    tag = "logspiral_printpath";
    
    const Real Lx = 0.5 * (0.25 * std::pow(53, 1.4));
    const Real Ly = 0.5 * (0.75 * 7);
    
    const Real relArea = parser.parse<Real>("-res", 0.0005);
    
    // initialize geometry
    RectangularPlate geometry(Lx, Ly, relArea, {false,false}, {false,false});
    mesh.init(geometry, false);
    
    // get densities
    const Real filament_thickness = parser.parse<Real>("-filthickness",0.400); // email sabetta Aug 11, 2016 -- 400µm nozzle diameter
    const Real max_filament_spacing = 0.25*(std::pow(53, 1.4) - std::pow(52, 1.4)); // use distance between last and second-last filament
    
    const Real sigma_prefac = parser.parse<Real>("-sigmafac", 1.25);
    const Real sigma = parser.parse<Real>("-sigma", sigma_prefac*max_filament_spacing); // use sigma = 2 * distance
    const Real dpath = parser.parse<Real>("-dpath", 0.5);
    
    const std::string filename_bot = parser.parse<std::string>("-fname_bot", "/Users/wvanrees/Documents/project_2/plates/4D_printing/path2rho/logspiral_vec_1.dat");
    const std::string filename_top = parser.parse<std::string>("-fname_top", "/Users/wvanrees/Documents/project_2/plates/4D_printing/path2rho/logspiral_vec_2.dat");
    
    DensityInputs densityInputs(filename_bot, filename_top, sigma, filament_thickness, dpath);
    
    // set up the specifics for this geometry that are non-default
    densityInputs.target_size_Y = parser.parse<Real>("-rescaleY", 5.25);
    densityInputs.dropSize = parser.parse<Real>("-dropsize", 1.0);
    densityInputs.dropFac = parser.parse<Real>("-dropfac", 1.0);
    
    Eigen::VectorXd density_bot, density_top;
    Eigen::MatrixXd dir_bot, dir_top;
    get_density_from_svg(densityInputs, density_bot, density_top, dir_bot, dir_top);
    
    setup_general_printpath(filament_thickness, density_bot, density_top, dir_bot, dir_top);
}


void Sim_Bilayer_4DFilaments::run_folding_flower_printpath()
{
    tag = "folding_flower_printpath";
    
    // initialize geometry
    const Real flower_h = 15;//13.5 // 20.0
    const Real flower_r0 = 2.5;
    const Real flower_dy = 1.5;
    const Real edgeLength = parser.parse<Real>("-res", 0.25);
    Geometry_FoldingFlower geometry(flower_h, flower_r0, flower_dy, edgeLength);
    mesh.init(geometry, false);
    
    // get densities
    const Real filament_thickness = parser.parse<Real>("-filthickness", 0.400); // email sabetta Oct 3rd, 2016 -- 400µm nozzle diameter
    const Real max_filament_spacing = parser.parse<Real>("-maxfilspacing", 0.7);//0.5
    
    const Real sigma_prefac = parser.parse<Real>("-sigmafac", 1.25);
    const Real sigma = parser.parse<Real>("-sigma", sigma_prefac*max_filament_spacing);
    const Real dpath = parser.parse<Real>("-dpath", 0.5);
    const std::string filename_bot = parser.parse<std::string>("-fname_bot", "/Users/wvanrees/Documents/project_2/plates/4D_printing/path2rho/folding_flower_vec_1.dat");
    const std::string filename_top = parser.parse<std::string>("-fname_top", "/Users/wvanrees/Documents/project_2/plates/4D_printing/path2rho/folding_flower_vec_2.dat");
    
    DensityInputs densityInputs(filename_bot, filename_top, sigma, filament_thickness, dpath);
    
    // set up the specifics for this geometry that are non-default
    densityInputs.target_size_Y = parser.parse<Real>("-rescaleY", (flower_h + flower_dy)*(1 + std::cos(M_PI/5)));
    densityInputs.origin_Y = parser.parse<Real>("-originY", 15);
    densityInputs.dropSize = parser.parse<Real>("-dropsize", 1.0);
    densityInputs.dropFac = parser.parse<Real>("-dropfac", 1.0);
    
    Eigen::VectorXd density_bot, density_top;
    Eigen::MatrixXd dir_bot, dir_top;
    get_density_from_svg(densityInputs, density_bot, density_top, dir_bot, dir_top);
    
    setup_general_printpath(filament_thickness, density_bot, density_top, dir_bot, dir_top);
}


void Sim_Bilayer_4DFilaments::run_orchid_printpath()
{
    tag = "orchid_printpath";
    
    // initialize geometry
    const Real edgeLength = parser.parse<Real>("-res", 0.1);
    Geometry_Orchid geometry(edgeLength);
    mesh.init(geometry, false);

    // set up the basic material properties
    const Real filament_thickness = parser.parse<Real>("-filthickness", 0.325);
    const Real max_filament_spacing = 1.20;
    
    const Real sigma_prefac = parser.parse<Real>("-sigmafac", 1.25);
    const Real sigma = parser.parse<Real>("-sigma", sigma_prefac*max_filament_spacing); // 1.2 is pretty large - make it 1.25x so that it becomes 1.5 -- just enough to cover the central hole in the top layer
    const Real dpath = parser.parse<Real>("-dpath", 0.5);
    
    const std::string filename_bot = parser.parse<std::string>("-fname_bot", "/Users/wvanrees/Documents/project_2/plates/4D_printing/path2rho/orchid_vec_flipshort_1.dat");
    const std::string filename_top = parser.parse<std::string>("-fname_top", "/Users/wvanrees/Documents/project_2/plates/4D_printing/path2rho/orchid_vec_flipshort_2.dat");
    
    DensityInputs densityInputs(filename_bot, filename_top, sigma, filament_thickness, dpath);
    
    // set up the specifics for this geometry that are non-default
    densityInputs.target_size_Y = parser.parse<Real>("-rescaleY", 34.3);
    densityInputs.origin_Y = parser.parse<Real>("-originY", 15.36);
    densityInputs.dropSize = parser.parse<Real>("-dropsize", 1.0);
    densityInputs.dropFac = parser.parse<Real>("-dropfac", 1.0);
    
    Eigen::VectorXd density_bot, density_top;
    Eigen::MatrixXd dir_bot, dir_top;
    get_density_from_svg(densityInputs, density_bot, density_top, dir_bot, dir_top);
    
    // add the plane
    const bool simulatePlane = parser.parse<bool>("-simulateplane", true);
    if(simulatePlane)
    {
        // set up boundary conditions : fix the closest point to the origin
        std::pair<Real, int> minDistVertex = {1e9, -1};
        const auto vertices = mesh.getRestConfiguration().getVertices();
        const int nVertices = mesh.getNumberOfVertices();
        for(int i=0;i<nVertices;++i)
        {
            const Real minDist = vertices.row(i).norm();
            if(minDist < minDistVertex.first)
                minDistVertex = std::make_pair(minDist, i);
        }
        
        auto vertexBoundaryConditions = mesh.getBoundaryConditions().getVertexBoundaryConditions();
        for(int d=0;d<3;++d)
            vertexBoundaryConditions(minDistVertex.second,d) = true;
    }
    
    setup_general_printpath(filament_thickness, density_bot, density_top, dir_bot, dir_top, simulatePlane);
}


void Sim_Bilayer_4DFilaments::setup_general_printpath(const Real filament_thickness, const Eigen::Ref<const Eigen::VectorXd> density_bot, const Eigen::Ref<const Eigen::VectorXd> density_top, const Eigen::Ref<const Eigen::MatrixXd> dir_bot, const Eigen::Ref<const Eigen::MatrixXd> dir_top, const bool simulatePlane)
{
    printf("Done with creating growth profile - setting up general simulation\n");
    
    const int nFaces = mesh.getNumberOfFaces();
    const Real Young = parser.parse<Real>("-E", 1e3);
    const Real poisson = parser.parse<Real>("-nu", 0.4);
    
    // bilayer thickness = 2xfilament thickness
    const Real botSwellFac = parser.parse<Real>("-botswellfac", 1.0);
    const Real topSwellFac = parser.parse<Real>("-topswellfac", 1.0);
    FilamentLayerConfig<bottom> botlayer(nFaces, Young, poisson, 2*filament_thickness, botSwellFac);
    FilamentLayerConfig<top> toplayer(nFaces, Young, poisson, 2*filament_thickness, topSwellFac);
    botlayer.setFilamentDensity(density_bot);
    toplayer.setFilamentDensity(density_top);
    botlayer.setFilamentDirectionFromMatrix(dir_bot);
    toplayer.setFilamentDirectionFromMatrix(dir_top);
    
    std::vector<Real> swellingRates = {0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    const bool swellThickness = parser.parse<bool>("-swellthickness", true);
    
    runSwelling(botlayer, toplayer, swellingRates, swellThickness, simulatePlane);
}


void Sim_Bilayer_4DFilaments::runSwelling(const FilamentLayerConfig<bottom> &botlayer, const FilamentLayerConfig<top> &toplayer, const std::vector<Real> &swellrates_in, const bool swellThickness, const bool simulatePlane)
{
    const std::string basetag = tag;

    const bool doOrtho = parser.parse<bool>("-ortho", false);
    Real ortho_fac = -1, ortho_G12_fac = -1;
    if(doOrtho)
    {
        ortho_fac = parser.parse<Real>("-orthofac", 0.5); // Young_transverse  = ortho_fac * Young_parallel = ortho_fac * material_*.Young;
        ortho_G12_fac = parser.parse<Real>("-orthoshearfac", 0.25); // G_12 = ortho_G12_fac * Young_parallel
    }
    
    const bool computeEnergyOnly = doOrtho ? false : parser.parse<bool>("-compute_energy", false); // cant do it if orthotropic elasticity
    
    const std::string subtag = computeEnergyOnly ? "energy" : "growth";
    
    // these guys are fixed
    const Real growthFac_orthogon_final = 0.26;
    const Real growthFac_parallel_final = 0.12;

    {
        FILE * f = fopen((basetag+"_"+subtag+"_"+"result.dat").c_str(),"w");
        fclose(f);
    }
    
    // dump initial
    dumpAll(botlayer, toplayer, basetag+"_"+subtag+"_"+helpers::ToString(0, 3));
    
    bool planeFlip = false;
    Real planePenalizationFac = -1;
    if(simulatePlane)
    {
        planeFlip = parser.parse<bool>("-planeflip", planeFlip);
        planePenalizationFac = 10*botlayer.Young;
    }
    
    std::vector<Real> swellrates = swellrates_in;
    if(computeEnergyOnly)
    {
        const std::string fname_energy = parser.parse<std::string>("-fname_energy");
        const Real swellrate_energy = parser.parse<Real>("-swellrate_energy");
        
        // override the vector with only one entry
        swellrates = {swellrate_energy};
        
        // initialize the vertices
        IOGeometry geometry(fname_energy);
        Eigen::MatrixXd vertices_energy;
        Eigen::MatrixXi face2vertices_energy;
        Eigen::MatrixXb vertices_bc_energy;
        geometry.get(vertices_energy, face2vertices_energy, vertices_bc_energy);
        auto vertices = mesh.getCurrentConfiguration().getVertices();
        vertices = vertices_energy;
        mesh.updateDeformedConfiguration();
    }
    
    const int nGrowth = (int)swellrates.size();
    for(int g=0;g<nGrowth;++g)
    {
        tag = basetag+"_"+subtag+"_"+helpers::ToString(g+1, 3);
        
        const Real growthFac_orthogon = swellrates[g] * growthFac_orthogon_final;
        const Real growthFac_parallel = swellrates[g] * growthFac_parallel_final;
        
        // apply growth to mesh
        botlayer.applyGrowthToMesh(mesh, growthFac_parallel, growthFac_orthogon);
        toplayer.applyGrowthToMesh(mesh, growthFac_parallel, growthFac_orthogon);
        
        // get material properties and energy operator
        Real eps = -1;
        if(not doOrtho)
        {
            MaterialProperties_Iso_Array material_bot = swellThickness ? botlayer.getMaterialProperties(botlayer.swellFac * growthFac_orthogon) : botlayer.getMaterialProperties();
            MaterialProperties_Iso_Array material_top = swellThickness ? toplayer.getMaterialProperties(toplayer.swellFac * growthFac_orthogon) : toplayer.getMaterialProperties();
            
            // create energy operator
            CombinedOperator_Parametric<tMesh, Material_Isotropic, bottom> combined_bot(material_bot);
            CombinedOperator_Parametric<tMesh, Material_Isotropic, top> combined_top(material_top);
            EnergyOperatorList<tMesh> engOp({&combined_bot, &combined_top});
            
            eps = minimizePlateEnergy(engOp, simulatePlane, planeFlip, planePenalizationFac);
            
            dumpBilayerEnergyDecomposition(g, swellrates[g], material_bot, material_top, combined_bot, combined_top, (g==0 ? true : false), basetag+"_"+subtag+"_energies_"+helpers::ToString(g+1, 3));
        }
        else
        {
            MaterialProperties_Ortho_Array material_bot = swellThickness ? botlayer.getMaterialPropertiesOrtho(mesh, ortho_fac, ortho_G12_fac, botlayer.swellFac * growthFac_orthogon) : botlayer.getMaterialPropertiesOrtho(mesh, ortho_fac, ortho_G12_fac);
            MaterialProperties_Ortho_Array material_top = swellThickness ? toplayer.getMaterialPropertiesOrtho(mesh, ortho_fac, ortho_G12_fac, toplayer.swellFac * growthFac_orthogon) : toplayer.getMaterialPropertiesOrtho(mesh, ortho_fac, ortho_G12_fac);
            
            // create energy operator
            CombinedOperator_Parametric<tMesh, Material_Orthotropic, bottom> combined_bot(material_bot);
            CombinedOperator_Parametric<tMesh, Material_Orthotropic, top> combined_top(material_top);
            EnergyOperatorList<tMesh> engOp({&combined_bot, &combined_top});
            
            eps = minimizePlateEnergy(engOp, simulatePlane, planeFlip, planePenalizationFac);
        }
        
        // dump result
        dumpAll(botlayer, toplayer, tag);
        {
            const Real thicknessFac_bot = swellThickness ? (1.0 + botlayer.swellFac * growthFac_orthogon)*botlayer.thickness : botlayer.thickness;
            const Real thicknessFac_top = swellThickness ? (1.0 + toplayer.swellFac * growthFac_orthogon)*toplayer.thickness : toplayer.thickness;
            
            FILE * f = fopen((basetag+"_"+subtag+"_"+"result.dat").c_str(),"a");
            fprintf(f, "%d \t %10.10e \t %10.10e \t %10.10e \t %10.10e \t %10.10e \t %10.10e\n",g, swellrates[g], thicknessFac_bot, thicknessFac_top, botlayer.swellFac, toplayer.swellFac, eps);
            fclose(f);
        }
    }
    
    // reset tag
    tag = basetag;
}


void Sim_Bilayer_4DFilaments::dumpBilayerEnergyDecomposition(const int idx, const Real swelling_fac, const MaterialProperties_Iso_Array & matprop_bot, const MaterialProperties_Iso_Array & matprop_top, const CombinedOperator_Parametric<tMesh, Material_Isotropic, bottom> & engOp_bot, const CombinedOperator_Parametric<tMesh, Material_Isotropic, top> & engOp_top, const bool overWrite, const std::string filename)
{
    const int nFaces = mesh.getNumberOfFaces();
    
    Eigen::MatrixXd energy_bi_bot(nFaces,3);
    Eigen::MatrixXd energy_bi_top(nFaces,3);
    energy_bi_bot.setZero();
    energy_bi_top.setZero();
    const Real eng_bi_bot = engOp_bot.computePerFaceEnergies(mesh, energy_bi_bot);
    const Real eng_bi_top = engOp_top.computePerFaceEnergies(mesh, energy_bi_top);
    
    
    Mesh mesh_tmp;
    Geometry_Dummy geometry_dummy(mesh.getCurrentConfiguration().getVertices(), mesh.getTopology().getFace2Vertices());
    mesh_tmp.init(geometry_dummy);
    // copy edge directors of current configuration
    mesh_tmp.getCurrentConfiguration().getEdgeDirectors() = mesh.getCurrentConfiguration().getEdgeDirectors();
    mesh_tmp.updateDeformedConfiguration();
    
    // create new material properties
    Eigen::VectorXd thickness_bot(nFaces), thickness_top(nFaces);
    Eigen::VectorXd thickness_mono(nFaces), Young_mono(nFaces);
    const Real nu = matprop_bot.getFaceMaterial(0).Poisson;
    for(int i=0;i<nFaces;++i)
    {
        Young_mono(i) = matprop_bot.getFaceMaterial(i).getYoung(); // should not matter bot or top -- both the same
        thickness_bot(i) = matprop_bot.getFaceMaterial(i).getThickness();
        thickness_top(i) = matprop_top.getFaceMaterial(i).getThickness();
        
        thickness_mono(i) = 0.5*(matprop_bot.getFaceMaterial(i).getThickness() + matprop_top.getFaceMaterial(i).getThickness());
    }
    MaterialProperties_Iso_Array matprop_mono(Young_mono, nu, thickness_mono);
    
    // get new enegy operator
    Eigen::MatrixXd energy_mono(nFaces,3);
    energy_mono.setZero();
    CombinedOperator_Parametric<Mesh, Material_Isotropic> engOp_mono(matprop_mono);
    
    // set metrics
    tVecMat2d & aform_tmp = mesh_tmp.getRestConfiguration().getFirstFundamentalForms();
    tVecMat2d & bform_tmp = mesh_tmp.getRestConfiguration().getSecondFundamentalForms();
    const tVecMat2d & aform_bot = mesh.getRestConfiguration().getFirstFundamentalForms<bottom>();
    const tVecMat2d & aform_top = mesh.getRestConfiguration().getFirstFundamentalForms<top>();
    
    for(int i=0;i<nFaces;++i)
    {
        const Real h1 = matprop_bot.getFaceMaterial(i).getThickness();
        const Real h2 = matprop_top.getFaceMaterial(i).getThickness();
        const Real abot_prefac = h1*(4*h2*h2 - h1*h2 + h1*h1)/std::pow(h1 + h2,3);
        const Real atop_prefac = h2*(4*h1*h1 - h1*h2 + h2*h2)/std::pow(h1 + h2,3);
        const Real b_prefac = 6*h1*h2/std::pow(h1 + h2,3);
        
        aform_tmp[i] = abot_prefac*aform_bot[i] + atop_prefac*aform_top[i];
        bform_tmp[i] = b_prefac * (aform_bot[i] - aform_top[i]);
    }
    const Real eng_mono = engOp_mono.computePerFaceEnergies(mesh_tmp, energy_mono);
    
    // compute bilayer constant term
    Eigen::VectorXd energy_constant(nFaces);
    energy_constant.setZero();
    Real eng_bi_const = 0.0;
    {
        for(int i=0;i<nFaces;++i)
        {
            const Real h1 = matprop_bot.getFaceMaterial(i).getThickness();
            const Real h2 = matprop_top.getFaceMaterial(i).getThickness();
            const Real area = 0.5*std::sqrt(aform_tmp[i].determinant());
            const Eigen::Matrix2d shapeOp = aform_tmp[i].inverse() * bform_tmp[i];
            const Real trace = shapeOp.trace();
            const Real trace_sq = (shapeOp * shapeOp).trace();
            
            const Real fac1 = matprop_mono.getFaceMaterial(i).getStVenantFactor1();
            const Real fac2 = matprop_mono.getFaceMaterial(i).getStVenantFactor2();
            
            const Real cfac = std::pow(h1 + h2,3) * (h2*h2 - h1*h2 + h1*h1)/(228*h1*h2);
            const Real eng = cfac * (fac1*trace*trace + fac2*trace_sq)*area; // factor one-half already included in StVenantFactors!!
            energy_constant(i) = eng;
            eng_bi_const += eng;
        }
    }
    std::cout << "ENERGY BI CONSTANT = " << eng_bi_const << "\t" << energy_constant.sum() << std::endl;
    for(int i=0;i<nFaces;++i)
    {
        const Real h1 = matprop_bot.getFaceMaterial(i).getThickness();
        const Real h2 = matprop_top.getFaceMaterial(i).getThickness();
        const Real abot_prefac = h1*(4*h2*h2 - h1*h2 + h1*h1)/std::pow(h1 + h2,3);
        const Real atop_prefac = h2*(4*h1*h1 - h1*h2 + h2*h2)/std::pow(h1 + h2,3);
        const Real b_prefac = 6*h1*h2/std::pow(h1 + h2,3);
        
        aform_tmp[i] = abot_prefac*aform_bot[i] + atop_prefac*aform_top[i];
        bform_tmp[i] = -b_prefac * (aform_bot[i] - aform_top[i]);
    }
    const Real eng_mono_flip = engOp_mono.compute(mesh_tmp);
    
    // compute bilayer constant term
    Real eng_bi_const_flip = 0.0;
    {
        for(int i=0;i<nFaces;++i)
        {
            const Real h1 = matprop_bot.getFaceMaterial(i).getThickness();
            const Real h2 = matprop_top.getFaceMaterial(i).getThickness();
            const Real area = 0.5*std::sqrt(aform_tmp[i].determinant());
            const Eigen::Matrix2d shapeOp = aform_tmp[i].inverse() * bform_tmp[i];
            const Real trace = shapeOp.trace();
            const Real trace_sq = (shapeOp * shapeOp).trace();
            
            const Real fac1 = matprop_mono.getFaceMaterial(i).getStVenantFactor1();
            const Real fac2 = matprop_mono.getFaceMaterial(i).getStVenantFactor2();
            
            const Real cfac = std::pow(h1 + h2,3) * (h2*h2 - h1*h2 + h1*h1)/(228*h1*h2);
            const Real eng = cfac * (fac1*trace*trace + fac2*trace_sq)*area; // factor one-half already included in StVenantFactors!!
            eng_bi_const_flip += eng;
        }
    }
    
    
    FILE * f = fopen((tag+"_bilayer_energies.dat").c_str(), (overWrite ? "w" : "a"));
    if(overWrite) fprintf(f, "# idx \t swelling_fac \t eng_bi_bot \t eng_bi_top \t eng_bi_total \t eng_bi_const \t eng_bi_const_flip \t eng_mono \t eng_mono_flip\n");
    fprintf(f, "%d \t %10.10e \t %10.10e \t %10.10e \t %10.10e \t %10.10e \t %10.10e \t %10.10e \t %10.10e\n", idx, swelling_fac, eng_bi_bot, eng_bi_top, eng_bi_bot + eng_bi_top, eng_bi_const, eng_bi_const_flip, eng_mono, eng_mono_flip);
    fclose(f);
    
    if(filename != "")
    {
        const auto cvertices = mesh.getCurrentConfiguration().getVertices();
        const auto cface2vertices = mesh.getTopology().getFace2Vertices();
        WriteVTK writer(cvertices, cface2vertices);
        writer.addScalarFieldToFaces(energy_bi_bot.col(0), "bi_bot_aa");
        writer.addScalarFieldToFaces(energy_bi_bot.col(1), "bi_bot_bb");
        writer.addScalarFieldToFaces(energy_bi_bot.col(2), "bi_bot_ab");
        writer.addScalarFieldToFaces(energy_bi_top.col(0), "bi_top_aa");
        writer.addScalarFieldToFaces(energy_bi_top.col(1), "bi_top_bb");
        writer.addScalarFieldToFaces(energy_bi_top.col(2), "bi_top_ab");
        writer.addScalarFieldToFaces(energy_mono.col(0), "mono_aa");
        writer.addScalarFieldToFaces(energy_mono.col(1), "mono_bb");
        writer.addScalarFieldToFaces(thickness_bot, "h_bot");
        writer.addScalarFieldToFaces(thickness_top, "h_top");
        writer.addScalarFieldToFaces(thickness_mono, "h_mono");
        //writer.addScalarFieldToFaces(energy_mono.col(2), "mono_ab"); // this will be zero
        writer.addScalarFieldToFaces(energy_constant, "bi_const");
        writer.write(filename);
    }
}

Real Sim_Bilayer_4DFilaments::minimizePlateEnergy(const EnergyOperator<tMesh> & engOp, const bool simulatePlane, const bool planeFlip, const Real planePenalizationFac)
{
    Real eps = 1e-2;
    if(simulatePlane)
    {
        Parametrizer_FixedZPlane<tMesh> parametrizer(mesh, planePenalizationFac , 0.1, (planeFlip ? +1 : -1)*0.1, planeFlip);
        
        // create the energy minimizer
        HLBFGS_Methods::HLBFGS_EnergyOp_Parametrized<tMesh, Parametrizer_FixedZPlane, true> hlbfgs_wrapper(mesh, engOp, parametrizer);
        
        // minimize
        const Real epsMin = std::numeric_limits<Real>::epsilon();
        hlbfgs_wrapper.minimize(tag+"_diagnostics.dat", epsMin);
        eps = hlbfgs_wrapper.get_lastnorm();
    }
    else
    {
        // minimize energy using default method
        minimizeEnergy(engOp, eps);
    }
    
    return eps;
}

void Sim_Bilayer_4DFilaments::vtu_to_stl()
{
    // user info
    // whatever the thicknessvalue : the integral goes from -thick_bot/2 to 0, and 0 to thick_top/2 --> total thickness = 0.5(thick_bot + thick_top)
//    const std::string filename = "filaments_catenoid_thick000_botFac000_growth010_final";
//    const Real thick_bot = 1e3 * 0.500; // 2*thickness of bottom layer (since we define bottomlayer thickness as hbot/2)
//    const Real thick_top = 1e3 * 0.500; // 2*thickness of top layer
    
    const std::string filename_raw = parser.parse<std::string>("-fname");
    const std::string filename = helpers::removeExtension(filename_raw, ".vtp");
    const Real filament_thickness = parser.parse<Real>("-filthickness");
    
    // read vtu file
    Eigen::MatrixXd vertices, attributes;
    Eigen::MatrixXi faces;
    std::vector<std::pair<std::string, int>> attribute_entries;
    ReadVTK::read(filename, vertices, faces, attributes, attribute_entries);
    
    // initialize mesh
    IOGeometry geometry(filename+".vtp");
    mesh.init(geometry);
    
    // extract the filament densities for each layer
    int idx_fildens_bot = -1, idx_fildens_top = -1;
    int abs_idx = 0;
    for(size_t i=0;i<attribute_entries.size();++i)
    {
        const std::string attributeName = attribute_entries[i].first;
        if(attributeName == "dens_bot")
            idx_fildens_bot = abs_idx;
        else if(attributeName == "dens_top")
            idx_fildens_top = abs_idx;
        
        abs_idx += attribute_entries[i].second;
    }
    if( (idx_fildens_bot < 0) or (idx_fildens_top < 0) )
    {
        std::cout << "Did not find one of the column indices : " << idx_fildens_bot << "\t" << idx_fildens_top << std::endl;
        std::exit(-1);
    }
    
    // create thickness arrays
    const int nFaces = mesh.getNumberOfFaces();
    Eigen::VectorXd h_bot(nFaces), h_top(nFaces);
    for(int i=0;i<nFaces;++i)
    {
        h_bot(i) = attributes(i, idx_fildens_bot) * filament_thickness * 2; // times 2 because we multiply by 0.5 in the STL routine
        h_top(i) = attributes(i, idx_fildens_top) * filament_thickness * 2;
    }
    
    WriteSTLVolumetricCurved::write(mesh.getTopology(), mesh.getCurrentConfiguration().getVertices(), mesh.getCurrentConfiguration().getTriangleInfos(), h_bot, h_top, filename);
}


void Sim_Bilayer_4DFilaments::dumpAll(const FilamentLayerConfig<bottom> & botlayer, const FilamentLayerConfig<top> & toplayer, const std::string filename)
{
    
    WriteVTK writer(mesh.getCurrentConfiguration().getVertices(), mesh.getTopology().getFace2Vertices());
    writer.addVectorFieldToFaces(botlayer.filamentDirection, "dir_bot");
    writer.addVectorFieldToFaces(toplayer.filamentDirection, "dir_top");
    writer.addScalarFieldToFaces(botlayer.relFilamentDensity, "dens_bot");
    writer.addScalarFieldToFaces(toplayer.relFilamentDensity, "dens_top");
    
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
