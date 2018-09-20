//
//  Sim_Bilayer_4DFilaments.hpp
//  Elasticity
//
//  Created by Wim van Rees on 8/10/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef Sim_Bilayer_4DFilaments_hpp
#define Sim_Bilayer_4DFilaments_hpp


#include "Sim.hpp"
#include "Mesh.hpp"
#include "MaterialProperties.hpp"
#include "GrowthHelper.hpp"
#include "CombinedOperator_Parametric.hpp"

template<MeshLayer layer>
struct FilamentLayerConfig
{
    static_assert((layer==top) or (layer==bottom), "FilamentLayerConfig<layer> only works on bilayer");
    
    typedef std::function<Eigen::Vector3d(const Real, const Real, const Real)> func3d;
    
    const int nFaces;
    const Real Young;
    const Real nu;
    const Real thickness;
    const Real swellFac;
    const bool adjustYoung;
    
    Eigen::VectorXd relFilamentDensity;
    Eigen::MatrixXd filamentDirection;
    
    FilamentLayerConfig(const int nFaces, const Real Y, const Real n, const Real h, const Real swellFac_in = 1.0, const bool adjustYoung = false):
    nFaces(nFaces),
    Young(Y),
    nu(n),
    thickness(h),
    swellFac(swellFac_in),
    adjustYoung(adjustYoung),
    relFilamentDensity(Eigen::VectorXd::Constant(nFaces, 1)),
    filamentDirection(Eigen::MatrixXd::Constant(nFaces, 3, 0))
    {}
    
    void setFilamentDensity(const Eigen::Ref<const Eigen::VectorXd> fildens)
    {
        assert(fildens.rows() == nFaces);
        relFilamentDensity = fildens;
    }
    
    void setFilamentDirectionFromMatrix(const Eigen::Ref<const Eigen::MatrixXd> fildir)
    {
        assert(fildir.rows() == nFaces);
        assert(fildir.cols() == 3);
        filamentDirection = fildir;
    }
    
    void setFilamentDirectionFromFunction(const BilayerMesh & mesh, const func3d parallel_dir_func)
    {
        filamentDirection.resize(nFaces, 3);
        for(int i=0;i<nFaces;++i)
        {
            const TriangleInfo info = mesh.getRestConfiguration().getTriangleInfoLite(mesh.getTopology(), i);
            const Eigen::Vector3d restfacecenter = info.computeFaceCenter();
            const Real myX = restfacecenter(0);
            const Real myY = restfacecenter(1);
            const Real myZ = restfacecenter(2);
            
            filamentDirection.row(i) = parallel_dir_func(myX, myY, myZ);
        }
    }
    
    MaterialProperties_Iso_Array getMaterialProperties(const Real swellThickness = 0.0) const
    {
        // create thickness and Young arrays
        Eigen::VectorXd thickness_array(nFaces), Young_array(nFaces);
        for(int i=0;i<nFaces;++i)
        {
            thickness_array(i) = (1.0 + swellThickness) * relFilamentDensity(i) * thickness;
            Young_array(i) = adjustYoung ? relFilamentDensity(i) * Young : Young;
        }
        MaterialProperties_Iso_Array retval(Young_array, nu, thickness_array);
        return retval;
    }
    
    MaterialProperties_Ortho_Array getMaterialPropertiesOrtho(const BilayerMesh & mesh, const Real orthoFac, const Real shearFac, const Real swellThickness = 0.0) const
    {
        // first create constant array
        const Real Y1 = Young;
        const Real Y2 = orthoFac * Young;
        const Real G12 = shearFac * Young;
        const Real nu1 = nu;
        const Real h = thickness;
        
        const int nFaces = mesh.getNumberOfFaces();
        
        MaterialProperties_Ortho_Array retval(nFaces, Y1, Y2, G12, nu1, h);
        
        // now we adjust
        for(int i=0;i<nFaces;++i)
        {
            // adjust thickness
            retval.materials[i].thickness = (1.0 + swellThickness) * relFilamentDensity(i) * thickness;
            
            // adjust Young if needed
            if(adjustYoung)
            {
                retval.materials[i].Young_d1 *= relFilamentDensity(i);
                retval.materials[i].Young_d2 *= relFilamentDensity(i);
                retval.materials[i].G_12 *= relFilamentDensity(i);
            }
            
            // now the direction
            const TriangleInfo info = mesh.getRestConfiguration().getTriangleInfoLite(mesh.getTopology(), i);
            
            // compute face normal
            const Eigen::Vector3d restface_normal = ((info.e2).cross(info.e0)).normalized();
            
            // compute orthogonal directions for this face
            const Eigen::Vector3d orthodir_d1 = filamentDirection.row(i); // parallel
            const Eigen::Vector3d a1 = (orthodir_d1 - (orthodir_d1.dot(restface_normal))*restface_normal).normalized(); // a1 mapped onto the triangle
            //                const Eigen::Vector3d a2 = (restface_normal.cross(a1)).normalized(); // a2 --> no need to do this
            
            // create the aform_bar
            // e1.dot(e1) = (e1.dot(a1))^2 + (e1.dot(a2))^2
            const Real a11_d1 = std::pow((info.e1).dot(a1), 2);
            const Real a12_d1 = (info.e1).dot(a1) * (info.e2).dot(a1);
            const Real a22_d1 = std::pow((info.e2).dot(a1), 2);
            
            retval.materials[i].aform_bar_d1 << a11_d1, a12_d1, a12_d1, a22_d1;
        }
        
        
        return retval;
    }
    
    void applyGrowthToMesh(BilayerMesh & mesh, const Real growthFac_p, const Real growthFac_o) const
    {
        // now we create the abar
        const int nFaces = mesh.getNumberOfFaces();
        const Eigen::VectorXd growthRates_parallel(Eigen::VectorXd::Constant(nFaces, swellFac*growthFac_p));
        const Eigen::VectorXd growthRates_orthogonal(Eigen::VectorXd::Constant(nFaces, swellFac*growthFac_o));
        Eigen::VectorXd angles(nFaces);
        
        for(int i=0;i<nFaces;++i)
        {
            const Eigen::Vector3d parallelDir = filamentDirection.row(i);
            angles(i) = std::atan2(parallelDir(1), parallelDir(0)); // dirX = cos(phi), dirY = sin(phi)
        }
        tVecMat2d & aforms = mesh.getRestConfiguration().getFirstFundamentalForms<layer>();
        
        GrowthHelper<BilayerMesh>::computeAbarsOrthoGrowth(mesh, angles, growthRates_parallel, growthRates_orthogonal, aforms);
    }
};

struct DensityInputs
{
    const std::string filename_bot, filename_top;
    const Real sigma;
    const Real filament_thickness;
    const Real dpath;

    Real target_size_X, target_size_Y;
    Real origin_X, origin_Y;
    Real min_extent_X, max_extent_X, min_extent_Y, max_extent_Y;
    Real dropSize, dropFac;
    
    DensityInputs(const std::string fname_bot, const std::string fname_top, const Real sigma, const Real d_f, const Real dpath):
    filename_bot(fname_bot),
    filename_top(fname_top),
    sigma(sigma),
    filament_thickness(d_f),
    dpath(dpath)
    {
        // set defaults
        target_size_X = -1;
        target_size_Y = -1;
        origin_X = 0.0;
        origin_Y = 0.0;
        min_extent_X = min_extent_Y = std::numeric_limits<Real>::lowest();
        max_extent_X = max_extent_Y = std::numeric_limits<Real>::max();
        dropSize = 1.0;
        dropFac = 1.0;
    }
};


class Sim_Bilayer_4DFilaments : public Sim<BilayerMesh>
{
    typedef BilayerMesh tMesh;
    
protected:
    
    void get_density_from_svg(const DensityInputs & densityInputs, Eigen::VectorXd & density_bot, Eigen::VectorXd & density_top, Eigen::MatrixXd & dir_bot, Eigen::MatrixXd & dir_top);

    void run_catenoid_printpath();
    void run_helicoid_printpath();
    void run_sombrero_printpath();
    void run_logspiral_printpath();
    void run_folding_flower_printpath();
    void run_orchid_printpath();
    void setup_general_printpath(const Real filament_thickness, const Eigen::Ref<const Eigen::VectorXd> density_bot, const Eigen::Ref<const Eigen::VectorXd> density_top, const Eigen::Ref<const Eigen::MatrixXd> dir_bot, const Eigen::Ref<const Eigen::MatrixXd> dir_top, const bool simulatePlane = false);
    
    void vtu_to_stl();
    
    void runSwelling(const FilamentLayerConfig<bottom> & botlayer, const FilamentLayerConfig<top> & toplayer, const std::vector<Real> & swellrates, const bool swellThickness = false, const bool simulatePlane = false);
    Real minimizePlateEnergy(const EnergyOperator<tMesh> & engOp, const bool simulatePlane, const bool planeFlip, const Real planePenalizationFac);
    void dumpAll(const FilamentLayerConfig<bottom> & botlayer, const FilamentLayerConfig<top> & toplayer, const std::string filename);
    
    void dumpBilayerEnergyDecomposition(const int idx, const Real swelling_fac, const MaterialProperties_Iso_Array & matprop_bot, const MaterialProperties_Iso_Array & matprop_top, const CombinedOperator_Parametric<tMesh, Material_Isotropic, bottom> & engOp_bot, const CombinedOperator_Parametric<tMesh, Material_Isotropic, top> & engOp_top, const bool overWrite, const std::string filename);
    
public:
    
    Sim_Bilayer_4DFilaments(ArgumentParser & parser):
    Sim<tMesh>(parser)
    {}
    
    void init() override;
    void run() override;
};

#endif /* Sim_Bilayer_4DFilaments_hpp */
