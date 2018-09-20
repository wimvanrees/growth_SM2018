//
//  HLBFGS_Wrapper_Parametrized.hpp
//  Elasticity
//
//  Created by Wim van Rees on 30/03/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef HLBFGS_Wrapper_Parametrized_hpp
#define HLBFGS_Wrapper_Parametrized_hpp

#ifdef USEHLBFGS

#include "HLBFGS_Wrapper.hpp"
#include "EnergyOperator.hpp"

namespace HLBFGS_Methods
{
    /**
     * @brief Wrapper around HBLFGS library that generalizes HLBFGS_Energy together with a parametrizer class
     * @tparam tMesh Type of the mesh
     * @tparam tMeshOperator Type of the energy operator
     * @tparam tParametrizer Type of the parametrization
     * @tparam verbose whether we should run the optimization silently or not
     *
     * @details
     * This class generalizes on HLBFGS_Energy in that it does not directly work on variables that are in the mesh, but rather on variables inside an instance of tParametrizer. These could be mesh variables inside a function, or derived variables from the mesh variables (eg volume). The gradients are handled through chain rules inside the tParametrizer class
     * @see Parametrizer     
     */
    
    template<typename tMesh, typename tMeshOperator, template<typename> class tParametrizer, bool verbose>
    class HLBFGS_Energy_Parametrized : public HLBFGS_Wrapper
    {
    protected:
        Profiler profiler;
        Profiler profiler_lbfgs;
        
        tMesh & mesh;
        const tMeshOperator & op;
        tParametrizer<tMesh> & parametrizer;
        const int nVariables_engop;
        const int nVariables_prms;
        Eigen::VectorXd gradient_container;
        
        std::string outFileName;
        Real last_gnorm;
        
    public:
        HLBFGS_Energy_Parametrized(tMesh & mesh, const tMeshOperator & engop, tParametrizer<tMesh> & parametrizer_in):
        mesh(mesh),
        op(engop),
        parametrizer(parametrizer_in),
        nVariables_engop(op.getNumberOfVariables(mesh)),
        nVariables_prms(parametrizer.getNumberOfVariables()),
        outFileName("hlbfgs_output.dat"),
        last_gnorm(-1)
        {
            gradient_container.resize(nVariables_engop);
        }
        
        int minimize(const std::string outFileName_in, const Real eps = 1e-5, const int Mval = 10)
        {
            this->outFileName = outFileName_in;
            return minimize(eps, Mval);
        }
        
        int minimize(const Real eps = 1e-5, const int Mval = 10)
        {
            mesh.updateDeformedConfiguration();
            const Real energy0 = op.compute(mesh) + parametrizer.computeEnergyContribution();
            
            // get variable pointer
            Real * x = parametrizer.getDataPointer();
            
            // init HLBFGS
            double parameter[20];
            int info[20];
            default_setup(parameter, info, eps, verbose);
            
            const int ret = HLBFGS(nVariables_prms, Mval, x, HLBFGS_Methods::evaluate, 0, HLBFGS_UPDATE_Hessian, HLBFGS_Methods::newiteration, this, parameter, info);
            
            if(verbose) std::cout << "HLBFGS return value = " << ret << std::endl;
            parametrizer.updateSolution();
            const Real energy1 = op.compute(mesh) + parametrizer.computeEnergyContribution();
            if(verbose) printf("Energy went from %10.10e to %10.10e, using epsilon = %e, final eps = %e \n", energy0, energy1, eps, last_gnorm);
            
            {
                FILE * ff = fopen(outFileName.c_str(), "a");
                fprintf(ff, "# ---- \n %d \t\t %d \t %d \t\t %10.10e \t %10.10e \t %10.10e\n# ---- \n",ret, info[2], info[1], energy0, energy1, last_gnorm);
                fclose(ff);
            }
            
            return (ret==2 or ret==3 or ret==4) ? 0 : 1;
        }
        
        
        // this function is called by the optimizer
        void evaluate(int /*N*/, double* /*x*/, double */*prev_x*/, double* f, double* g) override
        {
            //*f is the energy, *g is the gradient
            
            this->profiler.push_start("update mesh");
            // all the mesh updating is done in the parametrizer
            parametrizer.updateSolution();
            this->profiler.pop_stop();
                        
            // compute energy and gradient
            gradient_container.setZero();
            this->profiler.push_start("compute energy/gradient");
            const Real energy = op.compute(mesh, gradient_container) + parametrizer.computeEnergyContribution();
            this->profiler.pop_stop();

            // update the gradient
            parametrizer.updateGradient(nVariables_engop, gradient_container, g);
            
            *f = energy;
        }
        
        // this function is also called by the optimizer
        void newiteration(int iter, int call_iter, double* /* x */, double* f, double* /* g */,  double* gnorm) override
        {
            last_gnorm = *gnorm;
            
            if(iter%1000==0 and verbose)
            {
                printf("%d: \t %d \t %10.10e \t %10.10e\n",iter,call_iter, *f, *gnorm);
                op.printProfilerSummary();
            }
            
            if(iter%100==0)
            {
                FILE * ff = fopen(outFileName.c_str(), "a");
                fprintf(ff, "%d \t %d \t %10.10e \t %10.10e\n", iter, call_iter, *f, *gnorm);
                fclose(ff);
            }
        }
        
        /**
         * Exposes the last value of the residual
         *
         */
        Real get_lastnorm() const
        {
            return last_gnorm;
        }
        
    };
    
    template<typename tMesh, template<typename> class tParametrizer, bool verbose>
    class HLBFGS_EnergyOp_Parametrized : public HLBFGS_Energy_Parametrized<tMesh,EnergyOperator<tMesh>, tParametrizer, verbose>
    {
    public:
        HLBFGS_EnergyOp_Parametrized(tMesh & mesh, const EnergyOperator<tMesh> & engop, tParametrizer<tMesh> & parametrizer_in):
        HLBFGS_Energy_Parametrized<tMesh,EnergyOperator<tMesh>, tParametrizer, verbose>(mesh, engop, parametrizer_in)
        {}
        
    };
}
#endif

#endif /* HLBFGS_Wrapper_Parametrized_hpp */
