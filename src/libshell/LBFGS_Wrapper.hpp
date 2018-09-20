//
//  LBFGS_Wrapper.hpp
//  Elasticity
//
//  Created by Wim van Rees on 3/2/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef LBFGS_Wrapper_h
#define LBFGS_Wrapper_h

#ifdef USELIBLBFGS

#include "lbfgs.h"
#include "Profiler.hpp"
#include "Mesh.hpp"

namespace LBFGS
{
    template<typename T>
    class LBFGS_Wrapper
    {
    protected:
        
    public:
        
        virtual T evaluate(const int n, const T * const x, T * gradient) = 0;
        
        virtual int progress(const T * const /*x*/, const T * const /*gradient*/, const T value, const T xnorm, const T gnorm, const T step, const int /*n*/, const int k, const int /*ls*/)
        {
            if(k%1000==0)
            {
                printf("Iteration %d:\n", k);
                printf("  function value = %10.10e\n", value);
                printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
                printf("\n");
            }
            return 0;
        }
    };    
    
    static lbfgsfloatval_t evaluate(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t /*step*/)
    {
        LBFGS_Wrapper<lbfgsfloatval_t> * myLBFGS = reinterpret_cast<LBFGS_Wrapper<lbfgsfloatval_t>*>(instance);
        return myLBFGS->evaluate(n, x, g);
    }
    
    static int progress(void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls)
    {
        LBFGS_Wrapper<lbfgsfloatval_t> * myLBFGS = reinterpret_cast<LBFGS_Wrapper<lbfgsfloatval_t>*>(instance);
        return myLBFGS->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
    }


    template<typename T, typename tMesh, typename tMeshOperator, bool verbose>
    class LBFGS_Energy : public LBFGS_Wrapper<T>
    {
    protected:
        Profiler profiler;
        Profiler profiler_lbfgs;
        
        tMesh & mesh;
        const tMeshOperator & op;
        
        void reportResult(const int result)
        {
            std::cout << "L-BFGS optimization terminated : ";
            if(result==LBFGS_SUCCESS)
                std::cout << "successfully";
            else if(result==LBFGS_ALREADY_MINIMIZED)
                std::cout << "function was already minimized";
            else if(result==LBFGSERR_UNKNOWNERROR)
                std::cout << "unknown error";
            else if(result==LBFGSERR_LOGICERROR)
                std::cout << "logic error";
            else if(result==LBFGSERR_OUTOFMEMORY)
                std::cout << "insufficient memory";
            else if(result==LBFGSERR_CANCELED)
                std::cout << "canceled";
            else if(result==LBFGSERR_OUTOFINTERVAL)
                std::cout << "line-search step went out of interval of uncertainty";
            else if(result==LBFGSERR_INCORRECT_TMINMAX)
                std::cout << "logic error, or interval of uncertainty becamse too small";
            else if(result==LBFGSERR_ROUNDING_ERROR)
                std::cout << "rounding error occurred, or no line-search step satisfies the sufficient decrease and curvature conditions";
            else if(result==LBFGSERR_MINIMUMSTEP)
                std::cout << "line-search step became smaller than min_step";
            else if(result==LBFGSERR_MAXIMUMSTEP)
                std::cout << "line-search step became larger than max_step";
            else if(LBFGSERR_MAXIMUMLINESEARCH)
                std::cout << "line-search reached the maximum number of iterations";
            else if(LBFGSERR_MAXIMUMITERATION)
                std::cout << "algorithm reaches maximum number of iterations";
            else if(LBFGSERR_WIDTHTOOSMALL)
                std::cout << "relative width of interval of uncertainty is at most xtol";
            else if(LBFGSERR_INVALIDPARAMETERS)
                std::cout << "logic error - negative line-search step";
            else if(LBFGSERR_INCREASEGRADIENT)
                std::cout << "current search direction increases objective function value";

            std::cout << std::endl;
        }
    public:
        
        LBFGS_Energy(tMesh & mesh, const tMeshOperator & engop):
            LBFGS_Wrapper<T>(),
            mesh(mesh),
            op(engop)
        {}
        
        int minimize(const Real eps = 1e-5)
        {
            const int nVertices = mesh.getNumberOfVertices();
            const int nEdges = mesh.getNumberOfEdges();
            const int nVariables = 3*nVertices + nEdges;
            

            
            mesh.updateDeformedConfiguration();
            const Real energy0 = op.compute(mesh);
            
            // declare variables
            lbfgsfloatval_t fx; // function value
            lbfgsfloatval_t *x = mesh.getDataPointer(); // unknowns: we use the raw pointer from the mesh
            
            /* Initialize the parameters for the L-BFGS optimization. */
            lbfgs_parameter_t param;
            lbfgs_parameter_init(&param); // put default values

            param.m = 100;//12; // inverse hessian # of corrections (default = 6)
            param.epsilon = eps; // epsilon for convergence test (default = 1e-5)
            param.past = 0; // distance for delta-based convergence tests (default = 0)
            param.delta = 1e-5; // delta for convergence tests (default = 1e-5)
            param.max_iterations = 0; // maximum number of iterations (default = 0)
            param.linesearch = LBFGS_LINESEARCH_DEFAULT;//LBFGS_LINESEARCH_DEFAULT; // line search algorithm to be used (default = more-thuente)
            param.max_linesearch = 40; // maximum number of trials for the linesearch (default = 40)
            param.min_step = 1e-20;
            param.max_step = 1e20;
            param.ftol = 1e-4; // control the accuracy of the linesearch routine (>0, <0.5, default = 1e-4)
            param.wolfe = 0.9; // coefficient for wolfe routine (only with backtracking_wolfe or backtracking_strong_wolfe)
            param.gtol = 0.9; // control the accuracy of the linesearch routine (larger than ftol, smaller than 1, default = 0.9)
            param.xtol = std::numeric_limits<T>::epsilon(); // machine precision for floating point values (default = 1e-16)
            // rest is for Least Squares
//            param.orthantwise_c = ;
//            param.orthantwise_start = ;
//            param.orthantwise_end = ;
            
            /*param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;*/
            
            /*
             Start the L-BFGS optimization; this will invoke the callback functions
             evaluate() and progress() when necessary.
             */
            
            this->profiler_lbfgs.push_start("lbfgs");
            const int ret = lbfgs(nVariables, x, &fx, LBFGS::evaluate, LBFGS::progress, this, &param);
            this->profiler_lbfgs.pop_stop();
            
            /* Report the result. */
            if(verbose) reportResult(ret);
            if(verbose) printf("Energy went from %f to %f, using epsilon = %e \n", energy0, fx, eps);
//            {
//                FILE * f = fopen("LBFGS_res.dat","aw");
//                fprintf(f,"%10.10e \t %10.10e \t %10.10e\n",energy0,fx,eps);
//                fclose(f);
//            }

//            lbfgs_free(x);
            return (ret==LBFGS_SUCCESS or ret==LBFGS_ALREADY_MINIMIZED) ? 0 : 1;
        }
        
        virtual int progress(const T * const /*x*/, const T * const /*gradient*/, const T value, const T xnorm, const T gnorm, const T step, const int /*n*/, const int k, const int /*ls*/) override
        {
            if(k%1000==0 and verbose)
            {
                printf("Iteration %d:\n", k);
                printf("  function value = %10.10e\n", value);
                printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
                profiler_lbfgs.printSummary();
                profiler.printSummary();
                op.printProfilerSummary();
                printf("\n");
            }
            return 0;
        }
        
#ifndef NDEBUG
        virtual T evaluate(const int n, const T * const x, T * gradient) override
#else
        virtual T evaluate(const int /*n*/, const T * const /*x*/, T * gradient) override
#endif
        {
            this->profiler_lbfgs.pop_stop();
            
            const int nVertices = mesh.getNumberOfVertices();
            const int nEdges = mesh.getNumberOfEdges();
            
            assert(n==(3*nVertices + nEdges));

            this->profiler.push_start("update edges");
            mesh.updateDeformedConfiguration();
            this->profiler.pop_stop();

            // wrap Eigen structures around the gradient pointer
            Eigen::Map<Eigen::MatrixXd> eigenGrad_vertices(gradient, nVertices, 3);
            Eigen::Map<Eigen::VectorXd> eigenGrad_edges(gradient + 3*nVertices, nEdges);
            
            eigenGrad_vertices.setZero();
            eigenGrad_edges.setZero();
            
            // compute energy
            this->profiler.push_start("compute energy");
            const Real energy = op.compute(mesh, eigenGrad_vertices, eigenGrad_edges);
            this->profiler.pop_stop();
            
            this->profiler_lbfgs.push_start("lbfgs");
            return energy;
        }
        
    };
}
#endif

#endif /* LBFGS_Wrapper_h */
