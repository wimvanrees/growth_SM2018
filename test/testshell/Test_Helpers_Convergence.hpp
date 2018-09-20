//
//  Test_Helpers_Convergence.hpp
//  Elasticity
//
//  Created by Wim van Rees on 8/20/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef Test_Helpers_Convergence_hpp
#define Test_Helpers_Convergence_hpp

#include "common.hpp"

namespace Test_Helpers
{
    class Convergence
    {
    public:
        static std::vector<Real> computeConvergence(const std::vector<std::pair<int, Real>> & results, const Real exactVal, const bool verbose=false)
        {
            const size_t nResults = results.size();
            std::vector<Real> slopes(nResults-1);
            if(verbose)
            {
                const Real error = std::abs(results[0].second-exactVal);
                printf("[----------] nFaces \t computed Eng \t error \t errorRel \t order wrt prev \t (exact Energy = %10.10e )\n", exactVal);
                printf("[----------] -------------------------------------------\n");
                printf("[----------] %d \t %10.10e \t %10.10e \t %10.10e \t 0\n", results[0].first, results[0].second, error, error/exactVal);
            }
            for(size_t i=1;i<nResults;++i)
            {
                const int nFacesCoarse = results[i-1].first;
                const Real errorCoarse = std::abs(results[i-1].second - exactVal);
                
                const int nFacesFine = results[i].first;
                const Real errorFine = std::abs(results[i].second - exactVal);
                
                const Real slope = (errorCoarse/errorFine) * (1.0*nFacesCoarse/nFacesFine);
                slopes[i-1] = slope;
                if(verbose)
                {
                    printf("[----------] %d \t %10.10e \t %10.10e \t %10.10e \t %10.10e \n", results[i].first, results[i].second, errorFine, errorFine/exactVal, slope);
                }
            }
            
            return slopes;
        }
        
        static bool testConvergence(const std::vector<Real> & slopes, const Real targetSlope, const Real acceptableError)
        {
            std::vector<bool> isSlopeAccepted(slopes.size());
            for(size_t i=0;i<slopes.size();++i)
                isSlopeAccepted[i] = (std::abs(slopes[i] - targetSlope) <= acceptableError);
            
            // test if true
            const bool allAccepted = (std::all_of(std::begin(isSlopeAccepted),std::end(isSlopeAccepted),[](bool i){return i;}));
            return allAccepted;
        }
    };
}
#endif /* Test_Helpers_Convergence_hpp */
