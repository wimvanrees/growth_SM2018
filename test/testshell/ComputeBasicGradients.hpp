//
//  ComputeBasicGradients.hpp
//  Elasticity
//
//  Created by Wim van Rees on 4/28/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef ComputeBasicGradients_hpp
#define ComputeBasicGradients_hpp

#include "common.hpp"
#include "Mesh.hpp"
#include "ExtendedTriangleInfo.hpp"
#include "EnergyOperator.hpp"

namespace TestBasicGradients
{
    template<typename tMesh, int tIdx>
    class Test_GradHeight : public EnergyOperator_DCS<tMesh>
    {
        typedef typename tMesh::tCurrentConfigData tCurrentConfigData;
        typedef typename tMesh::tReferenceConfigData tReferenceConfigData;
        
    protected:
        Real computeAll(const tMesh & mesh, Eigen::Ref<Eigen::MatrixXd> gradient_vertices, Eigen::Ref<Eigen::VectorXd> , const bool computeGradient) const override
        {
            if(tIdx < 0 or tIdx > 2) return 0;
            
            const auto & currentState = mesh.getCurrentConfiguration();
            
            // do the calculations
            const int nFaces = mesh.getNumberOfFaces();
            Real energy = 0.0;
            for(int i=0;i<nFaces;++i)
            {
                const ExtendedTriangleInfo & info = currentState.getTriangleInfo(i);
                
                const Real h = info.height(tIdx);
                
                energy += h*h;
                
                if(not computeGradient) continue;
                
                const Eigen::Vector3d gradv0_h = info.get_gradv_height<0,tIdx>();
                const Eigen::Vector3d gradv1_h = info.get_gradv_height<1,tIdx>();
                const Eigen::Vector3d gradv2_h = info.get_gradv_height<2,tIdx>();
                
                for(int d=0;d<3;++d)
                {
                    gradient_vertices(info.idx_v0,d) += 2.0*h*gradv0_h(d);
                    gradient_vertices(info.idx_v1,d) += 2.0*h*gradv1_h(d);
                    gradient_vertices(info.idx_v2,d) += 2.0*h*gradv2_h(d);
                }
            }
            return energy;
        }
    };
    
    template<typename tMesh, int tIdx>
    class Test_GradTheta : public EnergyOperator_DCS<tMesh>
    {
        typedef typename tMesh::tCurrentConfigData tCurrentConfigData;
        typedef typename tMesh::tReferenceConfigData tReferenceConfigData;
        
    protected:
        Real computeAll(const tMesh & mesh, Eigen::Ref<Eigen::MatrixXd> gradient_vertices, Eigen::Ref<Eigen::VectorXd> , const bool computeGradient) const override
        {
            if(tIdx < 0 or tIdx > 2) return 0;
            
            const auto & currentState = mesh.getCurrentConfiguration();
            
            // do the calculations
            const int nFaces = mesh.getNumberOfFaces();
            Real energy = 0.0;
            
            for(int i=0;i<nFaces;++i)
            {
                const ExtendedTriangleInfo & info = currentState.getTriangleInfo(i);
                
                const Real theta = info.theta(tIdx);
                
                energy += theta*theta;
                
                if(not computeGradient) continue;
                
                // compute the gradient depending on tIdx = eidx_rel
                
                Eigen::Vector3d gradv0, gradv1, gradv2, gradv_other;
                int vidx_other;
                
                const bool isInterior = info.get_gradv_theta<tIdx>(gradv0, gradv1, gradv2, gradv_other, vidx_other);
                if(not isInterior) continue;
                
                for(int d=0;d<3;++d)
                {
                    gradient_vertices(info.idx_v0,d) += 2.0*theta*gradv0(d);
                    gradient_vertices(info.idx_v1,d) += 2.0*theta*gradv1(d);
                    gradient_vertices(info.idx_v2,d) += 2.0*theta*gradv2(d);
                    gradient_vertices(vidx_other ,d) += 2.0*theta*gradv_other(d);
                }
            }
            
            return energy;
        }
    };
    
}

#endif /* ComputeBasicGradients_hpp */
