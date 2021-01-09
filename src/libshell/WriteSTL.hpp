//
//  WriteSTL.hpp
//  Elasticity
//
//  Created by Wim van Rees on 4/12/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef WriteSTL_h
#define WriteSTL_h

#include <igl/writeSTL.h>
#include "TopologyData.hpp"
#include "DCSConfigurations.hpp"


class WriteSTL
{
    
public:
    
    static void write(const TopologyData & topo, const DCSBasicConfiguration & config, const std::string filename, const bool ascii=false)
    {
        const int nVertices = config.getNumberOfVertices();
        const int nFaces = topo.getNumberOfFaces();
        

        const auto vertices = config.getVertices();
        const auto face2vertices = topo.getFace2Vertices();
        
        Eigen::MatrixXd myVertices(nVertices,3);
        Eigen::MatrixXi myFace2Vertices(nFaces,3);
        Eigen::MatrixXd myFaceNormals = config.computeFaceNormals(topo);
        
        // copy over basic data
        for(int i=0;i<nVertices;++i)
        {
            myVertices(i,0) = vertices(i,0);
            myVertices(i,1) = vertices(i,1);
            myVertices(i,2) = vertices(i,2);
        }
        
        for(int i=0;i<nFaces;++i)
        {
            myFace2Vertices(i,0) = face2vertices(i,0);
            myFace2Vertices(i,1) = face2vertices(i,1);
            myFace2Vertices(i,2) = face2vertices(i,2);
        }
        
        //igl::writeSTL(filename+".stl", myVertices, myFace2Vertices, myFaceNormals, ascii); // no ascii (last argument)
        igl::FileEncoding encoding = (ascii ? igl::FileEncoding::Ascii : igl::FileEncoding::Binary);
        igl::writeSTL(filename+".stl", myVertices, myFace2Vertices, myFaceNormals, encoding);
    }
};
#endif /* WriteSTL_h */
