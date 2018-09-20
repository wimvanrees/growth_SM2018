//
//  Sim_IOops.cpp
//  Elasticity
//
//  Created by Wim van Rees on 12/19/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#include "Sim_IOops.hpp"
#include "Geometry.hpp"
#include "ReadVTK.hpp"
#include "WriteVTK.hpp"
#include "WriteSTL.hpp"
#include "WriteINP.hpp"
#include <igl/writeOBJ.h>
#include <ctime>

void Sim_IOops::convertCSV(const std::string & fname_in, const std::string & fname_out_raw)
{
    // read the VTP data
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi face2vertices;
    Eigen::MatrixXd attributes;
    std::vector<std::pair<std::string, int>> attribute_entries;
    
    ReadVTK::read(fname_in, vertices, face2vertices, attributes, attribute_entries);
    
    // write data to csv, incl attributes
    const std::string fname_out = helpers::removeExtension(fname_out_raw, ".csv")+".csv";
    
    const int nEntries = (int)attribute_entries.size();
    const int nCols = attributes.cols();
    
    std::cout << "FOUND " << nEntries << " entries" << std::endl;
    
    FILE * f = fopen(fname_out.c_str(), "w");
    // write header
    fprintf(f,"#");
    for(int n=0;n<nEntries;++n)
        fprintf(f,"\t %s", attribute_entries[n].first.c_str());
    
    const int nFaces = face2vertices.rows();
    for(int i=0;i<nFaces;++i)
    {
        for(int n=0;n<nCols;++n)
        {
            fprintf(f,"%10.10e", attributes(i,n));
            if(n < nCols - 1) fprintf(f,"\t");
            else fprintf(f,"\n");
        }
    }
    fclose(f);   
}

void Sim_IOops::convertSTL(const std::string & fname_out)
{
    const bool ascii = parser.parse<bool>("-ascii", false);
    const std::string fname_out_noext = helpers::removeExtension(fname_out, ".stl");
    WriteSTL::write(mesh.getTopology(), mesh.getCurrentConfiguration(), fname_out_noext, ascii);
}

void Sim_IOops::convertOBJ(const std::string & fname_out)
{
    const Eigen::MatrixXd vertices = mesh.getCurrentConfiguration().getVertices();
    const Eigen::MatrixXi face2vertices = mesh.getTopology().getFace2Vertices();
    igl::writeOBJ<Eigen::MatrixXd, Eigen::MatrixXi>(fname_out, vertices, face2vertices);
}

void Sim_IOops::convertINP(const std::string & fname_out)
{
    // include file for abaqus .inp format
    const std::string elType = parser.parse<std::string>("-element", "STRI3");
    const std::string fname_out_noext = helpers::removeExtension(fname_out, ".inp");
    WriteINP::write(mesh.getTopology(), mesh.getCurrentConfiguration(), fname_out_noext, elType);
}

void Sim_IOops::init()
{
    
}

void Sim_IOops::run()
{
    // read an input file and initialize the mesh
    const std::string fname_in = parser.parse<std::string>("-filename","");
    IOGeometry geometry(fname_in);
    mesh.init(geometry);

    // now we write the output data based on the extension
    const std::string fname_out = parser.parse<std::string>("-filename_out","");    
    const std::string::size_type idx = fname_out.rfind('.');
    std::string extension = "";
    if(idx != std::string::npos)
    {
        extension = fname_out.substr(idx+1);
        std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);        
    }
    
    if(extension=="stl")
    {
        convertSTL(fname_out);
    }
    else if(extension=="obj")
    {
        convertOBJ(fname_out);
    }
    else if(extension=="inp")
    {
        convertINP(fname_out);
    }    
    else
    {
        convertCSV(fname_in, fname_out);
    }
}
