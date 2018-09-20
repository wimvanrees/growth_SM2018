//
//  WriteVTK_Quads.hpp
//  Elasticity
//
//  Created by Wim van Rees on 7/9/17.
//  Copyright © 2017 Wim van Rees. All rights reserved.
//

#ifndef WriteVTK_Quads_hpp
#define WriteVTK_Quads_hpp

#include "common.hpp"

#ifdef USEVTK
#include <vtkVersion.h>
#include <vtkPointData.h>
//#include <vtkCellData.h>
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
//#include <vtkDoubleArray.h>
#include <vtkPolyLine.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkCleanPolyData.h>
#include <vtkAppendPolyData.h>
#endif

#ifdef USEVTK
class WriteVTK_Quads
{
protected:
    vtkSmartPointer<vtkPolyData> polydata;
    
    void add_line_to_polydata(const Eigen::Ref<const Eigen::MatrixXd> pts)
    {
        const int nPoints = pts.rows();
        assert(nPoints > 0);
        assert(pts.cols() == 3);
        
        // create points
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        points->SetNumberOfPoints(nPoints);
        for(int i=0;i<nPoints;++i)
        {
            points->SetPoint(i, pts(i,0), pts(i,1), pts(i,2) );
        }
        
        // create polyline
        vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New();
        polyLine->GetPointIds()->SetNumberOfIds(nPoints);
        for(int i=0; i<nPoints; i++)
        {
            polyLine->GetPointIds()->SetId(i,i);
        }
        // create cell array to store lines in
        vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
        cells->InsertNextCell(polyLine);
        
        // Create our polydata object and add the points and linesto it.
        vtkSmartPointer<vtkPolyData> polydata_tmp = vtkSmartPointer<vtkPolyData>::New();
        polydata_tmp->SetPoints(points);
        polydata_tmp->SetLines(cells);
        
        // append it to our current polydata
        vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
        appendFilter->AddInputData(polydata);
        appendFilter->AddInputData(polydata_tmp);
        appendFilter->Update();
        
        // Remove any duplicate points.
        vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
        cleanFilter->SetInputConnection(appendFilter->GetOutputPort());
        cleanFilter->Update();
        
        polydata = cleanFilter->GetOutput();
    }
    
public:
    WriteVTK_Quads()
    {}
    
    WriteVTK_Quads(const Eigen::Ref<const Eigen::MatrixXd> line_pts)
    {
        add_line(line_pts);
    }
    
    void add_line(const Eigen::Ref<const Eigen::MatrixXd> line_pts)
    {
        add_line_to_polydata(line_pts);
    }

    vtkSmartPointer<vtkDoubleArray> convertScalarFieldToVTK(const Eigen::Ref<const Eigen::VectorXd> field, const std::string fieldname)
    {
        const int N = field.rows();
        vtkSmartPointer<vtkDoubleArray> scalar_array = vtkSmartPointer<vtkDoubleArray>::New();
        scalar_array->SetNumberOfComponents(1);
        scalar_array->SetNumberOfTuples(N);
        scalar_array->SetName(fieldname.c_str());
        for(int i=0;i<N;++i)
            scalar_array->SetTuple1(i, field(i));
        
        return scalar_array;
    }
    
    void addScalarFieldToLines(const Eigen::Ref<const Eigen::VectorXd> field, const std::string fieldname)
    {
        assert(field.rows() == polydata->GetNumberOfCells());
        
        vtkSmartPointer<vtkDoubleArray> scalar_array = convertScalarFieldToVTK(field, fieldname);
        polydata->GetCellData()->AddArray(scalar_array);
    }
    
    int write(const std::string filename)
    {
        vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
#if VTK_MAJOR_VERSION <= 5
        writer->SetInput(polydata);
#else
        writer->SetInputData(polydata);
#endif
        writer->SetFileName( (filename+".vtp").c_str());
        const int retval = writer->Write(); // return value 1 for success, 0 for failure
        
        return retval;
    }
};
#else
class WriteVTK_Quads
{
protected:
    
public:
    
    WriteVTK_Quads()
    {
        std::cout << "NO VTK -- WriteVTK_Quads not supported " << std::endl;
    }
    
    WriteVTK_Quads(const Eigen::Ref<const Eigen::MatrixXd> )
    {
        std::cout << "NO VTK -- WriteVTK_Quads not supported " << std::endl;
    }
    
    void add_line(const Eigen::Ref<const Eigen::MatrixXd> )
    {
        std::cout << "NO VTK -- WriteVTK_Quads not supported " << std::endl;
    }
    
    int write(const std::string)
    {
        return -1;
    }
};
#endif


#endif /* WriteVTK_Quads_h */
