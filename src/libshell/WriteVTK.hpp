//
//  WriteVTK.hpp
//  Elasticity
//
//  Created by Wim van Rees on 2/17/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef WriteVTK_hpp
#define WriteVTK_hpp

#include "common.hpp"

#ifdef USEVTK
#include <vtkVersion.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkTriangle.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>
#endif

/*! \class WriteVTK
 * \brief Write a triangle mesh in VTK output to be read by Paraview or Visit
 *
 * If USEVTK is defined, this method will use VTK libraries to write output in binary format using VTP (polygonal) format. If USEVTK is not defined, it will write plain ASCII using inefficient legacy VTK format.
 */

#ifdef USEVTK
class WriteVTK
{
protected:
    const int nVertices;
    const int nFaces;
    vtkSmartPointer<vtkPolyData> polydata;
    
public:
    WriteVTK(const Eigen::Ref<const Eigen::MatrixXd> verts, const Eigen::Ref<const Eigen::MatrixXi> faces):
    nVertices(verts.rows()),
    nFaces(faces.rows())
    {
        assert(nVertices > 0);
        assert(nFaces > 0);
        assert(verts.cols() == 3);
        assert(faces.cols() == 3);
        
        // create points
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        points->SetNumberOfPoints(nVertices);
        for(int i=0;i<nVertices;++i)
        {
            points->SetPoint(i, verts(i,0), verts(i,1), verts(i,2) );
        }
        
        // Create our polydata object and add the points to it.
        polydata = vtkSmartPointer<vtkPolyData>::New();
        polydata->SetPoints(points);
        
        // create faces
        vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
        for(int i=0;i<nFaces;++i)
        {
            // Create a triangle
            vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
            
            triangle->GetPointIds()->SetId(0,faces(i,0));
            triangle->GetPointIds()->SetId(1,faces(i,1));
            triangle->GetPointIds()->SetId(2,faces(i,2));
            
            triangles->InsertNextCell(triangle);
        }
        
        polydata->SetPolys(triangles);
    }
    
    WriteVTK(const Eigen::Ref<const Eigen::MatrixXd> verts):
    nVertices(verts.rows()),
    nFaces(-1)
    {
        assert(nVertices > 0);
        assert(verts.cols() == 3);
        
        // create points
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        points->SetNumberOfPoints(nVertices);
        for(int i=0;i<nVertices;++i)
        {
            points->SetPoint(i, verts(i,0), verts(i,1), verts(i,2) );
        }
        
        // Create our polydata object and add the points to it.
        polydata = vtkSmartPointer<vtkPolyData>::New();
        polydata->SetPoints(points);
        
        // no faces
    }
    
    vtkSmartPointer<vtkDoubleArray> convertVectorFieldToVTK(const Eigen::Ref<const Eigen::MatrixXd> field, const std::string fieldname)
    {
        const int N = field.rows();
        assert(field.cols() == 3);
        vtkSmartPointer<vtkDoubleArray> vector_array = vtkSmartPointer<vtkDoubleArray>::New();
        vector_array->SetNumberOfComponents(3);
        vector_array->SetNumberOfTuples(N);
        vector_array->SetName(fieldname.c_str());
        for(int i=0;i<N;++i)
            vector_array->SetTuple3(i, field(i,0),field(i,1),field(i,2));
        
        return vector_array;
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
    
    void addTensorFieldToFaces(const tVecMat3d & field, const std::string fieldname)
    {
        const int N = (int)field.size();
        assert(N == nFaces);
        // assuming the tensor is symmetric
        
        vtkSmartPointer<vtkDoubleArray> tensor_array = vtkSmartPointer<vtkDoubleArray>::New();
        tensor_array->SetNumberOfComponents(9);
        tensor_array->SetNumberOfTuples(N);
        tensor_array->SetName(fieldname.c_str());
        for(int i=0;i<N;++i)
        {
            const Eigen::Matrix3d & tensor = field[i];
            tensor_array->SetTuple9(i, tensor(0,0), tensor(0,1), tensor(0,2), tensor(1,0), tensor(1,1), tensor(1,2), tensor(2,0), tensor(2,1), tensor(2,2));
        }
        polydata->GetCellData()->AddArray(tensor_array);
    }
    
    
    void addVectorFieldToFaces(const Eigen::Ref<const Eigen::MatrixXd> field, const std::string fieldname)
    {
        assert(field.rows() == nFaces);
        assert(field.cols() == 3);
        
        vtkSmartPointer<vtkDoubleArray> vector_array = convertVectorFieldToVTK(field, fieldname);
        polydata->GetCellData()->AddArray(vector_array);
    }
    
    void addScalarFieldToFaces(const Eigen::Ref<const Eigen::VectorXd> field, const std::string fieldname)
    {
        assert(field.rows() == nFaces);
        
        vtkSmartPointer<vtkDoubleArray> scalar_array = convertScalarFieldToVTK(field, fieldname);
        polydata->GetCellData()->AddArray(scalar_array);
    }
    
    
    void addVectorFieldToVertices(const Eigen::Ref<const Eigen::MatrixXd> field, const std::string fieldname)
    {
        assert(field.rows() == nVertices);
        assert(field.cols() == 3);
        
        vtkSmartPointer<vtkDoubleArray> vector_array = convertVectorFieldToVTK(field, fieldname);
        polydata->GetPointData()->AddArray(vector_array);
    }
    
    void addScalarFieldToVertices(const Eigen::Ref<const Eigen::VectorXd> field, const std::string fieldname)
    {
        assert(field.rows() == nVertices);
        
        vtkSmartPointer<vtkDoubleArray> scalar_array = convertScalarFieldToVTK(field, fieldname);
        polydata->GetPointData()->AddArray(scalar_array);
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
#else /* USEVTK */

    // ASCII! not for larger domains

class WriteVTK
{
protected:
    const int nVertices;
    const int nFaces;
    const Eigen::Ref<const Eigen::MatrixXd> verts;
    const Eigen::Ref<const Eigen::MatrixXi> faces;
    
public:
    WriteVTK(const Eigen::Ref<const Eigen::MatrixXd> verts, const Eigen::Ref<const Eigen::MatrixXi> faces):
    nVertices(verts.rows()),
    nFaces(faces.rows()),
    verts(verts),
    faces(faces)
    {
        assert(nVertices > 0);
        assert(nFaces > 0);
        assert(verts.cols() == 3);
        assert(faces.cols() == 3);

    }
    
    void addVectorFieldToFaces(const Eigen::Ref<const Eigen::MatrixXd> , const std::string )
    {
        std::cout << "addVectorFieldToFaces not implemented for WriteVTK if not compiled with VTK support\n";
    }
    
    void addScalarFieldToFaces(const Eigen::Ref<const Eigen::VectorXd> , const std::string )
    {
        std::cout << "addScalarFieldToFaces not implemented for WriteVTK if not compiled with VTK support\n";
    }
    
    void addVectorFieldToVertices(const Eigen::Ref<const Eigen::MatrixXd> , const std::string )
    {
        std::cout << "addVectorFieldToVertices not implemented for WriteVTK if not compiled with VTK support\n";
    }
    
    void addScalarFieldToVertices(const Eigen::Ref<const Eigen::VectorXd> , const std::string )
    {
        std::cout << "addScalarFieldToVertices not implemented for WriteVTK if not compiled with VTK support\n";
    }
    
    int write(const std::string filename)
    {
        std::ofstream ofs(filename+".vtk");
        if (not ofs.is_open()) return 0;
        
        ofs << "# vtk DataFile Version 1.0" << std::endl;
        ofs << "2D Unstructured Grid of Linear Triangles" << std::endl;
        ofs << "ASCII" << std::endl;
        ofs << std::endl;
        ofs << "DATASET UNSTRUCTURED_GRID" << std::endl;
        ofs << "POINTS "+std::to_string(nVertices)+" float" << std::endl;
        
        for(int i=0;i<nVertices;++i)
            ofs << verts(i,0) << " " << verts(i,1) << " " << verts(i,2) << std::endl;
        
        ofs << std::endl;
        ofs << "CELLS "+std::to_string(nFaces)+" "+std::to_string(nFaces*4) << std::endl;
        for(int i=0;i<nFaces;++i)
            ofs << 3 << " " << faces(i,0) << " " << faces(i,1) << " " << faces(i,2) << std::endl;
        
        ofs<<std::endl;
        ofs << "CELL_TYPES "+std::to_string(nFaces) << std::endl;
        for(int i=0;i<nFaces;++i)
            ofs << "5" << std::endl;
        
        ofs.close();
        
        return 1;
    }
};
#endif /* USEVTK */

#endif /* WriteVTK_hpp */

