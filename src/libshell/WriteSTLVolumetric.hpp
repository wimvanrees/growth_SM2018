//
//  WriteSTLVolumetric.h
//  Elasticity
//
//  Created by Wim van Rees on 8/2/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef WriteSTLVolumetric_h
#define WriteSTLVolumetric_h

#include "common.hpp"
#include "IndexMapper.hpp"
#include <igl/boundary_facets.h>
#include <igl/unique_edge_map.h>
#include <igl/extract_manifold_patches.h>
#include <igl/writeSTL.h>

class WriteSTLVolumetric
{
public:
    static void write(const Eigen::MatrixXi & face2vertices, const Eigen::MatrixXd & vertices, const Eigen::VectorXi & faces_STL, const Real thickness, const std::string filename)
    {
        // face2vertices gives the vertex indices for each face for the entire mesh (nFaces x 3)
        // vertices gives the vertex locations for the entire mesh (nVertices x 3)
        // faces_STL is a vector with the face indices that need to be saved out of all faces in the mesh (to store a subset)
        
        
        std::cout << " === Writing volumetric STL file \t " << filename << std::endl;
        const int nFaces_STL = faces_STL.rows();
        std::cout << "Midsurface selection has \t " << nFaces_STL << " \t faces" << std::endl;
        Eigen::MatrixXi face2vertices_STL(nFaces_STL,3); // this will be our face2vertices list for the midsurface faces that need to be stored
        
        // if we extract the faces that we want to store out of the global array the corresponding vertex indices will generally be sparse and scattered, rather than a contiguous dense list. Therefore we remap the vertex indices that we need to a unique, compact list
        IndexMapper indexmapper;
        for(int i=0;i<nFaces_STL;++i)
        {
            const int faceidx = faces_STL(i);
            // global vertex indices
            const int vidx_0 = face2vertices(faceidx,0);
            const int vidx_1 = face2vertices(faceidx,1);
            const int vidx_2 = face2vertices(faceidx,2);
            // remapped vertex indices
            face2vertices_STL(i,0) = indexmapper.mapIndex(vidx_0);
            face2vertices_STL(i,1) = indexmapper.mapIndex(vidx_1);
            face2vertices_STL(i,2) = indexmapper.mapIndex(vidx_2);
        }
        const int vertices_cnt = indexmapper.getTotalIndices(); // this is the total number of vertices we mapped
        std::cout << "Midsurface selection has \t " << vertices_cnt << "\t unique vertices " << std::endl;
        
        // we have to create the manifold patches: a list of (list of connecting faces that form a single volume)
        Eigen::MatrixXi E, uE;
        Eigen::VectorXi EMAP;
        std::vector<std::vector<size_t> > uE2E;
        igl::unique_edge_map(face2vertices_STL, E, uE, EMAP, uE2E);
        Eigen::VectorXi face2patch;
        const size_t npatches = igl::extract_manifold_patches(face2vertices_STL, EMAP, uE2E, face2patch);
        std::cout << "Midsurface selection has \t " << npatches << " \t individual manifold patches" << std::endl;
        
        // face2patch contains for each face the patch it belongs to, but we want to have for each patch the list of faces belonging to it : reverse the storage
        std::vector<std::vector<int>> patch2faces(npatches);
        for(int i=0;i<face2vertices_STL.rows();++i)
            patch2faces[face2patch(i)].push_back(i);
        
        // now we can create the data structures for the 3D STL file : a face2vertices instance that contains all the faces, and a vertices instance that contains all the vertices
        // first we create all the vertices
        Eigen::MatrixXd vertices_STL_3D(2*vertices_cnt,3);// twice the vertices because we have two layers
        for(int i=0;i<vertices_cnt;++i)
        {
            // to get the vertex coordinates we have to look in the global array, which means we have to unmap the local index to retrieve the global one
            const int unmap_idx = indexmapper.unMapIndex(i);
            
            // layer one
            for(int d=0;d<3;++d)
                vertices_STL_3D(i,d) = vertices(unmap_idx,d) - (d==2 ? 0.5*thickness : 0.0);
            // layer two
            for(int d=0;d<3;++d)
                vertices_STL_3D(i+vertices_cnt,d) = vertices(unmap_idx,d) + (d==2 ? 0.5*thickness : 0.0);
        }
        
        // now we create the faces. We dont know yet how many faces there will be since we dont know the boundaries for each of the patches. Threfore we create a temporary dynamic structures to store the faces for each patch, and then later combine everything into a big array
        std::vector<std::vector<std::array<int, 3>>> face2vertices_STL_3D_container(npatches);
        for(size_t p=0;p<npatches;++p)
        {
            const int nFacesInPatch = patch2faces[p].size();
            
            // want to know the boundary edges for this patch : use igl::boundary_facets. However that function needs again an ordered list of faces, so we do the mapping trick once more
            Eigen::MatrixXi face2vertices_patch(nFacesInPatch, 3);
            IndexMapper idxmapper_vertices_patch;
            for(int i=0;i<nFacesInPatch;++i)
            {
                const int face_idx = patch2faces[p][i];
                const int vidx_0 = face2vertices_STL(face_idx,0);
                const int vidx_1 = face2vertices_STL(face_idx,1);
                const int vidx_2 = face2vertices_STL(face_idx,2);
                
                face2vertices_patch(i,0) = idxmapper_vertices_patch.mapIndex(vidx_0);
                face2vertices_patch(i,1) = idxmapper_vertices_patch.mapIndex(vidx_1);
                face2vertices_patch(i,2) = idxmapper_vertices_patch.mapIndex(vidx_2);
            }
            
            // now we can identify the boundary edges
            Eigen::MatrixXi boundary_edges;
            igl::boundary_facets(face2vertices_patch, boundary_edges);
            const int nBoundaries = boundary_edges.rows();
            
            // the total number of faces is twice the midsurface faces, plus twice the number of boundary edges (each pair of boundary edges will be triangulated using 2 right-angled triangles)
            const int nFacesInPatch_3D = 2*nFacesInPatch + 2*nBoundaries;
            face2vertices_STL_3D_container[p].resize(nFacesInPatch_3D);
            
            // now we store the faces, first for the bottom and top layers
            for(int i=0;i<nFacesInPatch;++i)
            {
                // we need to unmap the vertices again to be consistent with vertices_STL_3D
                const int vidx_0 = idxmapper_vertices_patch.unMapIndex(face2vertices_patch(i,0));
                const int vidx_1 = idxmapper_vertices_patch.unMapIndex(face2vertices_patch(i,1));
                const int vidx_2 = idxmapper_vertices_patch.unMapIndex(face2vertices_patch(i,2));
                
                // create the two layers
                std::array<int, 3> face_bot = {vidx_0, vidx_1, vidx_2};
                std::array<int, 3> face_top = {vidx_0 + vertices_cnt, vidx_1 + vertices_cnt, vidx_2 + vertices_cnt}; // since the top layer vertices are the bottom one shifted by vertices_cnt
                
                face2vertices_STL_3D_container[p][i] = face_bot;
                face2vertices_STL_3D_container[p][i + nFacesInPatch] = face_top;
            }
            // now we have filled face2vertices_STL_3D_container[p] from 0 to 2*nFacesInPatch with the bottom/top layers : need to add the edges next
            for(int i=0;i<nBoundaries;++i)
            {
                // unmap to be consistent with vertices_STL_3D
                const int vidx_0 = idxmapper_vertices_patch.unMapIndex(boundary_edges(i,0));
                const int vidx_1 = idxmapper_vertices_patch.unMapIndex(boundary_edges(i,1));
                
                // create the boundary triangles
                std::array<int, 3> face_boundary_1 = {vidx_0, vidx_1, vidx_1 + vertices_cnt};
                std::array<int, 3> face_boundary_2 = {vidx_0, vidx_1 + vertices_cnt, vidx_0 + vertices_cnt};
                
                // store the boundary faces after the surface faces
                face2vertices_STL_3D_container[p][2*nFacesInPatch + 2*i + 0] = face_boundary_1;
                face2vertices_STL_3D_container[p][2*nFacesInPatch + 2*i + 1] = face_boundary_2;
            }
        }
        
        // count the total number of faces we have created
        int nFaces_STL_3D = 0;
        for(size_t p=0;p<npatches;++p)
            nFaces_STL_3D += face2vertices_STL_3D_container[p].size();
        
        // create finally the global face2vertices_STL_3D structure that contains all the faces
        Eigen::MatrixXi face2vertices_STL_3D(nFaces_STL_3D,3);
        int face2vertices_STL_3D_idx = 0;
        for(size_t p=0;p<npatches;++p)
        {
            const size_t nFacesInPatch_3D = face2vertices_STL_3D_container[p].size();
            for(size_t i=0;i<nFacesInPatch_3D;++i,++face2vertices_STL_3D_idx)
                for(int d=0;d<3;++d)
                    face2vertices_STL_3D(face2vertices_STL_3D_idx,d) = face2vertices_STL_3D_container[p][i][d];
        }
        
        // fill the normals array : not sure why I need this or why this should be correct - does not seem to matter
        Eigen::MatrixXd facenormals(nFaces_STL_3D,3);
        facenormals.setZero();
        for(int i=0;i<nFaces_STL_3D;++i)
            facenormals(i,2) = 1.0;
        
        // write
        //igl::writeSTL(filename+".stl", vertices_STL_3D, face2vertices_STL_3D, facenormals, false); // no ascii (last argument)        
        igl::writeSTL(filename+".stl", vertices_STL_3D, face2vertices_STL_3D, facenormals, igl::FileEncoding::Binary);        
    }
};

#endif /* WriteSTLVolumetric_h */
