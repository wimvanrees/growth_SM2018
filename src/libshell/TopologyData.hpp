//
//  TopologyData.hpp
//  Elasticity
//
//  Created by Wim van Rees on 23/04/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef TopologyData_hpp
#define TopologyData_hpp

#include <igl/edge_topology.h>
#include <igl/vertex_triangle_adjacency.h>
#include "common.hpp"

struct TopologyData
{
    // logical information
    Eigen::MatrixXi face2vertices;
    Eigen::MatrixXi edge2vertices;
    Eigen::MatrixXi edge2oppositevertices;
    Eigen::MatrixXi edge2faces;
    Eigen::MatrixXi face2edges;
    std::vector<std::vector<int> > vertex2faces;
    
    void clear()
    {
        face2vertices.resize(0,0);
        edge2vertices.resize(0,0);
        edge2oppositevertices.resize(0,0);
        edge2faces.resize(0,0);
        face2edges.resize(0,0);
        vertex2faces.clear();
    }
    
    int getNumberOfFaces() const
    {
        return face2vertices.rows();
    }

    int getNumberOfEdges() const
    {
        return edge2vertices.rows();
    }
    
    void init(const Eigen::MatrixXd & vertices_in, const Eigen::MatrixXi & faces_in)
    {
        // directly set face2vertices
        face2vertices = faces_in;
        
        // compute edge topology
        igl::edge_topology(vertices_in,face2vertices,edge2vertices,face2edges,edge2faces);
        
        // compute triangly adjacency
        std::vector<std::vector<int>> dummy; // vertex2faces_indices --> not needed
        igl::vertex_triangle_adjacency(vertices_in.rows(), face2vertices, vertex2faces, dummy);
        
        updateOppositeVertices();
    }
    
    void updateOppositeVertices()
    {
        if(edge2oppositevertices.rows() != edge2vertices.rows())
            edge2oppositevertices.resize(edge2vertices.rows(),2);
        
        // find the vertex opposite to each edge (or put -1 if it does not exist)
        
        for(int i=0;i<edge2oppositevertices.rows();++i)
        {
            // these are the two vertex indices on either side of the edge
            const int idx1 = edge2vertices(i,0);
            const int idx2 = edge2vertices(i,1);
            
            // do the first face
            for(int j=0;j<2;++j)
            {
                const int face_idx = edge2faces(i,j);
                if(face_idx < 0)
                    edge2oppositevertices(i,j) = -1; // face does not exist
                else
                {
                    // get the three face vertices
                    const int v_idx1 = face2vertices(face_idx,0);
                    const int v_idx2 = face2vertices(face_idx,1);
                    const int v_idx3 = face2vertices(face_idx,2);
                    
                    // find the one face vertex that is not part of the edge
                    if( (v_idx1 == idx1 && v_idx2 == idx2) or (v_idx1 == idx2 && v_idx2 == idx1) )
                        edge2oppositevertices(i,j) = v_idx3;
                    else if( (v_idx1 == idx1 && v_idx3 == idx2) or (v_idx1 == idx2 && v_idx3 == idx1) )
                        edge2oppositevertices(i,j) = v_idx2;
                    else
                    {
                        assert( (v_idx2 == idx1 && v_idx3 == idx2) or (v_idx2 == idx2 && v_idx3 == idx1) );
                        edge2oppositevertices(i,j) = v_idx1;
                    }
                }
            }
        }
    }
    
    const Eigen::Ref<const Eigen::MatrixXi> getFace2Vertices() const
    {
        return face2vertices;
    }
    
    const Eigen::Ref<const Eigen::MatrixXi> getEdge2Vertices() const
    {
        return edge2vertices;
    }
    
    const Eigen::Ref<const Eigen::MatrixXi> getEdge2OppositeVertices() const
    {
        return edge2oppositevertices;
    }
    
    const Eigen::Ref<const Eigen::MatrixXi> getEdge2Faces() const
    {
        return edge2faces;
    }
    
    const Eigen::Ref<const Eigen::MatrixXi> getFace2Edges() const
    {
        return face2edges;
    }
    
    const std::vector<std::vector<int>> & getVertex2Faces() const
    {
        return vertex2faces;
    }
};


#endif /* TopologyData_hpp */
