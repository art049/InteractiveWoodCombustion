// --------------------------------------------------------------------------
// Copyright(C) 2009-2016
// Tamy Boubekeur
// 
// Permission granted to use this code only for teaching projects and 
// private practice.
//
// Do not distribute this code outside the teaching assignements.                                                                           
// All rights reserved.                                                       
// --------------------------------------------------------------------------

#include "TreeGraph.h"
#include "Tree.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>

using namespace std;

Tree::Tree(TreeGraph * graph) {
    const int circular_subdiv = 3;
    const float radius = .08f;
    m_positions.resize(2 * circular_subdiv * graph->edges_count);
    m_triangles.resize(2 * circular_subdiv * graph->edges_count);
    for(int i = 0; i < graph->edges_count; i++){
        Edge * e = &graph->edges[i];
        Vec3f axis = e->n_target->position - e->n_source->position;
        Vec3f orth_ref = axis.orthogonal();
        orth_ref.normalize();
        std::cout<< orth_ref << std::endl;
        orth_ref *= radius;
        unsigned int index = 2 * i * circular_subdiv;
        for(int j = 0; j < circular_subdiv; j++){
            int angle = j * (360 / circular_subdiv);
            Vec3f current = orth_ref.rotated(axis, angle);
            m_positions[index + j] = current + e->n_source->position;
            m_positions[index + circular_subdiv + j] = current + e->n_target->position;
            m_triangles[index + j] = Triangle(index+j, 
                                              index+(j+1)%circular_subdiv, 
                                              index+j+circular_subdiv);
            m_triangles[index + circular_subdiv + j] = Triangle(index+(j+1)%circular_subdiv,
                                                                index+circular_subdiv+(j+2)%circular_subdiv,
                                                                index+j+circular_subdiv);
        }        
    }
    centerAndScaleToUnit ();
    recomputeNormals ();
}


void Tree::clear () {
    m_positions.clear ();
    m_normals.clear ();
    m_triangles.clear ();
}
void Tree::recomputeNormals () {
    m_normals.clear ();
    m_normals.resize (m_positions.size (), Vec3f (1.f, 0.f, 0.f));
    for (unsigned int i = 0; i < m_triangles.size (); i++) {
        Vec3f e01 = m_positions[m_triangles[i][1]] -  m_positions[m_triangles[i][0]];
        Vec3f e02 = m_positions[m_triangles[i][2]] -  m_positions[m_triangles[i][0]];
        Vec3f n = cross (e01, e02);
        n.normalize ();
        for (unsigned int j = 0; j < 3; j++)
            m_normals[m_triangles[i][j]] += n;
    }
    for (unsigned int i = 0; i < m_normals.size (); i++)
        m_normals[i].normalize ();
}

void Tree::centerAndScaleToUnit () {
    Vec3f c;
    for  (unsigned int i = 0; i < m_positions.size (); i++)
        c += m_positions[i];
    c /= m_positions.size ();
    float maxD = dist (m_positions[0], c);
    for (unsigned int i = 0; i < m_positions.size (); i++){
        float m = dist (m_positions[i], c);
        if (m > maxD)
            maxD = m;
    }
    for  (unsigned int i = 0; i < m_positions.size (); i++)
        m_positions[i] = (m_positions[i] - c) / maxD;
}
