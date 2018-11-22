#pragma once
#include <vector>
#include "Vec3.h"

struct Node
{
    Vec3f position;
    Edge *previous_edge, *next_edge;
};
struct Edge
{
    Node *n_source, *n_target;
    Edge *previous_edge, *next_edge;
    float length;
};

class TreeGraph
{
public:
    Node * root_node;
    int nodes_count, edges_count;
    Node * nodes;
    Edge * edges;
    inline TreeGraph(){
      // Sample tree g
      this->nodes_count = 3;
      this->edges_count = 2;
      this->nodes = (Node*) malloc(nodes_count*sizeof(Node));
      this->edges = (Edge*) malloc(edges_count*sizeof(Edge));

      this->nodes[0].position = Vec3f();
      this->nodes[1].position = Vec3f(5.f,0.f,0.f);
      this->nodes[2].position = Vec3f(1.f,1.f,0.f);

      this->edges[0].n_source = &nodes[0];
      this->edges[0].n_target = &nodes[1];
      this->edges[0].length = 1.f;
      this->edges[0].previous_edge = NULL;
      this->edges[0].next_edge = NULL;

      this->edges[1].n_source = &nodes[0];
      this->edges[1].n_target = &nodes[2];
      this->edges[1].length = 0.5f;
      this->edges[1].previous_edge = NULL;
      this->edges[1].next_edge = NULL;

      this->root_node = &nodes[0];
    }
    inline TreeGraph(Node * root_node, 
      int nodes_count, 
      Node * nodes,
      int edges_count,
      Edge * edges): 
        root_node(root_node),
        nodes(nodes), 
        nodes_count(nodes_count),
        edges(edges), 
        edges_count(edges_count) {};
    inline ~TreeGraph() {};
};

