#pragma once
#include <vector>
#include "Vec3.h"

struct Edge;
struct Node
{
    Vec3f position;
    Edge *previous_edge, *next_edge;
};
struct Edge
{
    Node *n_source, *n_target;
    Edge *previous_edge, *next_edge;
};

class TreeGraph
{
public:
    int nodes_count;
    int edges_count;
    Node * root_node;
    Node * nodes;
    Edge * edges;
    inline TreeGraph(){
      // Sample tree graph
      this->nodes_count =10;
      this->edges_count = 5;
      this->nodes = (Node*) malloc(nodes_count*sizeof(Node));
      this->edges = (Edge*) malloc(edges_count*sizeof(Edge));

      this->nodes[0].position = Vec3f();
      this->nodes[1].position = Vec3f(0.f,5.f,0.f);
      
      this->nodes[2].position = Vec3f(0.f,1.f,0.f);
      this->nodes[3].position = Vec3f(1.f,2.f,0.f);
      
      this->nodes[4].position = Vec3f(0.f,2.f,0.f);
      this->nodes[5].position = Vec3f(-1.f,3.f,-1.f);
      
      this->nodes[6].position = Vec3f(0.f,3.f,0.f);
      this->nodes[7].position = Vec3f(1.5f,4.f,1.f);

      this->nodes[8].position = Vec3f(0.f,4.f,0.f);
      this->nodes[9].position = Vec3f(-.5f,5.f,-1.f);

      this->edges[0].n_source = &nodes[0];
      this->edges[0].n_target = &nodes[1];
      this->edges[1].n_source = &nodes[2];
      this->edges[1].n_target = &nodes[3];
      this->edges[2].n_source = &nodes[4];
      this->edges[2].n_target = &nodes[5];
      this->edges[3].n_source = &nodes[6];
      this->edges[3].n_target = &nodes[7];
      this->edges[4].n_source = &nodes[8];
      this->edges[4].n_target = &nodes[9];

      this->root_node = &nodes[0];
    }
    inline TreeGraph(Node * root_node, 
      int nodes_count, 
      Node * nodes,
      int edges_count,
      Edge * edges): 
        edges_count(edges_count),
        nodes_count(nodes_count),
        root_node(root_node),
        nodes(nodes), 
        edges(edges) {};
    inline ~TreeGraph() {};
};

