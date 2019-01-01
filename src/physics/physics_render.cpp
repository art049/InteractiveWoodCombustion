#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "physics.h"
#include <iostream>
#include "../Vec3.h"

void Physics::initSmokeQuads(){
    smokeQuadsPositions = new float[4*3*grid3d->NFLAT()];
    smokeQuadsColors = new float[4*4*grid3d->NFLAT()];
    smokeIndexes = new uint[4*grid3d->NFLAT()];
    float hx = grid3d->hd[0];
    float hy = grid3d->hd[1];
    float hz = grid3d->hd[2];
    int sparsity = 1;
    for(uint x = 0; x < GRID_WIDTH; x+=sparsity)
    {
        for(uint y = 0; y < GRID_HEIGHT; y+=sparsity)
        {
            for(uint z = 0; z < GRID_DEPTH; z+=sparsity)
            {
                std::array<float,12> vertexes = {
                    x*hx, y*hy, z*hz,
                    (x+sparsity)*hx, y*hy, z*hz,
                    (x+sparsity)*hx, (y+sparsity)*hy, z*hz,
                    x*hx, (y+sparsity)*hy, z*hz
                };
                std::copy(vertexes.begin(), vertexes.end(), smokeQuadsPositions + 12*grid3d->flatten(x,y,z));

                glm::vec3 center = glm::make_vec3(grid3d->ld.begin());
                glm::vec3 pos((x+0.5)*hx, (y+0.5)*hy, (z+0.5)*hz); 
                center = center * 0.5f;
                float alpha = expf( (-10.f)* glm::length(center - pos));
                std::array<float,4> color = {1.f,1.f,1.f,alpha};
                for(uint i =0; i< 4; i++){
                    std::copy(color.begin(), color.end(), smokeQuadsColors + (4*i+16*grid3d->flatten(x,y,z)));
                }
            }
        }
    }
    for(int i=0; i< 4*grid3d->NFLAT(); i++) smokeIndexes[i] = i;
}

void Physics::renderSmokeQuads(){
    glDisable(GL_CULL_FACE);
    glClear(GL_DEPTH_BUFFER_BIT);
    glDisable(GL_DEPTH_TEST);
    glVertexPointer (3, GL_FLOAT, 3*sizeof(float), (GLvoid*)(smokeQuadsPositions));
    glColorPointer (4, GL_FLOAT, 0, (GLvoid*)(smokeQuadsColors));
    glDrawElements (GL_QUADS,4*grid3d->NFLAT(), GL_UNSIGNED_INT, smokeIndexes);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
}

void Physics::renderGrid(){
    glColor3f(.9,.9,.5);
    Vec3f gridSize(1,1,1);
    float hx = gridSize[0]/GRID_WIDTH;
    float hy = gridSize[1]/GRID_HEIGHT;
    float hz = gridSize[0]/GRID_DEPTH;
    int sparsity = 20;
    glBegin(GL_LINES);
    for(uint x = 0; x <= GRID_WIDTH; x+=sparsity)
    {
        for(uint y = 0; y <= GRID_HEIGHT; y+=sparsity)
        {
            glVertex3f(x*hx, y*hy, 0);
            glVertex3f(x*hx, y*hy, gridSize[2]);
        }        
    }
    for(uint z = 0; z <= GRID_DEPTH; z+=sparsity)
    {
        for(uint y = 0; y <= GRID_HEIGHT; y+=sparsity)
        {
            glVertex3f(0, y*hy, z*hz);
            glVertex3f(gridSize[0], y*hy, z*hz);
        }        
    }
    for(uint z = 0; z <= GRID_DEPTH; z+=sparsity)
    {
        for(uint x = 0; x <= GRID_WIDTH; x+=sparsity)
        {
            glVertex3f(x*hx, 0, z*hz);
            glVertex3f(x*hx, gridSize[1], z*hz);
        }        
    }
    glEnd();
}

void Physics::render() {
    if(gridEnabled) renderGrid();
    renderSmokeQuads();
}
