#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "physics.h"
#include <iostream>
#include "../Vec3.h"

void Physics::initSmokeQuads(){
    smokeQuadsPositions = new float[4*3*grid3d->NFLAT()];
    smokeQuadsColors = new float[4*4*grid3d->NFLAT()];
    smokeIndexes = new uint[4*grid3d->NFLAT()];
    int sparsity = 1;
    for(uint x = 0; x < GRID_COUNT; x+=sparsity)
    {
        for(uint y = 0; y < GRID_COUNT; y+=sparsity)
        {
            for(uint z = 0; z < GRID_COUNT; z+=sparsity)
            {
                std::array<float,12> vertexes = {
                    x*BLOCK_SIZE, y*BLOCK_SIZE, z*BLOCK_SIZE,
                    (x+sparsity)*BLOCK_SIZE, y*BLOCK_SIZE, z*BLOCK_SIZE,
                    (x+sparsity)*BLOCK_SIZE, (y+sparsity)*BLOCK_SIZE, z*BLOCK_SIZE,
                    x*BLOCK_SIZE, (y+sparsity)*BLOCK_SIZE, z*BLOCK_SIZE
                };
                std::copy(vertexes.begin(), vertexes.end(), smokeQuadsPositions + 12*grid3d->flatten(x,y,z));

                glm::vec3 center = glm::make_vec3(grid3d->ld.begin());
                glm::vec3 pos((x+0.5)*BLOCK_SIZE, (y+0.5)*BLOCK_SIZE, (z+0.5)*BLOCK_SIZE); 
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
    int sparsity = 20;
    glBegin(GL_LINES);
    for(uint x = 0; x <= GRID_COUNT; x+=sparsity)
    {
        for(uint y = 0; y <= GRID_COUNT; y+=sparsity)
        {
            glVertex3f(x*BLOCK_SIZE, y*BLOCK_SIZE, 0);
            glVertex3f(x*BLOCK_SIZE, y*BLOCK_SIZE, GRID_SIZE);
        }        
    }
    for(uint z = 0; z <= GRID_COUNT; z+=sparsity)
    {
        for(uint y = 0; y <= GRID_COUNT; y+=sparsity)
        {
            glVertex3f(0, y*BLOCK_SIZE, z*BLOCK_SIZE);
            glVertex3f(GRID_SIZE, y*BLOCK_SIZE, z*BLOCK_SIZE);
        }        
    }
    for(uint z = 0; z <= GRID_COUNT; z+=sparsity)
    {
        for(uint x = 0; x <= GRID_COUNT; x+=sparsity)
        {
            glVertex3f(x*BLOCK_SIZE, 0, z*BLOCK_SIZE);
            glVertex3f(x*BLOCK_SIZE, GRID_SIZE, z*BLOCK_SIZE);
        }        
    }
    glEnd();
}

void Physics::render() {
    if(gridEnabled) renderGrid();
    renderSmokeQuads();
}
