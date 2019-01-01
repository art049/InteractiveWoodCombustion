#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "physics.h"
#include <iostream>
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
                glm::vec3 pos(x*hx, y*hy, z*hz); 
                center = center * 0.5f;
                float alpha = expf( (-10.f)* glm::length(center - pos));
                std::array<float,4> color = {1.f,1.f,1.f,alpha};
                for(uint i =0; i< 4; i++){
                    std::copy(color.begin(), color.end(), smokeQuadsColors + (4*i+16*grid3d->flatten(x,y,z)));
                }
            }
        }
    }
    for(uint i=0; i< 4*grid3d->NFLAT(); i++) smokeIndexes[i] = i;
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

void Physics::render() {
    renderSmokeQuads();
}
