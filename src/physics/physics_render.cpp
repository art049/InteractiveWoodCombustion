#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "physics.h"
#include <iostream>
#include "../Vec3.h"

inline int flatten(int col, int row, int z) {
    return col + row*(GRID_COUNT+1) + z*(GRID_COUNT+1)*(GRID_COUNT+1);
}

void Physics::initSmokeQuads(){
    int nflat = grid3d->NFLAT();
    smokeQuadsPositions = new float[3* 4*3*nflat];
    smokeQuadsColors = new float[4*4*nflat];
    smokeIndexes = new uint[4*nflat];
    for(uint x = 0; x < GRID_COUNT; x++){
        for(uint y = 0; y < GRID_COUNT; y++){
            for(uint z = 0; z < GRID_COUNT; z++){
                std::array<float,12> vertexes = {
                    x*BLOCK_SIZE, y*BLOCK_SIZE, z*BLOCK_SIZE,
                    x*BLOCK_SIZE, y*BLOCK_SIZE, (z+1)*BLOCK_SIZE,
                    x*BLOCK_SIZE, (y+1)*BLOCK_SIZE, (z+1)*BLOCK_SIZE,
                    x*BLOCK_SIZE, (y+1)*BLOCK_SIZE, z*BLOCK_SIZE
                };
                std::copy(vertexes.begin(), vertexes.end(), 
                          smokeQuadsPositions +  12*grid3d->flatten(x,y,z));
            }
        }
    }
    for(uint x = 0; x < GRID_COUNT; x++){
        for(uint y = 0; y < GRID_COUNT; y++){
            for(uint z = 0; z < GRID_COUNT; z++){
                std::array<float,12> vertexes = {
                    x*BLOCK_SIZE, y*BLOCK_SIZE, z*BLOCK_SIZE,
                    (x+1)*BLOCK_SIZE, y*BLOCK_SIZE, z*BLOCK_SIZE,
                    (x+1)*BLOCK_SIZE, y*BLOCK_SIZE, (z+1)*BLOCK_SIZE,
                    x*BLOCK_SIZE, y*BLOCK_SIZE, (z+1)*BLOCK_SIZE
                };
                std::copy(vertexes.begin(), vertexes.end(), 
                          smokeQuadsPositions + 1* 4*3*nflat +  12*grid3d->flatten(x,y,z));
            }
        }
    }
    for(uint x = 0; x < GRID_COUNT; x++){
        for(uint y = 0; y < GRID_COUNT; y++){
            for(uint z = 0; z < GRID_COUNT; z++){
                std::array<float,12> vertexes = {
                    x*BLOCK_SIZE, y*BLOCK_SIZE, z*BLOCK_SIZE,
                    (x+1)*BLOCK_SIZE, y*BLOCK_SIZE, z*BLOCK_SIZE,
                    (x+1)*BLOCK_SIZE, (y+1)*BLOCK_SIZE, z*BLOCK_SIZE,
                    x*BLOCK_SIZE, (y+1)*BLOCK_SIZE, z*BLOCK_SIZE
                };
                std::copy(vertexes.begin(), vertexes.end(), 
                          smokeQuadsPositions + 2* 4*3*nflat +  12*grid3d->flatten(x,y,z));
            }
        }
    }
    

    glGenBuffers(1, &smokeQuadVBO);
    glBindBuffer(GL_ARRAY_BUFFER, smokeQuadVBO);
    glBufferData(GL_ARRAY_BUFFER, 3*nflat*4*3*sizeof(float), smokeQuadsPositions, GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glGenBuffers(1, &smokeQuadIndexBO);
    for(int i=0; i< 4*grid3d->NFLAT(); i++) smokeIndexes[i] = i;
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, smokeQuadIndexBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 4*nflat*sizeof(uint), smokeIndexes, GL_STATIC_DRAW);
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void Physics::renderSmokeQuads(uint cameraAxis){
    glDisable(GL_CULL_FACE);
    glClear(GL_DEPTH_BUFFER_BIT);
    glDisable(GL_DEPTH_TEST);
    glBindBuffer(GL_ARRAY_BUFFER, smokeQuadVBO);
    glVertexPointer (3, GL_FLOAT, 0, 4*3*grid3d->NFLAT() );
    glBindBuffer(GL_ARRAY_BUFFER, smokeColorBufferObj);
    glColorPointer (4, GL_UNSIGNED_BYTE, 0, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, smokeQuadIndexBO);
    glDrawElements (GL_QUADS,/*4**/grid3d->NFLAT(), GL_UNSIGNED_INT, (void*) 0);
    
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
}

void Physics::renderGrid(){
    int sparsity = 20;
    glBegin(GL_LINES);
    glColor3f(0,0,1);
    for(uint x = 0; x <= GRID_COUNT; x+=sparsity)
    {
        for(uint y = 0; y <= GRID_COUNT; y+=sparsity)
        {
            glVertex3f(x*BLOCK_SIZE, y*BLOCK_SIZE, 0);
            glVertex3f(x*BLOCK_SIZE, y*BLOCK_SIZE, GRID_SIZE);
        }        
    }
    glColor3f(1,0,0);
    for(uint z = 0; z <= GRID_COUNT; z+=sparsity)
    {
        for(uint y = 0; y <= GRID_COUNT; y+=sparsity)
        {
            glVertex3f(0, y*BLOCK_SIZE, z*BLOCK_SIZE);
            glVertex3f(GRID_SIZE, y*BLOCK_SIZE, z*BLOCK_SIZE);
        }        
    }
    glColor3f(0,1,0);
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

void Physics::renderLightRays(){
    const float tdefault = 0.3;
    float t;
    int3 voxel;
    glBegin(GL_LINES);
    for(int k_x = 0; k_x < SMOKE_RAY_SQRT_COUNT; k_x++){
        for(int k_y = 0; k_y < SMOKE_RAY_SQRT_COUNT; k_y++){
            vec3 ray_orig = vec3(SMOKE_LIGHT_POS) + vec3(0, (k_x+0.5) * BLOCK_SIZE, (k_y+0.5) * BLOCK_SIZE);
            vec3 rayDir = vec3(SMOKE_LIGHT_DIR);

            bool inter = rayGridIntersect(ray_orig, rayDir, &voxel, &t);
            if(!inter){
                glColor3f(1,0,0);
                glVertex3fv((GLfloat *) &ray_orig);
                glVertex3f(ray_orig.x() + tdefault * rayDir.x(),
                           ray_orig.y() + tdefault * rayDir.y(),
                           ray_orig.z() + tdefault * rayDir.z());
            }
            else {
                glColor3f(0,1,0);
                glVertex3fv((GLfloat *) &ray_orig);
                glVertex3f(ray_orig.x() + t * rayDir.x(),
                           ray_orig.y() + t * rayDir.y(),
                           ray_orig.z() + t * rayDir.z());
            }
        }
    }
    glEnd();
}

void Physics::renderExternalForce(){
    glBegin(GL_LINES);
        glColor3f(1,0,0);
        glVertex3f(0,0,0);
        glVertex3fv((GLfloat*) &externalForce);
    glEnd();
}

void Physics::render(uint cameraAxis) {
    if(gridEnabled) renderGrid();
    if(raysEnabled) renderLightRays();
    renderSmokeQuads(cameraAxis);
    renderExternalForce();
}
