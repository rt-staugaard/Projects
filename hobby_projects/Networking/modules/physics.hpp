#pragma once
#include <cmath>
#include <thread>
#include "config.hpp"

struct Grid{
    float Xs[gridLength * gridHeight];
    float Ys[gridLength * gridHeight];

    Grid(){
        for (int y = 0; y < gridHeight; y++) {
            for (int x = 0; x < gridLength; x++) {
                int i = x + y * gridLength;
            
                Xs[i] = 0.5f + (float)x; 
                Ys[i] = 0.5f + (float)y;
            }
        }
    }
};

int findIdx(int x, int y){
    return x + y * gridLength;
    //return (x + gridLength) % gridLength + ((y + gridHeight) % gridHeight) * gridLength;
}

void updateRows(int startRow, int endRow, Grid* curr, Grid* prev, Grid* buffer, float* forceX, float* forceY){
    float L = 1.0f; 
    float k = 0.1f; 
    float inertia = 1; 
    float dt = 0.1;

    for (int y = startRow; y < endRow; y++) {
        for (int x = 0; x < gridLength; x++) {
            int idx = x + y * gridLength;

            if (y == 0 || y == gridHeight - 1 || x == 0 || x == gridLength - 1){
                buffer->Xs[idx] = curr->Xs[idx];
                buffer->Ys[idx] = curr->Ys[idx];
                continue;
            }

            float totalFx = 0, totalFy = 0;

            int neighbors[4] = { findIdx(x, y-1), findIdx(x, y+1), findIdx(x-1, y), findIdx(x+1, y) };

            for (int nIdx : neighbors) {
                float dx = curr->Xs[nIdx] - curr->Xs[idx];

                //if (dx > gridLength * 0.5f)  dx -= gridLength;
                //if (dx < -gridLength * 0.5f) dx += gridLength;

                float dy = curr->Ys[nIdx] - curr->Ys[idx];
                float dist = sqrt(dx*dx + dy*dy);

                if (dist > 0.0001f) { 
                    float forceMag = k * (dist - L);
                    totalFx += (dx / dist) * forceMag;
                    totalFy += (dy / dist) * forceMag;
                }
            }
            forceX[idx] = totalFx;
            forceY[idx] = totalFy;

            buffer->Xs[idx] = curr->Xs[idx] + (curr->Xs[idx]- prev->Xs[idx]) * 0.99f + inertia * forceX[idx] * dt * dt; 
            buffer->Ys[idx] = curr->Ys[idx] + (curr->Ys[idx]- prev->Ys[idx]) * 0.99f + inertia * forceY[idx] * dt * dt; 
        }
    }
}


void physicsStep(Grid* current, Grid* previous, Grid* buffer, float* forceX, float* forceY, int numThreads){
    std::vector<std::thread> workers;
    int rowsPerThread = gridHeight / numThreads;

    for (int i = 0; i < numThreads; ++i){
        int start = i * rowsPerThread;
        int end = (i == numThreads - 1) ? gridHeight : (start + rowsPerThread);

        workers.emplace_back(updateRows, start, end, current, previous, buffer, forceX, forceY);
    }

    for (auto& t : workers){
        t.join();
    }
}