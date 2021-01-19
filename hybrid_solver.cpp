#define N 64
#define cellSize 1.0f/64.0f
#define IX(i,j) ((i)+(N+2)*(j))
#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;}
#define FOR_EACH_CELL for ( i=1 ; i<=N ; i++ ) { for ( j=1 ; j<=N ; j++ ) {
#define END_FOR }}

#include <vector>
#include "hybrid_solver.h"
#include <cstdlib>
#include <bits/stdc++.h>

using namespace std;


void seedingParticles(vector<Particula> &particulas){
    int i,j;
    FOR_EACH_CELL
        for(int k=0; k<4; k++)
        {
            Particula part;
            part.x = i*(cellSize) + ((double) rand() / ((RAND_MAX)*N));
            part.y = j*(cellSize) + ((double) rand() / ((RAND_MAX)*N));
            part.u = 10;//(double) rand()*50;
            part.v = 10;//(double) rand()*50;
            part.col = part.x/(cellSize) - fmod(part.x, cellSize);
            part.lin = part.y/(cellSize) - fmod(part.y, cellSize);
            particulas.push_back(part);
        }
    END_FOR

}

void iniciaCelulasAux(vector<Celula> &celulas){
    for(int i=0; i<N+2; i++)
    for(int j=0; j<N+2; j++)
    {
        Celula cel;
        cel.k=0;
        cel.u=10;
        cel.v=10;
        celulas.push_back(cel);
    }
}

float H(float r){
    if(r<=1 && r>0)
        return 1 - r;
    else if(r>=-1 && r<=0)
        return 1 + r;

    return 0;
}

float K(float x,float y){
    return H(x/(cellSize))*H(y/(cellSize));
}

/*float newCellVelocityX(float gridCenterX, float gridCenterY, vector<Particula> particulas){
    float numerador;
    float denominador;

    for(int i=0; i<particulas.size(); i++)
    {
        numerador += particulas[i].u*K(particulas[i].x - gridCenterX,particulas[i].y - gridCenterY);
        denominador += K(particulas[i].x - gridCenterX,particulas[i].y - gridCenterY);
    }

    return numerador/denominador;
}

float newCellVelocityY(float gridCenterX, float gridCenterY, vector<Particula> particulas){
    float numerador = 0;
    float denominador = 0;

    for(int i=0; i<particulas.size(); i++)
    {
        numerador += particulas[i].v*K(particulas[i].x - gridCenterX,particulas[i].y - gridCenterY);
        denominador += K(particulas[i].x - gridCenterX,particulas[i].y - gridCenterY);
    }

    return numerador/denominador;
}*/

void particleToGridTransfer(float *u, float *v, vector<Celula> &celulas, vector<Particula> &particulas){


   int cellX, cellY;

    for(int k=0; k<N+2; k++){
     for(int l=0; l<N+2; l++){
        celulas[IX(k,l)].k = 0;
        celulas[IX(k,l)].u = 0;
        celulas[IX(k,l)].v = 0;
     }
    }
   for(int k=0; k<particulas.size(); k++){
        cellX = particulas[k].col*(cellSize) + (cellSize)/2;
        cellY = particulas[k].lin*(cellSize) + (cellSize)/2;

        celulas[IX(particulas[k].col,particulas[k].lin)].k += K(particulas[k].x - cellX,particulas[k].y - cellY);
        celulas[IX(particulas[k].col,particulas[k].lin)].u += particulas[k].u*K(particulas[k].x - cellX,particulas[k].y - cellY);
        celulas[IX(particulas[k].col,particulas[k].lin)].v += particulas[k].v*K(particulas[k].x - cellX,particulas[k].y - cellY);
    }

    for(int i=0; i<N+2; i++){
     for(int j=0; j<N+2; j++){

      if(celulas[IX(i,j)].k !=0)
      {
        u[IX(i,j)] = celulas[IX(i,j)].u/celulas[IX(i,j)].k;
        v[IX(i,j)] = celulas[IX(i,j)].v/celulas[IX(i,j)].k;
      }else{
        u[IX(i,j)] = celulas[IX(i,j)].u;
        v[IX(i,j)] = celulas[IX(i,j)].v;
      }

     }
     }

}

void gridToParticleTransfer(float *u, float *v, vector<Particula> &particulas){

    float cellX, cellY;
    float distanciasX, distanciasY;

    for(int i=0; i<particulas.size(); i++){

        distanciasX = 0;
        distanciasY = 0;
        cellX = particulas[i].col*(cellSize) + cellSize/2;
        cellY = particulas[i].lin*(cellSize) + cellSize/2;

        particulas[i].u = 0.0;
        particulas[i].v = 0.0;
        particulas[i].u += u[IX(particulas[i].col,particulas[i].lin)]*K(particulas[i].x - cellX,particulas[i].y - cellY);
        particulas[i].v += v[IX(particulas[i].col,particulas[i].lin)]*K(particulas[i].x - cellX,particulas[i].y - cellY);
        distanciasX += abs(particulas[i].x - cellX);
        distanciasY += abs(particulas[i].y - cellY);

        if(particulas[i].col>0)
        {
            particulas[i].u += u[IX(particulas[i].col-1,particulas[i].lin)]*K(particulas[i].x - cellX - (cellSize), particulas[i].y - cellY);
            particulas[i].v += v[IX(particulas[i].col-1,particulas[i].lin)]*K(particulas[i].x - cellX - (cellSize), particulas[i].y - cellY);
            distanciasX += abs(particulas[i].x - cellX - (cellSize));
            distanciasY += abs(particulas[i].y - cellY);
        }

        if(particulas[i].col<N+2)
        {
            particulas[i].u += u[IX(particulas[i].col+1,particulas[i].lin)]*K(particulas[i].x - cellX + (cellSize), particulas[i].y - cellY);
            particulas[i].v += v[IX(particulas[i].col+1,particulas[i].lin)]*K(particulas[i].x - cellX + (cellSize), particulas[i].y - cellY);
            distanciasX += abs(particulas[i].x - cellX + (cellSize));
            distanciasY += abs(particulas[i].y - cellY);
        }

        if(particulas[i].lin>0)
        {
            particulas[i].u += u[IX(particulas[i].col-1,particulas[i].lin)]*K(particulas[i].x - cellX, particulas[i].y - cellY - (cellSize));
            particulas[i].v += v[IX(particulas[i].col-1,particulas[i].lin)]*K(particulas[i].x - cellX, particulas[i].y - cellY - (cellSize));
            distanciasX += abs(particulas[i].x - cellX);
            distanciasY += abs(particulas[i].y - cellY - (cellSize));
        }

        if(particulas[i].lin<N+2)
        {
            particulas[i].u += u[IX(particulas[i].col+1,particulas[i].lin)]*K(particulas[i].x - cellX, particulas[i].y - cellY - (cellSize));
            particulas[i].v += v[IX(particulas[i].col+1,particulas[i].lin)]*K(particulas[i].x - cellX, particulas[i].y - cellY - (cellSize));
            distanciasX += abs(particulas[i].x - cellX);
            distanciasY += abs(particulas[i].y - cellY + (cellSize));
        }

        if(distanciasX == 0)
        {
            cout<<"teste"<<endl;
            distanciasX = 1;
        }
        if(distanciasY == 0)
        {
            cout<<"teste"<<endl;
            distanciasY = 1;
        }

        if(particulas[i].u>0 || particulas[i].v>0)
        {
            cout<<"teste"<<endl;
        }

        particulas[i].u = particulas[i].u/distanciasX;
        particulas[i].v = particulas[i].v/distanciasY;
    }
}

float rungeKuta(float velocidade, float dt){

    float k1, k2, k3;

    k1 = velocidade;
    k2 = velocidade + 0.5f*dt*k1;
    k3 = velocidade + 0.75f*dt*k2;

    return (2.0f/9.0f)*k1 + (3.0f/9.0f)*k2 + (4.0f/9.0f)*k3;
}

void advectParticles(float dt, vector<Particula> &particulas){


    for(int i=0; i<particulas.size(); i++){
        particulas[i].x += rungeKuta(particulas[i].u,dt);
        particulas[i].y += rungeKuta(particulas[i].v,dt);
        particulas[i].col = particulas[i].x/(cellSize) - fmod(particulas[i].x, cellSize);
        particulas[i].lin = particulas[i].y/(cellSize) - fmod(particulas[i].y, cellSize);


        if(particulas[i].col > N)
        {
            particulas[i].x = cellSize*(N-1) + cellSize/2;
            particulas[i].u*=-1;
            particulas[i].col = particulas[i].x/(cellSize) - fmod(particulas[i].x, cellSize);
        }else if(particulas[i].col < 0)
        {
            particulas[i].x = cellSize/2;
            particulas[i].u*=-1;
            particulas[i].col = particulas[i].x/(cellSize) - fmod(particulas[i].x, cellSize);
        }

        if(particulas[i].lin > N)
        {
            particulas[i].y = cellSize*(N-1) + cellSize/2;
            particulas[i].lin = particulas[i].y/(cellSize) - fmod(particulas[i].y, cellSize);
        }else if(particulas[i].lin < 0)
        {
            particulas[i].y = cellSize + 2*cellSize/2;
            particulas[i].lin = particulas[i].y/(cellSize) - fmod(particulas[i].y, cellSize);
        }
    }
}

