#ifndef HYBRID_SOLVER_H_INCLUDED
#define HYBRID_SOLVER_H_INCLUDED

#include <vector>
using namespace std;

struct Particula{
    double x;
    double y;

    double u;
    double v;

    int col;
    int lin;
};
struct Celula{
    float k;
    float u;
    float v;
};


void seedingParticles(vector<Particula> &particulas);
void iniciaCelulasAux(vector<Celula> &celulas);
float H(float r);
float K(float x,float y);
float newCellVelocityX(float gridCenterX, float gridCenterY, vector<Particula> &particulas);
float newCellVelocityY(float gridCenterX, float gridCenterY, vector<Particula> &particulas);
void particleToGridTransfer(float *u, float *v, vector<Celula> &celulas, vector<Particula> &particulas);
void gridToParticleTransfer(float *u, float *v, vector<Particula> &particulas);
float rungeKuta(float velocidade, float dt);
void advectParticles(float dt, vector<Particula> &particulas);

#endif // HYBRID_SOLVER_H_INCLUDED
