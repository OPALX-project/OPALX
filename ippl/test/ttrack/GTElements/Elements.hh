#ifndef ELEMENTS_H
#define ELEMENTS_H
#include "Physics.hh"

void dipole(double t, double *y, double *yp, double *coeff)
{
    /*
              coeff[0] == mass   / kg
              coeff[1] == k0     / T
    */
    const double m = coeff[0];
    const double k0 = coeff[1];
    const double psq=pow(y[3],2)+pow(y[4],2)+pow(y[5],2);
    const double gamma=sqrt(1.0+psq/pow(m*Physics::c,2));
    const double vx = y[3]/(m*gamma);
    const double vy = y[4]/(m*gamma);
    const double vz = y[5]/(m*gamma);
    yp[3] = -Physics::q_e*vz*k0;
    yp[4] = 0.0; 
    yp[5] =  Physics::q_e*vx*k0;
    yp[0] = vx;
    yp[1] = vy;
    yp[2] = vz;
}

void quadrupole(double t, double *y, double *yp, double *coeff)
{
    /*
      coeff[0] == mass   / kg
      coeff[2] == k1     / T/m
    */
    const double m = coeff[0];
    const double k1 = coeff[2];
    const double psq=pow(y[3],2)+pow(y[4],2)+pow(y[5],2);
    const double gamma=sqrt(1.0+psq/pow(m*Physics::c,2));
    const double vx = y[3]/(m*gamma);
    const double vy = y[4]/(m*gamma);
    const double vz = y[5]/(m*gamma);
    
    yp[3] = -Physics::q_e*vz*k1*y[0];
    yp[4] =  Physics::q_e*vz*k1*y[1];
    yp[5] =  Physics::q_e*k1*(vx*y[0]-vy*y[1]);
    yp[0] = vx;
    yp[1] = vy;
    yp[2] = vz;
}
#endif
