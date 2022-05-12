#ifndef _MATRIX_RHS_H

#include <string>
#include <vector>
#include <cmath>
#include "fonction.h"

    /* double _dt ;
    int _Nx , _Ny, _i1, _iN, _me , _Np,_n ;
    fonctions* _fct;
    charge_* _ch; */
    std::vector<double> matvec(std::vector<double> x, int i1, int iN, int me,int Nx,int Ny,int Np, double dt);
    std::vector<double> RHS( std::vector<double> u, double t,  int Nx, int Ny, double Lx, double Ly,double dt );

#define _MATRIX_RHS_H
#endif
