#ifndef _MATRIX_RHS_H

#include <string>
#include <vector>
#include <cmath>
#include "fonction.h"
class matrix_RHS
{
private :
    double _dt ;
    int _Nx , _Ny, _i1, _iN, _me , _Np,_n ;
    fonctions* _fct;
    charge_* _ch;
public :
    matrix_RHS(double dt, int Nx, int Ny, fonctions* fct, charge_* ch);
    std::vector<double> matvec(std::vector<double> x);
    std::vector<double> RHS(std::vector<double> u, double t);
};
#define _MATRIX_RHS_H
#endif
