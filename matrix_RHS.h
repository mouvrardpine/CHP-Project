#ifndef _MATRIX_RHS_H

#include <string>
#include <vector>
#include <cmath>
#include "fonction.h"

std::vector<double> matvec(std::vector<double> x, int i1, int iN, int me,int Nx,int Ny,int Np, double dt);
std::vector<double> RHS( std::vector<double> u, double t,  int Nx, int Ny, int i1, int iN, double Lx, double Ly,double dt, int pb );

#define _MATRIX_RHS_H
#endif
