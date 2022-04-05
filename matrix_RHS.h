#ifndef _MATRIX_RHS_H

#include <string>
#include <vector>
#include <cmath>

std::vector<double> matvec(double dt,int Nx,int Ny,std::vector<double> x);
std::vector<double> RHS( double dt, int Nx, int Ny, std::vector<double> u, double t);

#define _MATRIX_RHS_H
#endif
