#ifndef _SOLV_LIN_H

#include "matrix_RHS.h"

using namespace std;

double ps(vector <double> x1, vector <double> x2);
vector<double> substract(vector <double> x1, vector <double> x2);
vector<double> add(vector <double> x1, vector <double> x2);
double norm(vector <double> x);
vector <double> mult(double a , vector<double> x);
vector<double> GC(std::vector <double> x0 , std::vector <double> b, int i1, int iN, int kmax,int me, int Nx, int Ny, int Np, double dt, double eps );
double normL2_2D(vector<double> x, double dx, double dy);

#define _SOLV_LIN_H
#endif
