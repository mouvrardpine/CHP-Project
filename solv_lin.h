#ifndef _SOLV_LIN_H

#include "matrix_RHS.h"

using namespace std;
double ps(vector <double> x1, vector <double> x2);
vector<double> substract(vector <double> x1, vector <double> x2);
vector<double> add(vector <double> x1, vector <double> x2);
double norm(vector <double> x);
vector <double> mult(double a , vector<double> x);
vector<double> GC(vector <double> x0 , vector <double> b , double eps , int kmax,int Nx, int Ny,double dt);
double normL2_2D(vector<double> x, double dx, double dy);

#define _SOLV_LIN_H
#endif
