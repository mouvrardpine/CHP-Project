#ifndef _SOLV_LIN_H

#include "matrix_RHS.h"

using namespace std;


class solv_lin 
{
private :
    int _kmax,_Nx, _Ny,_i1,_iN,_me, _n, _Np;
    double _eps,_dt;
    matrix_RHS* _mRHS;
    charge_* _ch;

public :
    solv_lin(int kmax, int Nx , int Ny , double eps , double dt, matrix_RHS* mrhs,charge_* ch);
    double ps(vector <double> x1, vector <double> x2);
    vector<double> substract(vector <double> x1, vector <double> x2);
    vector<double> add(vector <double> x1, vector <double> x2);
    double norm(vector <double> x);
    vector <double> mult(double a , vector<double> x);
    vector<double> GC(vector <double> x0 , vector <double> b );
    double normL2_2D(vector<double> x, double dx, double dy);


};
#define _SOLV_LIN_H
#endif
