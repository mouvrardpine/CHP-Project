#pragma once


#include <stdio.h>
#include <vector>
#include <math.h>
#include <mpi.h>

class fonctions
{
private:
    int _pb;
    double _Lx, _Ly;
public:
    fonctions(int pb, double Lx, double Ly);
    ~fonctions();

    double f1(double x, double y, double t); 
    double g1(double x, double y, double t); 
    double h1(double x, double y, double t);
    double uex(double x, double y);

};


class charge_
{
private :
    int _i1 , _iN , _me, _n, _Np ;

public:
    charge_(int n);
    void charge(int *i1,int *iN, int me, int n, int Np);
    int whichMe (int i);
    
    int Geti1() {return _i1;};
    int GetiN() {return _iN;};
    int Getn() {return _n;};
};

