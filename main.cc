#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>
#include "fonction.h"
#include "matrix_RHS.h"
#include "solv_lin.h"


int main()
{
    int Nx(3), Ny(4) ; 
    double dt(0.1),dx,dy,D, Lx(1), Ly(1) ; 
    std::vector<double> v(Nx*Ny,1), u(Nx*Ny);
    u=matvec(dt, Nx, Ny, v);
    for (int k(0); k<Nx*Ny; k++)
    {
    std::cout<< u[k]<<std :: endl;
    }



    return 0;

}