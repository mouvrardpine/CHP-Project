#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>
#include "fonction.h"
#include "matrix_RHS.h"
#include "solv_lin.h"


int main()
{
    int Nx(3), Ny(4),k(0),kmax(1000); 
    double dt(0.1),dx,dy,D, Lx(1), Ly(1), tmax(1),eps(0.0001); 
    std::vector<double> u(Nx*Ny,0),b(Nx*Ny,0),x0(Nx*Ny);



    b=RHS(dt,Nx,Ny,u,dt*k); 
    u=GC(x0,b,eps,kmax,Nx,Ny,dt); 

    while (dt*k<tmax)
    {
       b=RHS(dt,Nx,Ny,u,dt*k); 
       u=GC(x0,b,eps,kmax,Nx,Ny,dt); 
       k++;
       
    }
    

    for (int k(0); k<Nx*Ny; k++)
    {
        std::cout<< u[k]<<std :: endl;
    }
    return 0;

}
