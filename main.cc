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
    double dt(1.0),dx,dy,D, Lx(1), Ly(1), tmax(50.0),eps(0.0001), t(0.);
    std::vector<double> u(Nx*Ny,0),b(Nx*Ny,0),x0(Nx*Ny);



  //b=RHS(dt,Nx,Ny,u,dt*k);
  //u=GC(x0,b,eps,kmax,Nx,Ny,dt);

    while (t<tmax)
     {
        b=RHS(dt,Nx,Ny,u,t);
        u=GC(x0,b,eps,kmax,Nx,Ny,dt);
        //k++;

        cout << "t =" << t << endl;
              for (int k(0); k<Nx*Ny; k++)
             {
                 std::cout<< u[k]<<std :: endl;
             }

      t += dt;
     }


    return 0;

}
