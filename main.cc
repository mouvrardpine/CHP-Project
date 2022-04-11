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
    std::vector<double> u(Nx*Ny,0),b(Nx*Ny,0),x0(Nx*Ny), uex(Nx*Ny,0), test(Nx*Ny,1), err(Nx*Ny,0);
    /* u=matvec(dt,Nx, Ny,test);
    for (int k(0); k<Nx*Ny; k++)
        {
          std::cout<< u[k]<<std :: endl;
        } */

  dx=1/(Nx+1);
  dy=1/(Ny+1);

  //b=RHS(dt,Nx,Ny,u,dt*k);
  //u=GC(x0,b,eps,kmax,Nx,Ny,dt);
  for (int j(0); j<Ny; j++)
  {
    for (int i(0);i<Nx; i++ )
    {
      uex[j*Nx+i]=i*dx*(1-i*dx)*j*dy*(1-j*dy); 
      
    }
  }
    b=RHS(dt,Nx,Ny,u,t);
    u=GC(x0,b,eps,kmax,Nx,Ny,dt);
    for (int k(0); k<Nx*Ny; k++)
        {
          std::cout<< u[k]<<std :: endl;
        }
    //while (t<tmax)
     {
        t += dt;
        b=RHS(dt,Nx,Ny,u,t);
        u=GC(u,b,eps,kmax,Nx,Ny,dt);
        cout << "t =" << t << endl;
        for (int k(0); k<Nx*Ny; k++)
        {
          std::cout<< u[k]<<std :: endl;
        }
        err= substract(u,uex);
        cout << "erruer"<< endl;
        for (int k(0); k<Nx*Ny; k++)
        {
          std::cout<< err[k]<<std :: endl;
        }
        cout << "err =" << norm(substract(u,uex)) << endl;
          
        

            
      }
     


    return 0;

}
