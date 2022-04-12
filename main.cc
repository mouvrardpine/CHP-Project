#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string>
#include "fonction.h"
#include "matrix_RHS.h"
#include "solv_lin.h"


int main()
{
    int Nx(12), Ny(16),k(0),kmax(1000);
    double dt(0.1),dx,dy,D, Lx(1), Ly(1), tmax(5.0),eps(pow(10,-10)), t(0.);
    std::vector<double> u(Nx*Ny,0),b(Nx*Ny,0),x0(Nx*Ny), uex(Nx*Ny,0), test(Nx*Ny,1), err(Nx*Ny,0);

  dx=1./(Nx+1);
  dy=1./(Ny+1);

  //b=RHS(dt,Nx,Ny,u,dt*k);
  //u=GC(x0,b,eps,kmax,Nx,Ny,dt);
  for (int j(0); j<Ny +0; j++)
  {
    for (int i(0);i<Nx+0; i++ )
    {
      //uex[j*Nx+i]=(i+1)*dx*(1-(i+1)*dx)*(j+1)*dy*(1-(j+1)*dy);
      uex[j*Nx+i]=sin((i+1)*dx) + cos((j+1)*dy);

    }
  }
    b=RHS(dt,Nx,Ny,u,t);
    u=GC(u,b,eps,kmax,Nx,Ny,dt);

    while (t<tmax)
     {
        t += dt;
        b=RHS(dt,Nx,Ny,u,t);
        u=GC(u,b,eps,kmax,Nx,Ny,dt);
        cout << "t =" << t << endl;
        /*for (int k(0); k<Nx*Ny; k++)
         {
           std::cout<< u[k]<<std :: endl;
        }
*/        // err= substract(u,uex);

        string filename("Output/sol_" +to_string(int(t/dt))+".dat");
        fstream file_out;

        file_out.open(filename, std::ios_base::out);
        if (file_out) {
        
        for (int j(0); j<Ny +0; j++)
        {
          for (int i(0);i<Nx+0; i++ )
          {
            file_out<< (i+1)*dx << " " << (j+1)*dy << " " << u[j*Nx+i] <<endl;

          }
        }



        for (int k(0); k<Nx*Ny; k++)
         {
           file_out<< u[k]<<endl;
        }
        //file_out << "Some random text to write." << endl;
        cout << "Done Writing!" << endl;
    } else {
          cout << "failed to open " << filename << '\n';
      
    }

        //cout << "err =" << normL2_2D(substract(u,uex), dx, dy) << endl;
      }



    return 0;

}
