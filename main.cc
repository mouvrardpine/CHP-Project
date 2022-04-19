#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string>
#include <mpi.h>
#include "fonction.h"
#include "matrix_RHS.h"
#include "solv_lin.h"






int main(int argc, char ** argv)
{
  int Nx(3), Ny(4),k(0),kmax(1000), i1,iN,me,np,pb(1);
  double dt(0.1),dx,dy,D, Lx(1), Ly(1), tmax(5.0),eps(pow(10,-10)), t(0.) ;
  MPI_Status status ;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&me); 
  MPI_Comm_size(MPI_COMM_WORLD,&np);


  //-----initialisation objets-----
  fonctions* fct;
  fct =new fonctions(pb,Lx,Ly);

  charge_* ch;
  ch = new charge_(Nx*Ny);

  matrix_RHS* mRHS;
  mRHS = new matrix_RHS(dt,Nx,Ny,fct,ch);

  solv_lin* sl;
  sl= new solv_lin(kmax,Nx , Ny ,eps ,dt,mRHS,ch);

  iN=ch->GetiN();
  i1=ch->Geti1();
  int size = iN - i1 + 1;

  std::vector<double> u(size,0),b(size,0),x0(size), uex(size,0), test(size,1), err(size,0);

  dx=1./(Nx+1);
  dy=1./(Ny+1);

  //b=RHS(dt,Nx,Ny,u,dt*k);
  //u=GC(x0,b,eps,kmax,Nx,Ny,dt);
 /*  for (int j(0); j<Ny +0; j++)
  {
    for (int i(0);i<Nx+0; i++ )
    {
      uex[j*Nx+i]=(i+1)*dx*(1-(i+1)*dx)*(j+1)*dy*(1-(j+1)*dy);
      //uex[j*Nx+i]=sin((i+1)*dx) + cos((j+1)*dy);

    }
  } */
  
    b=mRHS->RHS(u,t);
    
    u=sl->GC(u,b);
    
    while (t<tmax)
     {
        t += dt;
        cout<<"----------------------------------------------------------------- t= "<< t<< "me="<<me<<" --------------------------"<<endl;
        
        b=mRHS->RHS(u,t);
        

        //pb ici
        u=sl->GC(u,b);
        
        //cout << "t =" << t << endl;

        string filename("Output/sol" +to_string(int(t/dt))+ to_string(me) +".dat");
        fstream file_out;

        file_out.open(filename, std::ios_base::out);
        if (file_out) {
        
          for (int i(i1);i<=iN; i++ )
          {
            file_out<< (i%Nx+1)*dx << " " << (i/Nx+1)*dy << " " << u[i] <<endl;
            //cout << i/Nx+1 << endl;
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
       /*  if (me == 0)
        {
        cout << "err =" << sl->normL2_2D(sl->substract(u,uex), dx, dy) << endl;
        } */
      
    }



    MPI_Finalize(); 
    return 0;

}
