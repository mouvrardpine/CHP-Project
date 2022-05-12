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
  int Nx(3), Ny(4),k(0),kmax(500), i1,iN,me,Np,pb(1);
  double dt(0.1),dx,dy,D, Lx(1), Ly(1), tmax(5.0),eps(pow(10,-10)), t(0.),err(0);
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  MPI_Comm_size(MPI_COMM_WORLD,&Np);
  MPI_Status status ;

   //-----initialisation objets-----

   //solv_lin* sl;
   //sl= new solv_lin(kmax,Nx , Ny ,eps ,dt,mRHS,ch);
   int n = Nx*Ny;
   //cout << n << endl;
   charge(&i1,&iN, me,n,Np);
 
   int size = iN - i1 + 1;
  
   std::vector<double> u(size,0),b(size,0),x0(size,1), uex_vec(size,0), test(size,1), utest(n),u1(size,0);

   dx=1./(Nx+1);
   dy=1./(Ny+1);
   //x0=mRHS->matvec(x0, i1, iN, me, Nx,Ny,Np,dt);

   // //cout<< "me=   "<<me<<endl;
   // MPI_Status Status;
   // int i1b , iNb;
   // for (int k(0);k<x0.size();k++)
   // {
   //
   //   utest[k+i1]=x0[k];
   //
   // }
   // if (me==0)
   // {
   //   for (int l(1);l<np;l++)
   //   {
   //     ch->charge(&i1b,&iNb,l,Nx*Ny,np);
   //     MPI_Recv(&utest[i1b] , iN-i1+1, MPI_DOUBLE, l, 10, MPI_COMM_WORLD, &Status);
   //   }
   //   for (int f(0);f<Nx*Ny;f++)
   //   {
   //     cout<<utest[f]<<endl;
   //   }
   // }
   // else
   // {
   //   MPI_Send(&x0[0], iN-i1+1, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
   // }

   for (int i = 0; i < size; i++) {
     if (i%2 ==0)
     {
       x0[i] = 2.;
     }
   }
 
   x0=matvec(x0, i1, iN, me, Nx,Ny,Np,dt);
   //cout<< "me=   "<<me<<endl;
    MPI_Status Status;
    int i1b , iNb;

    for (int k(0);k<x0.size();k++)
    {
      utest[k+i1]=x0[k];
      //cout << "k = " << k + i1 << endl;
    }
  
    if (me==0)
    {
      for (int l(1);l<Np;l++)
      {
        charge(&i1b,&iNb,l,n,Np);
        MPI_Recv(&utest[i1b] , iN-i1+1, MPI_DOUBLE, l, 10, MPI_COMM_WORLD, &Status);
      }
      for (int f(0);f<Nx*Ny;f++)
      {
        cout<<utest[f]<<endl;
      }
    }
    else
    {
      MPI_Send(&x0[0], iN-i1+1, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
    }

      // MATHIS BOUH BOUH BOUH
    /* for (int j(0); j<Ny +0; j++)
    {
      for (int i(0);i<Nx+0; i++ )
      {
        uex[j*Nx+i]=(i+1)*dx*(1-(i+1)*dx)*(j+1)*dy*(1-(j+1)*dy);
        //uex[j*Nx+i]=sin((i+1)*dx) + cos((j+1)*dy);

      }
    } */

   b=RHS(u,t,Nx,Ny,i1, iN,Lx,Ly,dt,pb);
// cout << "me = " << me << " size = " << b.size() << endl;

     u1=GC(u,b, i1, iN, kmax, me, Nx, Ny, Np, dt, eps);
     u=u1;


  while (t<tmax)
     {
        t += dt;
        cout<<"-------------------- t= "<< t<< "me="<<me<<" --------------------------"<<endl;

        b=RHS(u,t,Nx,Ny,i1, iN,Lx,Ly,dt,pb);

        cout<<"RHS_"<< me<<endl;

        u=GC(u,b, i1, iN, kmax, me, Nx, Ny, Np, dt, eps);

        //cout << "t =" << t << endl;

        string filename("Output/sol" +to_string(me)+ to_string(int(t/dt)) +".dat");
        fstream file_out;

        file_out.open(filename, std::ios_base::out);
        if (file_out) {

          for (int i(i1);i<=iN; i++ )
          {
            file_out<< (i%Nx+1)*dx << " " << (i/Nx+1)*dy << " " << u[i-i1] <<endl;
            //cout << i/Nx+1 << endl;
          }


        // for (int k(0); k<size; k++)return 0;
        //  {
        //    file_out<< u[k]<<endl;
        // }
        //cout<<"saucisse3 "<< me<<endl;
        //file_out << "Some random text to write." << endl;
        cout << "Done Writing!" << me<<endl;
    } else {
          cout << "failed to open " << filename << '\n';

    }
    for (int i(i1); i<=iN ; i++)
    {
        uex_vec[i-i1]=uex((i%Nx+1)*dx, (i/Nx+1)*dy, pb);
        //uex_vec[i-i1]=0;
    } 
    
    err = normL2_2D(substract(u,uex_vec), dx, dy);
    //err = norm(substract(u,uex_vec));
    cout<<"me = " << me << " err=  "<<err<<endl;

    }

    MPI_Finalize();
    return 0;

}
