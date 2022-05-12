#ifndef _MATRIX_RHS_CPP
#include "matrix_RHS.h"
#include "fonction.h"
#include <cmath>
#include <iostream>

using namespace std;


std::vector<double>   matvec(std::vector<double> x, int i1, int iN, int me,int Nx,int Ny,int Np, double dt)
{

	int he, size(iN-i1+1);
	double dx = 1./(Nx+1) , dy = 1./(Ny+1), msg;
	double alpha = 2*dt*(1/(pow(dx,2)) +1 / pow(dy,2)) + 1 ;
	double beta_x = - dt/pow(dx,2) , beta_y = - dt/pow(dy,2) ;
	std::vector<double> Ax(size,0), grand_x(size + 2*Nx);
	//cout << "alpha= "<< alpha << "	beta_x="<<beta_x<< "	beta_y="<<beta_y<<endl;

	// Envoi des Nx premiers éléments de me à me-1 et réception des Nx derniers éléments de me envoyés par me+1
	for (int i = 0; i < Nx; i++) {
		MPI_Status Status;
		int tag = i +1;
		//cout<<" me "<<me << "  j'envoie à "<<he <<" saucisse1"<<endl;
		if (me !=0) MPI_Send(&x[i], 1, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD);
		//cout<<" me "<<me << "  j'attends à "<<he <<" saucisse1"<<endl;
		if (me != Np -1) MPI_Recv(&grand_x[grand_x.size()-1-(Nx-i-1)], 1, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD, &Status);
	}

	// Envoi des Nx derniers éléments de me à me+1 et réception des Nx derniers éléments de me envoyés par me-1
	for (int i = 0; i < Nx; i++) {
		MPI_Status Status;
		int tag = 10*i +1;
		//cout<<" me "<<me << "  j'envoie à "<<he <<" saucisse1"<<endl;
		if (me != Np-1) MPI_Send(&x[size- 1-(Nx-i-1)], 1, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD);
		//cout<<" me "<<me << "  j'attends à "<<he <<" saucisse1"<<endl;
		if (me != 0) MPI_Recv(&grand_x[i], 1, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD, &Status);
	}

	//remplissage des autres éléments (éléments de me)
	for (int i = Nx; i < Nx + size ; i++)
	{
		grand_x[i] = x[i - Nx];
	}

	 // for (int i = 0; i < grand_x.size(); i++)
	 // {
	 // 	cout << "test" << grand_x[i] << endl ;
	 // }

	int iloc;
	for (int iglob = i1; iglob <= iN; ++iglob)
	{
	 	iloc = iglob - i1;
		Ax[iloc] = x[iloc]*alpha;
		//cout<< "i= "<< i << " me "<< me<<endl;
    if (iglob < Nx*Ny - Nx)
      {
				//cout << "ok" << endl;
        Ax[iloc] += beta_y*grand_x[iloc+2*Nx];
      }
    if (iglob > Nx-1)
      {
        Ax[iloc] += beta_y*grand_x[iloc];
      }

    if ((iglob+1)%Nx != 0)
      {
        Ax[iloc] += beta_x*grand_x[iloc+1+Nx];
      }

		if (iglob%Nx !=0 && iglob != 0)
      {
        Ax[iloc] += beta_x*grand_x[iloc-1+Nx];
      }

	}

	return Ax;
}

std::vector<double> RHS( std::vector<double> u, double t, int Nx, int Ny, double Lx, double Ly, double dt )
{
	double dx = 1./(Nx+1) , dy = 1./(Ny+1);
	double alpha = -2*dt*(1/(pow(dx,2)) +1 / pow(dy,2)) + 1 ;
	double beta_x = dt/pow(dx,2) , beta_y = dt/pow(dy,2) ;
	std::vector<double> F(Nx*Ny,0.);
	int pb(2);

	for (int i = 0; i < Nx*Ny ; ++i)
	{
		double x = (i%Nx+1)*dx, y = (i/Nx+1)*dy;

		F[i] = dt*f1(x, y, t + dt, pb, Lx, Ly) + u[i];

		if ((i+1)%Nx == 1)
		{
			F[i] += beta_x* h1(0,y,pb);
		}

		if ((i+1)%Nx == 0)
		{
			F[i] += beta_x* h1(1.,y,pb);
		}

		if (i/Nx == 0)
		{
			F[i] += beta_y* g1(x,0,pb);
		}

		if (i/Nx + 1 == Ny)
		{
			F[i] += beta_y* g1(x,1.,pb);
		}

		/*if ( (x == dx) || (x == 1.- dx))
		{
			F[i] += beta_x*h1(x,y,t+dt, 1);
		}

		if ( (y == dy) || (y == 1.-dy) )
		{
			F[i] += beta_y*g1(x,y,t+dt, 1);
		}*/

	}

	return F;
}

#define _MATRIX_RHS_CPP
#endif
