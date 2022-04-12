#ifndef _MATRIX_RHS_CPP

#include "matrix_RHS.h"
#include "fonction.h"
#include <cmath>
#include <iostream>

using namespace std;

std::vector<double> matvec( double dt,  int Nx,  int Ny,  std::vector<double> x)
{

	double dx = 1./(Nx+1) , dy = 1./(Ny+1);
	double alpha = 2*dt*(1/(pow(dx,2)) +1 / pow(dy,2)) + 1 ;
	double beta_x = - dt/pow(dx,2) , beta_y = - dt/pow(dy,2) ;
	//cout << "alpha= "<< alpha << "	beta_x="<<beta_x<< "	beta_y="<<beta_y<<endl;
	std::vector<double> Ax(Nx*Ny,0);

	for (int i = 0; i < Nx*Ny; ++i)
	{
		Ax[i] = x[i]*alpha;

		if (i < Nx*Ny - Nx)
		{
			Ax[i] += beta_y*x[i+Nx];
		}
		if (i > Nx-1)
		{
			Ax[i] += beta_y*x[i-Nx];
		}

		if ((i+1)%Nx != 0)
		{
			Ax[i] += beta_x*x[i+1];
		}
		if (i%Nx !=0 && i != 0)
			{
				Ax[i] += beta_x*x[i-1];
			}

	}

	return Ax;
}

std::vector<double> RHS( double dt, int Nx, int Ny, std::vector<double> u, double t)
{
	double dx = 1./(Nx+1) , dy = 1./(Ny+1);
	double alpha = -2*dt*(1/(pow(dx,2)) +1 / pow(dy,2)) + 1 ;
	double beta_x = - dt/pow(dx,2) , beta_y = - dt/pow(dy,2) ;
	std::vector<double> F(Nx*Ny,0.);


	//for (int i = 0; i < Nx*Ny; i++)
	//{

		//F[i] = 0;
	//}

	for (int i = 0; i < Nx*Ny ; ++i)
	{
		double x = (i%Nx+1)*dx, y = (i/Nx+1)*dy; // +1 a vérifier

		F[i] = dt*f1(x, y, t + dt, 1, 1. ,1.) + u[i];

		if (x == dx || x == 1.-dx)
		{
			F[i] += beta_x*h1(x,y,t+dt, 1);    // attention à refaire LES DEUX !!!!
		}

		if (y == dy || y == 1.-dy)
		{
			F[i] += beta_y*g1(x,y,t+dt, 1);
		}
		//cout<< "x="<< x << "  y=  "<< y<< "   F="<< F[i]<< endl;
	}

	return F;
}

#define _MATRIX_RHS_CPP
#endif
