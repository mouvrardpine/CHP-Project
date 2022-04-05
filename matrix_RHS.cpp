#ifndef _MATRIX_RHS_CPP

#include "matrix_RHS.h"
#include <cmath>
#include <iostream>

using namespace std;

std::vector<double> matvec( double dt,  int Nx,  int Ny,  std::vector<double> x)
{

	double dx = 1./(Nx+1) , dy = 1./(Ny+1);
	double alpha = 2*dt/(pow(dx,2) + pow(dy,2)) - 1 ;
	double beta_x = - dt/pow(dx,2) , beta_y = - dt/pow(dy,2) ; 
	cout << "alpha= "<< alpha << "	beta_x="<<beta_x<< "	beta_y="<<beta_y<<endl;
	std::vector<double> Ax(Nx*Ny,0);

	for (int i = 0; i < Nx*Ny; ++i)
	{
		Ax[i] = x[i]*alpha;

		if (i < Nx*Ny - Nx)
		{
			Ax[i] += beta_y*x[i+Nx];
		}
		if (i > Nx)
		{
			Ax[i] += beta_y*x[i-Nx];
		}

		if ((i+1)%Nx != 0)
		{
			if (i < Nx*Ny)
			{
				Ax[i] += beta_x*x[i+1];
			}
			if (i > 0)
			{
				Ax[i] += beta_x*x[i-1];
			}
			
		}
	}

	return Ax;
}





#define _MATRIX_RHS_CPP
#endif