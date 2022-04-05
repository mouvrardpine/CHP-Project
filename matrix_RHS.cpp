#ifndef _MATRIX_RHS_CPP

#include "matrix_RHS.h"
#include <cmath>
#include <iostream>

using namespace std;

std::vector<double> matvec(const double dt, const int Nx, const int Ny, const std::vector<double> x) const
{

	double dx = 1./(Nx+1) ; dy = 1./(Ny+1);
	double alpha = 2*dt/(pow(dx,2) + pow(dy,2)) - 1 ;
	double beta_x = - dt/pow(dx,2) ;  beta_y = - dt/pow(dy,2) ; 

	std::vector<double> Ax;
	Ax.resize(Nx*Ny);

	for (int i = 0; i < Nx*Ny; ++i)
	{
		Ax(i) = x(i*(Nx*Ny+1))*alpha;

		if (i < Nx*Ny - 3)
		{
			Ax(i) += beta_y*x(i(Nx*Ny + 1) + Nx);
		}
		if (i > 2)
		{
			Ax(i) += beta_y*x((i+Nx)*(Nx*Ny) + i);
		}

		if ((i+1)%Nx /= 0)
		{
			Ax(i) += beta_x*x(i*(Nx*Ny + 1) + (Nx-2));
			Ax(i) += beta_x*x((i+Nx-2)*Nx*Ny + i);
		}
	}

	return Ax;
}





#define _MATRIX_RHS_CPP
#endif