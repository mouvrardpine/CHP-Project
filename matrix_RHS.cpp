#ifndef _MATRIX_RHS_CPP
#include "matrix_RHS.h"
#include "fonction.h"
#include <cmath>
#include <iostream>

using namespace std;

// Produit matrice vecteur pour le problème posé
std::vector<double>   matvec(std::vector<double> x, int i1, int iN, int me,int Nx,int Ny,int Np, double dt)
{
	//Taille du vecteur local
	int size(iN-i1+1);

	//Pas d'espace du maillage
	double dx = 1./(Nx+1) , dy = 1./(Ny+1);

	//Stockage des coefficients de la matrice A
	double alpha = 2*dt*(1/(pow(dx,2)) +1 / pow(dy,2)) + 1 ;
	double beta_x = - dt/pow(dx,2) , beta_y = - dt/pow(dy,2) ;

	//Initialisation du vecteur résultat (Ax) et de grand_x, vecteur contenant toutes les informations nécessaires au caclul de Ax
	std::vector<double> Ax(size,0.), grand_x(size + 2*Nx,0.);

	// Envoi des Nx premiers éléments de me à me-1 et réception des Nx derniers éléments de me envoyés par me+1
	for (int i = 0; i < Nx; i++) {
		MPI_Status Status;
		int tag = i +1;
		if (me !=0) MPI_Send(&x[i], 1, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD);
		if (me != Np -1) MPI_Recv(&grand_x[grand_x.size()-1-(Nx-1-i)], 1, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD, &Status);
	}

	// Envoi des Nx derniers éléments de me à me+1 et réception des Nx derniers éléments de me envoyés par me-1
	for (int i = 0; i < Nx; i++) {
		MPI_Status Status;
		int tag = 10*i +1;
		if (me != Np-1) MPI_Send(&x[size- 1-(Nx-i-1)], 1, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD);
		if (me != 0) MPI_Recv(&grand_x[i], 1, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD, &Status);
	}

	//remplissage des autres éléments (éléments de me)
	for (int i = Nx; i < Nx + size ; i++)
	{
		grand_x[i] = x[i - Nx];
	}

	//indice local
	int iloc;

	//Calcul de Ax
	for (int iglob = i1; iglob <= iN; ++iglob)
	{
	 	iloc = iglob - i1;
		Ax[iloc] = x[iloc]*alpha;
    if (iglob < Nx*Ny - Nx)
    {
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

// Calcul du second membre
std::vector<double> RHS( std::vector<double> u, double t, int Nx, int Ny, int i1, int iN, double Lx, double Ly, double dt, int pb )
{
	//pas d'espace
	double dx = 1./(Nx+1) , dy = 1./(Ny+1);

	//coefficients de A
	double alpha = -2*dt*(1/(pow(dx,2)) +1 / pow(dy,2)) + 1 ;
	double beta_x = dt/pow(dx,2) , beta_y = dt/pow(dy,2) ;

	//Initialisation du second membre
	std::vector<double> F(u.size(),0.);

	//Calcul du second membre
	for (int i = i1; i <= iN ; ++i)
	{
		double x = (i%Nx+1)*dx, y = (i/Nx+1)*dy;

		F[i-i1] = dt*f1(x, y, t + dt, pb, Lx, Ly) + u[i-i1];

		if ((i+1)%Nx == 1)
		{
			F[i-i1] += beta_x* h1(0,y,pb);
		}

		if ((i+1)%Nx == 0)
		{
			F[i-i1] += beta_x* h1(1.,y,pb);
		}

		if (i/Nx == 0)
		{
			F[i-i1] += beta_y*g1(x,0,pb);
		}

		if (i/Nx + 1 == Ny)
		{
			F[i-i1] += beta_y* g1(x,1.,pb);
		}
	}

	return F;
}

#define _MATRIX_RHS_CPP
#endif
