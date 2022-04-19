#ifndef _MATRIX_RHS_CPP

#include "matrix_RHS.h"
#include "fonction.h"
#include <cmath>
#include <iostream>

using namespace std;

matrix_RHS :: matrix_RHS(double dt , int Nx, int Ny, fonctions* fct,charge_* ch) : _dt(dt), _Nx(Nx) , _Ny(Ny), _fct(fct) , _ch(ch)
{
	MPI_Comm_rank(MPI_COMM_WORLD,&_me); 
  	MPI_Comm_size(MPI_COMM_WORLD,&_Np);
	_i1= _ch->Geti1();
	_n=ch->Getn();
	_iN=ch->GetiN();
}

std::vector<double>  matrix_RHS :: matvec(std::vector<double> x)
{

	int he;
	double dx = 1./(_Nx+1) , dy = 1./(_Ny+1), msg;
	double alpha = 2*_dt*(1/(pow(dx,2)) +1 / pow(dy,2)) + 1 ;
	double beta_x = - _dt/pow(dx,2) , beta_y = - _dt/pow(dy,2) ;
	std::vector<double> Ax(_iN-_i1+1,0);

	for (int i = _i1; i <= _iN; ++i)
	{
		Ax[i-_i1] = x[i-_i1]*alpha;
		cout<< "i= "<< i << " me "<< _me<<endl;
		if (i < _Nx*_Ny - _Nx)
		{
			if (i+_Nx > _iN)
			{
				he = _ch->whichMe(i+ _Nx);
				MPI_Status Status;
				//cout<<" me "<<_me << "  j'envoie à "<<he <<endl; 
				MPI_Send(&x[i-_i1], 1, MPI_DOUBLE, he, 101, MPI_COMM_WORLD);
				//cout<<" me "<<_me << "  j'attends à "<<he <<endl; 
				MPI_Recv(&msg, 1, MPI_DOUBLE, he, 102, MPI_COMM_WORLD, &Status);
				
				
				Ax[i-_i1] += beta_y*msg;
			}
			else
			{
				Ax[i-_i1] += beta_y*x[i-_i1+_Nx];
			}
		}
		if (i > _Nx-1)
		{
			if (i-_Nx < _i1)
			{
				he = _ch->whichMe(i- _Nx);
				MPI_Status Status;
				//cout<<"me "<<_me << "ok2"<< " j'envoie à "<<he <<endl; 
				MPI_Send(&x[i-_i1], 1, MPI_DOUBLE, he, 102, MPI_COMM_WORLD);
				//cout<<"me "<<_me << "ok2"<< " j'attends à "<<he <<endl; 
				MPI_Recv(&msg , 1, MPI_DOUBLE, he, 101, MPI_COMM_WORLD, &Status);
				
				Ax[i-_i1] += beta_y*msg;
			}
			else
			{
				Ax[i-_i1] += beta_y*x[i-_i1-_Nx];
			}
			
		}

		if ((i+1)%_Nx != 0)
		{
			if (i+1> _iN)
			{
				he = _ch->whichMe(i+1);
				MPI_Status Status;
				//cout<<"me "<<_me << "ok3"<< " j'envoie à "<<he <<endl; 
				MPI_Send(&x[i-_i1], 1, MPI_DOUBLE, he, 301, MPI_COMM_WORLD);
				//cout<<"me "<<_me << "ok3"<< " j'attends à "<<he <<endl; 
				MPI_Recv(&msg , 1, MPI_DOUBLE, he, 302, MPI_COMM_WORLD, &Status);
				
				Ax[i-_i1] += beta_x*msg;
			}
			else
			{
				Ax[i-_i1] += beta_x*x[i-_i1+1];
			}
			
			
		}
		if (i%_Nx !=0 && i != 0)
			if (i-1< _i1)
			{
				he = _ch->whichMe(i-1);
				MPI_Status Status;
				//cout<<"me "<<_me << "ok4"<< " j'envoie à "<<he <<endl; 
				MPI_Send(&x[i-_i1], 1, MPI_DOUBLE, he, 302, MPI_COMM_WORLD);

				//cout<<"me "<<_me << "ok4"<< " j'attends à "<<he <<endl; 
				MPI_Recv(&msg , 1, MPI_DOUBLE, he, 301, MPI_COMM_WORLD, &Status);
				
				Ax[i-_i1] += beta_x*msg;
			}
			else
			{
				Ax[i-_i1] += beta_x*x[i-_i1-1];
			}
	}

	return Ax;
}

std::vector<double> matrix_RHS :: RHS( std::vector<double> u, double t)
{
	double dx = 1./(_Nx+1) , dy = 1./(_Ny+1);
	double alpha = -2*_dt*(1/(pow(dx,2)) +1 / pow(dy,2)) + 1 ;
	double beta_x = _dt/pow(dx,2) , beta_y = _dt/pow(dy,2) ;
	std::vector<double> F(_Nx*_Ny,0.);
	int pb(2);

	for (int i = 0; i < _Nx*_Ny ; ++i)
	{
		double x = (i%_Nx+1)*dx, y = (i/_Nx+1)*dy; 

		F[i] = _dt*_fct->f1(x, y, t + _dt) + u[i];

		if ((i+1)%_Nx == 1)
		{
			F[i] += beta_x*_fct->h1(0,y,t+_dt);
		}

		if ((i+1)%_Nx == 0)
		{
			F[i] += beta_x*_fct->h1(1.,y,t+_dt);
		}

		if (i/_Nx == 0)
		{
			F[i] += beta_y*_fct->g1(x,0,t+_dt);
		}

		if (i/_Nx + 1 == _Ny)
		{
			F[i] += beta_y*_fct->g1(x,1.,t+_dt);
		}

		/*if ( (x == dx) || (x == 1.- dx))
		{
			F[i] += beta_x*h1(x,y,t+_dt, 1);   
		}

		if ( (y == dy) || (y == 1.-dy) )
		{
			F[i] += beta_y*g1(x,y,t+_dt, 1);
		}*/

	}

	return F;
}

#define _MATRIX_RHS_CPP
#endif
