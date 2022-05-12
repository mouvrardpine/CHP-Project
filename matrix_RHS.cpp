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

	int he, size(_iN-_i1+1);
	double dx = 1./(_Nx+1) , dy = 1./(_Ny+1), msg;
	double alpha = 2*_dt*(1/(pow(dx,2)) +1 / pow(dy,2)) + 1 ;
	double beta_x = - _dt/pow(dx,2) , beta_y = - _dt/pow(dy,2) ;
	std::vector<double> Ax(size,0), grand_x(size + 2*_Nx);
	//cout << "alpha= "<< alpha << "	beta_x="<<beta_x<< "	beta_y="<<beta_y<<endl;

	// Envoi des _Nx premiers éléments de me à me-1 et réception des _Nx derniers éléments de me envoyés par me+1
	for (int i = 0; i < _Nx; i++) {
		MPI_Status Status;
		int tag = i +1;
		//cout<<" me "<<_me << "  j'envoie à "<<he <<" saucisse1"<<endl;
		if (_me !=0) MPI_Send(&x[i], 1, MPI_DOUBLE, _me-1, tag, MPI_COMM_WORLD);
		//cout<<" me "<<_me << "  j'attends à "<<he <<" saucisse1"<<endl;
		if (_me != _Np -1) MPI_Recv(&grand_x[grand_x.size()-1-(_Nx-i-1)], 1, MPI_DOUBLE, _me+1, tag, MPI_COMM_WORLD, &Status);
	}

	// Envoi des _Nx derniers éléments de me à me+1 et réception des _Nx derniers éléments de me envoyés par me-1
	for (int i = 0; i < _Nx; i++) {
		MPI_Status Status;
		int tag = 10*i +1;
		//cout<<" me "<<_me << "  j'envoie à "<<he <<" saucisse1"<<endl;
		if (_me != _Np-1) MPI_Send(&x[size- 1-(_Nx-i-1)], 1, MPI_DOUBLE, _me+1, tag, MPI_COMM_WORLD);
		//cout<<" me "<<_me << "  j'attends à "<<he <<" saucisse1"<<endl;
		if (_me != 0) MPI_Recv(&grand_x[i], 1, MPI_DOUBLE, _me-1, tag, MPI_COMM_WORLD, &Status);
	}

	//remplissage des autres éléments (éléments de me)
	for (int i = _Nx; i < _Nx + size ; i++)
	{
		grand_x[i] = x[i - _Nx];
	}

	 // for (int i = 0; i < grand_x.size(); i++)
	 // {
	 // 	cout << "test" << grand_x[i] << endl ;
	 // }

	int iloc;
	for (int iglob = _i1; iglob <= _iN; ++iglob)
	{
	 	iloc = iglob - _i1;
		Ax[iloc] = x[iloc]*alpha;
		//cout<< "i= "<< i << " me "<< _me<<endl;
    if (iglob < _Nx*_Ny - _Nx)
      {
				//cout << "ok" << endl;
        Ax[iloc] += beta_y*grand_x[iloc+2*_Nx];
      }
    if (iglob > _Nx-1)
      {
        Ax[iloc] += beta_y*grand_x[iloc];
      }

    if ((iglob+1)%_Nx != 0)
      {
        Ax[iloc] += beta_x*grand_x[iloc+1+_Nx];
      }

		if (iglob%_Nx !=0 && iglob != 0)
      {
        Ax[iloc] += beta_x*grand_x[iloc-1+_Nx];
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
		//	F[i] += beta_y*_fct->g1(x,1.,t+_dt);
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
