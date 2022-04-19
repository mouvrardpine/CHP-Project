#ifndef _SOLV_LIN_CPP

#include "solv_lin.h"
#include "math.h"
#include <iostream>

using namespace std;

solv_lin :: solv_lin(int kmax, int Nx , int Ny , double eps , double dt, matrix_RHS* mrhs,charge_* ch) : _kmax(kmax), _Nx(Nx) , _Ny(Ny), _eps(eps), _dt(dt) ,_mRHS(mrhs),_ch(ch)
{
    MPI_Comm_rank(MPI_COMM_WORLD,&_me); 
  	MPI_Comm_size(MPI_COMM_WORLD,&_Np);
	_i1= _ch->Geti1();
	_n= _ch->Getn();
	_iN= _ch->GetiN();
}

 double solv_lin :: ps(vector <double> x1, vector <double> x2)
{
    double res(0),sum(0) ;
    int n(x1.size());
    for (int i(0); i< n; i++)
    {
        res += x1[i]*x2[i];
    }
    MPI_Reduce(&res,&sum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    return sum;

}

 vector <double> solv_lin :: substract(vector <double> x1, vector <double> x2)
{
    int n(x1.size());
    vector<double> res(n,0) ;

    for (int i(0); i< n; i++)
    {
        res[i] =x1[i]-x2[i];
    }
    return res;

}


 vector<double> solv_lin :: add(vector <double> x1, vector <double> x2)
{
int n(x1.size());
vector<double> res(n,0) ;

    for (int i(0); i< n; i++)
    {
        res[i] =x1[i]+x2[i];
    }
    return res;

}


 double solv_lin :: norm(vector<double> x)
{
    return sqrt(ps(x,x));
}

vector <double> solv_lin :: mult(double a , vector<double> x)
{
    int n(x.size());
    vector<double> res(n,0);

    for (int i(0); i< n;i++)
    {
       res[i]=a*x[i];
    }
    return res ;
}

double solv_lin :: normL2_2D(vector<double> x, double dx, double dy)
{
  double res(0),sum(0);
    for (int i = 0; i < x.size(); i++) {
      res += pow(x[i],2);
    };
    MPI_Reduce(&res,&sum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    sum = sqrt(sum*dx*dy);
    
    return sum;
}

std::vector<double> solv_lin :: GC(std::vector <double> x0 , std::vector <double> b )
{
    int k(0), n(x0.size()), test1(2),test2(3); //test1  et test2 a supprimer : seulement pour esquiver bug de compilation
    vector<double> r(n,0),x(n,0), d(n,0), z(n,0), rp(n,0);
    double beta, gamma,alpha ;

    x=x0;
    r=substract(b,_mRHS->matvec(x));
    d=r;
    beta= norm(r); //reduction pour la norme de r
    while ((beta>_eps)&&(k<_kmax))
    {
        z=_mRHS->matvec(d);
        gamma=beta*beta; 
        alpha=gamma/ps(z,d); //reduction pour le produit scalaire  !!!!!!!!!!!!!!!!!!!!!! vérification en print sur le gamma!!!!!!!!!!!!!!!!
        x=add(x,mult(alpha,d)); 
        rp=substract(r,mult(alpha,z));
        d= add(rp,mult((ps(rp,rp)/pow(beta,2)),d)); // reduction pour le ps rp 
        r=rp;
        beta = norm(r);  //racine de ps(rp,rp)
        //cout << beta << endl;
        k=k+1;
        //cout<<"k= "<<k<< "  beta= " << beta<< "  beta reel =" << norm(substract(b,matvec(_dt,_Nx, _Ny,x))) << endl;
    }

    if (k>_kmax)
    {
        printf("tolérance non atteinte");
    }
    return x;

}



#define _SOLV_LIN_CPP
#endif
