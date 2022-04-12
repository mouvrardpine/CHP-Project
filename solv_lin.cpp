#ifndef _SOLV_LIN_CPP

#include "solv_lin.h"
#include "math.h"
#include <iostream>

using namespace std;
double ps(vector <double> x1, vector <double> x2)
{
    double res(0) ;
    int n(x1.size());
    for (int i(0); i< n; i++)
    {
        res += x1[i]*x2[i];
    }
    return res;

}

vector <double> substract(vector <double> x1, vector <double> x2)
{
    int n(x1.size());
    vector<double> res(n,0) ;

    for (int i(0); i< n; i++)
    {
        res[i] =x1[i]-x2[i];
    }
    return res;

}


vector<double> add(vector <double> x1, vector <double> x2)
{
int n(x1.size());
vector<double> res(n,0) ;

    for (int i(0); i< n; i++)
    {
        res[i] =x1[i]+x2[i];
    }
    return res;

}


double norm(vector<double> x)
{
    return sqrt(ps(x,x));
}

vector <double> mult(double a , vector<double> x)
{
    int n(x.size());
    vector<double> res(n,0);

    for (int i(0); i< n;i++)
    {
       res[i]=a*x[i];
    }
    return res ;
}

double normL2_2D(vector<double> x, double dx, double dy)
{
  double res(0);
    for (int i = 0; i < x.size(); i++) {
      res += pow(x[i],2);
    };

    res = sqrt(res*dx*dy);
    return res;
}

std::vector<double> GC(std::vector <double> x0 , std::vector <double> b , double eps , int kmax,int Nx, int Ny,double dt)
{
    int k(0), n(x0.size()), test1(2),test2(3); //test1  et test2 a supprimer : seulement pour esquiver bug de compilation
    vector<double> r(n,0),x(n,0), d(n,0), z(n,0), rp(n,0);
    double beta, gamma,alpha ;

    x=x0;
    r=substract(b,matvec(dt,Nx, Ny,x,test1,test2));
    d=r;
    beta= norm(r); //reduction pour la norme de r
    while ((beta>eps)&&(k<kmax))
    {
        z=matvec(dt,Nx,Ny,d,test1,test2);
        gamma=beta*beta; 
        alpha=gamma/ps(z,d); //reduction pour le produit scalaire  !!!!!!!!!!!!!!!!!!!!!! vérification en print sur le gamma!!!!!!!!!!!!!!!!
        x=add(x,mult(alpha,d)); 
        rp=substract(r,mult(alpha,z));
        d= add(rp,mult((ps(rp,rp)/pow(beta,2)),d)); // reduction pour le ps rp 
        r=rp;
        beta = norm(r);  //racine de ps(rp,rp)
        //cout << beta << endl;
        k=k+1;
        //cout<<"k= "<<k<< "  beta= " << beta<< "  beta reel =" << norm(substract(b,matvec(dt,Nx, Ny,x))) << endl;
    }

    if (k>kmax)
    {
        printf("tolérance non atteinte");
    }
    return x;

}



#define _SOLV_LIN_CPP
#endif
