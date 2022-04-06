#include "solv_lin.h"

using namespace std;
double ps(vector <double> x1, vector <double> x2)
{
    double res(0) ;
    for (int i(0); i< x1.size(); i++)
    {
        res += x1[i]*x2[i];
    }
    return res; 

}

vector <double> substract(vector <double> x1, vector <double> x2)
{
    vector<double> res(x1.size(),0) ;
    for (int i(0); i< x1.size(); i++)
    {
        res[i] =x1[i]-x2[i];
    }
    return res; 

}


vector<double> add(vector <double> x1, vector <double> x2)
{

vector<double> res(x1.size(),0) ;
    for (int i(0); i< x1.size(); i++)
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
    vector<double> res(x.size(),0);
    for (int i(0); i< x.size();i++)
    {
       res[i]=a*x[i];
    }
    return res ;
}

vector<double> GC(vector <double> x0 , vector <double> b , double eps , int kmax,int Nx, int Ny,double dt)
{
    
   /* vector<double> r(x0.size(),0),x(x0.size(),0), d(x0.size(),0), z(x0.size(),0), rp(x0.size(),0);
    double beta, gamma,alpha ; 
    int k ;


    x=x0; 
    /* r=b-matvec(dt,Nx, Ny,x);  */
    
    /*r=substract(b,matvec(dt,Nx, Ny,x));
    d=r; 
    beta= norm(r);
    while ((beta>eps)&&(k<kmax))
    {
        z=matvec(dt,Nx,Ny,d);
        gamma=ps(r,r);
        alpha=gamma/ps(z,d); 
        x=add(x,mult(alpha,d));
        rp=substract(r,mult(alpha,z));
        d= add(rp,mult(ps(rp,rp),d));
        r=rp; 
        beta = norm(r);
        k=k+1; 
    }
     
    if (k>kmax)
    {
        printf("tolérance non atteinte");
    }*/
    return x; 

}



