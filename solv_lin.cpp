#include "solv_lin.h"

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

vector<double> GC(vector <double> x0 , vector <double> b , double eps , int kmax,int Nx, int Ny,double dt)
{
    int k(0), n(x0.size());
    vector<double> r(n,0),x(n,0), d(n,0), z(n,0), rp(n,0);
    double beta, gamma,alpha ; 
    


    x=x0; 
    r=substract(b,matvec(dt,Nx, Ny,x));
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
    }
    return x; 

}



