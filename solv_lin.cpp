#ifndef _SOLV_LIN_CPP

#include "solv_lin.h"
#include "math.h"
#include <iostream>

using namespace std;

double ps(vector <double> x1, vector <double> x2)
{
    double res(0),sum(0) ;
    int n(x1.size());
    for (int i(0); i< n; i++)
    {
        res += x1[i]*x2[i];
    }
    //MPI_Allreduce(&res,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    //return sum;
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
  double res(0),sum(0);
    for (int i = 0; i < x.size(); i++) {
      res += pow(x[i],2);
    };
    MPI_Allreduce(&res,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    sum = sqrt(sum*dx*dy);
    
    return sum;
}

std::vector<double> GC(std::vector <double> x0 , std::vector <double> b )
{
    int k(0), n(x0.size()), test1(2),test2(3),cond(0),cond1; //test1  et test2 a supprimer : seulement pour esquiver bug de compilation
    vector<double> r(n,0),x(n,0),x1(n,0), d(n,0), z(n,0), rp(n,0),d1(n,0);
    double beta, gamma,alpha, beta_part, gamma_part, alpha_part, alpha_inter_part, alpha_inter, tau_part, tau ;

    x=x0;
    r=substract(b,matvec(x));
    d=r;

    /* for (int k(0);k<r.size();k++)
  {
  
   cout << "d["<<k<<"] = "<<d[k]<<endl;
  
} */
    beta_part= ps(r,r); //reduction pour la norme de r
    MPI_Allreduce(&beta_part,&beta,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    beta = sqrt(beta);
  
   // cout<<"beta = " << beta<<endl;
    while (cond1<1)
    {
        //cout<<"saucisse1"<<_me<<endl;
        z=matvec(d);
        //cout<<"saucisse2"<<_me<<endl;
        gamma_part=ps(r,r); 
        cout<<"k= "<<k <<"me = "<<_me<< "gamma_part = " << gamma_part << endl;
        MPI_Allreduce(&gamma_part,&gamma,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        cout<<"k= "<<k<<"me = "<<_me<< "gamma = " << gamma << endl;
        alpha_inter_part=ps(z,d); 
        MPI_Allreduce(&alpha_inter_part,&alpha_inter,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        
        alpha=gamma/alpha_inter; //reduction pour le produit scalaire  !!!!!!!!!!!!!!!!!!!!!! vérification en print sur le gamma!!!!!!!!!!!!!!!!
       
        x1=add(x,mult(alpha,d));
       // cout<<"saucisse4"<<_me<<endl;
        x=x1; 
        //cout<<"saucisse5"<<_me<<endl;
        rp=substract(r,mult(alpha,z));
        //cout<<"saucisse6"<<_me<<endl;
        tau_part= ps(rp,rp); //reduction pour la norme de r
        MPI_Allreduce(&tau_part,&tau,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        d1= add(rp,mult((tau/pow(beta,2)),d));
        //cout<< "saucisse7 "<< _me <<" k= "<< k<<endl;
        d=d1; // reduction pour le ps rp 
        r=rp;
        beta_part= ps(r,r); //reduction pour la norme de r
        MPI_Allreduce(&beta_part,&beta,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        beta = sqrt(beta);
        //cout<<"saucisse7"<<_me<<endl;
          //racine de ps(rp,rp)
        //cout << beta << endl;
        k=k+1;
          //racine de ps(rp,rp)
        //cout<<"k= "<<k<< "  beta= " << beta<< "  beta reel =" << norm(substract(b,matvec(_dt,_Nx, _Ny,x))) << endl;
        if(!((beta>_eps)&&(k<_kmax)))
        {
            cond=1;
        }
        MPI_Allreduce(&cond,&cond1,1,MPI_INT,MPI_PROD,MPI_COMM_WORLD);
        //cout<<"k = " << k <<"beta = " << beta << "me" << _me <<endl;
        //cout<< "k = " << k << " me = " << _me <<" gamma = " << gamma << " alpha = "<<alpha<<"beta = " << beta <<endl;

    }
    
    //cout<<"saucisse10"<<endl;
    if (k>_kmax)
    {
        printf("tolérance non atteinte");
    }
    //cout<<"beta = " << beta<<endl;
    return x;

}

#define _SOLV_LIN_CPP
#endif
