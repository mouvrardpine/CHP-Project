#ifndef _SOLV_LIN_CPP

#include "solv_lin.h"
#include "math.h"
#include <iostream>

using namespace std;


//Produit scalaire
double ps(vector <double> x1, vector <double> x2)
{
    double res(0),sum(0) ;
    int n(x1.size());
    for (int i(0); i< n; i++)
    {
      res += x1[i]*x2[i];
    }
    return res;

}

//Soustraction de deux vecteurs
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

//Addition de vecteurs
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

//Norme euclidienne d'un vecteur
double norm(vector<double> x)
{
    return sqrt(ps(x,x));
}

//Multiplication d'un vecteur par un réel
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

//Norme L2 (pour le calcul de l'erreur quadratique)
double normL2_2D(vector<double> x, double dx, double dy)
{
  double res(0),sum(0);
  for (int i = 0; i < x.size(); i++) {
    res += pow(x[i],2);
  };
  sum = sqrt(res*dx*dy);

  return sum;
}

//Gradient conjugué
std::vector<double> GC(std::vector <double> x0 , std::vector <double> b, int i1, int iN, int kmax,int me, int Nx, int Ny, int Np, double dt, double eps )
{/* code */
    //Déclaration/initisalisation des variables
    int k(0), n(x0.size());
    vector<double> r(n,0),x(n,0), d(n,0), z(n,0), rp(n,0),d1(n,0);
    double beta, gamma,alpha, beta_part, gamma_part, alpha_part, alpha_inter_part, alpha_inter, tau_part, tau ;

    //Initialisation de x
    x=x0;

    //Calcul du résidu initial
    r=substract(b,matvec(x, i1, iN, me, Nx, Ny, Np, dt));
    d=r;

    //Calcul de la norme du résidu
    beta_part= ps(r,r);
    MPI_Allreduce(&beta_part,&beta,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    beta = sqrt(beta);


    //Boucle sur k
    while ((beta>eps)&&(k<kmax))
    {
      z=matvec(d, i1, iN,me, Nx,Ny,Np,dt);

      gamma_part=ps(r,r);
      MPI_Allreduce(&gamma_part,&gamma,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

      alpha_inter_part=ps(z,d);
      MPI_Allreduce(&alpha_inter_part,&alpha_inter,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      alpha=gamma/alpha_inter;

      //Mise à jour de x
      x=add(x,mult(alpha,d));

      //Calcul du nouveau résidu
      rp=substract(r,mult(alpha,z));

      //Calcul et mise à jour de d
      tau_part= ps(rp,rp);
      MPI_Allreduce(&tau_part,&tau,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      d= add(rp,mult((tau/pow(beta,2)),d));

      //Mise à jour du résidu
      r=rp;
      beta_part= ps(r,r);
      MPI_Allreduce(&beta_part,&beta,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      beta = sqrt(beta);

      //Mise à jour de k
      k=k+1;
    }

    if (k>kmax)
    {
        printf("tolérance non atteinte");
    }
    return x;

}

#define _SOLV_LIN_CPP
#endif
