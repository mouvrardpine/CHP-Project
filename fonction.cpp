#include "fonction.h"
#include "math.h"


double f1(double x, double y, double t,int pb, double Lx, double Ly)
{
    if (pb==1)
    {
        return 2*(y-y*y+x-x*x);
    }
    else if (pb==2)
    {
        return sin(x)+cos(y);
    }
    else
    {
        return exp(-pow(x-(Lx/2),2))*exp(-pow(y-(Ly/2),2))*cos((M_PI/2)*t);
    }
}


double g1(double x, double y, int pb)
{
    if (pb==1)
    {
        return 0;
    }
    else if (pb==2)
    {
        return sin(x)+cos(y);
    }
    else
    {
        return 0;
    }
}

double h1(double x, double y, int pb)
{

    if (pb==1)
    {
        return 0;
    }
    else if (pb==2)
    {
        return sin(x)+cos(y);
    }
    else
    {
        return 1;
    }
}

double uex(double x, double y, int pb)
{
    if (pb==1)
    {

        return x*(1-x)*y*(1-y);
    }
    else if (pb==2)
    {
        return sin(x)+cos(y);
    }
    else
    {
        return 1;
    }

}

void charge(int i1,int iN, int me, int n, int Np)
{
  int r = n%Np;
  if (me<r)
  {
    i1 = me*(n/Np+1);
    iN = i1 + (n/Np+1)-1;

  }
  else
  {
    i1 = r +me*(n/Np);
    iN = i1+(n/Np)-1;
  }

}

// int charge_ :: whichMe(int i)
// {
//     int i1, iN;
//     for (int k(0); k<_Np;k++)
//     {
//         charge(&i1,&iN,k,_n,_Np);
//         if (( i<=iN )&&(i>= i1))
//         {
//             return k;
//         }
//     }
//     std::cout<<"erreur dans WhichME"<<std::endl;
// }
