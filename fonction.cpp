#include "fonction.h"


fonctions::fonctions(int pb, double Lx, double Ly) : _pb(pb) , _Lx(Lx) , _Ly(Ly)
{
}

fonctions::~fonctions()
{
}



double fonctions :: f1(double x, double y, double t)
{
    if (_pb==1)
    {
        return 2*(y-y*y+x-x*x); 
    }
    else if (_pb==2)
    {
        return sin(x)+cos(y);
    }
    else
    {
        return exp(-pow(x-(_Lx/2),2))*exp(-pow(y-(_Ly/2),2))*cos((3.14/2)*t);
    }
}


double fonctions :: g1(double x, double y, double t)
{
    if (_pb==1)
    {
        return 0; 
    }
    else if (_pb==2)
    {
        return sin(x)+cos(y);
    }
    else
    {
        return 0;
    }
}

double fonctions :: h1(double x, double y, double t)
{

    if (_pb==1)
    {
        return 0; 
    }
    else if (_pb==2)
    {
        return sin(x)+cos(y);
    }
    else
    {
        return 1;
    }
}

double fonctions :: uex(double x, double y)
{
    if (_pb==1)
    {
        
        return x*(1-x)*y*(1-y); 
    }
    else if (_pb==2)
    {
        return sin(x)+cos(y);
    }
    else
    {
        return 1;
    }

}

charge_ :: charge_ (int n) : _n(n)
{
    MPI_Comm_rank(MPI_COMM_WORLD,&_me);
    MPI_Comm_size(MPI_COMM_WORLD,&_Np);
    charge( &_i1, &_iN, _me, _n, _Np);
    

    /* int r = _n%_Np;
    if (_me<r)
    {
        _i1=_me*(_n/_Np+1);
        _iN= _i1 + (_n/_Np+1)-1;

    }
    else 
    {
        _i1=r+_me*(_n/_Np);
        _iN=_i1+(_n/_Np)-1;
    }     */
}

void charge_ :: charge(int *i1,int *iN, int me, int n, int Np)
{
  int r = n%Np;
  if (me<r)
  {
    *i1=me*(n/Np+1);
    *iN= *i1 + (n/Np+1)-1;

  }
  else 
  {
    *i1=r+me*(n/Np);
    *iN=*i1+(n/Np)-1;
  }    
}

int charge_ :: whichMe(int i)
{
    int i1, iN;
    for (int k(0); k<_Np;k++)
    {
        charge(&i1,&iN,k,_n,_Np);
        if (( _n<=iN )&&(_n>= i1))
        {
            return k;   
        }
    }
    std::cout<<"erreur dans WhichME"<<std::endl;
}