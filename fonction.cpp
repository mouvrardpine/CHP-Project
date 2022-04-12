#include "fonction.h"


double f1(double x, double y, double t, int pb, double Lx, double Ly)
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
        return exp(-pow(x-(Lx/2),2))*exp(-pow(y-(Ly/2),2))*cos((3.14/2)*t);
    }
}


double g1(double x, double y, double t, int pb)
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

double h1(double x, double y, double t, int pb)
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