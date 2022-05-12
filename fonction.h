#pragma once

#include <stdio.h>
#include <vector>
#include <math.h>
#include <mpi.h>


double f1(double x, double y, double t,int pb, double Lx, double Ly); 
double g1(double x, double y, int pb); 
double h1(double x, double y, int pb);
double uex(double x, double y, int pb);
void charge(int i1,int iN, int me, int n, int Np);

