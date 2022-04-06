#include "matrix_RHS.h"

using namespace std;
double ps(vector <double> x1, vector <double> x2);
vector<double> substract(vector <double> x1, vector <double> x2);
vector<double> add(vector <double> x1, vector <double> x2);
double norm(vector <double> x);
vector <double> mult(double a , vector<double> x);
vector<double> GC(vector <double> x0 , vector <double> b , double eps , int kmax,int Nx, int Ny,double dt);


