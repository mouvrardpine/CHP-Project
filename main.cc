#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string>
#include <mpi.h>
#include "fonction.h"
#include "matrix_RHS.h"
#include "solv_lin.h"
#include <ctime>

int main(int argc, char ** argv)
{


  // Lecture des paramètres

  //parametres
  int Nx, Ny,kmax,pb;
  double dt,dx,dy, Lx, Ly, tmax;
  ifstream fichier("parametres.txt", ios::in);

  if(fichier)
  {

  fichier >> Nx >> Ny >> kmax  >> pb >> dt >> Lx >> Ly >> tmax;  /*on lit jusqu'à l'espace et on stocke ce qui est lu dans la variable indiquée */

  fichier.close();
  }
  else
  cerr << "Impossible d'ouvrir le fichier !" << endl;

  int k(0),i1,iN,me,Np;
  double eps(pow(10,-16)),err(0), t(0.);

//démarrage du chrono
  std::clock_t c_start = std::clock();
  //your_algorithm



  // Initialiser l'environnement
  MPI_Init(&argc, &argv);

  // Numéro du processeur actif
  MPI_Comm_rank(MPI_COMM_WORLD,&me);

  // Nombre de processeurs
  MPI_Comm_size(MPI_COMM_WORLD,&Np);

  MPI_Status status ;

  // Taille globale du vecteur
  int n = Nx*Ny;

  // Calcul de la charge du processeur me
  charge(&i1,&iN, me,n,Np);

  // Taille locale du vecteur
  int size = iN - i1 + 1;

  //calcul du max des tailles

  int max_size;
  MPI_Allreduce(&size,&max_size,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);



  //Initialisation des vecteurs
  std::vector<double> u(size,0),b(size,0),x0(size,1), uex_vec(size,0);

  // Calcul des pas du maillage
  dx=1./(Nx+1);
  dy=1./(Ny+1);

  //Initialisation de la solution (premier pas de temps)
  b=RHS(u,t,Nx,Ny,i1, iN,Lx,Ly,dt,pb);
  u=GC(u,b, i1, iN, kmax, me, Nx, Ny, Np, dt, eps);


  // Calcul de la solution exacte (solutions stationnaires)
  for (int i(i1); i<=iN ; i++)
  {
      uex_vec[i-i1]=uex((i%Nx+1)*dx, (i/Nx+1)*dy, pb);
  }

  // Boucle en temps
  while (t<tmax)
  {
    //Mise à jour du temps t
    t += dt;

    //Calcul de la solution
    b=RHS(u,t,Nx,Ny,i1, iN,Lx,Ly,dt,pb);
    u=GC(u,b, i1, iN, kmax, me, Nx, Ny, Np, dt, eps);

    //écriture du fichier solution
    string filename("Output/sol" +to_string(me)+ to_string(int(t/dt)) +".dat");
    fstream file_out;

    file_out.open(filename, std::ios_base::out);
    if (file_out) {

      for (int i(i1);i<=iN; i++ )
      {
        file_out<< (i%Nx+1)*dx << " " << (i/Nx+1)*dy << " " << u[i-i1] <<endl;
      }
    }
    else
    {
    cout << "failed to open " << filename << '\n';
    }

    //Calcul des erreurs
    err = normL2_2D(substract(u,uex_vec), dx, dy);

    //Affichage des erreurs pour chaque processeur
    cout<< "t = " << t << "me = " << me << " err=  "<<err<< endl;
    }


    // Fin du chrono
    std::clock_t c_end = std::clock();

    double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;

    //affichage

    double total_time_ms;

    MPI_Allreduce(&time_elapsed_ms,&total_time_ms,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

    if (me ==0) {
      std::cout << "CPU time used: " << total_time_ms / 1000.0 << " s\n";
      //efficacité

      double E = float(n)/float(Np*max_size);

      cout << "Efficacité :" << E << endl;
    }

    MPI_Finalize();

    return 0;

}
