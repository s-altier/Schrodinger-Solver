#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "matrixmath.h"
#include <cstring>


using namespace std;

double pureQ(int n){
  return 0.867145*pow(((double)n+0.5), 4/3.)*pow(4., 1/3.);
}

int main(int argc, char *argv[]) {
  
  ofstream fout("pure_quartic.dat");

  fout << 26 << endl;
  for(int i = 0; i < 26; i++){
    fout << i << "\t" << pureQ(i) << endl;
  }
  fout.close();


  FILE *gnuplot = popen("gnuplot", "w");
  fprintf(gnuplot, "set out 'pure_quartic_eigenvalues.ps'\n");
  fprintf(gnuplot, "set term post land\n");
  fprintf(gnuplot, "set xlabel 'n'\n");
  fprintf(gnuplot, "set ylabel 'Eigenvalue'\n");
  

  fprintf(gnuplot, "plot 'pure_quartic.dat' u 1:2\n");
  pclose(gnuplot);

  system( (char *) "ps2pdf pure_quartic_eigenvalues.ps");
  
  return 0;
}

