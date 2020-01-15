#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <string>

using namespace std;

double quadV(double x){
  return 0.5*x*x;
}

double quartV(double x){
  return x*x*x*x + 0.5*x*x;
}

double sixV(double x){
  return x*x*x*x*x*x + x*x*x*x + 0.5*x*x;
}

int main(int argc, char *argv[]){

  int N;
  double dx = 0.1, x;
  
  if(argc < 2){
    cerr << "Usage: " << argv[0] << " numpts [dx]" << endl;
    cerr << "Where numpts is the number of points to generate in the potential with spacing 0.001" << endl;
    cerr << "If numpts is odd, it will be incremented by 1 to maintain a symmetric potential about the origin" << endl;
    cerr << "dx is the optional spacing of the data, default is 0.1" << endl;
    return -1;
  } else {
    N = atoi(argv[1]);
    if(N%2 == 0)
      N++;
    x = -N/2*dx;
    if(argc == 3){
      dx = atof(argv[2]);
    }
  }
  
  // generate file names
  string quadfile = "quadratic_" + to_string(N) + ".dat";
  string quartfile = "quartic_" + to_string(N) + ".dat";
  string sixfile = "6tic_" + to_string(N) + ".dat";
  
  // generate quadratic and quartic potentials

  // quadratic (harmonic oscillator) and quartic (oscillator with quartic anharmonicity)
  ofstream fqd(quadfile.c_str()), fqt(quartfile.c_str()), fsx(sixfile.c_str());
  fqd << N << endl;
  fqt << N << endl;
  fsx << N << endl;
  
  for(int i = 0; i < N; i++){
    fqd << x << "\t" << quadV(x) << endl;
    fqt << x << "\t" << quartV(x) << endl;
    fsx << x << "\t" << sixV(x) << endl;
    x += dx;
  }

  fqd.close();
  fqt.close();
  fsx.close();
 

  return 0;
  
}
