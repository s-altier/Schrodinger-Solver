#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <string>
#include "matrixmath.h"

using namespace std;

double Off(const vector< vector<double> >& B);

double indMaxEl(const vector< vector<double> > &B, int& i, int& j);

double getTan(const vector< vector<double> >& S);

int jacobi(const vector< vector<double> > &S, double tol, vector<double> &E, int M);

void ascendingSort(vector<double> &E);

int main(int argc, char* argv[]){

  // NOTE: we are assuming hbar = 1, we take the product m*omega^2 to be 1, and we want the

  
  string fname;
  int N, M;
  double tol = 0.0;
  double omega = 1.0;

  if(argc < 5){
    cerr << "Usage: " << argv[0] << " fname M tol omegaFactor [all]" << endl;
    cerr << "Where fname is the name of the file containing the potential in the format specified in the readme.txt" << endl;
    cerr << "M is the maximum number of Jacobi rotations to execute, tol is the tolerance to rotate to" << endl;
    cerr << "omegaFactor is the factor by which to multiply the discretization spacing when choosing a fundamental frequency "
         << "for the potential" << endl;
    cerr << "[all] is an optional command that will output all of the errors in eigenvalue, instead of capping it at 50,"
         << " any argument passed here will trigger this." << endl;
    return -1;
  } else {
    fname = argv[1];
    M = atoi(argv[2]);
    tol = atof(argv[3]);
    omega = atof(argv[4]);
    if(argc == 6)
      doAll = true;
  }
  cout << "Running Schrodinger Eigen-Energy Solver with Dirichlet Boundary Conditions" << endl;
  
  // start reading from file
  string line;
  ifstream fin(fname.c_str());
  getline(fin, line);
  stringstream ss(line);
  ss >> N;

  vector<double> V(N, 0.0), x(N, 0.0);

  int j = 0;
  while(getline(fin, line)){
    stringstream iss(line);
    iss >> x[j] >> V[j];
    j++;
  } // end while
  
  // construct our matrix of form:
  /*      [  2 + V1/t0     -1     ... 0 ...     0    ]
   *      [     -1     1 + V2/t0  ...   ...     0    ]
   * t0 * [      0         ...          ...    ...   ]
   *      [           ...  ...   ...    ...    ...   ]
   *      [      0    ...   0    ... ...   2 + Vn/t0 ]
   *
   */

  double h = 0.001;
  double t0 = omega*omega/2.;
  vector<double> rows(N, 0.0);
  vector< vector<double> > H(N, rows);

  for(int i = 0; i < N; i++){
    H[i][i] = 2.*t0 + V[i];
    if(i != 0){
      H[i][i-1] = -t0;
      H[i-1][i] = -t0;
    } // endif
  } // end i-loop


  // done constructing matrix

  // solve
  vector<double> E(N, 0.0);
  int reps = 0;
  reps = jacobi(H, tol, E, M);

  cout << "Took " << reps << " iterations" << endl;

  
  ascendingSort(E);
  print(E);

  // calculate theoretical eigenvalues
  int size;
  if(doAll){
    size = N;
  }
  else if(N<150){
    size = N;
  }
  else
    size = 150;
  vector<double> theoryE(size, 0.0);
  for(int i = 0; i < size; i++){
    theoryE[i] = 0.1*omega*(i+0.5);
  }

  string efile = "errors_dirichlet_omega" + to_string((int)omega) + "_N" + to_string(H.size()) + ".dat";
  ofstream fe(efile.c_str());

  fe << size << endl;
  for(int i = 0; i < size; i++){
    fe << i << "\t" << abs(E[i]-theoryE[i]) << "\t" << 1.0 << "\t" << E[i] << "\t" << theoryE[i] << endl; 
  }
  fe.close();
  
  return 0;
}

double Off(const vector< vector<double> >& B){

  double off = 0.0;
  int nrow = B.size();

  for(int i = 0; i < nrow-1; i++){
    for(int j = i + 1; j < nrow; j++){
      off += 2 * abs(B[i][j])*abs(B[i][j]);
    }
  }
  
  return off;
}

double indMaxEl(const vector< vector<double> > &B, int& i, int& j){

  int nrow = B.size();

  double max = 0.0;

  int maxi, maxj;
  
  for(int m = 0; m < nrow-1; m++){
    for(int n = m+1; n < nrow; n++){
      if(abs(B[m][n]) > max){
	max = abs(B[m][n]);
	maxi = m;
	maxj = n;
      }
    }
  }

  i = maxi;
  j = maxj;

  return max;

}


// determines the rotation angle to eliminate an off-diagonal element, using a 2x2 matrix
double getTan(const vector< vector<double> >& S, const int& i, const int& j){

  double tan = 0.0, tau = 0.0;

  double t = 0.0;
  double signBeta = 0.0;
  tau = (S[j][j] - S[i][i])/2/S[i][j];

  double signTau = (tau<0) ? -1 : ((tau>0) ? 1 : 0);

  tan = signTau / (abs(tau) + sqrt(tau*tau + 1));

  return tan;
  
}


// jacobi algorithm
int jacobi(const vector< vector<double> >& S, double tol, vector<double>& E, int M){

  int N = S.size();
  vector< vector<double> > B = S;

  int i, j;
  double t = 0.0, c, s;
  double offB = Off(B);
  int numiter = 0;


  indMaxEl(B, i, j);
  
    
  while(offB > tol && numiter < M){

    // storage variables for before matrix update
    double bii = 0., bjj = 0., bij = 0.;
    vector<double> bik(N, 0.0), bjk(N, 0.0);

    
    indMaxEl(B, i, j);

    t = getTan(B, i, j);

    // rotation elements
    c = 1 / sqrt(t*t + 1);
    s = c*t;

    // partial "diagonalization"

    bii = c*c*B[i][i] - 2*s*c*B[i][j] + s*s*B[j][j];
    bjj = s*s*B[i][i] + 2*s*c*B[i][j] + c*c*B[j][j];
    
    bij = (c*c-s*s)*B[i][j] + s*c*(B[i][i]-B[j][j]);

    // bik, bjk values
    for(int k = 0; k < N; k++){
      if(k != i && k != j){
	bik[k] = c*B[i][k] - s*B[j][k];
	bjk[k] = s*B[i][k] + c*B[j][k];
      } else
	bjk[k] = 0.0;
    }
    
    // perform update

    B[i][i] = bii;
    B[j][j] = bjj;
    B[i][j] = bij;
    B[j][i] = bij;

    for(int k = 0; k < N; k++){
      if(k != i && k != j){
	B[i][k] = bik[k];
	B[k][i] = bik[k];
	B[j][k] = bjk[k];
	B[k][j] = bjk[k];
      }
    }      

    numiter++;
    offB = Off(B);
  }

  

  // records diagonal eigenvalues
  for(int n = 0; n < N; n++){
    E[n] = B[n][n];
  }

  return numiter;
}

void ascendingSort(vector<double> &E){

  int N = E.size();

  for(int i = 0; i < N-1; i++){
    int min = i;
    for(int j = i+1; j < N; j++){
      if(abs(E[j]) < abs(E[min]))
	min = j;
    }

    if(min != i){
      double temp = E[min];
      E[min] = E[i];
      E[i] = temp;
    }

  }
}
