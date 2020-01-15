#include <iostream>
#include <fstream>
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

int main(int argc, char* argv[]){

  

  vector< vector<double> > A(3, vector<double>(3, 0.0));
  A[0][0] = 8.;
  A[0][1] = -2.;
  A[0][2] = -20.;
  A[1][0] = -2.0;
  A[1][1] = -3.;
  A[1][2] = -9.;
  A[2][0] = -20.;
  A[2][1] = -9.;
  A[2][2] = -3.;

  vector<double> E(3, 0.0);
  
  cout << Off(A) << endl;
  int reps = jacobi(A, 1e-8, E, 10000);

  print(E);
  cout << "Num. reps: " << reps << endl;
  
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

/*

// jacobi algorithm
int jacobi(const vector< vector<double> >& S, double tol, vector<double>& E, int M){

  int N = S.size();
  vector< vector<double> > B = S;

  int i, j;
  double t = 0.0, c, s;
  
  int numiter = 0;

  indMaxEl(B, i, j);
  
  cout << B[i][j] << endl;
  
  print(B);
    
  while(Off(B) > tol && numiter < N){

    // storage variables for before matrix update
    double bii = 0., bjj = 0., bij = 0.;
    vector<double> bik(N/2, 0.0), bjk(N/2, 0.0);

    
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
  }

  print(B);

  // records diagonal eigenvalues
  for(int n = 0; n < N; n++){
    E[n] = B[n][n];
  }

  return numiter;
}

*/
