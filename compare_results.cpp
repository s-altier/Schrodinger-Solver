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

double poly(vector <double> a, double x);
double chi2(int N, vector <double> a, vector <double>&x, vector <double> &y, vector <double> &err);

double poly(vector <double> a, double x){
  int M = a.size();
  double value = 0.0;
  for(int i = 0; i < M; i++){
    value += a[i]*pow(x,i);
  }

  return value;
}

double chi2(int N, vector <double> a, vector <double>&x, vector <double> &y, vector <double> &err){

  double value = 0.0;

  for (int i = 0; i < N; i++){
    value += pow((poly(a,x[i]) - y[i])/err[i] , 2);
  }

  return value; 
}

int main(int argc, char *argv[]) {
  int numdata, numparams;
  int i=0, j=0;
  
  string fname1, fname2;
  if(argc<3){
  	cerr << "Usage: "<<argv[0]<<" file1 file2"<<endl;
	cerr << "file1 and file2 are the two files to be compared, they must have the same number of points and follow the format from readme.txt" << endl;
  	return -1;
  }
  else{
    fname1 = argv[1];
    fname2 = argv[2]; 
  }

  int N1, N2;

  string line;
  ifstream file1(fname1.c_str()), file2(fname2.c_str());

  getline(file1, line);
  istringstream iss(line);
  iss >> N1;
  getline(file2, line);
  iss = istringstream(line);
  iss >> N2;

  if(N1!=N2){
    cerr << "Can only compare files of equal length" << endl;
    return -1;
  }

  vector<double> values1(N1, 0.0);
  vector<double> values2(N1, 0.0);
  vector<double> values(N1, 0.0);

  int k = 0;
  double temp;
  while(getline(file1, line)){
    iss = istringstream(line);
    iss >> temp >> values1[k];
    k++;
  }
  file1.close();

  k = 0;
  while(getline(file2, line)){
    iss = istringstream(line);
    iss >> temp >> values2[k];
    k++;
  }
  file2.close();
  

  // generate difference data
  for(int i = 0; i < N1; i++){
    values[i] = values1[i] - values2[i];
  }

  ofstream fp("difference.dat");

  fp << N1 << endl;
  for(int i = 0; i < N1; i++){
    fp << i << "\t" << values[i] << endl;
  }
  fp.close();

   
  return 0;
}

