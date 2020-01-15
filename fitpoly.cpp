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
  
  string fname;
  if(argc<3){
  	cout << "Usage: "<<argv[0]<<" M filename"<<endl;
	cout << "M: number of fit parameters for an (M-1)th-order polynomial fit." << endl; 
	cout << "The file should be in the form: " << endl;
	cout << "numdata\nx y err\n x y err\n ..." << endl;
  	return -1;
  }
  else{
    numparams = atoi(argv[1]);
    fname = argv[2]; 
  }

  string line;
  ifstream myfile (fname.c_str()); //input filestream

  vector <double> a(numparams, 0);
  vector <double> aerrors(numparams, 0);
  
  i = 0;
  
  vector < double > x(numdata,0);
  vector < double > y(numdata,0); 
  vector < double > sigma(numdata,1); 
      
  while ( getline(myfile,line) ){
    istringstream iss(line); 
    //    double sub;
    
    if(i==0){
      iss >> numdata;
      //numdata = (int) sub; 
      x.reserve(numdata);
      y.reserve(numdata);
      sigma.reserve(numdata);
    } else {
      j = i-1;
      iss >> x[j] >> y[j] >> sigma[j];
    }
    i++;
  }
  myfile.close();

  cout << "Read in data successfully..." << endl;

  vector < double > rows(numparams, 0); //M-dim vector
  vector < vector <double> > A(numdata, rows); //NxM matrix
  vector < double > b(numdata,0);
  vector < vector < double > > AT(numparams,b); //MxN matrix
  vector < vector < double > > C(numparams, rows); // MxM matrix
  vector < vector < double > > ATA(numparams, rows); // MxM matrix
  vector < double > ATb(numparams, 0);

  for(i = 0; i < numdata; i++){
    for(j = 0; j < numparams; j++){
      A[i][j] = pow(x[i], j)/sigma[i];
      AT[j][i] = pow(x[i], j)/sigma[i];
    }
    b[i] = y[i]/sigma[i];
  }
  //  AT = transpose(A);

  int check = matprod(AT, A, ATA);
  if(check == -1){
    cerr << "There's a problem with multiplication!" << endl;
    //cerr: acts like cout here; can have cerr go to a different file, like error log, rather than outside
  }

  double detA = inv(ATA, C);
  if(detA < 1e-10){
    cout << "The determinant is " << detA << ", be careful!" << endl;
  }

  check = matprod(AT, b, ATb);
  if(check ==-1){
    cerr << "There's a problem with multiplication!" << endl;
  }

  check = matprod(C, ATb, a);
  if(check==-1){
    cerr << "There's a problem with multiplication!" << endl;
  }

  for(j = 0; j < numparams; j++){
    aerrors[j] = sqrt(C[j][j]);
  }
  
  cout << "Our fit results are: " << endl;
  for(j = 0; j < numparams; j++){
    cout << "a[ " << j << "] = " << a[j] << " +/- " << aerrors[j] << endl;
  }    
  cout << "chi^2 / dof = " << chi2(numdata,a,x,y,sigma) << " / " << numdata - numparams << endl;
  
  //Plotting time


  // generate polystring
  string temp;
  string polystring = "f(x) = ";
  for(int i = numparams - 1; i >= 0; i--){
    temp = to_string(a[i]);
    polystring += temp;
    temp = "*x^" + to_string(i);
    polystring += temp;
    if(i!=0)
      polystring += "+";
  }
  cout << polystring << endl;
  //  double fitfunc = 0.0;
  double t = 0.0, dt = 0.01; 
  
  t = x[0];
  double tmax = x[numdata-1];

  size_t found = fname.find(".dat");
  string outplot = fname;
  outplot.replace(found, 4, ".ps");
  string fitfilename = "fit_poly_" + to_string(numparams-1) + "_" + fname;
  
ofstream fout(fitfilename.c_str());

  while(t <= tmax){
    fout << t << "\t" << poly(a,t) << endl;
    t+=dt;
  }

  fout.close();


  FILE *gnuplot = popen("gnuplot", "w");
  fprintf(gnuplot, "set out '%s'\n", outplot.c_str()); //%s == char string
  fprintf(gnuplot, "set term post land\n");

  string pltcmd = "plot '" + fname + "' u 1:2, " + "'" + fitfilename + "' u 1:2 w l title 'f(x)'\n";
  fprintf(gnuplot,"%s", pltcmd.c_str());
  pclose(gnuplot);

  system( (char *)( "ps2pdf "+outplot).c_str());
  
  return 0;
}

