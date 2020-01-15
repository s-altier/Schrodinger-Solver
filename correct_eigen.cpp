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

int main(int argc, char *argv[]) {
  int numdata, numparams;
  int i=0, j=0;
  
  string fname;
  if(argc<2){
  	cout << "Usage: "<<argv[0]<<" filename"<<endl;
	cout << "filename is the file to be corrected, using the smooth error function in smooth_errors_501.dat" << endl;
  	return -1;
  }
  else{
    numparams = 4;
    fname = argv[1]; 
  }

  string line;
  string smootherr = "smooth_errors_501.dat";
  ifstream myfile (smootherr.c_str()); //input filestream

  vector <double> a(numparams, 0);
  vector <double> aerrors(numparams, 0);
  
  i = 0;
  
  vector < double > x;
  vector < double > y;
  vector < double > sigma;
      
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
      iss >> x[j] >> y[j];
    }
    i++;
  }
  myfile.close();

  cout << "Read in data successfully..." << endl;


  // read in data to be corrected
  myfile.open(fname.c_str());
  int size;

  string line2;

  
  getline(myfile, line2);
  istringstream iss2(line2);
  iss2 >> size;
  

  size = 26; // degeneracy occurs after around this point for higher order potentials
  vector<double> values(size, 0.0);

  double temp;
  int k = 0;
  
  while(getline(myfile, line2) && k < size){

    
    istringstream iss3(line2);

    iss3 >> temp >> values[k];

    k++;
  
  }
  myfile.close();

  
  
  // data correction:
  for(k = 0; k < size; k++){
    values[k] += y[k];
  }
  

  // output data to corrected file
  string efile = "corrected_" + fname;
  
  ofstream fp(efile.c_str());

  fp << size << endl;
  for(k = 0; k < size; k++){
    fp << k << "\t" << values[k] << endl;
  }
  fp.close();

    
  return 0;
}

