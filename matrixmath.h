#ifndef MATRICES_H
#define MATRICES_H

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
using namespace std;


void print(vector< double > A);

void print(vector< vector<double> > A);

double det(vector< vector<double> > A);

vector < vector < double > > matscale(vector < vector<double> > X, double c);

int matprod(vector< vector<double> > A,
	    vector< vector<double> > B,
	    vector< vector<double> >& C);

int matprod(vector< vector<double> > A,
            vector<double> B,
            vector<double> & C);

int matadd(vector< vector<double> > A,
	   vector< vector<double> > B,
	   vector< vector<double> >& C);

vector<double> gauss(vector< vector<double> > A);

double inv(vector< vector<double> > A,vector< vector<double> >& Ainv) ;

vector < vector <double> > transpose(vector < vector <double> > A);

double max_element(vector< vector<double> > V);

double tridiaginv(vector< vector<double> > T,
		  vector< vector<double> >& Tinv);

double aveg(vector<double> sample);

double stdev(vector<double> sample);

double minvalue(vector<double> sample);

double maxvalue(vector<double> sample);

#endif
