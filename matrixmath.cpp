#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>

using namespace std;

void print(vector< vector<double> > A) {
  int n1 = A.size();
  int n2 = A[0].size();
  for (int i=0; i<n1; i++) {
    for (int j=0; j<n2; j++) {
      cout << A[i][j] << "\t";
    }
    cout << "\n";
  }
  cout << endl;
}

void print(vector< double > x) {
  int n1 = x.size();
  for (int i=0; i<n1; i++) {
    cout << x[i];
    cout << "\n";
  }
  cout << endl;
}

double det(vector< vector<double> > A) {
  int n = A.size();
  int piv = 1; // to determine pivots
  
  for (int i=0; i<n; i++) {
    // Search for maximum in this column
    double maxl = abs(A[i][i]);
    int maxr = i;
    for (int k=i+1; k<n; k++) {
      if (abs(A[k][i]) > maxl) {
	maxl = abs(A[k][i]);
	maxr = k;
      }
    }
	
    // Swap maximum row with current row (column by column)
    for (int k=i; k<n+1;k++) {
      double tmp = A[maxr][k];
      A[maxr][k] = A[i][k];
      A[i][k] = tmp;
      piv *= -1;
    }
	
    // Make all rows below this one 0 in current column
    for (int k=i+1; k<n; k++) {
      double c = -A[k][i]/A[i][i];
      for (int j=i; j<n+1; j++) {
	if (i==j) {
	  A[k][j] = 0;
	} else {
	  A[k][j] += c * A[i][j];
	}
      }
    }
  }
  
  // evaluate product of diag elements.
  double x = 1.0;
  for (int i=0; i<n; i++) {
    x*= A[i][i];
  }
  return piv*x;
}


vector < vector <double> > transpose(vector < vector <double> > A){
  int n1 = A.size();
  int n2 = A[0].size();
  // NOTE: the row/column number values were flipped
  vector<double> ll(n1,0);
  vector< vector<double> > AT(n2,ll);
  for (int i=0; i<n1; i++) {
    for (int j=0; j<n2; j++) {
      AT[j][i] = A[i][j];
    }
  }
  return AT;
}

vector < vector < double > > matscale(vector < vector<double> > X, double c){

  vector < vector<double> > Y(X.size(), vector<double>(X[0].size(), 0));

  int nrow = X.size();
  int ncol = X[0].size();

  for(int i = 0; i < nrow; i++){
    for(int j = 0; j < ncol; j++){
      Y[i][j] = c * X[i][j];
    }
  }
  
  return Y;
}

int matprod(vector< vector<double> > A,
	    vector< vector<double> > B,
	    vector< vector<double> >& C) {
  int na1 = A.size();
  int na2 = A[0].size();
  int nb1 = B.size();
  int nb2 = B[0].size();

  if(na2!=nb1){
    cout<<"Matrices do not match!"<<endl;
    return -1;
  }
  
  for (int i=0; i<na1; i++) {
    for (int k=0; k<nb2; k++) {
      double sum =0.0;
      for (int j=0; j<na2; j++) {
	sum += A[i][j]*B[j][k];
      }
      C[i][k] = sum;
    }
  }
  return 0;

}

int matprod(vector< vector<double> > A,
	    vector<double> B,
	    vector<double> & C) {
  int na1 = A.size();
  int na2 = A[0].size();
  int nb1 = B.size();

  if(na2!=nb1){
    cout<<"Matrices do not match!"<<endl;
    return -1;
  }
  
  for (int i=0; i<na1; i++) {
    
    double sum =0.0;
    for (int j=0; j<na2; j++) {
      sum += A[i][j]*B[j];
    }
    C[i] = sum;
    
  }
  return 0;

}

int matadd(vector< vector<double> > A,
	    vector< vector<double> > B,
	    vector< vector<double> >& C) {
  int na1 = A.size();
  int na2 = A[0].size();
  int nb1 = B.size();
  int nb2 = B[0].size();

  if(na1 != nb1 || na2 != nb2){
    cout<<"Matrices do not match!"<<endl;
    return -1;
  }
  
  for (int i=0; i<na1; i++) {
    for (int j=0; j<na2; j++) {
      C[i][j] = A[i][j] + B[i][j];
    }
  }
  return 0;

}

vector<double> gauss(vector< vector<double> > A) {
  int n = A.size();
  
  for (int i=0; i<n; i++) {
    // Search for maximum in this column
    double maxl = abs(A[i][i]);
    int maxr = i;
    for (int k=i+1; k<n; k++) {
      if (abs(A[k][i]) > maxl) {
	maxl = abs(A[k][i]);
	maxr = k;
      }
    }
	
    // Swap maximum row with current row (column by column)
    for (int k=i; k<n+1;k++) {
      double tmp = A[maxr][k];
      A[maxr][k] = A[i][k];
      A[i][k] = tmp;
    }
	
    // Make all rows below this one 0 in current column
    for (int k=i+1; k<n; k++) {
      double c = -A[k][i]/A[i][i];
      for (int j=i; j<n+1; j++) {
	if (i==j) {
	  A[k][j] = 0;
	} else {
	  A[k][j] += c * A[i][j];
	}
      }
    }
  }
  
  // Solve equation Ax=b for an upper triangular matrix A
  vector<double> x(n);
  for (int i=n-1; i>=0; i--) {
    x[i] = A[i][n]/A[i][i];
    for (int k=i-1;k>=0; k--) {
      A[k][n] -= A[k][i] * x[i];
    }
  }
  return x;
}

double inv(vector< vector<double> > A, vector< vector<double> >& Ainv) {
  int i,j,k;
  int n = A.size();
  int pivotsign = 1;
  
  int *index; // to mark pivot
  index = new int[n];
  
  vector<double> scale(n,0); //for scale factors
  Ainv = A; // Copy A to make sure Ainv is the right size.

  vector<double> ll(n,0);
  vector< vector<double> > b(n,ll);
  // Set diag elements to 1.0
  for(i = 0 ; i < n ; i++){
    b[i][i] = 1.0;
  }
  
  for (int i=0; i<n; i++) {
    index[i] = i; //initialize index list
    double scalemax = 0.0;
    for(j=0;j<n;j++){
      if(scalemax > abs(A[i][j]) )
	scalemax = scalemax;
      else
	scalemax = abs(A[i][j]);
    }
    scale[i] = scalemax;
  }  
  
  for (int i=0; i<n-1; i++) {
    // Select pivot row
    double ratiomax = 0.0;
    int jpivot = i;
    for(k=i;k<n;k++){
      double ratio = abs(A[index[k]][i]/scale[index[k]]);
      if(ratio>ratiomax){
	jpivot = k;
	ratiomax = ratio;
      }
    }
    // Perform pivoting using row index list
    int indexj = index[i];
    if(jpivot!=i){//PIVOT!
      indexj = index[jpivot];
      index[jpivot] = index[i];//Swap index jpivot and k
      index[i] = indexj;
      pivotsign *= -1;//flip sign of det
    }
    // Forward elimination
    for(k=i+1;k<n;k++){
      double c = A[index[k]][i]/A[indexj][i];
      for(j=i+1;j<n;j++){
	A[index[k]][j] -= c*A[indexj][j];
      }
      A[index[k]][i] = c;
      for(j=0;j<n;j++){
	b[index[k]][j] -= A[index[k]][i]*b[indexj][j];
      }
    }
  }
  double det = pivotsign;
  for(i=0;i<n;i++){
    det *= A[index[i]][i];
  }

  // Backsubstitution
  for(k=0;k<n;k++){
    Ainv[n-1][k] = b[index[n-1]][k]/A[index[n-1]][n-1];
    for(i=n-2;i>=0;i--){
      double sum = b[index[i]][k];
      for(j=i+1;j<n;j++){
	sum -= A[index[i]][j]*Ainv[j][k];
      }
      Ainv[i][k] = sum/A[index[i]][i];
    }
  }
  delete [] index;
  return det;
}

double max_element(vector< vector<double> > V){
  // this might be a problem, as we saw in class with using .size()
  int n = V.size();
  int m = V[0].size();

  double max = 0.0;

  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++){
      if(abs(V[i][j]) > max){
        max = abs(V[i][j]);
      }
    }
  }

  return max;
}

double tridiaginv(vector< vector<double> > T,
        vector< vector<double> >& Tinv){
        /*

        Here in real matrix form:
        T[i,i] = a_{i+1}
        T[i,i+1] = c_{i+1}
        T[i+1,i] = b_{i+2}

        T is actually going to be a 3 x n matrix:
        T[0][i] = a_{i+1}
        T[1][i] = c_{i+1}
        T[2][i] = b_{i+2}

        the elements that we don't want will just be
        zero and ignored like the 0th components
        of the vectors.
        */

        int n = T[0].size();
        vector < double > a(n+1,0);
        vector < double > b(n+1,0);
        vector < double > c(n+1,0);

        for(int i = 1 ; i < n+1 ; i++){
                a[i] =   T[0][i-1];
                c[i] =   T[2][i-1];
                if(i<n)
                  b[i+1] = T[1][i-1];
        }
        c[n] = 1.0;
        vector < double > A(n+1,0);

        A[0] = 1.0;
        A[1] = -a[1]*A[0]/c[1];
        for(int i=1; i<=n-1;i++){
                A[i+1] = -( a[i+1]*A[i] + b[i+1]*A[i-1] )/c[i+1];
        }

        vector < double > lll(n+1,0);
        vector < vector < double > > C(n+1,lll);
        vector < vector < double > > E(n+1,lll);
        for(int i=1; i<n+1;i++){
                E[i][i] = 1.0;
        }
        for(int i=0; i<n;i++){
                C[n][i+1] = -A[i]/A[n];
                C[n-1][i+1] = (E[n][i+1] - a[n]*C[n][i+1]  ) / c[n-1];
                for(int j = n-1 ; j >= 2 ; j--){
                        C[j-1][i+1] = (E[j][i+1]
                                - a[j]*C[j][i+1]
                                - b[j+1]*C[j+1][i+1]
                                ) / c[j-1] ;
                }
        }
        for(int i=1; i <= n ; i++){
                for(int j=1; j <= n ; j++){
                        Tinv[i-1][j-1] = C[i][j];
                }
        }
        double det=A[n];
        for(int k = 1; k <=n-2;k++)
                det *= -c[k];
        return det;
}

double aveg(vector<double> sample){
  int n = sample.size();

  double sum = 0.0;

  for(int i = 0; i < n; i++){
    sum += sample[i];
  }

  return sum / ((double) n);
}

double stdev(vector<double> sample){
  int n = sample.size();

  double avg = aveg(sample);

  double sum = 0.0;

  for(int i = 0; i < n; i++){
    sum += pow(avg - sample[i], 2.0);
  }

  sum /= ((double) (n - 1));
  return sqrt(sum);

}

double minvalue(vector<double> sample){
  double min = sample[0];
  int n = sample.size();

  for(int i = 0; i < n; i++){
    if(sample[i] < min)
      min = sample[i];
  }

  return min;
}

double maxvalue(vector<double> sample){
  double max = sample[0];
  int n = sample.size();

  for(int i = 0; i < n; i++){
    if(sample[i] > max)
      max = sample[i];
  }

  return max;
}

