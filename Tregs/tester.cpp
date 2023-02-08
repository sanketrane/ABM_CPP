#include <iostream>
#include "custom_functions.cpp"

using namespace std;


int main () {
  int i=0, j=2, k;
  double l=0.0, m=0.1, n;
  k=runif<double>(i,j);
  n=runif<double>(l,m);
  cout << k << endl;
  cout << n << endl;
  return 0;
}