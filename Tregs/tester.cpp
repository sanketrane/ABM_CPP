#include <iostream>
#include "custom_functions.cpp"
#include <filesystem>

using namespace std;


int main () {
  int i=0, j=2, k;
  double l=0.0, m=0.1, n;
  k=randunif<double>(i,j);
  n=randunif<double>(l,m);
  cout << k << endl;
  cout << n << endl;

  cout << "Current path is " << std::filesystem::current_path() << '\n'; 
  return 0;
}