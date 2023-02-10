#include <iostream>
#include "custom_functions.cpp"
#include <filesystem>

using namespace std;


int main () {
  int x=4, y=1, k;
  k=randnorm<double>(x,y);
  cout << k << endl;

  double xd, yd, n;
  cout << "Enter Mean and SD: " << endl;
  cin >> xd >> yd;
  n=randnorm<double>(xd,yd);
  cout << "Your sample from mean " << xd << " and sd " << yd << " is: " << n << endl;

  cout << "Current path is " << std::filesystem::current_path() << '\n'; 
  return 0;
}