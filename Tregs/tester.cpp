#include <iostream>
//#include "custom_functions.cpp"
//#include <filesystem>
#include <armadillo>

int main () {
  arma::vec M(5, arma::fill::randu);
  
  arma::mat B(5, 5, arma::fill::randu);
  arma::mat C = B.t() * B;
  
  arma::mat X = mvnrnd(M, C, 5);


  //cout << k << endl;

  //double xd, yd, n;
  //cout << "Enter Mean and SD: " << endl;
  //cin >> xd >> yd;
  //n=randnorm<double>(xd,yd);
  //cout << "Your sample from mean " << xd << " and sd " << yd << " is: " << n << endl;
//
  //cout << "Current path is " << std::filesystem::current_path() << '\n'; 
  return 0;
}