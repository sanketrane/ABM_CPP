#include <iostream>
//#include "custom_functions.cpp"
//#include <filesystem>
#include <armadillo>

int main () {
  
  // set random seed
  arma::arma_rng::set_seed_random();
  arma::vec M = {3, 5, 1, 10, 6};
  
  //arma::mat B(5, 5, arma::fill::randu);
  //arma::mat C = B.t() * B;

  arma::vec b = {0.2, 0.5, 1, 1.5, 0.6};
  arma::mat B = arma::diagmat(b);
  arma::mat C = B.t() * B;

  std::cout << C << std::endl;
  
  arma::colvec X = mvnrnd(M, C);

  std::cout << X << std::endl;


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