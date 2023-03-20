#include <iostream>
//#include "custom_functions.cpp"
//#include <filesystem>
#include <armadillo>

using namespace std;

//int main () {
//  
//  // set random seed
//  arma::arma_rng::set_seed_random();
//  arma::vec M = {0.3, 0.5};
//  
//  //arma::mat B(5, 5, arma::fill::randu);
//  //arma::mat C = B.t() * B;
//
//  arma::vec b = {0.1, 0.1};
//  arma::mat B = arma::diagmat(b);
//  arma::mat C = B.t() * B;
//
//  std::cout << C << std::endl;
//  
//  arma::vec X = mvnrnd(M, C);
//
//  std::cout << X << std::endl;
//  return 0;
//}

// spline1 --
float sp_numbers(int Time, float params[]){ 
    // these are true numbers - need to scale down for simulation (done elsewhere)
    float sp_numbers, theta0=params[1], nu=params[2], psi=params[3];

    // thymic FoxP3 negative SP4 T cell numbers - spline defs
    sp_numbers = psi * pow(10., theta0) * exp(-1. * nu * Time); 
    return sp_numbers;
}

int main () {

  float parms[3] = {6.4, 0.0024, 0.3};
  
  cout << "pars " << parms << endl;
  //std::vector<int> Time_pred(10);
  //std::iota (std::begin(Time_pred), std::end(Time_pred), 0); // Fill with 0, 1, ..., 99.
  float sp_cal = sp_numbers(10, parms);
  float sp_num = 0.3 * pow(10., 6.4) * exp(-1. * 0.0024 * 10); 
   
  std::cout << "sp_cal = "<< sp_cal << " sp_num = " << sp_num << std::endl;
  return 0;

}