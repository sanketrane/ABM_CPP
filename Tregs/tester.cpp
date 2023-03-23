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
    float sp_numbers, theta0=params[0], nu=params[1], psi=params[2];

    // thymic FoxP3 negative SP4 T cell numbers - spline defs
    sp_numbers = psi * pow(10., theta0) * exp(-1. * nu * Time); 
    return sp_numbers;
}

int main (int argc, char * const argv[]) {
  int i;
  float sp_cal[10];
  float parms[3] = {6.4, 0.0024, 0.3};

  std::vector<int> Time_pred(10); 
  float startNum=40, step=5;
  std::generate(Time_pred.begin(), Time_pred.end(), [&startNum, &step]{ return startNum+=step;});


  for (i=0; i<10; i++){
    std::cout << "Params " << parms[i] << '\n';
    std::cout << "Time_pred " << Time_pred[i] << '\n';
    sp_cal[i] = sp_numbers(Time_pred[i], parms);
    std::cout << "Sp_cal " << sp_cal[i] << '\n';
  }
  
  return 0;
}