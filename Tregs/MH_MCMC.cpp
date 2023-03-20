#include <iostream>
//#include "custom_functions.cpp"
//#include <filesystem>
#include <armadillo>

// spline1 --
float sp_numbers(float Time, float params[]){ 
    // these are true numbers - need to scale down for simulation (done elsewhere)
    float sp_numbers, theta0=params[1], nu=params[2], psi=params[3];

    // thymic FoxP3 negative SP4 T cell numbers - spline defs
    sp_numbers = psi * pow(10., theta0) * exp(-1. * nu * Time); 
    return sp_numbers;
}

void mhsampler (float initial_guess[], )
{
  // set random seed
  arma::arma_rng::set_seed_random();
  arma::vec M = initial_guess;

  arma::vec b = {0.2, 0.5};
  arma::mat B = arma::diagmat(b);
  arma::mat C = B.t() * B;

  std::cout << C << std::endl;
  
  arma::vec X = mvnrnd(M, C);

  std::cout << X << std::endl;

  return 0;
}

int main (int argc, char * const argv[]) {


    std::vector<float> pars = {6.4, 0.0024, 0.3};
    std::vector<int> Time_pred(10);
    std::iota (std::begin(Time_pred), std::end(Time_pred), 0); // Fill with 0, 1, ..., 99.

    float sp_cal = sp_numbers(Time_pred, pars);

    std::cout << "sp_num"<< sp_cal << std::endl;

    return 0;

}