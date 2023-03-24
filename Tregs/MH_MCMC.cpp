#include <iostream>
//#include "custom_functions.cpp"
//#include <filesystem>
#include <armadillo>

// spline1 --
float sp_numbers(float Time, float params[]){ 
    // these are true numbers - need to scale down for simulation (done elsewhere)
    float sp_numbers, theta0=params[0], nu=params[1], psi=params[2];

    // thymic FoxP3 negative SP4 T cell numbers - spline defs
    sp_numbers = psi * pow(10., theta0) * exp(-1. * nu * Time); 
    return sp_numbers;
}




void mhsampler (arma::vec initial_guess)
{
  // set random seed
  arma::arma_rng::set_seed_random();
  arma::vec mu_current = initial_guess;                     // initial guess for the current mu

  arma::vec sigma_proposal = {1, 0.5};                      // fixed width (sigmas) for the proposal dstribution
  arma::mat B = arma::diagmat(sigma_proposal);     
  arma::mat SIGMA_proposal = B.t() * B;                     // covariance matrix -- diagonal
  arma::vec X = mvnrnd(mu_current, SIGMA_proposal);         // proposal distribution -- MVN 
  std::cout << X << '\n';                                   // New sampled mu
}

int main (int argc, char * const argv[]) {
  
  int i;
  float sp_cal[10], parms[3] = {6.4, 0.0024, 0.3};
  arma::vec mu_vec = {4, 0.3};

  std::vector<int> Time_pred(100); 
  float startNum=40, step=1;
  std::generate(Time_pred.begin(), Time_pred.end(), [&startNum, &step]{ return startNum+=step;});


  //for (i=0; i<10; i++){
  //  sp_cal[i] = sp_numbers(Time_pred[i], parms);
  //  std::cout << "Sp_cal " << sp_cal[i] << '\n';
  //}

  mhsampler(mu_vec);
  
  return 0;
}