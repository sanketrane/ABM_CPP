#include <iostream>
#include "custom_functions.cpp"

// spline1 --
float sp_numbers(float Time, float params[]){ 
    // these are true numbers - need to scale down for simulation (done elsewhere)
    float sp_numbers, theta0=params[0], nu=params[1], psi=params[2];

    // thymic FoxP3 negative SP4 T cell numbers - spline defs
    sp_numbers = psi * pow(10., theta0) * exp(-1. * nu * Time); 
    return sp_numbers;
}

void mhsampler(
  arma::vec proposal_width, arma::vec initial_guess,
  arma::vec prior_mu, arma::vec prior_sd
  )
{
  int k = prior_mu.size();
  arma::arma_rng::set_seed_random();                        // set random seed
  arma::vec mu_current = initial_guess;                     // initial guess for the current mu
  
  arma::mat sigma_diag = arma::diagmat(proposal_width);     
  arma::mat Sigma_proposal = sigma_diag.t() * sigma_diag;                     // covariance matrix -- diagonal
  
  arma::vec mu_proposal = mvnrnd(mu_current, Sigma_proposal);         // proposal distribution -- MVN 
  std::cout << "mu_proposal: " << mu_proposal << '\n';                                   // New sampled mu

  arma::rowvec Mu_prior = arma::conv_to<arma::rowvec>::from(prior_mu);
  arma::mat sd = arma::diagmat(prior_sd);     
  arma::mat SD = sd.t() * sd;                     // covariance matrix -- diagonal
  
  // prior of the current mu
  arma::rowvec Mu_C(1, k), Mu_P(1, k);
  Mu_C = arma::conv_to<arma::rowvec>::from(mu_current); 
  Mu_P = arma::conv_to<arma::rowvec>::from(mu_proposal); 
  arma::vec p_current = dmvnrm_arma(Mu_C, Mu_prior, SD, true);  
  arma::vec p_proposal = dmvnrm_arma(Mu_P, Mu_prior, SD, true);  // New sampled mu

  std::cout << "p_current: " << p_current << "p_proposal: " << p_proposal << '\n';  
  
}

int main (int argc, char * const argv[]) {
  
  int i;
  float sp_cal[10], parms[3] = {6.4, 0.0024, 0.3};
  arma::vec init_guess = {4, 0.3, 0.2};
  arma::vec proposal_width = {0.5, 0.1, 0.1};
  arma::vec prior_mu = {6, 0.5, 0.3};
  arma::vec prior_sd = {1, 0.1, 0.02};

  std::vector<int> Time_pred(100); 
  float startNum=40, step=1;
  std::generate(Time_pred.begin(), Time_pred.end(), [&startNum, &step]{ return startNum+=step;});


  //for (i=0; i<10; i++){
  //  sp_cal[i] = sp_numbers(Time_pred[i], parms);
  //  std::cout << "Sp_cal " << sp_cal[i] << '\n';
  //}

  //read data
  std::vector<std::vector<std::string>> datatofit;
	std::vector<std::string> row;
	std::string line, word;
 
	std::ifstream Data_FILE(fulldatfile_path);
	if(Data_FILE.is_open())
	{
		while(getline(Data_FILE, line))
		{
			row.clear();
 
			std::stringstream str(line);
 
			while(getline(str, word, ','))
				row.push_back(word);
			datatofit.push_back(row);
		}
	}
	else
		std::cout<<"Could not open the file\n";
 
	for(int i=0;i<datatofit.size();i++)
	{
		for(int j=0;j<datatofit[i].size();j++)
		{
			std::cout<<datatofit[i][j]<<" ";
		}
		std::cout<<"\n";
	}
 


  arma::vec A1(5), A2(5), A3(5);

  A1.imbue( [&]() { return randnorm<float>(5, 0.5); } );
  A2.imbue( [&]() { return randnorm<float>(0.6, 0.3); } );
  A3.imbue( [&]() { return randnorm<float>(0.3, 0.1); } );

  arma::mat A = join_rows(A1, A2, A3);

  for(i=0; i<5; i++){
    arma::vec myval = arma::conv_to<arma::vec>::from(A.row(i));
    mhsampler(myval, proposal_width, prior_mu, prior_sd);
  }
  
   
  return 0;
}