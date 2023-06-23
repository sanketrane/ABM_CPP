#include <iostream>
#include "../Models/custom_functions.cpp"

// spline1 --
float sp_numbers(int Time, float params[]){ 

    // these are true numbers - need to scale down for simulation (done elsewhere)
    float tot_numbers, theta0=params[0], nu=params[1], b0=params[2];
    // thymic FoxP3 negative SP4 T cell numbers - spline defs
    tot_numbers = pow(10., b0) + pow(10., theta0) * exp(-1. * nu * Time);
    return tot_numbers;
}

float llfunction(
  std::vector<int>x, float params[], float sd, std::function<float(float, float[])> func,
  float y_obs[]){

    int i, size_dat = x.size();
    float y_pred[size_dat], ll_ind[size_dat], ll_est;

    for (i=0; i<size_dat; i++){
      y_pred[i] = log(func(x[i], params));
      std::cout << "y_pred " << y_pred[i] << " y_obs " << y_obs[i] << '\n';
      ll_ind[i] = normal_pdf(y_pred[i], y_obs[i], sd);
      std::cout << "ll_ind " << log(ll_ind[i]) << '\n';
      ll_est += log(ll_ind[i]);
    }
    return ll_est;
}

void mhsampler(
  arma::vec proposal_width, arma::vec initial_guess,                     // latent variables 
  arma::vec prior_mu, arma::vec prior_sd                                 // priors defined for each model
  )             
{             
  int k = prior_mu.size();             
  arma::arma_rng::set_seed_random();                                     // set random seed
  arma::vec mu_current = initial_guess;                                  // initial guess is used to set the current mu for sampler initialization
               
  arma::mat sigma_diag = arma::diagmat(proposal_width);                  // define a diagonal matrix for sigma of mvn proposal distribution
  arma::mat Sigma_proposal = sigma_diag.t() * sigma_diag;                // covariance matrix -- diagonal
  
  // propose a new mu using proposal distribution
  arma::vec mu_proposal = mvnrnd(mu_current, Sigma_proposal);            // proposal distribution -- MVN 
  std::cout << "mu_proposal: " << mu_proposal << '\n';                   // New sampled mu

  // definition of the prior distribution
  arma::rowvec Mu_prior = arma::conv_to<arma::rowvec>::from(prior_mu);   // conversion to a row vector to feed into dmvnorm function
  arma::mat sd_diag = arma::diagmat(prior_sd);     
  arma::mat SD_prior = sd_diag.t() * sd_diag;                                            // covariance matrix -- diagonal
  
  // defining prior probabilities of mu_current and mu_proposed
  arma::rowvec Mu_C(1, k), Mu_P(1, k), Sig_C(1, k), Sig_P(1, k);         // defining vectors with k elements for old and new mu
  Mu_C = arma::conv_to<arma::rowvec>::from(mu_current);                  // conversion to a row vector to feed into dmvnorm function
  Mu_P = arma::conv_to<arma::rowvec>::from(mu_proposal);                 
  arma::vec prior_current = dmvnrm_arma(Mu_C, Mu_prior, SD_prior, true);           // probability of mu_current based of prior distribution
  arma::vec prior_proposal = dmvnrm_arma(Mu_P, Mu_prior, SD_prior, true);          // probability of mu_proposed based of prior distribution

  std::cout << "p_current: " << prior_current << "p_proposal: " << prior_proposal << '\n';  
}


int main (int argc, char * const argv[]) {
  
  int i;
  float parms[3] = {6.5, 0.1, 5},  sig_sp = 1.5, sp_cal[10], sp_dat[10];

  arma::vec init_guess = {4, 0.03, 5}, proposal_width = {0.5, 0.5, 0.5}, prior_mu = {5, 0.5, 5}, prior_sd = {1, 0.1, 1};

  std::vector<int> Time_pred(10); 
  float startNum=40, step=1;
  std::generate(Time_pred.begin(), Time_pred.end(), [&startNum, &step]{ return startNum+=step;});

  //print_vector(Time_pred);
  //sp_dat = randnorm<float>(sp_cal, sig_sp);
  
  for (i=0; i<10; i++){
    sp_cal[i] = sp_numbers(Time_pred[i], parms);
    sp_dat[i] = log(randnorm(sp_cal[i], sig_sp));
    std::cout << " sp_dat " << sp_dat[i] << '\n';
  }

  float ll_val = llfunction(Time_pred, parms, sig_sp, &sp_numbers, sp_dat);
  std::cout << "myval " << ll_val << '\n';


//
  ////read data
  //std::vector<std::vector<std::string>> datatofit;
	//std::vector<std::string> row;
	//std::string line, word;
 //
	//std::ifstream Data_FILE(fulldatfile_path);
	//if(Data_FILE.is_open())
	//{
	//	while(getline(Data_FILE, line))
	//	{
	//		row.clear();
 //
	//		std::stringstream str(line);
 //
	//		while(getline(str, word, ','))
	//			row.push_back(word);
	//		datatofit.push_back(row);
	//	}
	//}
	//else
	//	std::cout<<"Could not open the file\n";
 //
	//for(int i=0;i<datatofit.size();i++)
	//{
	//	for(int j=0;j<datatofit[i].size();j++)
	//	{
	//		std::cout<<datatofit[i][j]<<" ";
	//	}
	//	std::cout<<"\n";
	//}
 //
//
//
//  arma::vec A1(5), A2(5), A3(5);
//
//  A1.imbue( [&]() { return randnorm<float>(5, 0.5); } );
//  A2.imbue( [&]() { return randnorm<float>(0.6, 0.3); } );
//  A3.imbue( [&]() { return randnorm<float>(0.3, 0.1); } );
//
//  arma::mat A = join_rows(A1, A2, A3);
//
//  for(i=0; i<5; i++){
//    arma::vec myval = arma::conv_to<arma::vec>::from(A.row(i));
//    mhsampler(myval, proposal_width, prior_mu, prior_sd);
//  }
  
  std::cout << "DONE" << std::endl;

  return 0;
}