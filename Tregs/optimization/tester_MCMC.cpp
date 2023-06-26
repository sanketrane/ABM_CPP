#include <iostream>
#include "../optimization/MH_MCMC.cpp"

// model --
float sp_numbers(int Time, arma::vec params){ 
    // these are true numbers - need to scale down for simulation (done elsewhere)
    float tot_numbers, theta0=params[0], nu=exp(params[1]), b0=params[2];
    // thymic FoxP3 negative SP4 T cell numbers - spline defs
    tot_numbers = pow(10., b0) + pow(10., theta0) * exp(-1. * nu * Time);
    
    return tot_numbers;
}

// Simulator function 
std::vector<float> model_output (arma::vec params){
  int i, tStart=40, tEnd = 120;
  std::vector<int> Time_pred(100); std::vector<float> sp_cal;
  std::generate(Time_pred.begin(), Time_pred.end(), [&startNum, &step]{ return startNum+=step;});
  // simuations
  for (i=tStart; i<tEnd; i++){
    sp_cal.push_back(log(sp_numbers(Time_pred[i], params)));
  }
  
  return sp_cal;
}

// Log likelihood (LL) -- to calculate LL of model predictions conditional on observed data
float llfunction(arma::vec params, std::vector<float> y_obs){
    int i, size_dat = y_obs.size();
    float  sd = params.at(3), ll_ind[size_dat], ll_est;
    std::vector<float> y_pred = model_output(params);       // model simulation

    for (i=0; i<size_dat; i++){
      ll_ind[i] = normal_pdf(y_pred[i], y_obs[i], sd);
      ll_est += log(ll_ind[i]);
    }
    return ll_est;
}

// Metropolis-Hastings sampler to generate posterior distribution of model parameters
void mhsampler(
  arma::vec proposal_width, arma::vec initial_guess,                            // Hyper parameters
  arma::vec prior_mu, arma::vec prior_sd,                                       // priors defined for each model
  std::vector<float> y_obs,                                                     // data as a vector
  int Num_samples
  )             
{             
  int k = prior_mu.size(),  s=0;          
  arma::arma_rng::set_seed_random();                                           // set random seed
  arma::vec current_mu = initial_guess;                                        // initial guess is used to set the current mu for sampler initialization
               
  arma::mat proposal_diag = arma::diagmat(proposal_width);                     // define a diagonal matrix for sigma of mvn proposal distribution
  arma::mat Sigma_proposal = proposal_diag.t() * proposal_diag;                // covariance matrix -- diagona
  
  arma::mat Posterior(1, k);
  Posterior = arma::conv_to<arma::rowvec>::from(current_mu);
  std::cout << " POST " << Posterior << '\n';

   while (s < Num_samples)
   {
      // propose a new mu using proposal distribution
      arma::vec new_sampled_mu = mvnrnd(current_mu, Sigma_proposal);                 // proposal distribution -- MVN 
      //std::cout << "new_sampled_mu: " << new_sampled_mu << '\n';                   // New sampled mu
    
      // definition of the prior distribution
      arma::rowvec Mu_prior = arma::conv_to<arma::rowvec>::from(prior_mu);          // conversion to a row vector to feed into dmvnorm function
      arma::mat sd_diag = arma::diagmat(prior_sd);     
      arma::mat SD_prior = sd_diag.t() * sd_diag;                                   // covariance matrix -- diagonal
      
      // defining prior probabilities of mu_current and mu_proposed
      arma::rowvec Current_Mu(1, k), New_Mu(1, k), Sig_C(1, k), Sig_P(1, k);         // defining vectors with k elements for old and new mu
      Current_Mu = arma::conv_to<arma::rowvec>::from(current_mu);                    // conversion to a row vector to feed into dmvnorm function
      New_Mu = arma::conv_to<arma::rowvec>::from(new_sampled_mu);                 
      arma::vec PriorProb_current = dmvnrm_arma(Current_Mu, Mu_prior, SD_prior, true);  // probability of mu_current based of prior distribution
      arma::vec PriorProb_new = dmvnrm_arma(New_Mu, Mu_prior, SD_prior, true);          // probability of mu_proposed based of prior distribution
    
      float  LL_current = llfunction(current_mu, y_obs);                       // LL of current mu conditional on y_obs
      float  LL_new = llfunction(new_sampled_mu, y_obs);                       // LL of new sampled mu conditional on y_obs
    
      // std::cout << "prior_current: " << PriorProb_current << "prior_proposal: " << PriorProb_new << '\n';  
      // std::cout << "ll_current: " << LL_current << " ll_proposal: " << LL_new << '\n';  
    
      float PostProb_current = LL_current + PriorProb_current[0];
      float PostProb_new = LL_new + PriorProb_new[0];
    
      // std::cout << "PostProb_current: " << PostProb_current << '\n';                   // New sampled mu
      // std::cout << "PostProb_new: " << PostProb_new << '\n';                   // New sampled mu
    
    
      float Acceptance = PostProb_new - PostProb_current;
      float Random_threshold = log(randunif<float>(0, 1));
      //std::cout << "Acceptance: " << Acceptance << '\n';                   // New sampled mu
      //std::cout << "Random_threshold: " << Random_threshold << '\n';                   // New sampled mu
    
    
      //std::cout << "old_mu: " << current_mu << '\n';                   // New sampled mu
      if (Random_threshold < Acceptance)
      {
        //update mu
        current_mu = new_sampled_mu;
        s++; 
        Posterior.insert_rows(s, arma::conv_to<arma::rowvec>::from(current_mu));
        std::cout << " Posterior has " << Posterior.n_rows << " rows." << '\n';

      }                  // New accepted mu
    }
    std::cout << Posterior << '\n'; 
}

int main (int argc, char * const argv[]) {

  int i, tStart=40, tEnd = 120, num_obs = 50;
  arma::vec init_guess = {5.0, -0.2, 5.0, 1.5}, proposal_width = {0.5, 0.5, 0.5, 0.5}, prior_mu = {5, -0.1, 5, 1.2}, prior_sd = {1, 0.5, 1, 0.5}, parms = {6.5, -0.08, 4.8, 1.5};
  float sig_sp = parms[3], sp_noise[num_obs];

  // lambda expression to generate random sequence for time point
  auto rand_sequence = [] (int begin, int end, int length) {
    std::vector<int> seqvec;  srand(99); //random seed //srand(time(0)); 
    for(int i=0; i<length; i++) {seqvec.push_back(begin + rand() % (end - begin));} return seqvec;};
    
  std::vector<int> Time_pred = rand_sequence(tStart, tEnd, num_obs); std::sort(Time_pred.begin(), Time_pred.end());
 

  // generating data
  std::vector<float> sp_cal = model_output(parms);
  print_vector(sp_cal);
  for (i=0; i<num_obs; i++){
    sp_noise[i] = randnorm<float>(sp_cal[i], sig_sp);
  }
  std::vector<float> sp_dat(std::begin(sp_noise), std::end(sp_noise));

  //mhsampler(proposal_width, init_guess, prior_mu, prior_sd, sp_dat, 10);


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
  
  std::cout << "DONE" << std::endl;

  return 0;
}