#include "custom_functions.cpp"

// model --
float sp_numbers(int Time, arma::vec params)
{
  int t0 = 40;
  // these are true numbers - need to scale down for simulation (done elsewhere)
  float tot_numbers, theta0 = params[0], nu = params[1];
  // thymic FoxP3 negative SP4 T cell numbers - spline defs
  tot_numbers = pow(10., theta0) * exp(-1. * nu * (Time - t0));
  // return log transformed tot_numbers
  return log(tot_numbers);
}

float ki_frac(int Time, arma::vec params)
{
  float tStep = 0.1, alpha = params[2];
  // thymic FoxP3 negative SP4 T cell numbers - spline defs
  float tot_numbers = exp(sp_numbers(Time, params));
  // division -- cells that are ki67 positive
  double prob_loss = 1 - exp(-alpha * tStep);
  // calculate fraction of dividing cells
  float div_frac = randbinom(tot_numbers, prob_loss) / tot_numbers;
  // return logit transformed fraction
  return log(div_frac / (1 - div_frac));
}

// Simulator function
std::vector<float> model_output(float (*func)(int, arma::vec), arma::vec params, std::vector<int> predictor_variable)
{
  std::vector<float> sim_cal;
  // simuations
  for (int i = 0; i < predictor_variable.size(); i++)
  {
    sim_cal.push_back(func(predictor_variable[i], params));
  }
  return sim_cal;
}

// Log likelihood (LL) -- to calculate LL of model predictions conditional on observed data
template <typename T, typename U>
auto llfunction(arma::vec params, std::vector<T> response_variable, std::vector<U> predictor_variable)
{
  int i, size_dat = response_variable.size();
  float sd1 = params.at(3), sd2 = params.at(4), LL_est = 0.;
  std::vector<float> model_pred = model_output(ki_frac, params, predictor_variable); // model simulation

  for (i = 0; i < size_dat; i++)
  {
    LL_est += log(normal_pdf(model_pred[i], response_variable[i], sd2));
  }
  return LL_est;
}

// Metropolis-Hastings sampler to generate posterior distribution of model parameters
template <typename T, typename U>
void mhsampler(
    arma::vec proposal_width, arma::vec initial_guess,                   // Hyper parameters
    arma::vec prior_mu, arma::vec prior_sd,                              // priors defined for each model
    std::vector<T> response_variable, std::vector<U> predictor_variable, // data as a vector
    int Num_samples)
{
  int k = prior_mu.size(), s = 0;
  arma::arma_rng::set_seed_random();    // set random seed
  arma::vec current_mu = initial_guess; // initial guess is used to set the current mu for sampler initialization

  arma::mat proposal_diag = arma::diagmat(proposal_width);      // define a diagonal matrix for sigma of mvn proposal distribution
  arma::mat Sigma_proposal = proposal_diag.t() * proposal_diag; // covariance matrix -- diagona

  arma::mat Posterior(1, k);
  Posterior = arma::conv_to<arma::rowvec>::from(current_mu);

  while (s < Num_samples)
  {
    // propose a new mu using proposal distribution
    arma::vec new_sampled_mu = mvnrnd(current_mu, Sigma_proposal); // proposal distribution -- MVN
    // std::cout << "new_sampled_mu: " << new_sampled_mu << '\n';                   // New sampled mu

    // definition of the prior distribution
    arma::rowvec Mu_prior = arma::conv_to<arma::rowvec>::from(prior_mu); // conversion to a row vector to feed into dmvnorm function
    arma::mat sd_diag = arma::diagmat(prior_sd);
    arma::mat SD_prior = sd_diag.t() * sd_diag; // covariance matrix -- diagonal

    // defining prior probabilities of mu_current and mu_proposed
    arma::rowvec Current_Mu(1, k), New_Mu(1, k), Sig_C(1, k), Sig_P(1, k); // defining vectors with k elements for old and new mu
    Current_Mu = arma::conv_to<arma::rowvec>::from(current_mu);            // conversion to a row vector to feed into dmvnorm function
    New_Mu = arma::conv_to<arma::rowvec>::from(new_sampled_mu);
    arma::vec PriorProb_current = dmvnrm_arma(Current_Mu, Mu_prior, SD_prior, true); // probability of mu_current based of prior distribution
    arma::vec PriorProb_new = dmvnrm_arma(New_Mu, Mu_prior, SD_prior, true);         // probability of mu_proposed based of prior distribution

    // calculating the Log likelihood of the current and new mu
    float LL_current = llfunction(current_mu, response_variable, predictor_variable); // LL of current mu conditional on y_obs
    float LL_new = llfunction(new_sampled_mu, response_variable, predictor_variable); // LL of new sampled mu conditional on y_obs

    std::cout << "LL_current: " << LL_current << '\n';
    std::cout << "LL_new: " << LL_new << '\n';

    // calculating the posterior probability of the current and new mu
    float PostProb_current = LL_current + PriorProb_current[0];
    float PostProb_new = LL_new + PriorProb_new[0];

    std::cout << "PostProb_current: " << PostProb_current << '\n';
    std::cout << "PostProb_new: " << PostProb_new << '\n';

    // calculating the acceptance probability
    float Acceptance = PostProb_new - PostProb_current;
    float Random_threshold = log(randunif<float>(0, 1));

    // updating the current mu
    if (Random_threshold < Acceptance)
    {
      // update mu
      current_mu = new_sampled_mu;
      s++;
      Posterior.insert_rows(s, arma::conv_to<arma::rowvec>::from(current_mu));

    } // New accepted mu
  }
  std::cout << " POST\n"
            << Posterior << '\n';
}

int main(int argc, char *const argv[])
{

  int i, tStart = 40, tEnd = 120, num_obs = 50;
  arma::vec init_guess = {5.0, 0.2, 5.0, 1.0, 0.25}, proposal_width = {0.5, 0.5, 0.5, 0.5, 0.5}, prior_mu = {5, 0.1, 5, 1.2, 0.5},
            prior_sd = {1, 0.5, 1, 0.5, 0.1}, parms = {5, 0.05, 0.05, 1.5, 0.5};
  float sig_num = parms[3], sig_frac = parms[4];

  // lambda expression to generate random sequence for time point
  auto rand_sequence = [](int begin, int end, int length)
  {
    std::vector<int> seqvec;  srand(99); //random seed //srand(time(0)); 
    for(int i=0; i<length; i++) {seqvec.push_back(begin + rand() % (end - begin));} return seqvec; };

  std::vector<int> Time_obs = rand_sequence(tStart, tEnd, num_obs);
  std::sort(Time_obs.begin(), Time_obs.end());

  // generating ordered sequence
  std::vector<float> Time_pred(10);
  float tStep = (tEnd - tStart) / 10.0, startNum = tStart / 1.0;
  std::generate(Time_pred.begin(), Time_pred.end(), [&startNum, &tStep]
                { return startNum += tStep; });

  // generating data
  std::vector<float> sp_cal = model_output(sp_numbers, parms, Time_obs);
  std::vector<float> sp_frac = model_output(ki_frac, parms, Time_obs);

  std::vector<float> num_dat, frac_dat;
  for (i = 0; i < num_obs; i++)
  {
    // std::cout << "sp num "  << sp_numbers(Time_obs[i], parms) << '\n';
    num_dat.push_back(sp_cal[i] + randnorm<float>(0, sig_num));
    frac_dat.push_back(sp_frac[i] + randnorm<float>(0, sig_frac));
  }
  // print_vector(num_dat);
  print_vector(frac_dat);

  // float llval_counts  = llfunction(sp_numbers, parms, num_dat, Time_obs);
  float llval_frac = llfunction(parms, frac_dat, Time_obs);
  // std::cout << "llval_counts: " << llval_counts << '\n';
  std::cout << "llval_frac: " << llval_frac << '\n';

  // calling mhsampler
  // mhsampler(proposal_width, init_guess, prior_mu, prior_sd, frac_dat, Time_obs, 10);

  std::cout << "mhsampler done\n";

  return 0;
}