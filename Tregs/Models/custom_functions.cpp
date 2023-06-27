#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <armadillo>
#include <filesystem>
#include <cmath>
#include <string>


// generatting samples from a uniform distribution with between range_from and range_to
template<typename T>
T randunif(T range_from, T range_to) {
    std::random_device                  rand_dev;                       // stock random number libraries
    std::mt19937                        generator(rand_dev());          // use the mersenne twister as the underlying RN generator
    std::uniform_real_distribution<T>   distr(range_from, range_to);
    return distr(generator);
}

// generatting samples from a univariate normal distribution with mean=mean and sd=sd
template<typename T>
T randnorm(T mean, T sd) {
    std::random_device                  rand_dev;                       // stock random number libraries
    std::mt19937                        generator(rand_dev());          // use the mersenne twister as the underlying RN generator
    std::normal_distribution<T>         distr(mean, sd);
    return distr(generator);
}

// generatting samples from a univariate normal distribution with mean=mean and sd=sd
int randbinom(int nTrials , float prob) {
    std::random_device                  rand_dev;                       // stock random number libraries
    std::mt19937                        generator(rand_dev());          // use the mersenne twister as the underlying RN generator
    std::binomial_distribution<int>         distr(nTrials, prob);
    return distr(generator);
}

// PDF for the univariate normal distribution
template <typename T>
T normal_pdf(T x, T m, T s)
{
    static const T inv_sqrt_2pi = 0.3989422804014327;
    T a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-T(0.5) * a * a);
}

// PDF for the multivariate normal distribution
arma::vec dmvnrm_arma(arma::mat const &x,  
                      arma::rowvec const &mean,  
                      arma::mat const &sigma, 
                      bool const logd = false) { 
    static double const log2pi = std::log(2.0 * M_PI);
    using arma::uword;
    uword const n = x.n_rows, 
             xdim = x.n_cols;
    arma::vec out(n);
    arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
    double const rootisum = arma::sum(log(rooti.diag())), 
                constants = -(double)xdim/2.0 * log2pi, 
              other_terms = rootisum + constants;
    
    arma::rowvec z;
    for (uword i = 0; i < n; i++) {
        z      = (x.row(i) - mean) * rooti;    
        out(i) = other_terms - 0.5 * arma::dot(z, z);     
    }  
      
    if (logd)
      return out;
    return exp(out);
}

// printing elements of a vector
template <typename S>
void print_vector(const std::vector<S>& vector,
                 std::string sep = " ")
{
    // Iterating vector by using index
    for (int i = 0; i < vector.size(); i++) {
        // Printing the element at
        // index 'i' of vector
        std::cout << vector[i] << sep;
    }
    std::cout << '\n';
}

