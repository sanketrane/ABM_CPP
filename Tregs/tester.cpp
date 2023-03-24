#include <iostream>
//#include "custom_functions.cpp"
//#include <filesystem>
#include <armadillo>

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
float sp_numbers(float Time, float params[]){ 

    float sp_cal, theta0=params[0], nu=params[1], psi=params[2];
    // thymic FoxP3 negative SP4 T cell numbers - spline defs
    sp_cal = psi * pow(10., theta0) * exp(-1. * nu * Time); 
    return sp_cal;
}

float normal_pdf(float x, float m, float s)
{
    static const float inv_sqrt_2pi = 0.3989422804014327;
    float a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}

//std::vector<float> MVN_pdf(arma::vec x, arma::vec m, arma::mat s)
//{
//    int k = x.size();
//    double det_s = arma::det(s);
//    static const float inv_sqrt_2pi = sqrt(pow(k, 2* M_PI) * det_s) ;
//
//    std::cout << "inv_sqrt_2pi" << inv_sqrt_2pi << '\n';
//    std::vector<float> value(k);
//
//}

int main (int argc, char * const argv[]) {
  //int i;
  //float sp_cal[10];
  //float parms[3] = {6.4, 0.0024, 0.3};
//
  //std::vector<int> Time_pred(10); 
  //std::cout << "Time_size " << Time_pred.size() << '\n';
//
  //float startNum=40, step=5;
  //std::generate(Time_pred.begin(), Time_pred.end(), [&startNum, &step]{ return startNum+=step;});
//
//
  //for (i=0; i<10; i++){
  //  std::cout << "Time_pred " << Time_pred[i] << '\n';
  //  sp_cal[i] = sp_numbers(Time_pred[i], parms);
  //  std::cout << "Sp_cal " << sp_cal[i] << '\n';
  //}
  float x=4, m=3, s=0.6, val;
  val = normal_pdf(x,m,s);
  //std::cout << "norm_pdf: " << val <<'\n';

  arma::vec X = {4, 0.6, 0.3};
  arma::vec Mu = {5, 1, 0.5};
  arma::vec b = {0.5, 0.5, 0.1};
  arma::mat B = arma::diagmat(b);
  arma::mat C = B.t() * B;

  arma::mat a = X - Mu;

  int k = X.size();
    double det_SiG = arma::det(C);
    arma::mat inv_SiG = arma::inv(C);
    static const float inv_sqrt_2pi_detS = 1/sqrt(pow(2* M_PI, k) * det_SiG);

  double mvval;


  arma::mat H = a.t() * a;
  arma::mat A = a * a.t();

  
  //mvval = inv_sqrt_2pi_detS * arma::exp(-0.5f * A * inv_SiG);

  
  std::cout << "a: " << a <<'\n';
  std::cout << "inv_SiG: " << inv_SiG <<'\n';
  std::cout << "H: " << H <<'\n';
  std::cout << "A: " << A <<'\n';
  //std::cout << "MVN_pdf: " << mvval <<'\n';


  return 0;
}