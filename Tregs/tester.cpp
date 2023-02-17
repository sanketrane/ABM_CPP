#include <iostream>
#include "custom_functions.cpp"
#include <filesystem>
#include </opt/homebrew/Cellar/libomp/15.0.7/include/omp.h>


int main () {
  #pragma omp parallel num_threads(3)
  {
    int id = omp_get_thread_num();
    int x=4, y=1;
    int data=randnorm<double>(x,y);
    int total = omp_get_num_threads();
    printf("Greetings from process %d out of %d with Data %d\n", id, total, data);
  }
  printf("parallel for ends.\n");
    
  //cout << k << endl;

  //double xd, yd, n;
  //cout << "Enter Mean and SD: " << endl;
  //cin >> xd >> yd;
  //n=randnorm<double>(xd,yd);
  //cout << "Your sample from mean " << xd << " and sd " << yd << " is: " << n << endl;
//
  //cout << "Current path is " << std::filesystem::current_path() << '\n'; 
  return 0;
}