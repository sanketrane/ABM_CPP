#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <cmath>
#include <string>
#include <vector>

//#include "custom_functions.cpp"


namespace fs = std::filesystem;

//setting current WD and filepaths for input and output
auto wdir = fs::current_path();
fs::path datfile ("artf_data.csv");
fs::path fulldatfile_path = wdir / datfile;

// spline1 --
float sp_numbers(float Time, float params[]){ 
    float sp_cal, theta0=params[0], nu=params[1], psi=params[2];
    // thymic FoxP3 negative SP4 T cell numbers - spline defs
    sp_cal = psi * std::pow(10., theta0) * std::exp(-1. * nu * Time); 
    return sp_cal;
}

int main (int argc, char * const argv[]) {
  //arma::vec A1(4), A2(4), A3(4);
  //A1.imbue( [&]() { return randnorm<float>(4, 0.5); } );
  //A2.imbue( [&]() { return randnorm<float>(0.6, 0.3); } );
  //A3.imbue( [&]() { return randnorm<float>(0.3, 0.1); } );
  //arma::mat A = join_rows(A1, A2, A3);

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

  return 0;
}