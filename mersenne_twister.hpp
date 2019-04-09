#ifndef MERSENNETWISTERDEF
#define MERSENNETWISTERDEF

/* This is a simple class to interface with C++ built in mersenne twister
random number generator. */ 

#include<random>

class mersenne_twister{
	
	public:
	
	    int seed;
		double lower_limit, upper_limit;
		std::mt19937_64 generator;
		std::uniform_real_distribution<> distribution;
		
		mersenne_twister(){}
		
		void set_seed(unsigned sd){
			generator.seed(sd);	
		}
		
		void set_limits(double l_lim, double u_lim){
			std::uniform_real_distribution<> set(l_lim,u_lim);
			distribution.param(set.param());
		}
		
		mersenne_twister(int sd, double l_lim, double u_lim){
			generator.seed(sd);
			seed=sd;
			std::uniform_real_distribution<> set(l_lim,u_lim);
			distribution.param(set.param());
			lower_limit=l_lim;
			upper_limit=u_lim;
			}
			
		double generate(void){
			double output=distribution(generator);
			return output;
		}
};

#endif
