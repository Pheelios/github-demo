/**************************************************************************
** C++ header file for abstract Toda field class.                        **
** Copyright (C) <2019> <Philip Giokas> <philipgiokas@gmail.com>         **
***************************************************************************
** This program is free software: you can redistribute it and/or modify  **
** it under the terms of the GNU General Public License as published by  **
** the Free Software Foundation, either version 3 of the License, or     **
** (at your option) any later version.                                   **
***************************************************************************
** This program is distributed in the hope that it will be useful,       **
** but WITHOUT ANY WARRANTY; without even the implied warranty of        **
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         **
** GNU General Public License for more details.                          **
***************************************************************************
** You should have received a copy of the GNU General Public License     **
** along with this program.  If not, see <https://www.gnu.org/licenses/> **
**************************************************************************/

/*The Toda Field Abstract class defines all the variables and methods common
to all Toda field monte carlo simulations, the variables that define the
paticular type of Toda theory that is simulated are defined in inherited 
classes.*/

#ifndef TODAFIELDABSTRACTDEF
#define TODAFIELDABSTRACTDEF
#include <cmath>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <string>
#include "mersenne_twister.hpp"


class toda_field_abstract
{
	    
	public:
	    
        std::vector<double> 
		decay_length, mean_field, alpha_n;
	
	    std::vector<std::vector<double>> 
	    alpha, beta_alpha, auto_correlation_function, 
	    wall_correlation_function;
	
	    std::vector<std::vector<std::vector<double>>>
	    field_lattice, wall_correlation_function_buffer;

	    std::vector<std::vector<std::vector<std::vector<double>>>>
		field_lattice_auto_correlation_function_buffer;

		double 
		beta, mass, mass_over_beta_square, update_weight, acceptance_ratio;
		
		int 
		lattice_length, lattice_length_mm, number_of_samples, auto_length {100},
		accepted_updates, number_of_wall_correlation_function_buffer_bins,
		number_samples_per_wall_correlation_function_buffer_bin_log_2,
		number_samples_per_wall_correlation_function_buffer_bin,
		number_of_samples_log_2;
		
		mersenne_twister 
		rand;
		
		toda_field_abstract(double, double, int, int, int);
		double calculate_wall_correlation_function_decay_length(int, double);
		double lattice_site_energy(int, int, int);
		void initialise();
		void set_update_weight(int, double);
		void set_mass(double, int, double);
		void field_update();
		void thermalise(int);
		void sample_mean_field();

		void calc_mean_field();
		void sample_wall_correlation_function_buffer(int);
		void sample_wall_correlation_function();
		void shift_field_lattice_auto();
		void sample_auto_correlation_function(int);
		void write_wall_correlation_function(std::string);
		void write_auto_correlation_function(std::string);
		void write_wall_correlation_function_buffer(std::string);
		virtual void alpha_set() = 0;
		virtual void alpha_n_set() = 0;
};		

#endif