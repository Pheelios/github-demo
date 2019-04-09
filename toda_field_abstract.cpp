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

#include"toda_field_abstract.hpp"

/*************************************************************************/

toda_field_abstract::toda_field_abstract                                              
(double beta_in, double mass_in,                                      
 int lattice_length_in, int number_samples_log_2_in, 
 int number_samples_per_wall_correlation_function_buffer_bin_log_2_in)
{
	number_samples_per_wall_correlation_function_buffer_bin_log_2 =         
    number_samples_per_wall_correlation_function_buffer_bin_log_2_in;

    number_samples_per_wall_correlation_function_buffer_bin = 
    1 << number_samples_per_wall_correlation_function_buffer_bin_log_2;       
                                                                                                                                
    number_of_samples_log_2 = number_samples_log_2_in;  
    number_of_samples = 1 << number_of_samples_log_2;                          

	number_of_wall_correlation_function_buffer_bins =                 
	number_of_samples/number_samples_per_wall_correlation_function_buffer_bin;
    update_weight = 1.0;
	beta =  beta_in;

	lattice_length = lattice_length_in;

	lattice_length_mm = lattice_length - 1;

	mass = mass_in;

	mass_over_beta_square = (mass / beta) * (mass / beta);

    decay_length.resize(2,0.0);

    mean_field.resize(2,0.0);

    alpha_n.resize(2,0.0);

    alpha.resize(3,std::vector<double>(2,0.0));

    beta_alpha.resize(3,std::vector<double>(2,0.0));

	field_lattice.resize(lattice_length, 
	std::vector<std::vector<double>>
	(lattice_length,std::vector<double>(2,0.0)));
				
    wall_correlation_function.resize(4, std::vector<double>(lattice_length, 0.0));

    wall_correlation_function_buffer.resize
    (4, std::vector<std::vector<double>>(lattice_length,std::vector<double>(
    		                 number_of_wall_correlation_function_buffer_bins)));

	field_lattice_auto_correlation_function_buffer.resize
	(lattice_length, std::vector<std::vector<std::vector<double>>>(lattice_length
		, std::vector<std::vector<double>>(2,std::vector<double>(auto_length,0.0))));

    auto_correlation_function.resize
    (2, std::vector<double>(auto_length, 0.0));
	
}

/*************************************************************************/

void toda_field_abstract::initialise()
{
	for(int j = 0 ; j < lattice_length ; j++)
	{
		for(int i = 0 ; i < lattice_length ; i++)
		{
			field_lattice[i][j][0] = 1 - 2 * rand.generate();
			field_lattice[i][j][1] = 1 - 2 * rand.generate();
		}
	}
}

/*************************************************************************/

void toda_field_abstract::set_update_weight(int N, double tol)
{
	bool U_lim_set{false}, L_lim_set{false};
	double U_lim{update_weight}, L_lim{update_weight};
	
	do
	{
		std::cout << "update weight = " 
		        << update_weight << " " <<std::flush;
		thermalise(N);
		
		if(!U_lim_set || !L_lim_set)
		{
			if(acceptance_ratio > 0.5 + tol)
			{
				L_lim = update_weight;
				L_lim_set = true;
				if(!U_lim_set)
				{
					update_weight *= 2.0;
				}
				else
				{
					update_weight = 0.5 * (U_lim + L_lim);
				}
					
					
			}
			else if(acceptance_ratio < 0.5 - tol)
			{
				U_lim = update_weight;
				U_lim_set = true;
				if(!L_lim_set)
				{
					update_weight *= 0.5;
				}
				else
				{
					update_weight = 0.5 * (U_lim + L_lim);
				}
			}
		}
		else
		{
			if(acceptance_ratio > 0.5 + tol)
			{
				L_lim = update_weight;
				update_weight = 0.5 * (U_lim + L_lim);
			}
			else if(acceptance_ratio < 0.5 - tol)
			{
				U_lim = update_weight;
				update_weight = 0.5 * (U_lim + L_lim);
			}
		}
		std::cout << "acceptance ratio = " << acceptance_ratio <<std::endl;
				  
	}while(fabs(acceptance_ratio - 0.5) >= tol);
	std::cout << " End of Ratio Set " << std::endl;
}

/*************************************************************************/

double toda_field_abstract::lattice_site_energy(int i, int j, int k)
{
	return

	0.5 *
	((field_lattice[i][j][k] - field_lattice[(i + 1) % lattice_length][j][k]) *
	( field_lattice[i][j][k] - field_lattice[(i + 1) % lattice_length][j][k]) +
	( field_lattice[i][j][k] - field_lattice[(i + lattice_length_mm) % lattice_length][j][k]) *
	( field_lattice[i][j][k] - field_lattice[(i + lattice_length_mm ) % lattice_length][j][k]) +
	( field_lattice[i][j][k] - field_lattice[i][(j + 1) % lattice_length][k]) *
	( field_lattice[i][j][k] - field_lattice[i][(j + 1) % lattice_length][k]) +
	( field_lattice[i][j][k] - field_lattice[i][(j + lattice_length_mm) % lattice_length][k]) *
	( field_lattice[i][j][k] - field_lattice[i][(j + lattice_length_mm) % lattice_length][k])) +
	mass_over_beta_square *
	(alpha_n[0] * (exp(beta_alpha[0][0] * field_lattice[i][j][0] +
	           beta_alpha[0][1] * field_lattice[i][j][1]) - 1.0) +
	alpha_n[1]  * (exp(beta_alpha[1][0] * field_lattice[i][j][0] +
	           beta_alpha[1][1] * field_lattice[i][j][1]) - 1.0) +
	alpha_n[2] * (exp(beta_alpha[2][0] * field_lattice[i][j][0] +
	    beta_alpha[2][1] * field_lattice[i][j][1]) - 1.0));
		 				     	
}

/*************************************************************************/

void toda_field_abstract::field_update()
{
		
	double energy_before, energy_after, field_old;		

	for(int j = 0 ; j < lattice_length; j++)
	{
		for(int i = 0 ; i < lattice_length ; i++)
		{
			for(int k = 0; k < 2 ; k++)
			{
			    energy_before = lattice_site_energy(i, j, k);
			    field_old = field_lattice[i][j][k];
			    field_lattice[i][j][k] += 2 * update_weight * (rand.generate() - 0.5);
				energy_after = lattice_site_energy(i, j, k);
		
				if(exp(energy_before - energy_after) > rand.generate()) accepted_updates++;
				else field_lattice[i][j][k] = field_old;
			}
		}
	}
}		

/*************************************************************************/	

void toda_field_abstract::thermalise(int N)
{
	for(int n = 0 ; n < N ; n++)
	{
		field_update();
	}
	
	accepted_updates = 0;
	
	for(int n = 0 ; n < N ; n++)
	{
		field_update();
	}
	acceptance_ratio = 
	((double)accepted_updates)/( 2 * N * lattice_length * lattice_length);	
}

/*************************************************************************/	

void toda_field_abstract::sample_mean_field()
{
	for(int i = 0; i < lattice_length ; i++)
	{
		for(int j = 0 ; j < lattice_length ; j++)
		{
			mean_field[0] += field_lattice[i][j][0];
			mean_field[1] += field_lattice[i][j][1];
		}
	}
}

/*************************************************************************/

void toda_field_abstract::sample_wall_correlation_function()
{
	double temp_0_0, temp_1_0, temp_0_n,temp_1_n, temp_0, 
	       temp_1, temp_2, temp_3;

	for(int n = 0 ; n < lattice_length ; n++)
	{
	    temp_0 = 0.0;
	    temp_1 = 0.0;
		temp_2 = 0.0;
		temp_3 = 0.0;

		for(int i = 0 ; i < lattice_length ; i++)
		{
			temp_0_0 = 0.0;
			temp_1_0 = 0.0;
            temp_0_n = 0.0;
			temp_1_n = 0.0;

			for(int j = 0 ; j < lattice_length ; j++)
			{
				temp_0_0 += field_lattice[i][j][0];
				temp_1_0 += field_lattice[i][j][1];
				temp_0_n += field_lattice[(i + n) % lattice_length][j][0];
				temp_1_n += field_lattice[(i + n) % lattice_length][j][1];
			}

			temp_0 += temp_0_0 * temp_0_n;
			temp_1 += temp_1_0 * temp_1_n;
			temp_2 += temp_0_0 * temp_1_n;
			temp_3 += temp_1_0 * temp_0_n;

		}

		wall_correlation_function[0][n] += temp_0;
		wall_correlation_function[1][n] += temp_1;
		wall_correlation_function[2][n] += temp_2;
		wall_correlation_function[3][n] += temp_3;

	}

	for(int n = 0 ; n < lattice_length ; n++)
	{
		temp_0 = 0.0;
		temp_1 = 0.0;
		temp_2 = 0.0;
		temp_3 = 0.0;
		
		for(int j = 0 ; j < lattice_length ; j++)
		{
			temp_0_0 = 0.0;
			temp_1_0 = 0.0;
			temp_0_n = 0.0;
			temp_1_n = 0.0;
			
			for(int i = 0 ; i < lattice_length ; i++)
			{
				temp_0_0 += field_lattice[i][j][0];
				temp_1_0 += field_lattice[i][j][1];
				temp_0_n += field_lattice[i][(j + n) % lattice_length][0];
				temp_1_n += field_lattice[i][(j + n) % lattice_length][1];
			}
			
			temp_0 += temp_0_0 * temp_0_n;
			temp_1 += temp_1_0 * temp_1_n;
			temp_2 += temp_0_0 * temp_1_n;
			temp_3 += temp_1_0 * temp_0_n;
		}

		wall_correlation_function[0][n] += temp_0;
		wall_correlation_function[1][n] += temp_1;
		wall_correlation_function[2][n] += temp_2;
		wall_correlation_function[3][n] += temp_3;
	}
}

/*************************************************************************/

void toda_field_abstract::
     sample_wall_correlation_function_buffer(int sample_number)
{   
	int 
	bin_number{sample_number / number_samples_per_wall_correlation_function_buffer_bin};
	
	double
    temp_0_n, temp_1_n, temp_0, temp_1, temp_2, 
    temp_3, temp_0_0, temp_1_0, 
    div_const
	{2.0 * number_samples_per_wall_correlation_function_buffer_bin * pow(lattice_length,3)};
	
	for(int n = 0 ; n < lattice_length ; n++)
	{
		temp_0 = 0.0;
		temp_1 = 0.0;
		temp_2 = 0.0;
		temp_3 = 0.0;
		
		for(int i = 0 ; i < lattice_length ; i++)
		{
			temp_0_0 = 0.0;
			temp_1_0 = 0.0;
			temp_0_n = 0.0;
			temp_1_n = 0.0;
			
			for(int j = 0; j < lattice_length ; j++)
			{
				temp_0_0 += field_lattice[i][j][0];
				temp_1_0 += field_lattice[i][j][1];
				temp_0_n += field_lattice[(i + n) % lattice_length][j][0];
				temp_1_n += field_lattice[(i + n) % lattice_length][j][1];
			}
			
			temp_0 += temp_0_0 * temp_0_n;
			temp_1 += temp_1_0 * temp_1_n;
			temp_2 += temp_0_0 * temp_1_n;
			temp_3 += temp_1_0 * temp_0_n;
		}
		
		wall_correlation_function_buffer[0][n][bin_number] += temp_0;
		wall_correlation_function_buffer[1][n][bin_number] += temp_1;
		wall_correlation_function_buffer[2][n][bin_number] += temp_2;
		wall_correlation_function_buffer[3][n][bin_number] += temp_3;
	}
	
	for(int n = 0 ; n < lattice_length ; n++)
	{
		temp_0 = 0.0;
		temp_1 = 0.0;
		temp_2 = 0.0;
		temp_3 = 0.0;
		
		for(int j = 0 ; j < lattice_length ; j++)
		{
			temp_0_0 = 0.0;
			temp_1_0 = 0.0;
			temp_0_n = 0.0;
			temp_1_n = 0.0;
			
			for(int i = 0 ; i < lattice_length ; i++)
			{
				temp_0_0 += field_lattice[i][j][0];
				temp_1_0 += field_lattice[i][j][1];
				temp_0_n += field_lattice[i][(j+n)%lattice_length][0];
				temp_1_n += field_lattice[i][(j+n)%lattice_length][1];
			}
			
			temp_0 += temp_0_0 * temp_0_n;
			temp_1 += temp_1_0 * temp_1_n;
			temp_2 += temp_0_0 * temp_1_n;
			temp_3 += temp_1_0 * temp_0_n;
		}
		wall_correlation_function_buffer[0][n][bin_number] += temp_0;
		wall_correlation_function_buffer[1][n][bin_number] += temp_1;
		wall_correlation_function_buffer[2][n][bin_number] += temp_2;
		wall_correlation_function_buffer[3][n][bin_number] += temp_3;
	}
	
	if((sample_number % number_samples_per_wall_correlation_function_buffer_bin)
		          == (number_samples_per_wall_correlation_function_buffer_bin - 1))
	{
		for (int n = 0 ; n < lattice_length ; n++)
		{   
			wall_correlation_function_buffer[0][n][bin_number] /= div_const;
			wall_correlation_function_buffer[1][n][bin_number] /= div_const;
			wall_correlation_function_buffer[2][n][bin_number] /= div_const;
			wall_correlation_function_buffer[3][n][bin_number] /= div_const;
		}
	}
}

/*************************************************************************/

void toda_field_abstract::shift_field_lattice_auto()
{
	for( int n = 0 ; n < (auto_length - 1) ; n++)
	{
		for(int i = 0 ; i < lattice_length ; i++)
		{
			for(int j = 0 ; j < lattice_length ; j++)
			{
				field_lattice_auto_correlation_function_buffer[i][j][0][n] =
				field_lattice_auto_correlation_function_buffer[i][j][0][n + 1];
				field_lattice_auto_correlation_function_buffer[i][j][1][n] =
				field_lattice_auto_correlation_function_buffer[i][j][1][n + 1];
			}
		}
	}

	for(int i = 0 ; i < lattice_length ; i++)
	{
		for(int j = 0 ; j < lattice_length ; j++)
		{
			field_lattice_auto_correlation_function_buffer[i][j][0][auto_length - 1] =
			field_lattice[i][j][0];
			field_lattice_auto_correlation_function_buffer[i][j][1][auto_length - 1] = 
			field_lattice[i][j][1];
		}
	}	
	
}

/*************************************************************************/

void toda_field_abstract::sample_auto_correlation_function(int N)
{	
	shift_field_lattice_auto();
    if(N > auto_length)
    {
    	for(int k = 0 ; k < 2 ; k++)
		{
			for(int n = 0 ; n < auto_length ; n++)
			{
				for(int i = 0 ; i < lattice_length ; i++)
				{
					for(int j = 0 ; j < lattice_length ; j++)
					{
						auto_correlation_function[k][n] +=
						field_lattice_auto_correlation_function_buffer[i][j][k][0] *
						field_lattice_auto_correlation_function_buffer[i][j][k][n];
					}
				}
			}
		}
    }
}

/*************************************************************************/

void toda_field_abstract::
     set_mass(double desired_decay_length, int number_iterations, double tolerance)
{
    mass = 1.0;
	
	double lower_limit, upper_limit, min_decay_length, old_mass;

	bool continue_flag{true}, upper_limit_flag{false}, lower_limit_flag{false};

	while(continue_flag)
	{
        mass_over_beta_square = (mass / beta) * (mass / beta); //no problem
		set_update_weight(20000, 0.05);                  //no problem  
		                               
		for(int i = 0 ; i < 4 ; i++)                         
		{
			for(int j = 0 ; j < lattice_length ; j++)
			{
				wall_correlation_function[i][j] = 0.0;   //no problem
			}
		}

		for (int n = 0 ; n < number_iterations ; n++)
		{
			field_update();
			sample_wall_correlation_function();    //no problem
		}
		decay_length[0] = calculate_wall_correlation_function_decay_length(0, 0.01);   //no problem
	    decay_length[1] = calculate_wall_correlation_function_decay_length(1, 0.01);   //no problem

		if(decay_length[0] <= decay_length[1])
		{ 
			min_decay_length = decay_length[0];                                        //no problem
		}
		else
		{
			min_decay_length = decay_length[1];                                        //no problem
		}

        
		if( (fabs(min_decay_length - desired_decay_length) < tolerance) || 
			(fabs((old_mass - mass)/old_mass) < 0.00001) )
		{
			continue_flag = false;
		}

		else if(min_decay_length < desired_decay_length)
		{
			upper_limit = mass;
			old_mass = mass;
			upper_limit_flag = true;

			if(!lower_limit_flag)
			{
				mass *= 0.5;
			}
			else
			{
				mass = 0.5 * (upper_limit + lower_limit);
			}
			
		}
		else if(min_decay_length > desired_decay_length)
		{
			lower_limit = mass;
			old_mass = mass;
			lower_limit_flag = true;
			if(!upper_limit_flag)
			{
				mass *= 2.0;
			}
			else
			{
				mass = 0.5 * (upper_limit + lower_limit);
			}
		}
	}

}

/*************************************************************************/

double toda_field_abstract::
       calculate_wall_correlation_function_decay_length(int field_component, double ratio)
{
	double decay_length, a, b, i_ratio,
	       f_sq{wall_correlation_function[field_component][lattice_length/2]};
	      
	for(int i = 0 ; i < lattice_length ; i++)
	{
		i_ratio = (wall_correlation_function[field_component][i] - f_sq) / 
			   (wall_correlation_function[field_component][0] - f_sq);
        std::cout << i << " " << i_ratio << std::endl;
		
		if(i_ratio == ratio)
		{
			decay_length = (double) i;
			break;
		}
		else if(i_ratio > ratio)
		{
			a = i_ratio;
		}
		else
		{
			b = i_ratio;

            if(a > b)
            {
            	decay_length = ((ratio - a) / (b - a)) + (i - 1);
			    break;
            }
            else
            {
            	decay_length = 1.0;
            	break;
            }
		}
	}
	
	return decay_length;
}

/*************************************************************************/

void toda_field_abstract::write_wall_correlation_function(std::string output_file)
{
	std::ofstream write_output(output_file);
	assert(write_output.is_open());
	write_output.precision(10);
	
	write_output << " % <phi_1(0)|phi_1(x)>_wall;" << std::endl;
	write_output << " w = [ ";

	for (int j=0 ; j<lattice_length ; j++)
	{
		write_output << wall_correlation_function[0][j] << " ";
	}
	write_output << " ];" << std::endl;
	
	write_output << " % <phi_2(0)|phi_2(x)>_wall;" << std::endl;
	write_output << " x = [ ";

	for (int j = 0 ; j < lattice_length ; j++)
	{
		write_output << wall_correlation_function[1][j] << " ";
	}
	write_output << " ];" << std::endl;
	
	write_output << " % <phi_1(0)|phi_2(x)>_wall;" << std::endl;
	write_output << " y = [ ";
	for (int j = 0 ; j < lattice_length ; j++)
	{
		write_output << wall_correlation_function[2][j] << " ";
	}
	write_output << " ];" << std::endl;
	
	write_output << " % <phi_2(0)|phi_2(x)>_wall;" << std::endl;
	write_output << " z = [ ";
	for (int j = 0 ; j < lattice_length ; j++)
	{
		write_output << wall_correlation_function[3][j] << " ";
	}
	write_output << " ];" << std::endl;
	write_output << "mass = " << mass << std::endl;
	write_output << "decay length phi_0 = " << decay_length[0] << std::endl;
	write_output << "decay length phi_1 = " << decay_length[1] << std::endl;

	write_output.close();
}

/*************************************************************************/

void toda_field_abstract::
	 write_wall_correlation_function_buffer(std::string output_file)
{
	std::ofstream write_output(output_file);
	assert(write_output.is_open());
	write_output.precision(10);
    
    for(int n = 0 ; n < number_of_wall_correlation_function_buffer_bins ; n++)
    {
    	for(int i = 0 ; i < 4 ; i++)
    	{
    		for(int j = 0 ; j < lattice_length ; j++)
    		{
    			write_output << wall_correlation_function_buffer[i][j][n] << "  ";
    		}
    		write_output << std::endl;
    	}
        write_output << std::endl;
    }
  
	write_output.close();
}

/*************************************************************************/

void toda_field_abstract::write_auto_correlation_function(std::string output_file)
{
	std::ofstream write_output(output_file);
	assert(write_output.is_open());
	write_output.precision(10);
	write_output << "% <phi_1_auto>" << std::endl
	             <<  "phi_1_auto = [ ";
	for (int j=0 ; j < auto_length - 1 ; j++)
	{
		write_output << auto_correlation_function[0][j] / auto_correlation_function[0][0] 
		             << " ";
	}
	write_output << " ];" << std::endl;
	
	write_output << "% <phi_2_auto>"  << std::endl 
	             <<  "phi_2_auto = [ ";
	for (int j = 0 ; j < auto_length - 1 ; j++)
	{
		write_output << auto_correlation_function[1][j] / auto_correlation_function[1][0]
		             << " ";
	}
	write_output << " ];" << std::endl;
	write_output << "figure;" << std::endl << "plot(0:" 
	             << auto_length - 2 << ", phi_1_auto);" << std::endl;
	write_output << "figure;" << std::endl << "plot(0:" 
	             << auto_length - 2 << ", phi_2_auto);" << std::endl;
	write_output.close();
}

/*************************************************************************/

void toda_field_abstract::calc_mean_field()
{
	long long div = (long long) number_of_samples 
	              * (long long) lattice_length 
	              * (long long) lattice_length;
	mean_field[0] /= div;
	mean_field[1] /= div;
}