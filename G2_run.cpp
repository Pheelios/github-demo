/**************************************************************************
** C++ header file for toda_field_lambda_root class                      **                                     **
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

/*this class is inherited from the toda_field_abstract_class it defines
the necessary values of the alpha and alpha_n vectors for the lambda 
paramterized Toda Field Theories */

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <fstream>
#include "toda_field_G2_root.hpp"
#include <iomanip>
#include <cassert>
//#include <cmath>

std::vector<double>
create_init_fit_data(std::vector<std::vector<std::vector<double>>>, int);

double initial_fit(size_t n, double y[], double params[]);

int expb_f (const gsl_vector * x , void *data, gsl_vector * f);

int expb_df (const gsl_vector * x, void *params, gsl_matrix * J);

int expb_fdf (const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J);

struct data 
{
	size_t n;
	double * y;
	double * sigma;
};

int main(int argc, char *argv[])
{

    std ::cout << sin(90);
	double beta {atof(argv[1])};
	int length {atoi(argv[2])};
    double mass {atof(argv[3])};

    std::string str_beta(argv[1]);
    std::string str_L(argv[2]);
    std::string str_mass(argv[3]);

	std::string file_out_1 {"results/G2_B_" + str_beta + "_L_" + str_L + 
	                                      "_wall_corr_mass" + str_mass + ".dat"};
	std::cout << file_out_1 << std::endl;
	                        
	toda_field_G2_root G2_mass_experiment( beta, mass, length,atoi(argv[4]), atoi(argv[5]));
    
	G2_mass_experiment.initialise();
    
	for(int N = 0 ; N < (1 << atoi(argv[4])) ; N++)
	{
		G2_mass_experiment.field_update();
		G2_mass_experiment.sample_wall_correlation_function_buffer(N);
		G2_mass_experiment.sample_mean_field();
        if (N % 10000 == 0) std::cout << N << std::endl;
	}
	
	int fit_length = 16;
	std::vector<double> init_fit_data;
	init_fit_data = create_init_fit_data(G2_mass_experiment.wall_correlation_function_buffer, fit_length);
	double y[init_fit_data.size()];
	std::copy(init_fit_data.begin(), init_fit_data.end(), y);
	std::string file_out_2 {"results/G2_B_" + str_beta + "_L_" + str_L + 
	                          "_wall_corr_buff_mass" + str_mass + ".dat"};

    G2_mass_experiment.write_wall_correlation_function_buffer(file_out_2);
    G2_mass_experiment.calc_mean_field();
    double params[8] = {1, 1, 1, G2_mass_experiment.mean_field[0], 1, 1, 1, G2_mass_experiment.mean_field[1]};
  
   std::ofstream write_output(file_out_1);
    assert(write_output.is_open());
    write_output << initial_fit(48,y,params) << std::endl;
    for(int i = 0 ; i < 8 ; i++)
    {
    	write_output << params[i] << std::endl;
    }
    write_output << G2_mass_experiment.mean_field[0] << std::endl
                 << G2_mass_experiment.mean_field[1] << std::endl;
    write_output.close();

	return 0;
}

std::vector<double>
create_init_fit_data(std::vector<std::vector<std::vector<double>>> wall_buff, int fit_length)
{
	std::vector<double> 
	init_fit_data(3 * fit_length, 0.0);
    int number_of_wall_correlation_function_buffer_bins = (int) wall_buff[0][0].size();
	for (int j = 0 ; j < fit_length ; j++)
	{
		for (int k = 0 ; k < number_of_wall_correlation_function_buffer_bins ; k++)
		{
			init_fit_data[j] += wall_buff[0][j][k];
			init_fit_data[fit_length + j] += wall_buff[1][j][k];
			init_fit_data[2 * fit_length + j] += (wall_buff[2][j][k] + wall_buff[3][j][k]) / 2;
		}
	
}	
	for (int j = 0 ; j < 3 * fit_length ; j++)
	{
		init_fit_data[j] /= number_of_wall_correlation_function_buffer_bins;	
	}
	return init_fit_data;
}


double initial_fit(size_t n, double y[], double params[])
{
	const gsl_multifit_fdfsolver_type *T;
    gsl_multifit_fdfsolver *s;

    int status;
   
    const size_t p = 8;

    gsl_matrix *covar = gsl_matrix_alloc (p, p);
    double sigma[n];  
   
    for (int i = 0 ; i < n ; i++)
    {
        sigma[i] = 1.0;
    }
	

	struct data d = { n, y, sigma};

	gsl_multifit_function_fdf f;

	gsl_vector_view x = gsl_vector_view_array (params, p);

	f.f = &expb_f;
	f.df = &expb_df;
	f.fdf = &expb_fdf;
	f.n = n;
	f.p = p;
	f.params = &d;

	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc (T, n, p);
	gsl_multifit_fdfsolver_set (s, &f, &x.vector);

	size_t iter = 0;

	do
    {
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);

		if (status)
		break;

		status = gsl_multifit_test_delta (s->dx, s->x,
		                                1e-4, 1e-4);
    }while (status == GSL_CONTINUE && iter < 500);

    double ratio_1 = (gsl_vector_get(s->x, 2))/(gsl_vector_get(s->x, 6));
	double ratio_2 = (gsl_vector_get(s->x, 6))/(gsl_vector_get(s->x, 2));
    double mass_ratio;

	if (ratio_1 > ratio_2) 
	{
		mass_ratio = ratio_1;
	}
	else
	{
		mass_ratio = ratio_2;
	}
	for (int i = 0 ; i < 8 ; i++)
    {
    	params[i] = gsl_vector_get(s->x, i);
    }


	gsl_multifit_fdfsolver_free (s);
            
	return mass_ratio;  

}


/*************************************************************************************/

int expb_f (const gsl_vector * x , void *data,
                                gsl_vector * f)
{
    size_t n = ((struct data *)data)->n;
    double *y = ((struct data *)data)->y;

    double a = gsl_vector_get (x, 0);
    double theta_1 = gsl_vector_get (x, 1);
    double m_1 = gsl_vector_get (x, 2);
    double S_1 = gsl_vector_get (x, 3);
    double b = gsl_vector_get (x, 4);
    double theta_2 = gsl_vector_get (x, 5);
    double m_2 = gsl_vector_get (x, 6);
    double S_2 = gsl_vector_get (x, 7);

   
    double t[n];
    for (int i = 0; i < n; i++)
    {
        t[i] = (i%16)*1.0;
    }

    for (size_t i = 0; i < n; i++)
    {  
		double Yi =

		(i < 16) * (a * (cos(theta_1) * cos(theta_1)) * exp(-m_1 * t[i]) +
		b * (sin(theta_2) * sin(theta_2)) * exp(-m_2 * t[i]) + S_1 * S_1) +
		(i > 15 && i < 32) * (a * (sin(theta_1) * sin(theta_1)) * exp(-m_1 * t[i]) +
		b * (cos(theta_2) * cos(theta_2)) * exp(-m_2 * t[i]) + S_2 * S_2) +
		(i > 31 && i < 48) * (-a * (cos(theta_1) * sin(theta_1)) * exp(-m_1 * t[i])+
		b * (sin(theta_2) * cos(theta_2)) * exp(-m_2 * t[i]) + S_1 * S_2);
		gsl_vector_set (f, i, Yi - y[i]);
    }

    return GSL_SUCCESS;
}

/*************************************************************************************/

int expb_df (const gsl_vector * x, void *params, gsl_matrix * J)
{
    size_t n = ((struct data *)params)->n;
    
    double a = gsl_vector_get (x, 0);
    double theta_1 = gsl_vector_get (x, 1);
    double m_1 = gsl_vector_get (x, 2);
    double S_1 = gsl_vector_get (x, 3);
    double b = gsl_vector_get (x, 4);
    double theta_2 = gsl_vector_get (x, 5);
    double m_2 = gsl_vector_get (x, 6);
    double S_2 = gsl_vector_get (x, 7);
   
    double t[n];

    for (int i = 0; i < n; i++)
    {
        t[i] = (i%16)*1.0;
    }

    double J_0, J_1, J_2, J_3, J_4, J_5, J_6, J_7;
    
	for (size_t i = 0; i < n; i++)
	{
		J_0 = (i < 16) * ((cos(theta_1) * cos(theta_1)) * exp(-m_1 * t[i])) +
		      (i > 15 && i < 32) * ((sin(theta_1) * sin(theta_1)) * exp(-m_1 * t[i])) +
		      (i > 31 && i < 48) * (-(cos(theta_1) * sin(theta_1)) * exp(-m_1 * t[i]));

		J_1 = (i < 16) * ( -2 * a * (sin(theta_1) * cos(theta_1)) * exp( -m_1 * t[i])) +
		      (i > 15 && i < 32) * (2 * a * (cos(theta_1) * sin(theta_1)) * exp(-m_1 * t[i])) +
		      (i > 31 && i < 48) * (-a * ((cos(theta_1) * cos(theta_1) - 
		      sin(theta_1) * sin(theta_1)) * exp(-m_1 * t[i])));

		J_2 = (i < 16) * (a * (cos(theta_1) * cos(theta_1)) * -t[i] * exp(-m_1 * t[i])) +
		      (i > 15 && i < 32) * (a * (sin(theta_1) * sin(theta_1)) * -t[i] * exp(-m_1 * t[i])) +
		      (i > 31 && i < 48) * (-a * (cos(theta_1) * sin(theta_1)) *  -t[i] * exp(-m_1 * t[i]));

		J_3 = (i < 16) * (2 * S_1) + (i > 15 && i < 32) * (0) +
		      (i > 31 && i < 48) * (S_2);

		J_4 = (i < 16) * ((sin(theta_2) * sin(theta_2)) * exp(-m_2 * t[i])) +
		      (i > 15 && i < 32) * ((cos(theta_2) * cos(theta_2)) * exp(-m_2 * t[i])) +
		      (i > 31 && i < 48) * ((sin(theta_2) * cos(theta_2)) * exp(-m_2 * t[i]));

		J_5 = (i < 16) * (2 * b * (cos(theta_2) * sin(theta_2)) * exp(-m_2 * t[i])) +
		      (i > 15 && i < 32) * ( -2 * b * (sin(theta_2) * cos(theta_2)) * exp(-m_2 * t[i])) +
		      (i > 31 && i < 48) * ( b * (cos(theta_2) * cos(theta_2)
		                            - sin(theta_2) * sin(theta_2)) * exp(-m_2 * t[i]));

		J_6 = b * (sin(theta_2) * sin(theta_2)) * -t[i] * exp(-m_2 * t[i]) +
		      (i > 15 && i < 32) * (b * (cos(theta_2) * cos(theta_2)) * -t[i] * exp(-m_2 * t[i])) +
		      (i > 31 && i < 48) * (b * (sin(theta_2) * cos(theta_2)) * -t[i] * exp(-m_2 * t[i]));

		J_7 = (i < 16) * (0) + (i > 15 && i < 32) * (2 * S_2) +
		      (i > 31 && i < 48) * (S_1);

		gsl_matrix_set(J, i, 0, J_0);
		gsl_matrix_set(J, i, 1, J_1);
		gsl_matrix_set(J, i, 2, J_2);
		gsl_matrix_set(J, i, 3, J_3);
		gsl_matrix_set(J, i, 4, J_4);
		gsl_matrix_set(J, i, 5, J_5);
		gsl_matrix_set(J, i, 6, J_6);
		gsl_matrix_set(J, i, 7, J_7);
	}	
    return GSL_SUCCESS;
}

/*************************************************************************************/

int expb_fdf (const gsl_vector * x, void *params,
          gsl_vector * f, gsl_matrix * J)
{
	expb_f (x, params, f);
	expb_df (x, params, J);

	return GSL_SUCCESS;
}


