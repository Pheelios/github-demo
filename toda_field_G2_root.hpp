/**************************************************************************
** C++ header file for toda_field_lambda_root class                      **                                     
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
the necessary values of the alpha and alpha_n vectors for the G2 affine
Toda Field Theory */

#ifndef G2ROOTDEF
#define G2ROOTDEF
#include "toda_field_abstract.hpp"

class toda_field_G2_root : public toda_field_abstract
{
	
	public:
		
		toda_field_G2_root(double b, double m, int L,int n_b,
		int n_s_b): toda_field_abstract(b , m , L , n_b , n_s_b)
		{
			alpha_set();
			alpha_n_set();
		}
			
			
		void alpha_set(void)
		{
			alpha[0][0] = sqrt(2.0);
	        alpha[0][1] = 0.0;
			alpha[1][0] = -1 / sqrt(2.0);
			alpha[1][1] = - sqrt(1.5);
			alpha[2][0] = - 1 / sqrt(2.0);
			alpha[2][1] = 1 / sqrt(6.0);
			
			for(int i = 0 ; i < 3 ; i++)
			{
				beta_alpha[i][0] = beta * alpha[i][0];
				beta_alpha[i][1] = beta * alpha[i][1];
			}
		}
		
		void alpha_n_set(void)
		{
			alpha_n[0] = 2.0;
			alpha_n[1] = 1.0;
			alpha_n[2] = 3.0;
		}
		
};
		
#endif