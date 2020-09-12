/*----------------------------------------------------------------------
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details:
  <http://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#ifndef __BOND_ORDERS_H_
#define __BOND_ORDERS_H_

#include "rexpon_types.h"

typedef struct{
  double C1dbo, C2dbo, C3dbo;
  double C1dbopi, C2dbopi, C3dbopi, C4dbopi;
  double C1dbopi2, C2dbopi2, C3dbopi2, C4dbopi2;
  double C1dDelta, C2dDelta, C3dDelta;
}dbond_coefficients;

void Add_dBond_to_Forces( rexpon_system*, int, int, storage*, rexpon_list** );
void Add_dBond_to_Forces_NPT( int, int, simulation_data*,
                              storage*, rexpon_list** );
int BOp(storage*, rexpon_list*, double, int, int, far_neighbor_data*,
        single_body_parameters*, single_body_parameters*, two_body_parameters*);
void BO( rexpon_system*, control_params*, simulation_data*,
         storage*, rexpon_list**, output_controls* );
#endif
