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

#ifndef __NONBONDED_H_
#define __NONBONDED_H_

#include "rexpon_types.h"

void vdW_Coulomb_Energy( rexpon_system*, control_params*, simulation_data*,
                         storage*, rexpon_list**, output_controls* );

void Tabulated_vdW_Coulomb_Energy( rexpon_system*, control_params*,
                                   simulation_data*, storage*,
                                   rexpon_list**, output_controls* );

void Compute_Polarization_Energy( rexpon_system*, simulation_data* );

void LR_vdW_Coulomb( rexpon_system*, storage*, control_params*,
                int, int, double, LR_data* );
#endif
