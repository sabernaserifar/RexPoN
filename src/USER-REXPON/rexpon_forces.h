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

#ifndef __FORCES_H_
#define __FORCES_H_

#include "rexpon_types.h"
#include "rexpon_defs.h"

extern interaction_function Interaction_Functions[NUM_INTRS];

void Init_Force_Functions( control_params* );
void Compute_Forces( rexpon_system*, control_params*, simulation_data*,
                     storage*, rexpon_list**, output_controls*, mpi_datatypes* );
void Estimate_Storages( rexpon_system*, control_params*, rexpon_list**,
                        int*, int*, int*, int*, MPI_Comm );
#endif
