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

#ifndef __ALLOCATE_H_
#define __ALLOCATE_H_

#include "rexpon_types.h"
int  PreAllocate_Space( rexpon_system*, control_params*, storage*, MPI_Comm );

int  Allocate_System( rexpon_system*, int, int, char* );
void DeAllocate_System( rexpon_system* );

int  Allocate_Workspace( rexpon_system*, control_params*, storage*,
                         int, int, MPI_Comm, char* );
void DeAllocate_Workspace( control_params*, storage* );

void ReAllocate( rexpon_system*, control_params*, simulation_data*, storage*,
                 rexpon_list**, mpi_datatypes* );
#endif
