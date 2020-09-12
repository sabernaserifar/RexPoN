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

#ifndef __RESET_TOOLS_H_
#define __RESET_TOOLS_H_

#include "rexpon_types.h"

void Reset_Pressures( simulation_data* );
void Reset_Simulation_Data( simulation_data*, int );
void Reset_Timing( rexpon_timing* );
void Reset_Workspace( rexpon_system*, storage* );
void Reset_Neighbor_Lists( rexpon_system*, control_params*, storage*,
                           rexpon_list**, MPI_Comm );
void Reset( rexpon_system*, control_params*, simulation_data*, storage*,
            rexpon_list**, MPI_Comm );
#endif
