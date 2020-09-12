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

#ifndef __IO_TOOLS_H_
#define __IO_TOOLS_H_

#include "rexpon_types.h"

int Init_Output_Files( rexpon_system*, control_params*,
                       output_controls*, mpi_datatypes*, char* );
int Close_Output_Files( rexpon_system*, control_params*,
                        output_controls*, mpi_datatypes* );
void  Output_Results( rexpon_system*, control_params*, simulation_data*,
                      rexpon_list**, output_controls*, mpi_datatypes* );
#endif
