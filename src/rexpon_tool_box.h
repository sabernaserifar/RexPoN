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

#ifndef __TOOL_BOX_H_
#define __TOOL_BOX_H_

#include "rexpon_types.h"
#include "rexpon_defs.h"

/* from system_props.h */
double Get_Time( );

/* from io_tools.h */
int   Tokenize( char*, char*** );

/* from lammps */
void *smalloc( rc_bigint, const char*, MPI_Comm );
void *scalloc( rc_bigint, rc_bigint, const char*, MPI_Comm );
void sfree( void*, const char* );
#endif
