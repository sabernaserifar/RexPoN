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

#ifndef __LIST_H_
#define __LIST_H_

#include "rexpon_types.h"

int  Make_List( int, int, int, rexpon_list*, MPI_Comm );
void Delete_List( rexpon_list*, MPI_Comm );

inline int  Num_Entries(int,rexpon_list*);
inline int  Start_Index( int, rexpon_list* );
inline int  End_Index( int, rexpon_list* );
inline void Set_Start_Index(int,int,rexpon_list*);
inline void Set_End_Index(int,int,rexpon_list*);

#if defined(LAMMPS_REXPON)
inline int Num_Entries( int i, rexpon_list *l )
{
  return l->end_index[i] - l->index[i];
}

inline int Start_Index( int i, rexpon_list *l )
{
  return l->index[i];
}

inline int End_Index( int i, rexpon_list *l )
{
  return l->end_index[i];
}

inline void Set_Start_Index( int i, int val, rexpon_list *l )
{
  l->index[i] = val;
}

inline void Set_End_Index( int i, int val, rexpon_list *l )
{
  l->end_index[i] = val;
}
#endif // LAMMPS_REXPON

#endif
