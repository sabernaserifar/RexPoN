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

#ifndef __LOOKUP_H_
#define __LOOKUP_H_

#include "rexpon_types.h"

void Tridiagonal_Solve( const double *a, const double *b,
                        double *c, double *d, double *x, unsigned int n);

void Natural_Cubic_Spline( const double *h, const double *f,
                           cubic_spline_coef *coef, unsigned int n,
                           MPI_Comm comm );

void Complete_Cubic_Spline( const double *h, const double *f, double v0, double vlast,
                            cubic_spline_coef *coef, unsigned int n,
                            MPI_Comm comm );

int Init_Lookup_Tables( rexpon_system*, control_params*, storage*,
                        mpi_datatypes*, char* );

void Deallocate_Lookup_Tables( rexpon_system* );

#endif
