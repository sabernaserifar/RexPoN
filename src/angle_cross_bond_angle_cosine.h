/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ANGLE_CLASS

AngleStyle(cross/ba/cosine, AngleCrossBondAngleCosine)

#else

#ifndef ANGLE_CROSS_BOND_ANGLE_COSINE_H
#define ANGLE_CROSS_BOND_ANGLE_COSINE_H

#include "stdio.h"
#include "angle.h"

namespace LAMMPS_NS {

class AngleCrossBondAngleCosine : public Angle {
 public:
  AngleCrossBondAngleCosine(class LAMMPS *);
  ~AngleCrossBondAngleCosine();
  void compute(int, int);
  void coeff(int, char **);
  double equilibrium_angle(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double single(int, int, int, int);

 private:
  double *theta0,*k2,*k3,*k4;
  double *bb_k,*bb_r1,*bb_r2;
  double *ba_k1,*ba_k2,*ba_r1,*ba_r2;
  int *setflag_a,*setflag_bb,*setflag_ba;

  void allocate();
};

}

#endif
#endif
