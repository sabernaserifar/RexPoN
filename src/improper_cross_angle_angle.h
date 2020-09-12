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

#ifdef IMPROPER_CLASS

ImproperStyle(cross/aa, ImproperCrossAngleAngle)

#else

#ifndef LMP_IMPROPER_CROSS_ANGLE_ANGLE_H
#define LMP_IMPROPER_CROSS_ANGLE_ANGLE_H

#include "stdio.h"
#include "improper.h"

namespace LAMMPS_NS {

class ImproperCrossAngleAngle : public Improper {
 public:
  ImproperCrossAngleAngle(class LAMMPS *);
  ~ImproperCrossAngleAngle();
  void compute(int, int);
  void coeff(int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);

 private:
  double *aa_k1,*aa_k2,*aa_k3,*aa_theta0_1,*aa_theta0_2,*aa_theta0_3;
  double PI;

  void allocate();
};

}

#endif
#endif
