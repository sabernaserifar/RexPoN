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

#ifdef DIHEDRAL_CLASS

DihedralStyle(cross/13bb,DihedralCross13BondBond)

#else

#ifndef LMP_DIHEDRAL_CROSS_13_BOND_BOND_H
#define LMP_DIHEDRAL_CROSS_13_BOND_BOND_H

#include "stdio.h"
#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralCross13BondBond : public Dihedral {
 public:
  DihedralCross13BondBond(class LAMMPS *);
  ~DihedralCross13BondBond();
  void compute(int, int);
  void coeff(int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);

 private:
  double *bb13t_k,*bb13t_r10,*bb13t_r30;
  void allocate();
};

}

#endif
#endif
