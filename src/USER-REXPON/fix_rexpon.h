/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(REXPON,FixRexPoN)

#else

#ifndef LMP_FIX_REXPON_H
#define LMP_FIX_REXPON_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRexPoN : public Fix {
  friend class PairRexPoN;
  //friend class PairRexPoNOMP;

 public:
  FixRexPoN(class LAMMPS *,int, char **);
  ~FixRexPoN();
  int setmask();

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);

 private:
  int maxbonds;              // max # of bonds for any atom
  int maxhbonds;             // max # of Hbonds for any atom
  int *num_bonds;            // # of bonds for each atom
  int *num_hbonds;           // # of Hbonds for each atom
};

}

#endif
#endif
