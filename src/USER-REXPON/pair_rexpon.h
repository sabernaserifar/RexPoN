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

/* ----------------------------------------------------------------------
   Contributing author: Saber Naserifar (Caltech, naserifar.saber@gmail.com)
   The LAMMPS REAXC code was modified to include REXPON force field
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(rexpon,PairRexPoN)

#else

#ifndef LMP_PAIR_REXPON_H
#define LMP_PAIR_REXPON_H

#include "pair.h"
#include "rexpon_types.h"

namespace LAMMPS_NS {

class PairRexPoN : public Pair {
 public:
  PairRexPoN(class LAMMPS *);
  ~PairRexPoN();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  virtual void init_style();
  double init_one(int, int);
  void *extract(const char *, int &);
  int fixbond_flag, fixspecies_flag;
  int **tmpid;
  double **tmpbo,**tmpr;

  control_params *control;
  rexpon_system *system;
  output_controls *out_control;
  simulation_data *data;
  storage *workspace;
  rexpon_list *lists;
  mpi_datatypes *mpi_data;

  bigint ngroup;

 protected:
  double cutmax;
  int nelements;                // # of unique elements
  char **elements;              // names of unique elements
  int *map;
  class FixRexPoN *fix_rexpon;

  double *chi,*eta,*gamma;
  int qeqflag;
  int setup_flag;
  int firstwarn;

  void allocate();
  void setup();
  void create_compute();
  void create_fix();
  void write_rexpon_atoms();
  void get_distance(rvec, rvec, double *, rvec *);
  void set_far_nbr(far_neighbor_data *, int, double, rvec);
  int estimate_rexpon_lists();
  int write_rexpon_lists();
  void read_rexpon_forces(int);

  int nmax;
  void FindBond();
  double memory_usage();

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Too many ghost atoms

Number of ghost atoms has increased too much during simulation and has exceeded
the size of rexpon arrays.  Increase safe_zone and min_cap in pair_style rexpon
command

*/
