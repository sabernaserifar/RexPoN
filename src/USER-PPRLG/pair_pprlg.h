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

#ifdef PAIR_CLASS

PairStyle(pprlg,PairPPRLG)

#else

#ifndef LMP_PAIR_PPRLG_H
#define LMP_PAIR_PPRLG_H

#include "pair.h"

namespace LAMMPS_NS {

class PairPPRLG : public Pair {
 public:
  PairPPRLG(class LAMMPS *);
  virtual ~PairPPRLG();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  double single(int, int, int, int, double, double, double, double &);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  void *extract(const char *, int &);


 protected:
  double **d0,**alpha,**r0;
  double **beta,**gamma,**eta,**delta;
  double **c6lg,**r6lg;
  double cut_inner,cut_outer;   // read on the input
  double cut_inner2,cut_outer2,cut_inner3,cut_outer3;
  double tap[8];
  double dtap[8];

  void init_taper();
  void pair_interact(double, double, double, double, double, double, 
                     double, double, double, double, 
                     double *, double *);

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

*/
