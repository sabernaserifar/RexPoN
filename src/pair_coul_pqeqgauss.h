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

#ifdef PAIR_CLASS

PairStyle(coul/pqeqgauss,PairCoulPqeqgauss)

#else

#ifndef LMP_PAIR_COUL_PQEQGAUSS_H
#define LMP_PAIR_COUL_PQEQGAUSS_H

#include "pair.h"

namespace LAMMPS_NS {

class PairCoulPqeqgauss : public Pair {
 public:
  PairCoulPqeqgauss(class LAMMPS *);
  ~PairCoulPqeqgauss();
  void compute(int, int);
  void allocate();
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  virtual void *extract(const char *, int &);

  // ============================
  // parameters of the qeq model
  double *chi;       // electronegativity chi_i corresponds to  dE/dq_i  
  double *idem;      // idempotential     idem_i = J_i corresponds to 1/2*d^2E/dq_i^2
  double *rcore;     // width of the core charge
  // for polarizable atoms 
  int *polar;        // =1 for polarizable atom, =0 for non polarizable
  double *qcore;     // the core charge (~ Z - #valence electrons)
  double *rshell;    // width of the (movable) shell charge
  double *kstring2;  // coefficient at the 2nd power string that holds the shell near the core
  double *kstring4;  // coefficient at the 4th power string that holds the shell near the core
  // ============================
  // auxiliarry array with erf coeficients
  // alpha = 0.5 / R**2 
  // alpha_ij = sqrt( alpha_i * alpha_j / (alpha_i + alpha_j) )
  // we assume that all coefficients are diagonal for now
  double **alphass;  // shell - shell 
  double **alphasc;  // shell - core
  double **alphacc;  // core - core
  double *swab; // cut inner and outer 
  double **specialc; // 1 - 2, 1 - 3 and 1 - 4 scaling

 protected:
  // cut offs for the taper
  double swa,swb;   // read on the input
  double swa2,swb2,swa3,swb3;
  // taper coefficients
  double Tap[8];
  double dTap[8];

  void init_taper();
  void coulomb(double,double,double *,double *);
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

E: Pair style coul/pqeqgauss requires atom attribute q

The atom style defined does not have these attributes.

E: Pair style coul/pqeqgauss requires atom attribute pqeq

The atom style defined does not have these attributes. Atom style has to be pqeq.

E: Pair inner cutoff >= Pair outer cutoff

The specified cutoffs for the pair style are inconsistent.

*/
