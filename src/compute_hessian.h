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

#ifdef COMPUTE_CLASS

ComputeStyle(hessian, ComputeHessian)

#else

#ifndef LMP_COMPUTE_HESSIAN_H
#define LMP_COMPUTE_HESSIAN_H

#include "compute.h"

/* zero-based row- and column-major indexing macros for the hessian. */
#define idx2_r(i, j, ldj) ((i * ldj) + j)
#define idx2_c(i, j, ldi) ((j * ldi) + i)

namespace LAMMPS_NS {

  class ComputeHessian : public Compute {
  public:
    ComputeHessian(class LAMMPS *, int, char **);
    ~ComputeHessian();
    //void init() {}
    void compute_vector();
    void init();
    void init_list(int,class NeighList *);

  
  protected:
    int mylocalsize;
    int myglobalsize;
  
    double *fglobal_ref, *fglobal_new, *fglobal_copy;
    double *hessian;
  
    double epsilon, iepsilon, masselect;

    int pack_flag;
    int nlevels_respa;
    class NeighList *list;
    class PairCoulPqeqgauss *pqeq;

    double swa, swb;      // lower/upper Taper cutoff radius

    bigint ngroup;
    double **sf;
    double **alphass, **alphasc, **alphacc;
    double *chi, *idem, *rcore, *rshell, *qcore;
    double *kstring2, *kstring4, *Tap, *dTap;
    int *polar;
    double cut_inner, cut_outer;
    double cut_inner2, cut_outer2;
    double cut_inner3, cut_outer3;

    //int pair_compute_flag;
    //int kspace_compute_flag;
  
    void force_clear();
    void sforce();
    void coulomb(double, double, double *, double *);
    int pack_forward_comm(int, int *, double *, int, int *);
    void unpack_forward_comm(int, int, double *);
    int pack_reverse_comm(int, int, double *);
    void unpack_reverse_comm(int, int *, double *);

  };

}

#endif
#endif
