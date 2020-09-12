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
   Contributing author: Saber Naserifar, Caltech, naseri@caltech.edu
   (fix_qeq_reax.* was used as the starting point)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(pqeq,FixPqeq)

#else

#ifndef LMP_FIX_PQEQ_H
#define LMP_FIX_PQEQ_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPqeq : public Fix {
 public:
  FixPqeq(class LAMMPS *, int, char **);
  ~FixPqeq();
  int setmask();
  void init();
  void init_list(int,class NeighList *);
  void init_storage();
  void setup_pre_force(int);
  void pre_force(int);

  void setup_pre_force_respa(int, int);
  void pre_force_respa(int, int, int);

  void min_setup_pre_force(int);
  void min_pre_force(int);

  int matvecs;
  double pqeq_time;

 protected:
  int nevery, method;
  int n, N, m_fill;
  int n_cap, nmax, m_cap;
  int pack_flag;
  int nlevels_respa;
  class NeighList *list;
  class PairCoulPqeqgauss *pqeq;

  double swa, swb;      // lower/upper Taper cutoff radius
  double netchg;        // total charge of the system 
  double damp;          // damping factor for charges and shells  
  double tolerance;     // tolerance for the norm of the rel residual in CG
  double cqf;
  int ncycle;           // number of PQEq cycles
  int numCG;            // number of CG iter 

  bigint ngroup;

  double **alphass, **alphasc, **alphacc;
  double *chi, *idem, *rcore, *rshell, *qcore;
  double *kstring2, *kstring4, *Tap, *dTap;
  int *polar;
  double cut_inner, cut_outer;

  double cut_inner2, cut_outer2;
  double cut_inner3, cut_outer3;


  // fictitious charges

  double *s, *t;
  double **s_hist, **t_hist;
  int nprev;

  typedef struct{
    int n, m;
    int *firstnbr;
    int *numnbrs;
    int *jlist;
    double *val;
  } sparse_matrix;

  sparse_matrix H;
  double *Hdia_inv;
  double *b_s, *b_t;
  double *b_prc, *b_prm;
  
  double *qf, **sf;

  //CG storage
  double *p, *q, *r, *d;

  //GMRES storage
  //double *g,*y;
  //double **v;
  //double **h;
  //double *hc, *hs;

  void allocate_storage();
  void deallocate_storage();
  void reallocate_storage();
  void allocate_matrix();
  void deallocate_matrix();
  void reallocate_matrix();

  void init_matvec();
  void init_H();
  void compute_H();
  double calculate_H(double, double);
  void calculate_Q();
  void calculate_Shell();

  int CG(double*, double*);
  //int GMRES(double*,double*);
  void sparse_matvec(sparse_matrix*, double*, double*);

  void qforce();
  void sforce();
  void coulomb(double, double, double *, double *);


  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

  double parallel_norm( double*, int );
  double parallel_dot( double*, double*, int );
  double parallel_vector_acc( double*, int );

  void vector_sum(double*,double,double*,double,double*,int);
  void vector_add(double*, double, double*,int);
};

}

#endif
#endif
