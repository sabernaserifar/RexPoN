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

/* ----------------------------------------------------------------------
   Contributing author: Saber Naserifar, Caltech, naseri@caltech.edu
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fix_pqeq.h"
#include "pair_coul_pqeqgauss.h"
#include "pair_rexpon.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "force.h"
#include "group.h"
#include "pair.h"
#include "respa.h"
#include "memory.h"
#include "citeme.h"
#include "error.h"
#include "math_const.h"
#include "rexpon_defs.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define EV_TO_KCAL_PER_MOL 14.399645
#define DANGER_ZONE     0.95
#define LOOSE_ZONE      0.7
#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define MIN_NBRS 100
#define MIN_CAP        200
#define SAFE_ZONE      2.0

static const char cite_fix_pqeq[] =
  "fix pqeq command:\n\n"
  "@Article{NaserifarPQEq,\n"
  " author = {S. Naserifar, D. J. Brooks, W. A. Goddard III, V. Cvicek},\n"
  " title = {Polarizable charge equilibration model for predicting accurate electrostatic interactions in molecules and solids},\n"
  " journal = {The Journal of Chemical Physics},\n"
  " year = {2017},\n"
  " volume = {146},\n"
  " pages = {124117}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

FixPqeq::FixPqeq(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_pqeq);

  if (narg < 2) error->all(FLERR,"Illegal fix pqeq command");

  // default values 
  method = 2;
  nevery = 1;
  netchg = 0.0;
  tolerance = 1.0e-6;
  damp = 1.0;
  ncycle = 2; 
  numCG = 200;

  // process optional keywords

  int iarg = 3;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"method") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix pqeq command");
      method = force->inumeric(FLERR,arg[iarg+1]);
      if ( method != 0 && method != 1 && method != 2)
         error->all(FLERR,"Illegal fix pqeq method");
      iarg += 2;
    } else if (strcmp(arg[iarg],"nevery") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix pqeq command");
      nevery = force->inumeric(FLERR,arg[iarg+1]);
      if (nevery <= 0) error->all(FLERR,"Illegal fix pqeq command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"charge") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix pqeq command");
      netchg = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"tolerance") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix pqeq command");
      tolerance = force->numeric(FLERR,arg[iarg+1]);
      if (tolerance < 0) error->all(FLERR,"Illegal fix pqeq command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"damp") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix pqeq command");
      damp = force->numeric(FLERR,arg[iarg+1]);
      if ( damp > 1.0 || damp < 0) error->all(FLERR,"Illegal fix pqeq command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"cycle") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix pqeq command");
      ncycle = force->numeric(FLERR,arg[iarg+1]);
      if ( ncycle < 0) error->all(FLERR,"Illegal fix pqeq command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"iter") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix pqeq command");
      numCG = force->numeric(FLERR,arg[iarg+1]);
      if ( numCG < 0) error->all(FLERR,"Illegal fix pqeq command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix pqeq command");
  }


  n = n_cap = 0;
  N = nmax = 0;
  m_fill = m_cap = 0;
  pack_flag = 0;
  s = NULL;
  t = NULL;
  nprev = 5;

  Hdia_inv = NULL;
  b_s = NULL;
  b_t = NULL;
  b_prc = NULL;
  b_prm = NULL;

  // CG
  p = NULL;
  q = NULL;
  r = NULL;
  d = NULL;

  // H matrix
  H.firstnbr = NULL;
  H.numnbrs = NULL;
  H.jlist = NULL;
  H.val = NULL;

  // shell
  qf = NULL;
  sf = NULL;
  cqf = 1.0;

  comm_forward = comm_reverse = 3;

  // perform initial allocation of atom-based arrays
  // register with Atom class

  s_hist = t_hist = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  for( int i = 0; i < atom->nmax; i++ )
    for (int j = 0; j < nprev; ++j )
      s_hist[i][j] = t_hist[i][j] = 0;

  pqeq = NULL;
  pqeq = (PairCoulPqeqgauss *) force->pair_match("pqeq",1);
}

/* ---------------------------------------------------------------------- */

FixPqeq::~FixPqeq()
{
  // unregister callbacks to this fix from Atom class

  if (copymode) return;

  atom->delete_callback(id,0);

  memory->destroy(s_hist);
  memory->destroy(t_hist);

  deallocate_storage();
  deallocate_matrix();
}

/* ---------------------------------------------------------------------- */

int FixPqeq::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= PRE_FORCE_RESPA;
  mask |= MIN_PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPqeq::allocate_storage()
{
  nmax = atom->nmax;

  memory->create(s,nmax,"pqeq:s");
  memory->create(t,nmax,"pqeq:t");

  memory->create(Hdia_inv,nmax,"pqeq:Hdia_inv");
  memory->create(b_s,nmax,"pqeq:b_s");
  memory->create(b_t,nmax,"pqeq:b_t");
  memory->create(b_prc,nmax,"pqeq:b_prc");
  memory->create(b_prm,nmax,"pqeq:b_prm");

  memory->create(p,nmax,"pqeq:p");
  memory->create(q,nmax,"pqeq:q");
  memory->create(r,nmax,"pqeq:r");
  memory->create(d,nmax,"pqeq:d");

  memory->create(qf,nmax,"pqeq:qf");
  memory->create(sf,nmax,3,"pqeq:sf");

}

/* ---------------------------------------------------------------------- */

void FixPqeq::deallocate_storage()
{
  memory->destroy(s);
  memory->destroy(t);

  memory->destroy( Hdia_inv );
  memory->destroy( b_s );
  memory->destroy( b_t );
  memory->destroy( b_prc );
  memory->destroy( b_prm );

  memory->destroy( p );
  memory->destroy( q );
  memory->destroy( r );
  memory->destroy( d );

  memory->destroy(qf);
  memory->destroy(sf);

}

/* ---------------------------------------------------------------------- */

void FixPqeq::reallocate_storage()
{
  deallocate_storage();
  allocate_storage();
  init_storage();
}

/* ---------------------------------------------------------------------- */

void FixPqeq::allocate_matrix()
{
  int i,ii,inum,m;
  int *ilist, *numneigh;

  int mincap;
  double safezone;

  /*if( reaxflag ) {
    mincap = reaxc->system->mincap;
    safezone = reaxc->system->safezone;
  } else {
    mincap = MIN_CAP;
    safezone = SAFE_ZONE;
  }*/


  mincap = MIN_CAP;
  safezone = SAFE_ZONE;

  n = atom->nlocal;
  n_cap = MAX( (int)(n * safezone), mincap );

  // determine the total space for the H matrix

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;

  m = 0;
  for( ii = 0; ii < inum; ii++ ) {
    i = ilist[ii];
    m += numneigh[i];
  }
  m_cap = MAX( (int)(m * safezone), mincap * MIN_NBRS );

  H.n = n_cap;
  H.m = m_cap;
  memory->create(H.firstnbr,n_cap,"pqeq:H.firstnbr");
  memory->create(H.numnbrs,n_cap,"pqeq:H.numnbrs");
  memory->create(H.jlist,m_cap,"pqeq:H.jlist");
  memory->create(H.val,m_cap,"pqeq:H.val");
}

/* ---------------------------------------------------------------------- */

void FixPqeq::deallocate_matrix()
{
  memory->destroy( H.firstnbr );
  memory->destroy( H.numnbrs );
  memory->destroy( H.jlist );
  memory->destroy( H.val );
}

/* ---------------------------------------------------------------------- */

void FixPqeq::reallocate_matrix()
{
  deallocate_matrix();
  allocate_matrix();
}

/* ---------------------------------------------------------------------- */

void FixPqeq::init()
{
  if (!atom->q_flag) error->all(FLERR,"Fix pqeq requires atom attribute q");

  pqeq = (PairCoulPqeqgauss *) force->pair_match("coul/pqeqgauss",1);
  if (pqeq == NULL) error->all(FLERR,"Cannot use fix pqeq without "
                  "pair_style coul/pqeqgauss");

  ngroup = group->count(igroup);
  if (ngroup == 0) error->all(FLERR,"Fix pqeq group has no atoms");

  // need a half neighbor list w/ Newton off and ghost neighbors
  // built whenever re-neighboring occurs

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->newton = 2;
  neighbor->requests[irequest]->ghost = 1;

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  int itmp;
  if (force->pair == NULL)
      error->all(FLERR,"fix pqeq is incompatible with Pair style");
  alphass = (double **) force->pair->extract("alphass",itmp);
  alphasc = (double **) force->pair->extract("alphasc",itmp);
  alphacc = (double **) force->pair->extract("alphacc",itmp);
  chi = (double *) force->pair->extract("chi",itmp);
  idem = (double *) force->pair->extract("idem",itmp);
  rcore = (double *) force->pair->extract("rcore",itmp);
  rshell = (double *) force->pair->extract("rshell",itmp);
  qcore = (double *) force->pair->extract("qcore",itmp);
  kstring2 = (double *) force->pair->extract("kstring2",itmp);
  kstring4 = (double *) force->pair->extract("kstring4",itmp);
  polar = (int *) force->pair->extract("polar",itmp);
  Tap = (double *) force->pair->extract("Tap",itmp);
  dTap = (double *) force->pair->extract("dTap",itmp);
  double *ptr_in = (double *) force->pair->extract("swa",itmp);
  double *ptr_out = (double *) force->pair->extract("swb",itmp);
  if (!alphass || !alphasc || !alphacc || !chi || !idem || !rcore 
      || !rshell || !qcore || !kstring2 || !kstring4 || !polar 
      || !Tap || !dTap || !ptr_in || !ptr_out)
      error->all(FLERR,"Fix pqeq is incompatible with Pair style");
  swa = *ptr_in;
  swb = *ptr_out;
}

/* ---------------------------------------------------------------------- */

void FixPqeq::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixPqeq::setup_pre_force(int vflag)
{
  // should not be needed
  // neighbor->build_one(list);

  deallocate_storage();
  allocate_storage();

  init_storage();

  deallocate_matrix();
  allocate_matrix();

  pre_force(vflag);

}

/* ---------------------------------------------------------------------- */

void FixPqeq::setup_pre_force_respa(int vflag, int ilevel)
{
  if (ilevel < nlevels_respa-1) return;
  setup_pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixPqeq::min_setup_pre_force(int vflag)
{
  setup_pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixPqeq::init_storage()
{
  int NN;

  qforce();

  NN = list->inum + list->gnum;

  for( int i = 0; i < NN; i++ ) {
    Hdia_inv[i] = 1. / idem[atom->type[i]];
    b_s[i] = -chi[atom->type[i]] +  EV_TO_KCAL_PER_MOL * qf[i];
    b_t[i] = -1.0;
    b_prc[i] = 0;
    b_prm[i] = 0;
    s[i] = t[i] = 0;
  }
}

/* ---------------------------------------------------------------------- */

void FixPqeq::pre_force(int vflag)
{
  double t_start, t_end;

  if (update->ntimestep % nevery) return;
  if( comm->me == 0 ) t_start = MPI_Wtime();

  n = atom->nlocal;
  N = atom->nlocal + atom->nghost;

  // grow arrays if necessary
  // need to be atom->nmax in length

  if( atom->nmax > nmax ) reallocate_storage();
  if( n > n_cap*DANGER_ZONE || m_fill > m_cap*DANGER_ZONE )
    reallocate_matrix();


  if ( method == 0 ) { // QEq with charge updates every ntimesteps
    cqf = 0.0;
    init_matvec();
    matvecs = CG(b_s, s);         // CG on s - parallel
    matvecs += CG(b_t, t);        // CG on t - parallel
    calculate_Q();
  }else if ( method == 1 ) { // only shell updates  
    sforce();               
    calculate_Shell();
  }else if ( method == 2 ) { // full PQEq
    cqf = 1.0; 
    for( int i = 0; i < ncycle; i++ ) {
      qforce();
      init_matvec();
      matvecs = CG(b_s, s);         
      matvecs += CG(b_t, t);        
      calculate_Q();
      sforce();               
      calculate_Shell();
    }
  }
 
  if( comm->me == 0 ) {
    t_end = MPI_Wtime();
    pqeq_time = t_end - t_start;
  }
}

/* ---------------------------------------------------------------------- */

void FixPqeq::pre_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixPqeq::min_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixPqeq::init_matvec()
{
  /* fill-in H matrix */
  compute_H();

  int nn, ii, i;
  int *ilist;

  nn = list->inum;
  ilist = list->ilist;

  for( ii = 0; ii < nn; ++ii ) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) {
      /* init pre-conditioner for H and init solution vectors */
      Hdia_inv[i] = 1. / idem[ atom->type[i] ];
      b_s[i]      = -chi[ atom->type[i] ] + cqf * EV_TO_KCAL_PER_MOL * qf[i];
      b_t[i]      = -1.0;

      /* linear extrapolation for s & t from previous solutions */
      //s[i] = 2 * s_hist[i][0] - s_hist[i][1];
      //t[i] = 2 * t_hist[i][0] - t_hist[i][1];

      /* quadratic extrapolation for s & t from previous solutions */
      //s[i] = s_hist[i][2] + 3 * ( s_hist[i][0] - s_hist[i][1] );
      t[i] = t_hist[i][2] + 3 * ( t_hist[i][0] - t_hist[i][1] );

      /* cubic extrapolation for s & t from previous solutions */
      s[i] = 4*(s_hist[i][0]+s_hist[i][2])-(6*s_hist[i][1]+s_hist[i][3]);
      //t[i] = 4*(t_hist[i][0]+t_hist[i][2])-(6*t_hist[i][1]+t_hist[i][3]);
    }
  }

  pack_flag = 2;
  comm->forward_comm_fix(this); //Dist_vector( s );
  pack_flag = 3;
  comm->forward_comm_fix(this); //Dist_vector( t );
}

/* ---------------------------------------------------------------------- */

void FixPqeq::compute_H()
{
  int inum, jnum, *ilist, *jlist, *numneigh, **firstneigh;
  int i, j, ii, jj, flag;
  double **x, SMALL = 0.0001;
  double dx, dy, dz, r_sqr;

  int *type = atom->type;
  tagint *tag = atom->tag;
  x = atom->x;
  int *mask = atom->mask;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // fill in the H matrix
  m_fill = 0;
  r_sqr = 0;
  for( ii = 0; ii < inum; ii++ ) {
    i = ilist[ii];
    if (mask[i] & groupbit) {
      jlist = firstneigh[i];
      jnum = numneigh[i];
      H.firstnbr[i] = m_fill;

      for( jj = 0; jj < jnum; jj++ ) {
        j = jlist[jj];
        j &= NEIGHMASK;

        dx = x[j][0] - x[i][0];
        dy = x[j][1] - x[i][1];
        dz = x[j][2] - x[i][2];
        r_sqr = SQR(dx) + SQR(dy) + SQR(dz);

        flag = 0;
        if (r_sqr <= SQR(swb)) {
          if (j < n) flag = 1;
          else if (tag[i] < tag[j]) flag = 1;
          else if (tag[i] == tag[j]) {
            if (dz > SMALL) flag = 1;
            else if (fabs(dz) < SMALL) {
              if (dy > SMALL) flag = 1;
              else if (fabs(dy) < SMALL && dx > SMALL)
                flag = 1;
	    }
	  }
	}

        if( flag ) {
          H.jlist[m_fill] = j;
          H.val[m_fill] = calculate_H( sqrt(r_sqr), alphacc[type[i]][type[j]] );
          m_fill++;
        }
      }
      H.numnbrs[i] = m_fill - H.firstnbr[i];
    }
  }

  if (m_fill >= H.m) {
    char str[128];
    sprintf(str,"H matrix size has been exceeded: m_fill=%d H.m=%d\n",
             m_fill, H.m );
    error->warning(FLERR,str);
    error->all(FLERR,"Fix pqeq has insufficient PQEq matrix size");
  }
}

/* ---------------------------------------------------------------------- */

double FixPqeq::calculate_H( double r, double alphacc )
{
  double Taper, coul;

  if (r > swa) {
    Taper  =  Tap[7];
    for(int n=6; n>=0; n--){
      Taper  = Taper * r + Tap[n];
    }
  } else Taper = 1.0;

  coul = erf(alphacc * r) / r;

  return Taper * EV_TO_KCAL_PER_MOL * coul;
}

/* ---------------------------------------------------------------------- */

int FixPqeq::CG( double *b, double *x )
{
  int  i, j, imax;
  double tmp, alpha, beta, b_norm;
  double sig_old, sig_new;

  int nn, jj;
  int *ilist;

  nn = list->inum;
  ilist = list->ilist;

  imax = numCG;

  pack_flag = 1;
  sparse_matvec( &H, x, q );
  comm->reverse_comm_fix( this ); //Coll_Vector( q );

  vector_sum( r , 1.,  b, -1., q, nn );

  for( jj = 0; jj < nn; ++jj ) {
    j = ilist[jj];
    if (atom->mask[j] & groupbit)
      d[j] = r[j] * Hdia_inv[j]; //pre-condition
  }

  b_norm = parallel_norm( b, nn );
  sig_new = parallel_dot( r, d, nn);

  for( i = 1; i < imax && sqrt(sig_new) / b_norm > tolerance; ++i ) {
    comm->forward_comm_fix(this); //Dist_vector( d );
    sparse_matvec( &H, d, q );
    comm->reverse_comm_fix(this); //Coll_vector( q );

    tmp = parallel_dot( d, q, nn);
    alpha = sig_new / tmp;

    vector_add( x, alpha, d, nn );
    vector_add( r, -alpha, q, nn );

    // pre-conditioning
    for( jj = 0; jj < nn; ++jj ) {
      j = ilist[jj];
      if (atom->mask[j] & groupbit)
        p[j] = r[j] * Hdia_inv[j];
    }

    sig_old = sig_new;
    sig_new = parallel_dot( r, p, nn);

    beta = sig_new / sig_old;
    vector_sum( d, 1., p, beta, d, nn );

  }

  if (i >= imax && comm->me == 0) {
    char str[128];
    sprintf(str,"Fix pqeq CG convergence failed after %d iterations "
            "at " BIGINT_FORMAT " step",i,update->ntimestep);
    error->warning(FLERR,str);
  }

  return i;
}


/* ---------------------------------------------------------------------- */

void FixPqeq::sparse_matvec( sparse_matrix *A, double *x, double *b )
{
  int i, j, itr_j;
  int nn, NN, ii;
  int *ilist;

  nn = list->inum;
  NN = list->inum + list->gnum;
  ilist = list->ilist;
  

  for( ii = 0; ii < nn; ++ii ) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit)
      b[i] = idem[ atom->type[i] ] * x[i];
  }

  for( ii = nn; ii < NN; ++ii ) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit)
      b[i] = 0;
  }

  for( ii = 0; ii < nn; ++ii ) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) {
      for( itr_j=A->firstnbr[i]; itr_j<A->firstnbr[i]+A->numnbrs[i]; itr_j++) {
        j = A->jlist[itr_j];
        b[i] += A->val[itr_j] * x[j];
        b[j] += A->val[itr_j] * x[i];
      }
    }
  }

}

/* ---------------------------------------------------------------------- */

void FixPqeq::calculate_Q()
{
  int i, k;
  double u, s_sum, t_sum;
  double q_new, q_old;
  double *q = atom->q;

  int nn, ii;
  int *ilist;

  nn = list->inum;
  ilist = list->ilist;

  s_sum = parallel_vector_acc( s, nn );
  t_sum = parallel_vector_acc( t, nn);
  u = (-netchg + s_sum) / t_sum;

  for( ii = 0; ii < nn; ++ii ) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) {
        q_old = q[i];
        //q[i] = s[i] - u * t[i];
        q_new = s[i] - u * t[i];
        q[i] = damp*q_new + (1.0-damp)*q_old; 
      /* backup s & t */
      for( k = 4; k > 0; --k ) {
        s_hist[i][k] = s_hist[i][k-1];
        t_hist[i][k] = t_hist[i][k-1];
      }
      s_hist[i][0] = s[i];
      t_hist[i][0] = t[i];
    }
  }

  pack_flag = 4;
  comm->forward_comm_fix( this ); //Dist_vector( atom->q );
}

/* ---------------------------------------------------------------------- */
void FixPqeq::calculate_Shell()
{
  double **rsx = atom->rsx;
  int *type = atom->type;
  double k2,r0,r1,r_new,r_old,delx,dely,delz;
  double rfx,rfy,rfz,rf;
  double rsx_old[3], rsx_new[3];

  int i, j, itype;
  int *mask = atom->mask;

  for (i = 0; i < n; i++) {

    if (mask[i] & groupbit) {
      itype = type[i];
      if (polar[itype]) {
        k2 = kstring2[itype];
        rfx = sf[i][0]; rfy = sf[i][1]; rfz = sf[i][2];
        rf = sqrt(rfx*rfx + rfy*rfy + rfz*rfz);
        //current position of the shell
        //delx = rsx[i][0]; dely = rsx[i][1]; delz = rsx[i][2];
        //r_old = sqrt(delx*delx + dely*dely + delz*delz);
        // calculate shell dist from the core using |-k2*r| = |(q+2)/q*rf|
        r1 =  rf/k2;
        
        if (r_new > 0.25){
            error->warning(FLERR,"Warning: core-shell dist more than 0.250 A");
        }  
        
        rsx_old[0] = rsx[i][0];
        rsx_old[1] = rsx[i][1];
        rsx_old[2] = rsx[i][2];

        if (rf > 1.e-10) {
           rsx_new[0] = r1 * rfx / rf;
           rsx_new[1] = r1 * rfy / rf;
           rsx_new[2] = r1 * rfz / rf;
        } else {
           rsx_new[0] = 0.0;
           rsx_new[1] = 0.0;
           rsx_new[2] = 0.0;
        }

        rsx[i][0] = damp * rsx_new[0] + (1.0 - damp) * rsx_old[0];
        rsx[i][1] = damp * rsx_new[1] + (1.0 - damp) * rsx_old[1];
        rsx[i][2] = damp * rsx_new[2] + (1.0 - damp) * rsx_old[2];
      }
    }
  }

  pack_flag = 5;
  comm->forward_comm_fix(this);
}

/* ---------------------------------------------------------------------- */
void FixPqeq::qforce()
{

  int i,j,ii,jj,inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int itype,jtype,polari,polarj;
  double ecoul,fpair,rsq,r,r2,factor_coul;
  double alpha,chii,idemi,kstring2i,kstring4i;
  double delx,dely,delz;
  double xi,yi,zi,rsxi,rsyi,rszi;
  double xj,yj,zj,rsxj,rsyj,rszj;
  double qi,qj,q1,q2,qcorei,qcorej,qshelli,qshellj;
  double shell_eng, coulombr, dcoulombr;

  double **x = atom->x;
  double *q = atom->q;
  double **rsx = atom->rsx;

  int *type = atom->type;
  int *tag = atom->tag;
  int flag;
  double SMALL = 0.0001;
  int nlocal = atom->nlocal;

  //double *special_coul = force->special_coul;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;
  double qe2f = force->qe2f;
  int *mask = atom->mask;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      qf[i] = 0.0;
  }

  // charge and shell forces, first forward then reverse

  pack_flag = 6;
  comm->forward_comm_fix(this);

  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
        qi = q[i];
        xi = x[i][0];
        yi = x[i][1];
        zi = x[i][2];
        itype = type[i];

        jlist = firstneigh[i];
        jnum = numneigh[i];
        polari = polar[itype];

        if (polari) {
          qcorei = qcore[itype];
          qshelli = -qcorei;
          rsxi = rsx[i][0];
          rsyi = rsx[i][1];
          rszi = rsx[i][2];
        } else qcorei = 0.0; 

        q1 = qcorei + qi;

        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          j &= NEIGHMASK;
          if (!(mask[j] & groupbit)) continue;

          factor_coul = 1.0; // we don't weight the charge and shell forces
          jtype = type[j];

              xj = x[j][0];
              yj = x[j][1];
              zj = x[j][2];

              delx = xi - xj;
              dely = yi - yj;
              delz = zi - zj;

          r2 = delx*delx + dely*dely + delz*delz;

          flag = 0;
          if (r2 <= SQR(swb)) {
            if (j < nlocal) flag = 1;
            else if (tag[i] < tag[j]) flag = 1;
            else if (tag[i] == tag[j]) {
              if (-delz > SMALL) flag = 1;
              else if (fabs(delz) < SMALL) {
                if (-dely > SMALL) flag = 1;
                else if (fabs(dely) < SMALL && -delx > SMALL)
                flag = 1;
              }
            }
          }

          if (flag) {
                qj = q[j];
                jtype = type[j];
                polarj = polar[jtype];

                if (polarj) {
                 qcorej = qcore[jtype];
                 qshellj = -qcore[jtype];
                 rsxj = rsx[j][0];
                 rsyj = rsx[j][1];
                 rszj = rsx[j][2];
                } else qcorej = 0.0;

                q2 = qcorej + qj;

                // core - core
                alpha = alphacc[itype][jtype];
                r = sqrt(r2);
                coulomb(r, alpha , &coulombr, &dcoulombr);
                // force on qi, it doesn't include qcorei which is for dipole
                qf[i] -= coulombr * qcorej/qqrd2e;
                if (newton_pair || j < nlocal) {
                  qf[j] -= coulombr * qcorei/qqrd2e;
                }

                // shell - shell
                // nothing to store 

                // shell i - core j
                if (polari) {
                  alpha = alphasc[itype][jtype];
                  delx = (rsxi + xi) - xj;
                  dely = (rsyi + yi) - yj;
                  delz = (rszi + zi) - zj;
                  r2 = delx*delx + dely*dely + delz*delz;
                  r = sqrt(r2);
                  coulomb(r, alpha , &coulombr, &dcoulombr);
                  if (newton_pair || j < nlocal) {
                    qf[j] -= coulombr * qshelli/qqrd2e;
                  }
                }

                // core i - shell j
                if (polarj) {
                  alpha = alphasc[itype][jtype];
                  delx = xi - (rsxj + xj);
                  dely = yi - (rsyj + yj);
                  delz = zi - (rszj + zj);
                  r2 = delx*delx + dely*dely + delz*delz;
                  r = sqrt(r2);
                  coulomb(r, alpha , &coulombr, &dcoulombr);
                  qf[i] -= coulombr * qshellj/qqrd2e;
                }

              }
          
        }
    
  }

  pack_flag = 7;
  comm->reverse_comm_fix( this );
}
/* ---------------------------------------------------------------------- */
void FixPqeq::sforce()
{

  int i,j,ii,jj,inum,jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int itype,jtype,polari,polarj;
  double ecoul,fpair,rsq,r,r2,factor_coul;
  double alpha,chii,idemi,kstring2i,kstring4i;
  double delx,dely,delz;
  double xi,yi,zi,rsxi,rsyi,rszi;
  double xj,yj,zj,rsxj,rsyj,rszj;
  double qi,qj,q1,q2,qcorei,qcorej,qshelli,qshellj;
  double shell_eng, coulombr, dcoulombr;

  double **x = atom->x;
  double *q = atom->q;
  double **rsx = atom->rsx;

  int *type = atom->type;
  int *tag = atom->tag;
  int flag;
  double SMALL = 0.0001;
  int nlocal = atom->nlocal;

  //double *special_coul = force->special_coul;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;
  double qe2f = force->qe2f;
  int *mask = atom->mask;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      sf[i][0] = 0.0;
      sf[i][1] = 0.0;
      sf[i][2] = 0.0;
  }

  // charge and shell forces, first forward then reverse

  pack_flag = 8;
  comm->forward_comm_fix(this);

  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
        qi = q[i];
        xi = x[i][0];
        yi = x[i][1];
        zi = x[i][2];
        itype = type[i];

        jlist = firstneigh[i];
        jnum = numneigh[i];
        polari = polar[itype];

        if (polari) {
          qcorei = qcore[itype];
          qshelli = -qcorei;
          rsxi = rsx[i][0];
          rsyi = rsx[i][1];
          rszi = rsx[i][2];
        } else qcorei = 0.0; 

        q1 = qcorei + qi;

        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          j &= NEIGHMASK;
          if (!(mask[j] & groupbit)) continue;

          factor_coul = 1.0; // we don't weight the charge and shell forces
          jtype = type[j];

              xj = x[j][0];
              yj = x[j][1];
              zj = x[j][2];

              delx = xi - xj;
              dely = yi - yj;
              delz = zi - zj;

          r2 = delx*delx + dely*dely + delz*delz;

          flag = 0;
          if (r2 <= SQR(swb)) {
            if (j < nlocal) flag = 1;
            else if (tag[i] < tag[j]) flag = 1;
            else if (tag[i] == tag[j]) {
              if (-delz > SMALL) flag = 1;
              else if (fabs(delz) < SMALL) {
                if (-dely > SMALL) flag = 1;
                else if (fabs(dely) < SMALL && -delx > SMALL)
                flag = 1;
              }
            }
          }

          if (flag) {
                qj = q[j];
                jtype = type[j];
                polarj = polar[jtype];

                if (polarj) {
                 qcorej = qcore[jtype];
                 qshellj = -qcore[jtype];
                 rsxj = rsx[j][0];
                 rsyj = rsx[j][1];
                 rszj = rsx[j][2];
                } else qcorej = 0.0;

                q2 = qcorej + qj;

               // core - core
               // nothing to store 
                
               // shell - shell
                if (polari and polarj) {
                  alpha = alphass[itype][jtype];
                  delx = rsxi - rsxj + delx;
                  dely = rsyi - rsyj + dely;
                  delz = rszi - rszj + delz;
                  r2 = delx*delx + dely*dely + delz*delz;
                  r = sqrt(r2);
                  coulomb(r, alpha , &coulombr, &dcoulombr);
                  fpair = - dcoulombr * qshelli * qshellj / r;
                  sf[i][0] += factor_coul*delx*fpair;
                  sf[i][1] += factor_coul*dely*fpair;
                  sf[i][2] += factor_coul*delz*fpair;
                  if (newton_pair || j < nlocal) {
                    sf[j][0] -= factor_coul*delx*fpair;
                    sf[j][1] -= factor_coul*dely*fpair;
                    sf[j][2] -= factor_coul*delz*fpair;
                  }
                }

                // shell i - core j
                if (polari) {
                  alpha = alphasc[itype][jtype];
                  delx = (rsxi + xi) - xj;
                  dely = (rsyi + yi) - yj;
                  delz = (rszi + zi) - zj;
                  r2 = delx*delx + dely*dely + delz*delz;
                  r = sqrt(r2);
                  coulomb(r, alpha , &coulombr, &dcoulombr);
                  fpair = - dcoulombr * qshelli * q2 / r;
                  sf[i][0] += factor_coul*delx*fpair;
                  sf[i][1] += factor_coul*dely*fpair;
                  sf[i][2] += factor_coul*delz*fpair;
                }

                // core i - shell j
                if (polarj) {
                  alpha = alphasc[itype][jtype];
                  delx = xi - (rsxj + xj);
                  dely = yi - (rsyj + yj);
                  delz = zi - (rszj + zj);
                  r2 = delx*delx + dely*dely + delz*delz;
                  r = sqrt(r2);
                  coulomb(r, alpha , &coulombr, &dcoulombr);
                  fpair = - dcoulombr * q1 * qshellj / r;
                  if (newton_pair || j < nlocal) {
                    sf[j][0] -= factor_coul*delx*fpair;
                    sf[j][1] -= factor_coul*dely*fpair;
                    sf[j][2] -= factor_coul*delz*fpair;
                  }
                }

              }
          
        }
    
  }

  pack_flag = 9;
  comm->reverse_comm_fix( this );
}
/* -----------------------------------------------------------------------
   Coulomb 
-------------------------------------------------------------------------- */

void FixPqeq::coulomb(double r, double alpha, double *coulombr, double *dcoulombr)
{
  int n;
  double taper,dtaper,screening,dscreening,coul,dcoul,r2;
  double qqrd2e = force->qqrd2e;

  r2 = r * r;
  // taper to smoothy go to zero between swa and swb
  if (r > swa) {
    taper  =  Tap[7];
    dtaper = dTap[7];
    for(int n=6; n>=0; n--){
      taper  = taper * r + Tap[n];
      dtaper  = dtaper * r + dTap[n];
    }
  } else {
    taper = 1.0;
    dtaper = 0.0;
  }

  if ( r > swb) {
    taper = 0.0;
    dtaper = 0.0; 
  }

  coul = qqrd2e / r;
  dcoul = - qqrd2e / r2;

  screening = erf(alpha * r);
  dscreening = 2.0 * alpha * exp(- alpha * alpha * r2) / MY_PIS;
  *coulombr = coul * taper * screening;
  *dcoulombr = dcoul * taper * screening + coul * dtaper * screening +
          coul * taper * dscreening;
}
/* -----------------------------------------------------------------------*/

int FixPqeq::pack_forward_comm(int n, int *list, double *buf,
                          int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;

  if ( pack_flag == 1)
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = d[j];
    }
  else if ( pack_flag == 2)
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = s[j];
    }
  else if ( pack_flag == 3)
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = t[j];
    }
  else if ( pack_flag == 4)
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = atom->q[j];
    }
  else if ( pack_flag == 5)
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = atom->rsx[j][0];
      buf[m++] = atom->rsx[j][1];
      buf[m++] = atom->rsx[j][2];
    }
  else if ( pack_flag == 6)
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = qf[j];
    }
  else if ( pack_flag == 8)
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = sf[j][0];
      buf[m++] = sf[j][1];
      buf[m++] = sf[j][2];
    }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixPqeq::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n ;

  if( pack_flag == 1)
    for (i = first; i < last; i++) {
      d[i] = buf[m++];
    } 
  else if( pack_flag == 2)
    for (i = first; i < last; i++) {
      s[i] = buf[m++];
    }
  else if( pack_flag == 3)
    for (i = first; i < last; i++) {
      t[i] = buf[m++];
    }
  else if( pack_flag == 4)
    for (i = first; i < last; i++) {
      atom->q[i] = buf[m++];
    }
  else if( pack_flag == 5)
    for (i = first; i < last; i++) {
      atom->rsx[i][0] = buf[m++];
      atom->rsx[i][1] = buf[m++];
      atom->rsx[i][2] = buf[m++];
    }
  else if( pack_flag == 6)
    for (i = first; i < last; i++) {
      qf[i] = buf[m++];
    }
  else if( pack_flag == 8)
    for (i = first; i < last; i++) {
      sf[i][0] = buf[m++];
      sf[i][1] = buf[m++];
      sf[i][2] = buf[m++];
    }

}

/* ---------------------------------------------------------------------- */

int FixPqeq::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  
  if (pack_flag == 7) 
    for (i = first; i < last; i++) {
      buf[m++] = qf[i];
    }
  else if (pack_flag == 9) 
    for (i = first; i < last; i++) {
      buf[m++] = sf[i][0];
      buf[m++] = sf[i][1];
      buf[m++] = sf[i][2];
    }
  else 
    for (i = first; i < last; i++) 
      buf[m++] = q[i];

  return m;
}

/* ---------------------------------------------------------------------- */

void FixPqeq::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  if (pack_flag == 7) 
    for (i = 0; i < n; i++) {
      j = list[i];
      qf[j] += buf[m++];
    }
  else if (pack_flag == 9) 
    for (i = 0; i < n; i++) {
      j = list[i];
      sf[j][0] += buf[m++];
      sf[j][1] += buf[m++];
      sf[j][2] += buf[m++];
    }
  else
    for (i = 0; i < n; i++) {
      j = list[i];
      q[j] += buf[m++];
    }

}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixPqeq::memory_usage()
{
  double bytes;

  bytes = atom->nmax*nprev*2 * sizeof(double); // s_hist & t_hist
  bytes += atom->nmax*13 * sizeof(double); // storage
  bytes += n_cap*2 * sizeof(int); // matrix...
  bytes += m_cap * sizeof(int);
  bytes += m_cap * sizeof(double);

  return bytes;
}

/* ----------------------------------------------------------------------
   allocate fictitious charge arrays
------------------------------------------------------------------------- */

void FixPqeq::grow_arrays(int nmax)
{
  memory->grow(s_hist,nmax,nprev,"pqeq:s_hist");
  memory->grow(t_hist,nmax,nprev,"pqeq:t_hist");
}

/* ----------------------------------------------------------------------
   copy values within fictitious charge arrays
------------------------------------------------------------------------- */

void FixPqeq::copy_arrays(int i, int j, int delflag)
{
  for (int m = 0; m < nprev; m++) {
    s_hist[j][m] = s_hist[i][m];
    t_hist[j][m] = t_hist[i][m];
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixPqeq::pack_exchange(int i, double *buf)
{
  for (int m = 0; m < nprev; m++) buf[m] = s_hist[i][m];
  for (int m = 0; m < nprev; m++) buf[nprev+m] = t_hist[i][m];
  return nprev*2;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixPqeq::unpack_exchange(int nlocal, double *buf)
{
  for (int m = 0; m < nprev; m++) s_hist[nlocal][m] = buf[m];
  for (int m = 0; m < nprev; m++) t_hist[nlocal][m] = buf[nprev+m];
  return nprev*2;
}

/* ---------------------------------------------------------------------- */

double FixPqeq::parallel_norm( double *v, int n )
{
  int  i;
  double my_sum, norm_sqr;

  int ii;
  int *ilist;

  ilist = list->ilist;

  my_sum = 0.0;
  norm_sqr = 0.0;
  for( ii = 0; ii < n; ++ii ) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit)
      my_sum += SQR( v[i] );
  }

  MPI_Allreduce( &my_sum, &norm_sqr, 1, MPI_DOUBLE, MPI_SUM, world );

  return sqrt( norm_sqr );
}

/* ---------------------------------------------------------------------- */

double FixPqeq::parallel_dot( double *v1, double *v2, int n)
{
  int  i;
  double my_dot, res;

  int ii;
  int *ilist;

  ilist = list->ilist;

  my_dot = 0.0;
  res = 0.0;
  for( ii = 0; ii < n; ++ii ) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit)
      my_dot += v1[i] * v2[i];
  }

  MPI_Allreduce( &my_dot, &res, 1, MPI_DOUBLE, MPI_SUM, world );

  return res;
}

/* ---------------------------------------------------------------------- */

double FixPqeq::parallel_vector_acc( double *v, int n )
{
  int  i;
  double my_acc, res;

  int ii;
  int *ilist;

  ilist = list->ilist;

  my_acc = 0.0;
  res = 0.0;
  for( ii = 0; ii < n; ++ii ) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) { 
        my_acc += v[i];
    }
  }

  MPI_Allreduce( &my_acc, &res, 1, MPI_DOUBLE, MPI_SUM, world );

  return res;
}

/* ---------------------------------------------------------------------- */

void FixPqeq::vector_sum( double* dest, double c, double* v,
                                double d, double* y, int k )
{
  int kk;
  int *ilist;

  ilist = list->ilist;

  for( --k; k>=0; --k ) {
    kk = ilist[k];
    if (atom->mask[kk] & groupbit)
      dest[kk] = c * v[kk] + d * y[kk];
  }
}

/* ---------------------------------------------------------------------- */

void FixPqeq::vector_add( double* dest, double c, double* v, int k )
{
  int kk;
  int *ilist;

  ilist = list->ilist;

  for( --k; k>=0; --k ) {
    kk = ilist[k];
    if (atom->mask[kk] & groupbit)
      dest[kk] += c * v[kk];
  }

}
