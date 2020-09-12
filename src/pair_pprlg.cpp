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
   Contributing authors: Saber Naserifar, naserifar.saber@gmail.com,
                         (Caltech 2018)
 ------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_pprlg.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairPPRLG::PairPPRLG(LAMMPS *lmp) : Pair(lmp) 
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairPPRLG::~PairPPRLG()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(d0);
    memory->destroy(alpha);
    memory->destroy(r0);
    memory->destroy(beta);
    memory->destroy(gamma);
    memory->destroy(eta);
    memory->destroy(delta);
    memory->destroy(c6lg);
    memory->destroy(r6lg);
  }
}

/* ---------------------------------------------------------------------- */
void PairPPRLG::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,r_0,dr,dexp,factor_lj;
  double de0,alpha0,beta0,gamma0,eta0,delta0;
  double c6lg0,r6lg0;
  double epair,depair;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;


  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {

      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cut_outer2) {
        r = sqrt(rsq);
        r_0 = r0[itype][jtype];
        de0 =  d0[itype][jtype];
        alpha0 = alpha[itype][jtype];
        beta0 = beta[itype][jtype];
        gamma0 = gamma[itype][jtype];
        eta0 = eta[itype][jtype];
        delta0 = delta[itype][jtype];
        c6lg0 = c6lg[itype][jtype];
        r6lg0 = r6lg[itype][jtype];

        pair_interact(r,r_0,de0,alpha0,beta0,gamma0,
                      eta0,delta0,c6lg0,r6lg0,
                      &epair,&depair);
        fpair = -factor_lj * depair / r;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        //if ( r < 1.5) epair = 0.0;
        
        if (eflag) {
          evdwl = epair;
          evdwl *= factor_lj;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairPPRLG::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq"); // needed for pair.ccp

  memory->create(d0,n+1,n+1,"pair:d0");
  memory->create(alpha,n+1,n+1,"pair:alpha");
  memory->create(r0,n+1,n+1,"pair:r0");
  memory->create(beta,n+1,n+1,"pair:beta");
  memory->create(gamma,n+1,n+1,"pair:gamma");
  memory->create(eta,n+1,n+1,"pair:eta");
  memory->create(delta,n+1,n+1,"pair:delta");
  memory->create(c6lg,n+1,n+1,"pair:c6lg");
  memory->create(r6lg,n+1,n+1,"pair:r6lg");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairPPRLG::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal pair_style command");

  cut_inner = force->numeric(FLERR,arg[0]);
  cut_outer = force->numeric(FLERR,arg[1]);

  if (cut_inner < 0.0 )
    error->all(FLERR,"Illegal pair_style command");
  if (cut_inner > cut_outer)
    error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs

------------------------------------------------------------------------- */

void PairPPRLG::coeff(int narg, char **arg)
{
  double z_one, z_two;

  if (narg < 11 )
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double d0_one = force->numeric(FLERR,arg[2]);
  double alpha_one = force->numeric(FLERR,arg[3]);
  double r0_one = force->numeric(FLERR,arg[4]);
  double b0_one = force->numeric(FLERR,arg[5]);
  double n0_one = force->numeric(FLERR,arg[6]);
  double c0_one = force->numeric(FLERR,arg[7]);
  double del0_one = force->numeric(FLERR,arg[8]);
  double r6lg_one = force->numeric(FLERR,arg[9]);
  double c6lg_one = force->numeric(FLERR,arg[10]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      d0[i][j] = d0_one;
      alpha[i][j] = alpha_one;
      r0[i][j] = r0_one;
      beta[i][j] = b0_one;
      gamma[i][j] = n0_one;
      eta[i][j] = c0_one;
      delta[i][j] = del0_one;
      c6lg[i][j] = c6lg_one;
      r6lg[i][j] = 2*r6lg_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairPPRLG::init_style()
{

  neighbor->request(this,instance_me);

  cut_inner2 = cut_inner * cut_inner;
  cut_outer2 = cut_outer * cut_outer;
  cut_outer3 = cut_outer2 * cut_outer;
  cut_inner3 = cut_inner2 * cut_inner;

  init_taper();
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairPPRLG::init_one(int i, int j)
{
  /*if (setflag[i][j] == 0) { 
    d0[i][j] = mix_distance(d0[i][i],d0[j][j]);
    alpha[i][j] = mix_distance(alpha[i][i],alpha[j][j]);
    r0[i][j] = mix_distance(r0[i][i],r0[j][j]);
    beta[i][j] = mix_distance(beta[i][i],beta[j][j]);
    gamma[i][j] = mix_distance(gamma[i][i],gamma[j][j]);
    eta[i][j] = mix_distance(eta[i][i],eta[j][j]);
    delta[i][j] = mix_distance(delta[i][i],delta[j][j]);
    c6lg[i][j] = mix_distance(c6lg[i][i],c6lg[j][j]);
    r6lg[i][j] = mix_distance(r6lg[i][i],r6lg[j][j]);
  }*/

  //error->all(FLERR,"All pair coeffs are not set");

    //c6lg[i][j] = mix_distance(c6lg[i][i],c6lg[j][j]);
    //r6lg[i][j] = mix_distance(r6lg[i][i],r6lg[j][j]);

  d0[j][i] = d0[i][j];
  alpha[j][i] = alpha[i][j];
  r0[j][i] = r0[i][j];
  beta[j][i] = beta[i][j];
  gamma[j][i] = gamma[i][j];
  eta[j][i] = eta[i][j];
  delta[j][i] = delta[i][j];
  c6lg[j][i] = c6lg[i][j];
  r6lg[j][i] = r6lg[i][j];


  return cut_outer;
}

/* ---------------------------------------------------------------------- */
void PairPPRLG::init_taper()
{
  double d7;

  if (fabs(cut_inner) > 0.01 && comm->me == 0)
    error->warning(FLERR,"Pair ppr has non-zero lower Taper radius cutoff");
  if (cut_outer < 0)
    error->all(FLERR, "Pair ppr has negative upper Taper radius cutoff");
  else if (cut_outer < 5 && comm->me == 0)
    error->warning(FLERR,"Pair ppr has very low Taper radius cutoff");

  d7 = pow( cut_outer - cut_inner, 7.0 );

  //taper
  tap[7] =  20.0 / d7;
  tap[6] = -70.0 * (cut_inner + cut_outer) / d7;
  tap[5] =  84.0 * (cut_inner2 + 3.0*cut_inner*cut_outer + cut_outer2) / d7;
  tap[4] = -35.0 * (cut_inner3 + 9.0*cut_inner2*cut_outer +
           9.0*cut_inner*cut_outer2 + cut_outer3 ) / d7;
  tap[3] = 140.0 * (cut_inner3*cut_outer + 3.0*cut_inner2*cut_outer2 +
           cut_inner*cut_outer3 ) / d7;
  tap[2] =-210.0 * (cut_inner3*cut_outer2 + cut_inner2*cut_outer3) / d7;
  tap[1] = 140.0 * cut_inner3 * cut_outer3 / d7;
  tap[0] = (-35.0*cut_inner3*cut_outer2*cut_outer2 +
            21.0*cut_inner2*cut_outer3*cut_outer2 +
            7.0*cut_inner*cut_outer3*cut_outer3 +
            cut_outer3*cut_outer3*cut_outer ) / d7;

  //taper derivative
  dtap[7] = 0.0;
  dtap[6] = 7.0*tap[7];
  dtap[5] = 6.0*tap[6];
  dtap[4] = 5.0*tap[5];
  dtap[3] = 4.0*tap[4];
  dtap[2] = 3.0*tap[3];
  dtap[1] = 2.0*tap[2];
  dtap[0] = tap[1];
}

/* ---------------------------------------------------------------------- */
void PairPPRLG::pair_interact(double r, double re, double d0, double a0, 
                            double b0, double g0, double e0, double del0, 
                            double clg, double rlg, double *Epair, double *dEpair)
{
  int n;
  double powrg,taper,dtaper;
  double dr,exp1,dexp1,exp2,dexp2,e_ppr,de_ppr;
  double r5,r6,rlg6,e_lg,de_lg;
  // taper to smoothy go to zero between cut_inner and cut_outer
  if (r > cut_inner) {
    taper  =  tap[7];
    dtaper = dtap[7];
    for(int n=6; n>=0; n--){
      taper  = taper * r + tap[n];
      dtaper = dtaper * r + dtap[n];
    }
  } else {
    taper = 1.0;
    dtaper = 0.0;
  }

  if ( r > cut_outer) {
    taper = 0.0;
    dtaper = 0.0;
  }

  dr = r-re;
  powrg = b0 * pow( r, g0 );
  exp1 = d0 * exp( -a0 * dr);
  dexp1 = -a0 * exp1;
  exp2 = exp( powrg + e0*r + del0);
  dexp2 = (g0*powrg/r+e0) * exp2;

  e_ppr = exp1 * exp2;
  de_ppr = dexp1*exp2 + exp1*dexp2;

  r5 = pow( r, 5.0 );
  r6 = pow( r, 6.0 );
  rlg6 = pow( rlg, 6.0 );
  e_lg = -clg/( r6 + rlg6 );
  de_lg = -6.0 * e_lg *  r5 / ( r6 + rlg6 ) ;

  *Epair = taper * ( e_ppr + e_lg );
  *dEpair = taper * ( de_ppr + de_lg ) + dtaper * (e_ppr + e_lg);

}

/* ---------------------------------------------------------------------- */

double PairPPRLG::single(int i, int j, int itype, int jtype, double rsq,
                         double factor_coul, double factor_lj,
                         double &fforce)
{
  double r,r_0,dr,phi,alf;
  double de0,alpha0,beta0,gamma0,eta0,delta0;
  double c6lg0,r6lg0;
  double epair,depair;

  r = sqrt(rsq);
  r_0 = r0[itype][jtype];
  de0 = d0[itype][jtype];
  alpha0 = alpha[itype][jtype];
  beta0 = beta[itype][jtype];
  gamma0 = gamma[itype][jtype];
  eta0 = eta[itype][jtype];
  delta0 = delta[itype][jtype];
  c6lg0 = c6lg[itype][jtype];
  r6lg0 = r6lg[itype][jtype];

  pair_interact(r,r_0,de0,alpha0,beta0,gamma0,eta0,delta0,c6lg0,r6lg0,
                &epair,&depair);

  fforce = -factor_lj * depair / r;

  phi = epair;
  return factor_lj*phi;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairPPRLG::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&d0[i][j],sizeof(double),1,fp);
        fwrite(&alpha[i][j],sizeof(double),1,fp);
        fwrite(&r0[i][j],sizeof(double),1,fp);
        fwrite(&beta[i][j],sizeof(double),1,fp);
        fwrite(&gamma[i][j],sizeof(double),1,fp);
        fwrite(&eta[i][j],sizeof(double),1,fp);
        fwrite(&delta[i][j],sizeof(double),1,fp);
        fwrite(&c6lg[i][j],sizeof(double),1,fp);
        fwrite(&r6lg[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairPPRLG::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&d0[i][j],sizeof(double),1,fp);
          fread(&alpha[i][j],sizeof(double),1,fp);
          fread(&r0[i][j],sizeof(double),1,fp);
          fread(&beta[i][j],sizeof(double),1,fp);
          fread(&gamma[i][j],sizeof(double),1,fp);
          fread(&eta[i][j],sizeof(double),1,fp);
          fread(&delta[i][j],sizeof(double),1,fp);
          fread(&c6lg[i][j],sizeof(double),1,fp);
          fread(&r6lg[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&d0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&alpha[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&r0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&beta[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&gamma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&eta[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&delta[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&c6lg[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&r6lg[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairPPRLG::write_restart_settings(FILE *fp)
{
  fwrite(&cut_inner,sizeof(double),1,fp);
  fwrite(&cut_outer,sizeof(double),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairPPRLG::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_inner,sizeof(double),1,fp);
    fread(&cut_outer,sizeof(double),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_inner,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_outer,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairPPRLG::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g %g %g %g %g %g\n",i,d0[i][i],alpha[i][i],r0[i][i], 
                beta[i][i],gamma[i][i],eta[i][i],delta[i][i],c6lg[i][i],r6lg[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairPPRLG::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g %g %g %g %g\n",
              i,j,d0[i][j],alpha[i][j],r0[i][j],beta[i][j],
              gamma[i][j],eta[i][j],delta[i][j],c6lg[i][j],r6lg[i][j]);
}
/* ---------------------------------------------------------------------- */

void *PairPPRLG::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"d0") == 0) return (void *) d0;
  if (strcmp(str,"r0") == 0) return (void *) r0;
  if (strcmp(str,"alpha") == 0) return (void *) alpha;
  if (strcmp(str,"beta") == 0) return (void *) beta;
  if (strcmp(str,"gamma") == 0) return (void *) gamma;
  if (strcmp(str,"eta") == 0) return (void *) eta;
  if (strcmp(str,"delta") == 0) return (void *) delta;
  if (strcmp(str,"c6lg") == 0) return (void *) c6lg;
  if (strcmp(str,"r6lg") == 0) return (void *) r6lg;
  return NULL;
}
