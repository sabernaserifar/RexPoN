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
#include "pair_exrydc.h"
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

PairEXRYDC::PairEXRYDC(LAMMPS *lmp) : Pair(lmp) 
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairEXRYDC::~PairEXRYDC()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(d0);
    memory->destroy(alpha);
    memory->destroy(r0);
    memory->destroy(L);
    memory->destroy(a0);
    memory->destroy(a1);
    memory->destroy(a2);
    memory->destroy(a3);
    memory->destroy(a4);
    memory->destroy(a5);
    memory->destroy(d0in);
    memory->destroy(a0in);
    memory->destroy(a1in);
    memory->destroy(a2in);
    memory->destroy(a3in);
    memory->destroy(a4in);
    memory->destroy(a5in);
    memory->destroy(flagatom);

  }
}

/* ---------------------------------------------------------------------- */
void PairEXRYDC::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,dr,dexp,factor_lj;
  double d_0,r_0,L_0,alpha_0,a_0,a_1,a_2,a_3,a_4,a_5;
  double d_0in,a_0in,a_1in,a_2in,a_3in,a_4in,a_5in;
  double epair_vdW,epair_inner,depair_vdW,depair_inner;
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
        d_0 = d0[itype][jtype];
        L_0 = L[itype][jtype];

        alpha_0 = alpha[itype][jtype];
        a_0 = a0[itype][jtype];
        a_1 = a1[itype][jtype];
        a_2 = a2[itype][jtype];
        a_3 = a3[itype][jtype];
        a_4 = a4[itype][jtype];
        a_5 = a5[itype][jtype];

        d_0in = d0in[itype][jtype];
        a_0in = a0in[itype][jtype];
        a_1in = a1in[itype][jtype];
        a_2in = a2in[itype][jtype];
        a_3in = a3in[itype][jtype];
        a_4in = a4in[itype][jtype];
        a_5in = a5in[itype][jtype];

        if ( (flagatom[itype][itype] + flagatom[jtype][jtype]) < -0.01 ) {
            d_0 = 0.0;
            d_0in = 0.0;
        }

        pair_interact(r,r_0,d_0,L_0,alpha_0,a_0,a_1,a_2,a_3,a_4,a_5,
                      &epair_vdW,&depair_vdW);

        pair_inner(r,d_0in,a_0in,a_1in,a_2in,a_3in,a_4in,a_5in,
                      &epair_inner,&depair_inner);

        // the current version does not include inner wall
        epair_inner = 0.0;
        depair_inner = 0.0;
        epair = epair_vdW + epair_inner;
        depair = depair_vdW + depair_inner;

        fpair = -factor_lj * depair / r;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        //if ( r > 1.6) epair = 0.0; 
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

void PairEXRYDC::allocate()
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
  memory->create(L,n+1,n+1,"pair:L");
  memory->create(a0,n+1,n+1,"pair:a0");
  memory->create(a1,n+1,n+1,"pair:a1");
  memory->create(a2,n+1,n+1,"pair:a2");
  memory->create(a3,n+1,n+1,"pair:a3");
  memory->create(a4,n+1,n+1,"pair:a4");
  memory->create(a5,n+1,n+1,"pair:a5");

  memory->create(d0in,n+1,n+1,"pair:d0in");
  memory->create(a0in,n+1,n+1,"pair:a0in");
  memory->create(a1in,n+1,n+1,"pair:a1in");
  memory->create(a2in,n+1,n+1,"pair:a2in");
  memory->create(a3in,n+1,n+1,"pair:a3in");
  memory->create(a4in,n+1,n+1,"pair:a4in");
  memory->create(a5in,n+1,n+1,"pair:a5in");
  memory->create(flagatom,n+1,n+1,"pair:flagatom");

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairEXRYDC::settings(int narg, char **arg)
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

void PairEXRYDC::coeff(int narg, char **arg)
{
  double z_one, z_two;

  if (narg < 18 )
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double d0_one = force->numeric(FLERR,arg[2]);
  double r0_one = force->numeric(FLERR,arg[3]);
  double L_one = force->numeric(FLERR,arg[4]);

  double alpha_one = force->numeric(FLERR,arg[5]);
  double a0_one = force->numeric(FLERR,arg[6]);
  double a1_one = force->numeric(FLERR,arg[7]);
  double a2_one = force->numeric(FLERR,arg[8]);
  double a3_one = force->numeric(FLERR,arg[9]);
  double a4_one = force->numeric(FLERR,arg[10]);
  double a5_one = force->numeric(FLERR,arg[11]);

  double d0in_one = force->numeric(FLERR,arg[12]);
  double a0in_one = force->numeric(FLERR,arg[13]);
  double a1in_one = force->numeric(FLERR,arg[14]);
  double a2in_one = force->numeric(FLERR,arg[15]);
  double a3in_one = force->numeric(FLERR,arg[16]);
  double a4in_one = force->numeric(FLERR,arg[17]);
  double a5in_one = force->numeric(FLERR,arg[18]);
  double flagatom_one = force->numeric(FLERR,arg[19]);


  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      d0[i][j] = d0_one;
      r0[i][j] = r0_one;
      L[i][j] = L_one;
      alpha[i][j] = alpha_one;
      a0[i][j] = a0_one;
      a1[i][j] = a1_one;
      a2[i][j] = a2_one;
      a3[i][j] = a3_one;
      a4[i][j] = a4_one;
      a5[i][j] = a5_one;
      d0in[i][j] = d0in_one;
      a0in[i][j] = a0in_one;
      a1in[i][j] = a1in_one;
      a2in[i][j] = a2in_one;
      a3in[i][j] = a3in_one;
      a4in[i][j] = a4in_one;
      a5in[i][j] = a5in_one;
      flagatom[i][j] = flagatom_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairEXRYDC::init_style()
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

double PairEXRYDC::init_one(int i, int j)
{
  /*if (setflag[i][j] == 0) { 
    d0[i][j] = mix_distance(d0[i][i],d0[j][j]);
    alpha[i][j] = mix_distance(alpha[i][i],alpha[j][j]);
    r0[i][j] = mix_distance(r0[i][i],r0[j][j]);
    a0[i][j] = mix_distance(a0[i][i],a0[j][j]);
    a1[i][j] = mix_distance(a1[i][i],a1[j][j]);
    a2[i][j] = mix_distance(a2[i][i],a2[j][j]);
    a3[i][j] = mix_distance(a3[i][i],a3[j][j]);
    a4[i][j] = mix_distance(a4[i][i],a4[j][j]);
    a5[i][j] = mix_distance(a5[i][i],a5[j][j]);
  }*/

  // force geometric or arithmetic average, need to check setflag 
  d0[i][j] = mix_distance(d0[i][i],d0[j][j]);
  r0[i][j] = mix_distance(r0[i][i],r0[j][j]);
  //d0[i][j] = 0.5*(d0[i][i]+d0[j][j]);
  //d0[i][j] = pow((0.5 * (pow(d0[i][i],6.0) + pow(d0[j][j],6.0))),1.0/6.0);
  //r0[i][j] = 0.5 * (r0[i][i]+r0[j][j]);
  L[i][j] = mix_distance(L[i][i],L[j][j]);
  alpha[i][j] = alpha[i][i];
  a0[i][j] = a0[i][i];
  a1[i][j] = a1[i][i];
  a2[i][j] = a2[i][i];
  a3[i][j] = a3[i][i];
  a4[i][j] = a4[i][i];
  a5[i][j] = a5[i][i];

  d0in[i][j] = 0.0;
  a0in[i][j] = 0.0;
  a1in[i][j] = 0.0;
  a2in[i][j] = 0.0;
  a3in[i][j] = 0.0;
  a4in[i][j] = 0.0;
  a5in[i][j] = 0.0;

  //error->all(FLERR,"All pair coeffs are not set");
  d0[j][i] = d0[i][j];
  r0[j][i] = r0[i][j];
  L[j][i] = L[i][j];
  alpha[j][i] = alpha[i][j];
  a0[j][i] = a0[i][j];
  a1[j][i] = a1[i][j];
  a2[j][i] = a2[i][j];
  a3[j][i] = a3[i][j];
  a4[j][i] = a4[i][j];
  a5[j][i] = a5[i][j];

  d0in[j][i] = d0in[i][j];
  a0in[j][i] = a0in[i][j];
  a1in[j][i] = a1in[i][j];
  a2in[j][i] = a2in[i][j];
  a3in[j][i] = a3in[i][j];
  a4in[j][i] = a4in[i][j];
  a5in[j][i] = a5in[i][j];
  //flagatom[j][i] = flagatom[i][j];

  return cut_outer;
}

/* ---------------------------------------------------------------------- */
void PairEXRYDC::init_taper()
{
  double d7;

  if (fabs(cut_inner) > 0.01 && comm->me == 0)
    error->warning(FLERR,"Pair EXRYDC has non-zero lower Taper radius cutoff");
  if (cut_outer < 0)
    error->all(FLERR, "Pair EXRYDC has negative upper Taper radius cutoff");
  else if (cut_outer < 5 && comm->me == 0)
    error->warning(FLERR,"Pair EXRYDC has very low Taper radius cutoff");

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
void PairEXRYDC::pair_interact(double r, double re, double de, double le,
                            double alf, double s0, double s1, double s2, double s3, 
                            double s4, double s5, double *Epair, double *dEpair)
{
  int n;
  double taper,dtaper;
  double dr,dr2,dr3,dr4,dr5;
  double exp1,dexp1,exp2,dexp2;
  double e_exryd,de_exryd;
  double L2,L3,L4,L5;

  // taper to smoothy go to zero between cut_inner and cut_outer
  if (r > cut_inner) {
    taper  =  tap[7];
    dtaper = dtap[7];
    for(n=6; n>=0; n--){
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

  L2 = le*le;
  L3 = le*L2;
  L4 = L2*L2;
  L5 = le*L4;

  alf = alf/le;
  s1 = s1/le;
  s2 = s2/L2;
  s3 = s3/L3;
  s4 = s4/L4;
  s5 = s5/L5;

  dr = r-re;
  dr2 = dr*dr;
  dr3 = dr*dr2;
  dr4 = dr*dr3;
  dr5 = dr*dr4;

  exp1 = -de * exp( -alf * dr);
  dexp1 = -alf * exp1;

  exp2 = s0 + s1*dr + s2*dr2 + s3*dr3 + s4*dr4 + s5*dr5 ;
  dexp2 = s1 + 2*s2*dr + 3*s3*dr2 + 4*s4*dr3 + 5*s5*dr4;
  
  e_exryd = exp1 * exp2;
  de_exryd = dexp1*exp2 + exp1*dexp2;

  *Epair = taper * e_exryd;
  *dEpair = taper * de_exryd + dtaper * e_exryd;

}

/* ---------------------------------------------------------------------- */
void PairEXRYDC::pair_inner(double r, double de, 
                            double s0, double s1, double s2, double s3,
                            double s4, double s5, double *Epair, double *dEpair)
{
  int n;
  double taper,dtaper;
  double r2,r3,r4,r5;
  double exp1,dexp1,exp2,dexp2;
  double e_inner,de_inner;

  // taper to smoothy go to zero between cut_inner and cut_outer
  if (r > cut_inner) {
    taper  =  tap[7];
    dtaper = dtap[7];
    for(n=6; n>=0; n--){
      taper  = taper * r + tap[n];
      dtaper = dtaper * r + dtap[n];
    }
  } else {
    taper = 1.0;
    dtaper = 0.0;
  }

  r2 = r*r;
  r3 = r*r2;
  r4 = r*r3;
  r5 = r*r4;

  exp1 = de;
  dexp1 = 0.0;

  exp2 = exp(-s0 - s1*r - s2*r2 - s3*r3 - s4*r4 - s5*r5) ;
  dexp2 = -(s1 + 2*s2*r + 3*s3*r2 + 4*s4*r3 + 5*s5*r4) * exp2;

  e_inner = exp1 * exp2;
  de_inner = dexp1*exp2 + exp1*dexp2;

  *Epair = taper * e_inner;
  *dEpair = taper * de_inner + dtaper * e_inner;

}
/* ---------------------------------------------------------------------- */

double PairEXRYDC::single(int i, int j, int itype, int jtype, double rsq,
                         double factor_coul, double factor_lj,
                         double &fforce)
{
  double r,dr,phi;
  double d_0,r_0,L_0,alpha_0,a_0,a_1,a_2,a_3,a_4,a_5;
  double d_0in,a_0in,a_1in,a_2in,a_3in,a_4in,a_5in;
  double epair_vdW,depair_vdW,epair_inner,depair_inner;
  double epair,depair;

  r = sqrt(rsq);
  r_0 = r0[itype][jtype];
  d_0 =  d0[itype][jtype];
  L_0 =  L[itype][jtype];
  alpha_0 = alpha[itype][jtype];
  a_0 = a0[itype][jtype];
  a_1 = a1[itype][jtype];
  a_2 = a2[itype][jtype];
  a_3 = a3[itype][jtype];
  a_4 = a4[itype][jtype];
  a_5 = a5[itype][jtype];

  d_0in = d0in[itype][jtype];
  a_0in = a0in[itype][jtype];
  a_1in = a1in[itype][jtype];
  a_2in = a2in[itype][jtype];
  a_3in = a3in[itype][jtype];
  a_4in = a4in[itype][jtype];
  a_5in = a5in[itype][jtype];


  if ( (flagatom[itype][itype] + flagatom[jtype][jtype]) < -0.01 ) {

     d_0 = 0.0;
     d_0in = 0.0;
  }

  pair_interact(r,r_0,d_0,L_0,alpha_0,a_0,a_1,a_2,a_3,a_4,a_5,
                &epair_vdW,&depair_vdW);

  pair_inner(r,d_0in,a_0in,a_1in,a_2in,a_3in,a_4in,a_5in,
                &epair_inner,&depair_inner);

  epair = epair_vdW + epair_inner;
  depair = depair_vdW + depair_inner;

  fforce = -factor_lj * depair / r;
  phi = epair;

  return factor_lj*phi;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairEXRYDC::write_restart(FILE *fp)
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
        fwrite(&L[i][j],sizeof(double),1,fp);
        fwrite(&a0[i][j],sizeof(double),1,fp);
        fwrite(&a1[i][j],sizeof(double),1,fp);
        fwrite(&a2[i][j],sizeof(double),1,fp);
        fwrite(&a3[i][j],sizeof(double),1,fp);
        fwrite(&a4[i][j],sizeof(double),1,fp);
        fwrite(&a5[i][j],sizeof(double),1,fp);

        fwrite(&d0in[i][j],sizeof(double),1,fp);
        fwrite(&a0in[i][j],sizeof(double),1,fp);
        fwrite(&a1in[i][j],sizeof(double),1,fp);
        fwrite(&a2in[i][j],sizeof(double),1,fp);
        fwrite(&a3in[i][j],sizeof(double),1,fp);
        fwrite(&a4in[i][j],sizeof(double),1,fp);
        fwrite(&a5in[i][j],sizeof(double),1,fp);
        fwrite(&flagatom[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairEXRYDC::read_restart(FILE *fp)
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
          fread(&L[i][j],sizeof(double),1,fp);
          fread(&a0[i][j],sizeof(double),1,fp);
          fread(&a1[i][j],sizeof(double),1,fp);
          fread(&a2[i][j],sizeof(double),1,fp);
          fread(&a3[i][j],sizeof(double),1,fp);
          fread(&a4[i][j],sizeof(double),1,fp);
          fread(&a5[i][j],sizeof(double),1,fp);

          fread(&d0in[i][j],sizeof(double),1,fp);
          fread(&a0in[i][j],sizeof(double),1,fp);
          fread(&a1in[i][j],sizeof(double),1,fp);
          fread(&a2in[i][j],sizeof(double),1,fp);
          fread(&a3in[i][j],sizeof(double),1,fp);
          fread(&a4in[i][j],sizeof(double),1,fp);
          fread(&a5in[i][j],sizeof(double),1,fp);
          fread(&flagatom[i][j],sizeof(double),1,fp);

        }
        MPI_Bcast(&d0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&alpha[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&r0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&L[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&a0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&a1[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&a2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&a3[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&a4[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&a5[i][j],1,MPI_DOUBLE,0,world);

        MPI_Bcast(&d0in[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&a0in[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&a1in[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&a2in[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&a3in[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&a4in[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&a5in[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&flagatom[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairEXRYDC::write_restart_settings(FILE *fp)
{
  fwrite(&cut_inner,sizeof(double),1,fp);
  fwrite(&cut_outer,sizeof(double),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairEXRYDC::read_restart_settings(FILE *fp)
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

void PairEXRYDC::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g \n",
             i,d0[i][i],alpha[i][i],r0[i][i],L[i][i],a0[i][i],a1[i][i],a2[i][i],
             a3[i][i],a4[i][i],a5[i][i],
             d0in[i][i],a0in[i][i],a1in[i][i],a2in[i][i],a3in[i][i],
             a4in[i][i],a5in[i][i],flagatom[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairEXRYDC::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
              i,j,d0[i][j],alpha[i][j],r0[i][j],L[i][j],a0[i][j],
              a1[i][j],a2[i][j],a3[i][j],a4[i][j],a5[i][j],
              d0in[i][j],a0in[i][j],a1in[i][j],a2in[i][j],a3in[i][j],
              a4in[i][j],a5in[i][j],flagatom[i][j]);
}
/* ---------------------------------------------------------------------- */

void *PairEXRYDC::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"d0") == 0) return (void *) d0;
  if (strcmp(str,"r0") == 0) return (void *) r0;
  if (strcmp(str,"L") == 0) return (void *) L;
  if (strcmp(str,"alpha") == 0) return (void *) alpha;
  if (strcmp(str,"a0") == 0) return (void *) a0;
  if (strcmp(str,"a1") == 0) return (void *) a1;
  if (strcmp(str,"a2") == 0) return (void *) a2;
  if (strcmp(str,"a3") == 0) return (void *) a3;
  if (strcmp(str,"a4") == 0) return (void *) a4;
  if (strcmp(str,"a5") == 0) return (void *) a5;

  if (strcmp(str,"d0in") == 0) return (void *) d0in;
  if (strcmp(str,"a0in") == 0) return (void *) a0in;
  if (strcmp(str,"a1in") == 0) return (void *) a1in;
  if (strcmp(str,"a2in") == 0) return (void *) a2in;
  if (strcmp(str,"a3in") == 0) return (void *) a3in;
  if (strcmp(str,"a4in") == 0) return (void *) a4in;
  if (strcmp(str,"a5in") == 0) return (void *) a5in;
  if (strcmp(str,"flagatom") == 0) return (void *) flagatom;

  return NULL;
}
