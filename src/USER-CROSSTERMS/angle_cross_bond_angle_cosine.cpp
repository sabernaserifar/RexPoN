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
   Contributing author: Eric Simon (Cray)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "angle_cross_bond_angle_cosine.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

#include "math_const.h";

using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL 0.001

/* ---------------------------------------------------------------------- */

AngleCrossBondAngleCosine::AngleCrossBondAngleCosine(LAMMPS *lmp) : Angle(lmp) {}

/* ---------------------------------------------------------------------- */

AngleCrossBondAngleCosine::~AngleCrossBondAngleCosine()
{
  if (allocated) {
    memory->destroy(setflag);

    memory->destroy(theta0);
    memory->destroy(ba_k1);
    memory->destroy(ba_k2);
    memory->destroy(ba_r1);
    memory->destroy(ba_r2);
  }
}

/* ---------------------------------------------------------------------- */

void AngleCrossBondAngleCosine::compute(int eflag, int vflag)
{
  int i1,i2,i3,n,type;
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double eangle,f1[3],f3[3];
  double dtheta,dtheta2,dtheta3,dtheta4,de_angle;
  double dr1,dr2,tk1,tk2,aa1,aa2,aa11,aa12,aa21,aa22;
  double rsq1,rsq2,r1,r2,c,s,a,a11,a12,a22,b1,b2;
  double vx11,vx12,vy11,vy12,vz11,vz12,vx21,vx22,vy21,vy22,vz21,vz22;

  eangle = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nanglelist; n++) {
    i1 = anglelist[n][0];
    i2 = anglelist[n][1];
    i3 = anglelist[n][2];
    type = anglelist[n][3];

    // 1st bond

    delx1 = x[i1][0] - x[i2][0];
    dely1 = x[i1][1] - x[i2][1];
    delz1 = x[i1][2] - x[i2][2];
    domain->minimum_image(delx1,dely1,delz1);

    rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    r1 = sqrt(rsq1);

    // 2nd bond

    delx2 = x[i3][0] - x[i2][0];
    dely2 = x[i3][1] - x[i2][1];
    delz2 = x[i3][2] - x[i2][2];
    domain->minimum_image(delx2,dely2,delz2);

    rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    r2 = sqrt(rsq2);

    // angle (cos and sin)

    c = delx1*delx2 + dely1*dely2 + delz1*delz2;
    c /= r1*r2;
        
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
        
    s = sqrt(1.0 - c*c);
    if (s < SMALL) s = SMALL;
    s = 1.0/s;

    // force & energy for bond-angle term

    dr1 = r1 - ba_r1[type];
    dr2 = r2 - ba_r2[type];
    dtheta = c - theta0[type];
    aa1 = -1 * dr1 * ba_k1[type];
    aa2 = -1 * dr2 * ba_k2[type];

    aa11 = aa1 * c / rsq1;
    aa12 = -aa1 / (r1 * r2);
    aa21 = aa2 * c / rsq1;
    aa22 = -aa2 / (r1 * r2);

    vx11 = (aa11 * delx1) + (aa12 * delx2);
    vx12 = (aa21 * delx1) + (aa22 * delx2);
    vy11 = (aa11 * dely1) + (aa12 * dely2);
    vy12 = (aa21 * dely1) + (aa22 * dely2);
    vz11 = (aa11 * delz1) + (aa12 * delz2);
    vz12 = (aa21 * delz1) + (aa22 * delz2);

    aa11 = aa1 * c / rsq2;
    aa21 = aa2 * c / rsq2;

    vx21 = (aa11 * delx2) + (aa12 * delx1);
    vx22 = (aa21 * delx2) + (aa22 * delx1);
    vy21 = (aa11 * dely2) + (aa12 * dely1);
    vy22 = (aa21 * dely2) + (aa22 * dely1);
    vz21 = (aa11 * delz2) + (aa12 * delz1);
    vz22 = (aa21 * delz2) + (aa22 * delz1);

    b1 = ba_k1[type] * dtheta / r1;
    b2 = ba_k2[type] * dtheta / r2;

    f1[0] = -1 * (vx11 + b1*delx1 + vx12);
    f1[1] = -1 * (vy11 + b1*dely1 + vy12);
    f1[2] = -1 * (vz11 + b1*delz1 + vz12);

    f3[0] = -1 * (vx21 + b2*delx2 + vx22);
    f3[1] = -1 * (vy21 + b2*dely2 + vy22);
    f3[2] = -1 * (vz21 + b2*delz2 + vz22);

    if (eflag) eangle = ba_k1[type]*dr1*dtheta + ba_k2[type]*dr2*dtheta;
    //printf("c: %lf theta0 %lf dtheta: %lf, dr1: %lf dr2: %lf energy: %lf\n",c,theta0[type],dtheta,dr1,dr2,eangle);
    // apply force to each of 3 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= f1[0] + f3[0];
      f[i2][1] -= f1[1] + f3[1];
      f[i2][2] -= f1[2] + f3[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    if (evflag) ev_tally(i1,i2,i3,nlocal,newton_bond,eangle,f1,f3,
			 delx1,dely1,delz1,delx2,dely2,delz2);
  }
}

/* ---------------------------------------------------------------------- */

void AngleCrossBondAngleCosine::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;
  memory->create(setflag,n+1,"angle:setflag");

  memory->create(theta0,n+1,"angle:theta0");

  memory->create(ba_k1,n+1,"angle:ba_k1");
  memory->create(ba_k2,n+1,"angle:ba_k2");
  memory->create(ba_r1,n+1,"angle:ba_r1");
  memory->create(ba_r2,n+1,"angle:ba_r2");

  for (int i = 1; i <= n; i++)
    setflag[i] = 0;
}

/* ---------------------------------------------------------------------- */

void AngleCrossBondAngleCosine::coeff(int narg, char **arg)
{
  if (narg != 6) error->all(FLERR,"Incorrect args for angle coefficients");

  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->nangletypes,ilo,ihi);

  int count = 0;


  double ba_r1_one = atof(arg[1]);
  double ba_r2_one = atof(arg[2]);
  double theta0_one = atof(arg[3]) * MY_PI/180;
  double ba_k1_one = atof(arg[4]);
  double ba_k2_one = atof(arg[5]);
  
  for (int i = ilo; i <= ihi; i++) {
    ba_k1[i] = ba_k1_one;
    ba_k2[i] = ba_k2_one;
    ba_r1[i] = ba_r1_one;
    ba_r2[i] = ba_r2_one;
    theta0[i] = cos(theta0_one);
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for angle coefficients");

}

/* ---------------------------------------------------------------------- */

double AngleCrossBondAngleCosine::equilibrium_angle(int i)
{
  return theta0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleCrossBondAngleCosine::write_restart(FILE *fp)
{
  fwrite(&theta0[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&ba_k1[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&ba_k2[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&ba_r1[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&ba_r2[1],sizeof(double),atom->nangletypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them 
------------------------------------------------------------------------- */

void AngleCrossBondAngleCosine::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&theta0[1],sizeof(double),atom->nangletypes,fp);
    fread(&ba_k1[1],sizeof(double),atom->nangletypes,fp);
    fread(&ba_k2[1],sizeof(double),atom->nangletypes,fp);
    fread(&ba_r1[1],sizeof(double),atom->nangletypes,fp);
    fread(&ba_r2[1],sizeof(double),atom->nangletypes,fp);
  }

  MPI_Bcast(&theta0[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&ba_k1[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&ba_k2[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&ba_r1[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&ba_r2[1],atom->nangletypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ---------------------------------------------------------------------- */

double AngleCrossBondAngleCosine::single(int type, int i1, int i2, int i3)
{
  double **x = atom->x;

  double delx1 = x[i1][0] - x[i2][0];
  double dely1 = x[i1][1] - x[i2][1];
  double delz1 = x[i1][2] - x[i2][2];
  domain->minimum_image(delx1,dely1,delz1);
  double r1 = sqrt(delx1*delx1 + dely1*dely1 + delz1*delz1);

  double delx2 = x[i3][0] - x[i2][0];
  double dely2 = x[i3][1] - x[i2][1];
  double delz2 = x[i3][2] - x[i2][2];
  domain->minimum_image(delx2,dely2,delz2);
  double r2 = sqrt(delx2*delx2 + dely2*dely2 + delz2*delz2);

  double c = delx1*delx2 + dely1*dely2 + delz1*delz2;
  c /= r1*r2;
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;
        
  double s = sqrt(1.0 - c*c);
  if (s < SMALL) s = SMALL;
  s = 1.0/s;

  double dtheta = c - theta0[type];
  double dr1 = r1 - ba_r1[type];
  double dr2 = r2 - ba_r2[type];

  energy = ba_k1[type]*dr1*dtheta + ba_k2[type]*dr2*dtheta;
  return energy;
}
