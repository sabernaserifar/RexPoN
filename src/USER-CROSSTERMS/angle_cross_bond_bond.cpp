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
#include "angle_cross_bond_bond.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define SMALL 0.001

/* ---------------------------------------------------------------------- */

AngleCrossBondBond::AngleCrossBondBond(LAMMPS *lmp) : Angle(lmp) {}

/* ---------------------------------------------------------------------- */

AngleCrossBondBond::~AngleCrossBondBond()
{
  if (allocated) {
    memory->destroy(setflag);

    memory->destroy(bb_k);
    memory->destroy(bb_r1);
    memory->destroy(bb_r2);
  }
}

/* ---------------------------------------------------------------------- */

void AngleCrossBondBond::compute(int eflag, int vflag)
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

    // force & energy for bond-bond term

    dr1 = r1 - bb_r1[type];
    dr2 = r2 - bb_r2[type];
    tk1 = bb_k[type] * dr1;
    tk2 = bb_k[type] * dr2;

    f1[0] = -delx1*tk2/r1;
    f1[1] = -dely1*tk2/r1;
    f1[2] = -delz1*tk2/r1;

    f3[0] = -delx2*tk1/r2;
    f3[1] = -dely2*tk1/r2;
    f3[2] = -delz2*tk1/r2;

    if (eflag) eangle = bb_k[type]*dr1*dr2;
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

void AngleCrossBondBond::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;

  memory->create(bb_k,n+1,"angle:bb_k");
  memory->create(bb_r1,n+1,"angle:bb_r1");
  memory->create(bb_r2,n+1,"angle:bb_r2");

  memory->create(setflag,n+1,"angle:setflag");
  for (int i = 1; i <= n; i++)
    setflag[i] = 0;
}

/* ---------------------------------------------------------------------- */

void AngleCrossBondBond::coeff(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR,"Incorrect args for angle coefficients");

  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->nangletypes,ilo,ihi);

  int count = 0;


  double bb_r1_one = atof(arg[1]);
  double bb_r2_one = atof(arg[2]);
  double bb_k_one = atof(arg[3]);
  
  for (int i = ilo; i <= ihi; i++) {
    bb_k[i] = bb_k_one;
    bb_r1[i] = bb_r1_one;
    bb_r2[i] = bb_r2_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for angle coefficients");
}

/* ---------------------------------------------------------------------- */
                                                                                                                                              
double AngleCrossBondBond::equilibrium_angle(int i)
{
  return 0;
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleCrossBondBond::write_restart(FILE *fp)
{
  fwrite(&bb_k[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&bb_r1[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&bb_r2[1],sizeof(double),atom->nangletypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them 
------------------------------------------------------------------------- */

void AngleCrossBondBond::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&bb_k[1],sizeof(double),atom->nangletypes,fp);
    fread(&bb_r1[1],sizeof(double),atom->nangletypes,fp);
    fread(&bb_r2[1],sizeof(double),atom->nangletypes,fp);
  }

  MPI_Bcast(&bb_k[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&bb_r1[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&bb_r2[1],atom->nangletypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ---------------------------------------------------------------------- */

double AngleCrossBondBond::single(int type, int i1, int i2, int i3)
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

  double dr1 = r1 - bb_r1[type];
  double dr2 = r2 - bb_r2[type];
  energy = bb_k[type]*dr1*dr2;

  return energy;
}
