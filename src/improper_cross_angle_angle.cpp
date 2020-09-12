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
#include "improper_cross_angle_angle.h"
#include "atom.h"
#include "neighbor.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define SMALL 0.001

/* ---------------------------------------------------------------------- */

ImproperCrossAngleAngle::ImproperCrossAngleAngle(LAMMPS *lmp) : Improper(lmp)
{
  PI = 4.0*atan(1.0);
}

/* ---------------------------------------------------------------------- */

ImproperCrossAngleAngle::~ImproperCrossAngleAngle()
{
  if (allocated) {
    memory->destroy(setflag);

    memory->destroy(aa_k1);
    memory->destroy(aa_k2);
    memory->destroy(aa_k3);
    memory->destroy(aa_theta0_1);
    memory->destroy(aa_theta0_2);
    memory->destroy(aa_theta0_3);
  }
}

/* ---------------------------------------------------------------------- */

void ImproperCrossAngleAngle::allocate()
{
  allocated = 1;
  int n = atom->nimpropertypes;

  memory->create(aa_k1,n+1,"improper:aa_k1");
  memory->create(aa_k2,n+1,"improper:aa_k2");
  memory->create(aa_k3,n+1,"improper:aa_k3");
  memory->create(aa_theta0_1,n+1,"improper:aa_theta0_1");
  memory->create(aa_theta0_2,n+1,"improper:aa_theta0_2");
  memory->create(aa_theta0_3,n+1,"improper:aa_theta0_3");

  memory->create(setflag,n+1,"improper:setflag");
  for (int i = 1; i <= n; i++)
    setflag[i] = 0;
}

/* ---------------------------------------------------------------------- */

void ImproperCrossAngleAngle::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->nimpropertypes,ilo,ihi);

  int count = 0;

  if (narg != 7) error->all(FLERR,"Incorrect args for improper coefficients");

  double k1_one = atof(arg[4]);
  double k2_one = atof(arg[5]);
  double k3_one = atof(arg[6]);
  double theta0_1_one = atof(arg[1]);
  double theta0_2_one = atof(arg[2]);
  double theta0_3_one = atof(arg[3]);
    
  // convert theta0's from degrees to radians

  for (int i = ilo; i <= ihi; i++) {
    aa_k1[i] = k1_one;
    aa_k2[i] = k2_one;
    aa_k3[i] = k3_one;
    aa_theta0_1[i] = theta0_1_one/180.0 * PI;
    aa_theta0_2[i] = theta0_2_one/180.0 * PI;
    aa_theta0_3[i] = theta0_3_one/180.0 * PI;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for improper coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file 
------------------------------------------------------------------------- */

void ImproperCrossAngleAngle::write_restart(FILE *fp)
{
  fwrite(&aa_k1[1],sizeof(double),atom->nimpropertypes,fp);
  fwrite(&aa_k2[1],sizeof(double),atom->nimpropertypes,fp);
  fwrite(&aa_k3[1],sizeof(double),atom->nimpropertypes,fp);
  fwrite(&aa_theta0_1[1],sizeof(double),atom->nimpropertypes,fp);
  fwrite(&aa_theta0_2[1],sizeof(double),atom->nimpropertypes,fp);
  fwrite(&aa_theta0_3[1],sizeof(double),atom->nimpropertypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them 
------------------------------------------------------------------------- */

void ImproperCrossAngleAngle::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&aa_k1[1],sizeof(double),atom->nimpropertypes,fp);
    fread(&aa_k2[1],sizeof(double),atom->nimpropertypes,fp);
    fread(&aa_k3[1],sizeof(double),atom->nimpropertypes,fp);
    fread(&aa_theta0_1[1],sizeof(double),atom->nimpropertypes,fp);
    fread(&aa_theta0_2[1],sizeof(double),atom->nimpropertypes,fp);
    fread(&aa_theta0_3[1],sizeof(double),atom->nimpropertypes,fp);
  }
  MPI_Bcast(&aa_k1[1],atom->nimpropertypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&aa_k2[1],atom->nimpropertypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&aa_k3[1],atom->nimpropertypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&aa_theta0_1[1],atom->nimpropertypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&aa_theta0_2[1],atom->nimpropertypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&aa_theta0_3[1],atom->nimpropertypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nimpropertypes; i++) setflag[i] = 1;
}

/* ------------------------------------------------------------------------- */

void ImproperCrossAngleAngle::compute(int eflag, int vflag)
{
  int i1,i2,i3,i4,i,j,k,n,type;
  double eimproper;
  double delxAB,delyAB,delzAB,rABmag2,rAB;
  double delxBC,delyBC,delzBC,rBCmag2,rBC;
  double delxBD,delyBD,delzBD,rBDmag2,rBD;
  double costhABC,thetaABC,costhABD;
  double thetaABD,costhCBD,thetaCBD,dthABC,dthCBD,dthABD;
  double sc1,t1,t3,r12;
  double dthetadr[3][4][3],fabcd[4][3];

  double **x = atom->x;
  double **f = atom->f;
  int **improperlist = neighbor->improperlist;
  int nimproperlist = neighbor->nimproperlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nimproperlist; n++) {
    i1 = improperlist[n][0];
    i2 = improperlist[n][1];
    i3 = improperlist[n][2];
    i4 = improperlist[n][3];
    type = improperlist[n][4];

    // difference vectors

    delxAB = x[i1][0] - x[i2][0];
    delyAB = x[i1][1] - x[i2][1];
    delzAB = x[i1][2] - x[i2][2];
    domain->minimum_image(delxAB,delyAB,delzAB);

    delxBC = x[i3][0] - x[i2][0];
    delyBC = x[i3][1] - x[i2][1];
    delzBC = x[i3][2] - x[i2][2];
    domain->minimum_image(delxBC,delyBC,delzBC);

    delxBD = x[i4][0] - x[i2][0];
    delyBD = x[i4][1] - x[i2][1];
    delzBD = x[i4][2] - x[i2][2];
    domain->minimum_image(delxBD,delyBD,delzBD);

    // bond lengths

    rABmag2 = delxAB*delxAB + delyAB*delyAB + delzAB*delzAB;
    rAB = sqrt(rABmag2);
    rBCmag2 = delxBC*delxBC + delyBC*delyBC + delzBC*delzBC;
    rBC = sqrt(rBCmag2);
    rBDmag2 = delxBD*delxBD + delyBD*delyBD + delzBD*delzBD;
    rBD = sqrt(rBDmag2);
        
    // angle ABC, ABD, CBD

    costhABC = (delxAB*delxBC + delyAB*delyBC + delzAB*delzBC) / (rAB * rBC);
    if (costhABC > 1.0)  costhABC = 1.0;
    if (costhABC < -1.0) costhABC = -1.0;
    thetaABC = acos(costhABC);

    costhABD = (delxAB*delxBD + delyAB*delyBD + delzAB*delzBD) / (rAB * rBD);
    if (costhABD > 1.0)  costhABD = 1.0;
    if (costhABD < -1.0) costhABD = -1.0;
    thetaABD = acos(costhABD);

    costhCBD = (delxBC*delxBD + delyBC*delyBD + delzBC*delzBD) /(rBC * rBD);
    if (costhCBD > 1.0)  costhCBD = 1.0;
    if (costhCBD < -1.0) costhCBD = -1.0;
    thetaCBD = acos(costhCBD);

    dthABC = thetaABC - aa_theta0_1[type];
    dthCBD = thetaCBD - aa_theta0_2[type];
    dthABD = thetaABD - aa_theta0_3[type];

    // energy

    if (eflag) eimproper = aa_k2[type] * dthABC * dthABD + 
		 aa_k1[type] * dthABC * dthCBD +
		 aa_k3[type] * dthABD * dthCBD;

    // d(theta)/d(r) array
    // angle i, atom j, coordinate k

    for (i = 0; i < 3; i++)
      for (j = 0; j < 4; j++)
	for (k = 0; k < 3; k++)
	  dthetadr[i][j][k] = 0.0;

    // angle ABC

    sc1 = sqrt(1.0/(1.0 - costhABC*costhABC));
    t1 = costhABC / rABmag2;
    t3 = costhABC / rBCmag2;
    r12 = 1.0 / (rAB * rBC);

    dthetadr[0][0][0] = sc1 * ((t1 * delxAB) - (delxBC * r12));
    dthetadr[0][0][1] = sc1 * ((t1 * delyAB) - (delyBC * r12));
    dthetadr[0][0][2] = sc1 * ((t1 * delzAB) - (delzBC * r12));
    dthetadr[0][1][0] = sc1 * ((-t1 * delxAB) + (delxBC * r12) +
			       (-t3 * delxBC) + (delxAB * r12));
    dthetadr[0][1][1] = sc1 * ((-t1 * delyAB) + (delyBC * r12) +
			       (-t3 * delyBC) + (delyAB * r12));
    dthetadr[0][1][2] = sc1 * ((-t1 * delzAB) + (delzBC * r12) +
			       (-t3 * delzBC) + (delzAB * r12));
    dthetadr[0][2][0] = sc1 * ((t3 * delxBC) - (delxAB * r12));
    dthetadr[0][2][1] = sc1 * ((t3 * delyBC) - (delyAB * r12));
    dthetadr[0][2][2] = sc1 * ((t3 * delzBC) - (delzAB * r12));

    // angle CBD

    sc1 = sqrt(1.0/(1.0 - costhCBD*costhCBD));
    t1 = costhCBD / rBCmag2;
    t3 = costhCBD / rBDmag2;
    r12 = 1.0 / (rBC * rBD);

    dthetadr[1][2][0] = sc1 * ((t1 * delxBC) - (delxBD * r12));
    dthetadr[1][2][1] = sc1 * ((t1 * delyBC) - (delyBD * r12));
    dthetadr[1][2][2] = sc1 * ((t1 * delzBC) - (delzBD * r12));
    dthetadr[1][1][0] = sc1 * ((-t1 * delxBC) + (delxBD * r12) +
			       (-t3 * delxBD) + (delxBC * r12));
    dthetadr[1][1][1] = sc1 * ((-t1 * delyBC) + (delyBD * r12) +
			       (-t3 * delyBD) + (delyBC * r12));
    dthetadr[1][1][2] = sc1 * ((-t1 * delzBC) + (delzBD * r12) +
			       (-t3 * delzBD) + (delzBC * r12));
    dthetadr[1][3][0] = sc1 * ((t3 * delxBD) - (delxBC * r12));
    dthetadr[1][3][1] = sc1 * ((t3 * delyBD) - (delyBC * r12));
    dthetadr[1][3][2] = sc1 * ((t3 * delzBD) - (delzBC * r12));

    // angle ABD

    sc1 = sqrt(1.0/(1.0 - costhABD*costhABD));
    t1 = costhABD / rABmag2;
    t3 = costhABD / rBDmag2;
    r12 = 1.0 / (rAB * rBD);

    dthetadr[2][0][0] = sc1 * ((t1 * delxAB) - (delxBD * r12));
    dthetadr[2][0][1] = sc1 * ((t1 * delyAB) - (delyBD * r12));
    dthetadr[2][0][2] = sc1 * ((t1 * delzAB) - (delzBD * r12));
    dthetadr[2][1][0] = sc1 * ((-t1 * delxAB) + (delxBD * r12) +
			       (-t3 * delxBD) + (delxAB * r12));
    dthetadr[2][1][1] = sc1 * ((-t1 * delyAB) + (delyBD * r12) +
			       (-t3 * delyBD) + (delyAB * r12));
    dthetadr[2][1][2] = sc1 * ((-t1 * delzAB) + (delzBD * r12) +
			       (-t3 * delzBD) + (delzAB * r12));
    dthetadr[2][3][0] = sc1 * ((t3 * delxBD) - (delxAB * r12));
    dthetadr[2][3][1] = sc1 * ((t3 * delyBD) - (delyAB * r12));
    dthetadr[2][3][2] = sc1 * ((t3 * delzBD) - (delzAB * r12));

    // angleangle forces

    for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++)
	fabcd[i][j] = - 
	  ((aa_k1[type] * 
	    (dthABC*dthetadr[1][i][j] + dthCBD*dthetadr[0][i][j])) +
	   (aa_k2[type] * 
	    (dthABC*dthetadr[2][i][j] + dthABD*dthetadr[0][i][j])) +
	   (aa_k3[type] *
	    (dthABD*dthetadr[1][i][j] + dthCBD*dthetadr[2][i][j])));

    // apply force to each of 4 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += fabcd[0][0];
      f[i1][1] += fabcd[0][1];
      f[i1][2] += fabcd[0][2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] += fabcd[1][0];
      f[i2][1] += fabcd[1][1];
      f[i2][2] += fabcd[1][2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += fabcd[2][0];
      f[i3][1] += fabcd[2][1];
      f[i3][2] += fabcd[2][2];
    }

    if (newton_bond || i4 < nlocal) {
      f[i4][0] += fabcd[3][0];
      f[i4][1] += fabcd[3][1];
      f[i4][2] += fabcd[3][2];
    }

    if (evflag)
      ev_tally(i1,i2,i3,i4,nlocal,newton_bond,eimproper,
	       fabcd[0],fabcd[2],fabcd[3],
	       delxAB,delyAB,delzAB,delxBC,delyBC,delzBC,delxBD,delyBD,delzBD);
  }
}
