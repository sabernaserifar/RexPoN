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
#include "dihedral_cross_13_bond_bond.h"
#include "atom.h"
#include "neighbor.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define TOLERANCE 0.05
#define SMALL     0.0000001

/* ---------------------------------------------------------------------- */

DihedralCross13BondBond::DihedralCross13BondBond(LAMMPS *lmp) : Dihedral(lmp) {}

/* ---------------------------------------------------------------------- */

DihedralCross13BondBond::~DihedralCross13BondBond()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(bb13t_k);
    memory->destroy(bb13t_r10);
    memory->destroy(bb13t_r30);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralCross13BondBond::compute(int eflag, int vflag)
{
  int i1,i2,i3,i4,i,j,k,n,type;
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z;
  double edihedral;
  double r1mag2,r1,r3mag2,r3;
  double r1_0,r3_0,dr1,dr2,tk1,tk2;
  double fabcd[4][3];

  edihedral = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;
  int **dihedrallist = neighbor->dihedrallist;
  int ndihedrallist = neighbor->ndihedrallist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < ndihedrallist; n++) {
    i1 = dihedrallist[n][0];
    i2 = dihedrallist[n][1];
    i3 = dihedrallist[n][2];
    i4 = dihedrallist[n][3];
    type = dihedrallist[n][4];

    // 1st bond

    vb1x = x[i1][0] - x[i2][0];
    vb1y = x[i1][1] - x[i2][1];
    vb1z = x[i1][2] - x[i2][2];
    domain->minimum_image(vb1x,vb1y,vb1z);

    // 2nd bond

    vb2x = x[i3][0] - x[i2][0];
    vb2y = x[i3][1] - x[i2][1];
    vb2z = x[i3][2] - x[i2][2];
    domain->minimum_image(vb2x,vb2y,vb2z);

    // 3rd bond

    vb3x = x[i4][0] - x[i3][0];
    vb3y = x[i4][1] - x[i3][1];
    vb3z = x[i4][2] - x[i3][2];
    domain->minimum_image(vb3x,vb3y,vb3z);

    // distances

    r1mag2 = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z;
    r1 = sqrt(r1mag2);
    r3mag2 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z;
    r3 = sqrt(r3mag2);

    // bond1/bond3 coupling

    if (fabs(bb13t_k[type]) > SMALL) {

      r1_0 = bb13t_r10[type];
      r3_0 = bb13t_r30[type];
      dr1 = r1 - r1_0;
      dr2 = r3 - r3_0;
      tk1 = -bb13t_k[type] * dr1 / r3;
      tk2 = -bb13t_k[type] * dr2 / r1;

      if (eflag) edihedral = bb13t_k[type]*dr1*dr2;
        
      fabcd[0][0] = tk2 * vb1x;
      fabcd[0][1] = tk2 * vb1y;
      fabcd[0][2] = tk2 * vb1z;

      fabcd[1][0] = -1 * (tk2 * vb1x);
      fabcd[1][1] = -1 * (tk2 * vb1y);
      fabcd[1][2] = -1 * (tk2 * vb1z);
        
      fabcd[2][0] = -1 * (tk1 * vb3x);
      fabcd[2][1] = -1 * (tk1 * vb3y);
      fabcd[2][2] = -1 * (tk1 * vb3z);

      fabcd[3][0] = tk1 * vb3x;
      fabcd[3][1] = tk1 * vb3y;
      fabcd[3][2] = tk1 * vb3z;
    }

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
      ev_tally(i1,i2,i3,i4,nlocal,newton_bond,edihedral,
	       fabcd[0],fabcd[2],fabcd[3],
	       vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralCross13BondBond::allocate()
{
  allocated = 1;
  int n = atom->ndihedraltypes;

  memory->create(bb13t_k,n+1,"dihedral:bb13t_k");
  memory->create(bb13t_r10,n+1,"dihedral:bb13t_r10");
  memory->create(bb13t_r30,n+1,"dihedral:bb13t_r30");

  memory->create(setflag,n+1,"dihedral:setflag");
  for (int i = 1; i <= n; i++)
    setflag[i] = 0;
}

/* ---------------------------------------------------------------------- */

void DihedralCross13BondBond::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->ndihedraltypes,ilo,ihi);

  int count = 0;

  if (narg != 4) error->all(FLERR,"Incorrect args for dihedral coefficients");

  double k_one = atof(arg[3]);
  double r10_one = atof(arg[1]);
  double r30_one = atof(arg[2]);
    
  for (int i = ilo; i <= ihi; i++) {
    bb13t_k[i] = k_one;
    bb13t_r10[i] = r10_one;
    bb13t_r30[i] = r30_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for dihedral coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file 
------------------------------------------------------------------------- */

void DihedralCross13BondBond::write_restart(FILE *fp)
{
  fwrite(&bb13t_k[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&bb13t_r10[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&bb13t_r30[1],sizeof(double),atom->ndihedraltypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them 
------------------------------------------------------------------------- */

void DihedralCross13BondBond::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&bb13t_k[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&bb13t_r10[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&bb13t_r30[1],sizeof(double),atom->ndihedraltypes,fp);
  }

  MPI_Bcast(&bb13t_k[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&bb13t_r10[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&bb13t_r30[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->ndihedraltypes; i++) setflag[i] = 1;
}
