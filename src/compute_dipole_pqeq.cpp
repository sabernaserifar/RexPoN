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

#include "math.h"
#include "math_const.h"
#include <string.h>
#include "compute_dipole_pqeq.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "domain.h"
#include "modify.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeDipolePqeq::ComputeDipolePqeq(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal compute dipole/pqeq command");

  vector_flag = 1;
  size_vector = 4;
  extvector = 0;

  vector = new double[4];
}

/* ---------------------------------------------------------------------- */

ComputeDipolePqeq::~ComputeDipolePqeq()
{
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeDipolePqeq::init()
{
  if (!atom->q_flag)
    error->all(FLERR,"dielectric calculations needs charge models");
}

/* ---------------------------------------------------------------------- */

void ComputeDipolePqeq::compute_vector()
{

  invoked_vector = update->ntimestep;

  double unwrap[3],dipole[4];

  dipole[0] = dipole[1] = dipole[2] = dipole[3] = 0.0;

  double **x = atom->x;
  double **rsx = atom->rsx;
  int *mask = atom->mask;
  int *type = atom->type;
  imageint *image = atom->image;
  double *q = atom->q;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i],image[i],unwrap);
      dipole[0] += q[i]*(unwrap[0]);
      dipole[1] += q[i]*(unwrap[1]);
      dipole[2] += q[i]*(unwrap[2]);
      dipole[0] += (-rsx[i][0]);
      dipole[1] += (-rsx[i][1]);
      dipole[2] += (-rsx[i][2]);

    }

  MPI_Allreduce(dipole,vector,4,MPI_DOUBLE,MPI_SUM,world);

  // total dipole 
  vector[3] = sqrt(vector[0]*vector[0]+vector[1]*vector[1]+
                     vector[2]*vector[2]);
}
