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
   Contributing author: Saber Naserifar (Caltech, naserifar.saber@gmail.com)

------------------------------------------------------------------------- */

#include "fix_rexpon.h"
#include "atom.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define MAX_REXPON_BONDS      30
#define MIN_REXPON_BONDS      15
#define MIN_REXPON_HBONDS     25

/* ---------------------------------------------------------------------- */

FixRexPoN::FixRexPoN(LAMMPS *lmp,int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  // perform initial allocation of atom-based arrays
  // register with atom class

  num_bonds = NULL;
  num_hbonds = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);

  // initialize arrays to MIN so atom migration is OK the 1st time

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++)
    num_bonds[i] = num_hbonds[i] = MIN_REXPON_BONDS;

  // set comm sizes needed by this fix

  comm_forward = 1;
}

/* ---------------------------------------------------------------------- */

FixRexPoN::~FixRexPoN()
{
  // unregister this fix so atom class doesn't invoke it any more

  atom->delete_callback(id,0);

  // delete locally stored arrays

  memory->destroy(num_bonds);
  memory->destroy(num_hbonds);
}

/* ---------------------------------------------------------------------- */

int FixRexPoN::setmask()
{
  int mask = 0;
  return mask;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixRexPoN::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax * 2 * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixRexPoN::grow_arrays(int nmax)
{
  memory->grow(num_bonds,nmax,"rexpon:num_bonds");
  memory->grow(num_hbonds,nmax,"rexpon:num_hbonds");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixRexPoN::copy_arrays(int i, int j, int delflag)
{
  num_bonds[j] = num_bonds[i];
  num_hbonds[j] = num_hbonds[i];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixRexPoN::pack_exchange(int i, double *buf)
{
  buf[0] = num_bonds[i];
  buf[1] = num_hbonds[i];
  return 2;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixRexPoN::unpack_exchange(int nlocal, double *buf)
{
  num_bonds[nlocal] = static_cast<int> (buf[0]);
  num_hbonds[nlocal] = static_cast<int> (buf[1]);
  return 2;
}

/* ---------------------------------------------------------------------- */

int FixRexPoN::pack_forward_comm(int n, int *list, double *buf,
                                int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = num_bonds[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixRexPoN::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    num_bonds[i] = static_cast<int> (buf[m++]);
}
