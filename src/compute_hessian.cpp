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

#include "mpi.h"
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "stdio.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "compute_hessian.h"
#include "atom.h"
#include "atom_vec.h"
#include "error.h"
#include "update.h"
#include "memory.h"
#include "domain.h"
#include "modify.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "comm.h"
#include "fix_pqeq.h"
#include "pair_coul_pqeqgauss.h"
#include "group.h"
#include "math_const.h"
#include "iostream"


using namespace LAMMPS_NS;
using namespace MathConst;

#define SQR(x) ((x)*(x))


/* ---------------------------------------------------------------------- */

ComputeHessian::ComputeHessian(LAMMPS *lmp, int narg, char **arg)
    : Compute(lmp, narg, arg) {
  if (narg != 5)
    error->all(FLERR, "Illegal compute hessian command");

  epsilon = atof(arg[3]);
  iepsilon = 1 / epsilon;

  masselect = atof(arg[4]);

  /* even though this is a massive 2d array, return the a vector instead.
   * we will explicitly manage the addressing in each dimension with a
   * preprocessor index macro. */
  vector_flag = 1;
  extvector = 0;

  /* these values will change if the system size changes. */
  int ndofs = atom->natoms * 3;
  size_vector = ndofs * ndofs;

  mylocalsize = 0;
  myglobalsize = 0;

  fglobal_ref = fglobal_new = fglobal_copy = NULL;
  hessian = NULL;

  // PQEq shell force 
  pack_flag = 0;
  sf = NULL;
  comm_forward = comm_reverse = 3;

}

/* ---------------------------------------------------------------------- */

ComputeHessian::~ComputeHessian() {
  free(fglobal_ref);
  free(fglobal_new);
  free(fglobal_copy);
  free(hessian);
}
/* ---------------------------------------------------------------------- */

void ComputeHessian::init()
{
  if (!atom->q_flag) error->all(FLERR,"Fix pqeq requires atom attribute q");

  pqeq = (PairCoulPqeqgauss *) force->pair_match("coul/pqeqgauss",1);
  if (pqeq == NULL) error->all(FLERR,"Cannot use fix pqeq without "
                  "pair_style coul/pqeqgauss");

  ngroup = group->count(igroup);
  if (ngroup == 0) error->all(FLERR,"Fix pqeq group has no atoms");


  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->newton = 2;
  neighbor->requests[irequest]->ghost = 1;

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

void ComputeHessian::init_list(int id, NeighList *ptr)
{
  list = ptr;
}
/* ---------------------------------------------------------------------- */
void ComputeHessian::compute_vector(void) {
  invoked_vector = update->ntimestep;

  /* tags must be defined and consecutive. */
  if (atom->tag_enable == 0)
    error->all(FLERR,
               "Cannot use Hessian compute unless atoms have IDs");
  if (atom->tag_consecutive() == 0)
    error->all(FLERR,
               "Atom IDs must be consecutive for Hessian compute");

  /* get pointers to all the original data. */
  //double **x = atom->x;
  //double **f = atom->f;
  double **rsx = atom->rsx;


  memory->create(sf,atom->nmax,3,"pqeq:sf");

  /* the global force and hessian arrays must be explicitly the correct size. */
  int needglobalsize = atom->natoms;
  int ndofs = atom->natoms * 3;
  bigint nhessianelements = ndofs * ndofs;
  if (needglobalsize != myglobalsize) {
    free (fglobal_ref); 
    free (fglobal_new); 
    free (fglobal_copy);
    free (hessian);

    fglobal_ref = (double *) malloc (ndofs * sizeof (double));   
    fglobal_new = (double *) malloc (ndofs * sizeof (double));   
    fglobal_copy = (double *) malloc (ndofs * sizeof (double));   
    hessian = (double *) malloc (nhessianelements * sizeof (double));


    /* always be sure to set the output vector since the address keeps changing. */
    vector = hessian;

    myglobalsize = needglobalsize;
  }

  /* a lot of the hessian will be zero, so start there. */
  memset (hessian, 0, nhessianelements * sizeof(double));

  /* set up a map if none exists so we can incrementally loop through all dofs
   * regardless of the location of the atom data. */
  int mapflag = 0;
  if (atom->map_style == 0) {
    mapflag = 1;
    atom->map_init();
    atom->map_set();
  }

  /* no energy or virial updates. */
  int eflag = 0;
  int vflag = 0;

  /* allow pair and kspace compute to be turned off via modify flags. */
  /*if (force->pair && force->pair->compute_flag)
    pair_compute_flag = 1;
  else
    pair_compute_flag = 0;
  if (force->kspace && force->kspace->compute_flag)
    kspace_compute_flag = 1;
  else
    kspace_compute_flag = 0;*/

  /* do a standard force call to get the reference forces. */
  /*comm->forward_comm();
  force_clear();
  if (modify->n_pre_force) modify->pre_force(vflag);

  if (pair_compute_flag) force->pair->compute(eflag, vflag);

  if (atom->molecular) {
    if (force->bond) force->bond->compute(eflag, vflag);
    if (force->angle) force->angle->compute(eflag, vflag);
    if (force->dihedral) force->dihedral->compute(eflag, vflag);
    if (force->improper) force->improper->compute(eflag, vflag);
  }

  if (kspace_compute_flag) force->kspace->compute(eflag, vflag);
  if (force->newton) comm->reverse_comm();
  if (modify->n_post_force) modify->post_force(vflag);*/


  // shell force call 
  sforce(); 


  /* construct fglobal_ref by explicit scatter and reduce to preserve atom-id
   * ordering. */
  int m, reduce_m;
  memset (&fglobal_copy[0], 0, myglobalsize * 3 * sizeof (double));
  for (int i = 1; i <= atom->natoms; i++) {
    m = atom->map(i);
    if (atom->mask[m]) {
      reduce_m = atom->tag[m] - 1;
      for (int j = 0; j < domain->dimension; j++)
        fglobal_copy[idx2_c(reduce_m, j, atom->natoms)] = sf[m][j];
    }
  }
  MPI_Allreduce (fglobal_copy, fglobal_ref, ndofs, MPI_DOUBLE, MPI_SUM, world);

  /* do numerical hessian compute by forward differences. */
  int n, reduce_n, index_a, index_b, global_atom_a, global_atom_b;
  double mass, difference, mass_weight, xstore;
  for (int i = 1; i <= atom->natoms; i++) {

    m = atom->map(i);
    if (atom->mask[m]) {
      /* global ids in lammps are handled by 1-based indexing, while everything
       * local is 0-based. */
      global_atom_a = atom->tag[m] - 1;
      //MPI_Bcast(&global_atom_a, 1, MPI_INT, comm->me, world);

      if (atom->rmass) {
        //mass = atom->rmass[m];
        mass = masselect;
        //MPI_Bcast(&mass, 1, MPI_DOUBLE, comm->me, world);
      } else {
        //mass = atom->mass[atom->type[m]];
        mass = masselect;
        //MPI_Bcast(&mass, 1, MPI_DOUBLE, comm->me, world);
      }
    }

    for (int j = 0; j < domain->dimension; j++) {
      /* increment the dof by epsilon on the right task. */
      if (atom->mask[m]) {
        xstore = rsx[m][j];
        rsx[m][j] += epsilon;
      }

      /* standard force call. */
      /*comm->forward_comm();
      force_clear();
      if (modify->n_pre_force) modify->pre_force(vflag);

      if (pair_compute_flag) force->pair->compute(eflag, vflag);

      if (atom->molecular) {
        if (force->bond) force->bond->compute(eflag, vflag);
        if (force->angle) force->angle->compute(eflag, vflag);
        if (force->dihedral) force->dihedral->compute(eflag, vflag);
        if (force->improper) force->improper->compute(eflag, vflag);
      }

      if (kspace_compute_flag) force->kspace->compute(eflag, vflag);*/

      sforce();


      /* put the original position back. */
      if (atom->mask[m]) rsx[m][j] = xstore;

      //if (force->newton) comm->reverse_comm();
      //if (modify->n_post_force) modify->post_force(vflag);

      /* construct fglobal_new by explicit scatter and reduce to preserve
       * atom-id ordering. */
      memset (&fglobal_copy[0], 0, myglobalsize * 3 * sizeof (double));
      for (int k = 1; k <= atom->natoms; k++) {
        n = atom->map(k);
        if (atom->mask[n]) {
          reduce_n = atom->tag[n] - 1;
          for (int l = 0; l < domain->dimension; l++)
            fglobal_copy[idx2_c(reduce_n, l, atom->natoms)] = sf[n][l];
        }
      }
      MPI_Allreduce (fglobal_copy, fglobal_new, ndofs, MPI_DOUBLE, MPI_SUM, world);

      /* compute the difference (not using symmetry so we can do an in-place
       * reduciton). */
      index_a = j + 3 * global_atom_a;
      for (int k = 1; k <= atom->natoms; k++) {
        n = atom->map(k);
        if (atom->mask[n]) {
          global_atom_b = atom->tag[n] - 1;

          /* don't need to broadcast the second mass because it will only be used
           * on this rank. */
          if (atom->rmass)
            //mass_weight = 1 / sqrt(mass * atom->rmass[n]);
            mass_weight = 1 / sqrt(mass * mass);

          else
            //mass_weight = 1 / sqrt(mass * atom->mass[atom->type[n]]);
            mass_weight = 1 / sqrt(mass * mass);


          /* once again, global arrays use 1-based indexing, so have to rebase
           * them to 0. */
          for (int l = 0; l < domain->dimension; l++) {
            index_b = l + 3 * global_atom_b;
            difference =
                fglobal_ref[idx2_c(global_atom_b, l, atom->natoms)] - \
                fglobal_new[idx2_c(global_atom_b, l, atom->natoms)];

            hessian[idx2_c(index_a, index_b, ndofs)] =
                difference * iepsilon * mass_weight;
          }
        }
      }
    }
  }

  /* only reduce the hessian to the root task. */
  MPI_Reduce(MPI_IN_PLACE, hessian, nhessianelements, MPI_DOUBLE, MPI_SUM, 0, world);

  //MPI_Allreduce (hessian,vector, nhessianelements, MPI_DOUBLE, MPI_SUM, world);
  std::cout<<"hessian "<<hessian[0]<<std::endl;


  /* destroy the atom map. */
  if (mapflag) {
    atom->map_delete();
    atom->map_style = 0;
  }

  /* do a standard force call to get the original forces back. */
  /*comm->forward_comm();
  force_clear();
  if (modify->n_pre_force) modify->pre_force(vflag);

  if (pair_compute_flag) force->pair->compute(eflag, vflag);

  if (atom->molecular) {
    if (force->bond) force->bond->compute(eflag, vflag);
    if (force->angle) force->angle->compute(eflag, vflag);
    if (force->dihedral) force->dihedral->compute(eflag, vflag);
    if (force->improper) force->improper->compute(eflag, vflag);
  }

  if (kspace_compute_flag) force->kspace->compute(eflag, vflag);
  if (force->newton) comm->reverse_comm();
  if (modify->n_post_force) modify->post_force(vflag);*/
 
  sforce();


}
/* ---------------------------------------------------------------------- */

void ComputeHessian::force_clear() {
  size_t nbytes;
  int nlocal = atom->nlocal;

  nbytes = sizeof(double) * nlocal;
  if (force->newton) nbytes += sizeof(double) * atom->nghost;

  if (nbytes) memset (&atom->f[0][0], 0, 3 * nbytes);
}
/* ---------------------------------------------------------------------- */
void ComputeHessian::sforce()
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


  // shell positions 
  pack_flag = 1;
  comm->forward_comm_compute(this);


  // shell force 
  pack_flag = 2;
  comm->forward_comm_compute(this);

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

  comm->reverse_comm_compute( this );
}
/* -----------------------------------------------------------------------
   Coulomb 
-------------------------------------------------------------------------- */

void ComputeHessian::coulomb(double r, double alpha, double *coulombr, double *dcoulombr)
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

int ComputeHessian::pack_forward_comm(int n, int *list, double *buf,
                          int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;

  if ( pack_flag == 1)
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = atom->rsx[j][0];
      buf[m++] = atom->rsx[j][1];
      buf[m++] = atom->rsx[j][2];
    }
  else if ( pack_flag == 2)
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = sf[j][0];
      buf[m++] = sf[j][1];
      buf[m++] = sf[j][2];
    }

  return m;
}

/* ---------------------------------------------------------------------- */
void ComputeHessian::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n ;

  if( pack_flag == 1)
    for (i = first; i < last; i++) {
      atom->rsx[i][0] = buf[m++];
      atom->rsx[i][1] = buf[m++];
      atom->rsx[i][2] = buf[m++];
    }
  else if( pack_flag == 2)
    for (i = first; i < last; i++) {
      sf[i][0] = buf[m++];
      sf[i][1] = buf[m++];
      sf[i][2] = buf[m++];
    }

}

/* ---------------------------------------------------------------------- */
int ComputeHessian::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  
  for (i = first; i < last; i++) {
      buf[m++] = sf[i][0];
      buf[m++] = sf[i][1];
      buf[m++] = sf[i][2];
    }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeHessian::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
      j = list[i];
      sf[j][0] += buf[m++];
      sf[j][1] += buf[m++];
      sf[j][2] += buf[m++];
  } 

}
