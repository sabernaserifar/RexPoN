/* ----------------------------------------------------------------------
 *
 *                    *** Smooth Mach Dynamics ***
 *
 * This file is part of the USER-SMD package for LAMMPS.
 * Copyright (2014) Georg C. Ganzenmueller, georg.ganzenmueller@emi.fhg.de
 * Fraunhofer Ernst-Mach Institute for High-Speed Dynamics, EMI,
 * Eckerstrasse 4, D-79104 Freiburg i.Br, Germany.
 *
 * ----------------------------------------------------------------------- */

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fix_smd_integrate_tlsph.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "error.h"
#include "pair.h"
#include "neigh_list.h"
#include <Eigen/Eigen>
#include "domain.h"
#include "neighbor.h"
#include "comm.h"
#include "modify.h"
#include <stdio.h>
#include <iostream>

using namespace Eigen;
using namespace LAMMPS_NS;
using namespace FixConst;
using namespace std;

/* ---------------------------------------------------------------------- */

FixSMDIntegrateTlsph::FixSMDIntegrateTlsph(LAMMPS *lmp, int narg, char **arg) :
                Fix(lmp, narg, arg) {
        if (narg < 3) {
                printf("narg=%d\n", narg);
                error->all(FLERR, "Illegal fix smd/integrate_tlsph command");
        }

        xsphFlag = false;
        vlimit = -1.0;
        int iarg = 3;

        if (comm->me == 0) {
                printf("\n>>========>>========>>========>>========>>========>>========>>========>>========\n");
                printf("fix smd/integrate_tlsph is active for group: %s \n", arg[1]);
        }

        while (true) {

                if (iarg >= narg) {
                        break;
                }

                if (strcmp(arg[iarg], "xsph") == 0) {
                        xsphFlag = true;
                        if (comm->me == 0) {
                                error->one(FLERR, "XSPH is currently not available");
                                printf("... will use XSPH time integration\n");
                        }
                } else if (strcmp(arg[iarg], "limit_velocity") == 0) {
                        iarg++;
                        if (iarg == narg) {
                                error->all(FLERR, "expected number following limit_velocity");
                        }

                        vlimit = force->numeric(FLERR, arg[iarg]);
                        if (comm->me == 0) {
                                printf("... will limit velocities to <= %g\n", vlimit);
                        }
                } else {
                        char msg[128];
                        sprintf(msg, "Illegal keyword for smd/integrate_tlsph: %s\n", arg[iarg]);
                        error->all(FLERR, msg);
                }

                iarg++;
        }

        if (comm->me == 0) {
                printf(">>========>>========>>========>>========>>========>>========>>========>>========\n");
        }

        time_integrate = 1;

        // set comm sizes needed by this fix

        atom->add_callback(0);

}

/* ---------------------------------------------------------------------- */

int FixSMDIntegrateTlsph::setmask() {
        int mask = 0;
        mask |= INITIAL_INTEGRATE;
        mask |= FINAL_INTEGRATE;
        return mask;
}

/* ---------------------------------------------------------------------- */

void FixSMDIntegrateTlsph::init() {
        dtv = update->dt;
        dtf = 0.5 * update->dt * force->ftm2v;
        vlimitsq = vlimit * vlimit;
}

/* ----------------------------------------------------------------------
 ------------------------------------------------------------------------- */

void FixSMDIntegrateTlsph::initial_integrate(int vflag) {
        double dtfm, vsq, scale;

        // update v and x of atoms in group

        double **x = atom->x;
        double **v = atom->v;
        double **vest = atom->vest;
        double **f = atom->f;
        double *rmass = atom->rmass;
        int *mask = atom->mask;
        int nlocal = atom->nlocal;
        int itmp;
        double vxsph_x, vxsph_y, vxsph_z;
        if (igroup == atom->firstgroup)
                nlocal = atom->nfirst;

        Vector3d *smoothVelDifference = (Vector3d *) force->pair->extract("smd/tlsph/smoothVel_ptr", itmp);

        if (xsphFlag) {
                if (smoothVelDifference == NULL) {
                        error->one(FLERR,
                                        "fix smd/integrate_tlsph failed to access smoothVel array. Check if a pair style exist which calculates this quantity.");
                }
        }

        for (int i = 0; i < nlocal; i++) {
                if (mask[i] & groupbit) {
                        dtfm = dtf / rmass[i];

                        // 1st part of Velocity_Verlet: push velocties 1/2 time increment ahead
                        v[i][0] += dtfm * f[i][0];
                        v[i][1] += dtfm * f[i][1];
                        v[i][2] += dtfm * f[i][2];

                        if (vlimit > 0.0) {
                                vsq = v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];
                                if (vsq > vlimitsq) {
                                        scale = sqrt(vlimitsq / vsq);
                                        v[i][0] *= scale;
                                        v[i][1] *= scale;
                                        v[i][2] *= scale;
                                }
                        }

                        if (xsphFlag) {

                                // construct XSPH velocity
                                vxsph_x = v[i][0] + 0.5 * smoothVelDifference[i](0);
                                vxsph_y = v[i][1] + 0.5 * smoothVelDifference[i](1);
                                vxsph_z = v[i][2] + 0.5 * smoothVelDifference[i](2);

                                vest[i][0] = vxsph_x + dtfm * f[i][0];
                                vest[i][1] = vxsph_y + dtfm * f[i][1];
                                vest[i][2] = vxsph_z + dtfm * f[i][2];

                                x[i][0] += dtv * vxsph_x;
                                x[i][1] += dtv * vxsph_y;
                                x[i][2] += dtv * vxsph_z;
                        } else {

                                // extrapolate velocity from half- to full-step
                                vest[i][0] = v[i][0] + dtfm * f[i][0];
                                vest[i][1] = v[i][1] + dtfm * f[i][1];
                                vest[i][2] = v[i][2] + dtfm * f[i][2];

                                x[i][0] += dtv * v[i][0]; // 2nd part of Velocity-Verlet: push positions one full time increment ahead
                                x[i][1] += dtv * v[i][1];
                                x[i][2] += dtv * v[i][2];
                        }
                }
        }

}

/* ---------------------------------------------------------------------- */

void FixSMDIntegrateTlsph::final_integrate() {
        double dtfm, vsq, scale;

// update v of atoms in group

        double **v = atom->v;
        double **f = atom->f;
        double *e = atom->e;
        double *de = atom->de;
        double *rmass = atom->rmass;
        int *mask = atom->mask;
        int nlocal = atom->nlocal;
        if (igroup == atom->firstgroup)
                nlocal = atom->nfirst;
        int i;

        for (i = 0; i < nlocal; i++) {
                if (mask[i] & groupbit) {
                        dtfm = dtf / rmass[i];

                        v[i][0] += dtfm * f[i][0]; // 3rd part of Velocity-Verlet: push velocities another half time increment ahead
                        v[i][1] += dtfm * f[i][1]; // both positions and velocities are now defined at full time-steps.
                        v[i][2] += dtfm * f[i][2];

                        // limit velocity
                        if (vlimit > 0.0) {
                                vsq = v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];
                                if (vsq > vlimitsq) {
                                        scale = sqrt(vlimitsq / vsq);
                                        v[i][0] *= scale;
                                        v[i][1] *= scale;
                                        v[i][2] *= scale;
                                }
                        }

                        e[i] += dtv * de[i];
                }
        }
}

/* ---------------------------------------------------------------------- */

void FixSMDIntegrateTlsph::reset_dt() {
        dtv = update->dt;
        dtf = 0.5 * update->dt * force->ftm2v;
        vlimitsq = vlimit * vlimit;
}

/* ---------------------------------------------------------------------- */
