/*----------------------------------------------------------------------
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details:
  <http://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#include "pair_rexpon.h"
#include "rexpon_bonds.h"
#include "rexpon_bond_orders.h"
#include "rexpon_list.h"
#include "rexpon_tool_box.h"
#include "rexpon_vector.h"

void Bonds( rexpon_system *system, control_params *control,
            simulation_data *data, storage *workspace, rexpon_list **lists,
            output_controls *out_control )
{
  int i, j, pj, natoms;
  int start_i, end_i;
  int type_i, type_j;
  double ebond, pow_BOs_be2, exp_be12, CEbo;
  double gp3, gp4, gp7, gp10, gp37;
  double exphu, exphua1, exphub1, exphuov, hulpov, estriph;
  double decobdbo, decobdboua, decobdboub;
  double bos,bos2,bos3,bos4,bos5,bos6;
  double exp1,exp2,part1,part2;
  double r_ij,frac;
  double r0,a0,b0,c0,n0,scal0,D0;
  double e_vdW,e_lg,e_qm,e_b;
  double r0nb,a0nb,n0nb,scal0nb;
  double re6,re0,r_ij5,r_ij6,c6lg;
  double CEvd_qm,CEvd_vdW,CEvd_lg;
  double drdbo;

  int flag_water = 0;
  double dexp1,dexp2,de_vdW;
  double d_vdW,alf_vdW,re_vdW,bet_vdW,gam_vdW,eta_vdW,del_vdW;
  double le,s0,s1,s2,s3,s4,s5,L2,L3,L4,L5;
  double dr,dr2,dr3,dr4,dr5;
  double Tap=1.0, dTap=0.0; 


  single_body_parameters *sbp_i, *sbp_j;
  two_body_parameters *twbp;
  bond_order_data *bo_ij;
  rexpon_list *bonds;

  bonds = (*lists) + BONDS;
  gp3 = system->rexpon_param.gp.l[3];
  gp4 = system->rexpon_param.gp.l[4];
  gp7 = system->rexpon_param.gp.l[7];
  gp10 = system->rexpon_param.gp.l[10];
  gp37 = (int) system->rexpon_param.gp.l[37];
  natoms = system->n;

  for( i = 0; i < natoms; ++i ) {
    start_i = Start_Index(i, bonds);
    end_i = End_Index(i, bonds);

    for( pj = start_i; pj < end_i; ++pj ) {
      j = bonds->select.bond_list[pj].nbr;

      if( system->my_atoms[i].orig_id > system->my_atoms[j].orig_id )
	continue;
      if( system->my_atoms[i].orig_id == system->my_atoms[j].orig_id ) {
        if (system->my_atoms[j].x[2] <  system->my_atoms[i].x[2]) continue;
      	if (system->my_atoms[j].x[2] == system->my_atoms[i].x[2] &&
      	    system->my_atoms[j].x[1] <  system->my_atoms[i].x[1]) continue;
        if (system->my_atoms[j].x[2] == system->my_atoms[i].x[2] &&
      	    system->my_atoms[j].x[1] == system->my_atoms[i].x[1] &&
      	    system->my_atoms[j].x[0] <  system->my_atoms[i].x[0]) continue;
      }

      /* set the pointers */
      type_i = system->my_atoms[i].type;
      type_j = system->my_atoms[j].type;
      sbp_i = &( system->rexpon_param.sbp[type_i] );
      sbp_j = &( system->rexpon_param.sbp[type_j] );
      twbp = &( system->rexpon_param.tbp[type_i][type_j] );
      bo_ij = &( bonds->select.bond_list[pj].bo_data );


      flag_water = 0;
      // find if interaction is between water molecules
      if ( !strcmp( system->rexpon_param.sbp[type_i].name, "HW" ) &&
           !strcmp( system->rexpon_param.sbp[type_j].name, "OW" )  )
            flag_water = 1;
      if ( !strcmp( system->rexpon_param.sbp[type_i].name, "OW" ) &&
           !strcmp( system->rexpon_param.sbp[type_j].name, "HW" ) )
            flag_water = 1;
      if ( !strcmp( system->rexpon_param.sbp[type_i].name, "OW" ) &&
           !strcmp( system->rexpon_param.sbp[type_j].name, "OW" ) )
            flag_water = 1;
      if ( !strcmp( system->rexpon_param.sbp[type_i].name, "HW" ) &&
           !strcmp( system->rexpon_param.sbp[type_j].name, "HW" ) )
            flag_water = 1;

      /* calculate the constants */
      /*pow_BOs_be2 = pow( bo_ij->BO_s, twbp->p_be2 );
      exp_be12 = exp( twbp->p_be1 * ( 1.0 - pow_BOs_be2 ) );
      CEbo = -twbp->De_s * exp_be12 *
	( 1.0 - twbp->p_be1 * twbp->p_be2 * pow_BOs_be2 );*/

      bos = bo_ij->BO_s;
      D0 = twbp->De_s;
      if (bos > 0.001 and D0>0.0) {

          /*bos2 = bos*bos;
          bos3 = bos*bos2;
          bos4 = bos2*bos2;
          bos5 = bos4*bos;
          bos6 = bos4*bos2;*/
          frac = (twbp->r_s-twbp->p_bo2) / twbp->p_bo1;
          r_ij = twbp->p_bo2 - frac * log((1.1-bos)/bos);
          drdbo = 1.1*frac/(bos*(1.1-bos));
          //D0 = twbp->De_s;
          a0 = twbp->p_be1;
          b0 = twbp->p_be2;
          c0 = twbp->bom;
          r0 = twbp->r_s;

          part1 = -D0 * exp(-b0*a0*(r_ij-r0));
          part2 = c0 + b0 * a0 * (r_ij-r0) * (c0 + a0*(r_ij-r0) * (c0 + a0*(r_ij-r0)) );
          e_qm = part1 * part2;
          CEvd_qm = -b0*a0 * part1 * part2 + part1 * a0* b0* ( 3*a0*a0*(r0-r_ij)*(r0-r_ij) + c0*(1.0-2.0*a0*r0+2.0*a0*r_ij) );
          CEvd_qm =  drdbo * CEvd_qm;

          // compute vdW 
          if (flag_water) {
              d_vdW = twbp->D;
              alf_vdW = twbp->beta;
              re_vdW = twbp->alpha;

              bet_vdW = twbp->repa0;
              gam_vdW = twbp->repr0;
              eta_vdW = twbp->repn;
              del_vdW = twbp->repscal;

              exp1 = d_vdW * exp(-alf_vdW*r_ij+re_vdW);
              exp2 = exp(bet_vdW*pow(r_ij,gam_vdW)+eta_vdW*r_ij+del_vdW);

              dexp1 = -alf_vdW * exp1;
              dexp2 = (bet_vdW*gam_vdW*pow(r_ij,(gam_vdW-1.0))+eta_vdW)*exp2;

              e_vdW = Tap*exp1*exp2;
              de_vdW = drdbo * (dexp1*exp2 + exp1*dexp2);

              CEvd_vdW = Tap * de_vdW + dTap * e_vdW;

          } else {
              d_vdW = twbp->D;
              re_vdW = twbp->r_vdW;
              le = twbp->alpha;

              alf_vdW = twbp->repa0;
              s0 = twbp->repr0;
              s1 = twbp->repn;
              s2 = twbp->repscal;
              s3 = twbp->reps;
              s4 = twbp->lgre/2; //initially multiplied by 2
              s5 = twbp->lgcij;

              L2 = le*le;
              L3 = le*L2;
              L4 = L2*L2;
              L5 = le*L4;

              alf_vdW = alf_vdW/le;
              s1 = s1/le;
              s2 = s2/L2;
              s3 = s3/L3;
              s4 = s4/L4;
              s5 = s5/L5;

              dr = r_ij-re_vdW;
              dr2 = dr*dr;
              dr3 = dr*dr2;
              dr4 = dr*dr3;
              dr5 = dr*dr4;

              exp1 = -d_vdW * exp( -alf_vdW * dr);
              exp2 = s0 + s1*dr + s2*dr2 + s3*dr3 + s4*dr4 + s5*dr5 ;

              dexp1 = -alf_vdW * exp1;
              dexp2 = s1 + 2*s2*dr + 3*s3*dr2 + 4*s4*dr3 + 5*s5*dr4;


              e_vdW = Tap*exp1 * exp2;
              de_vdW = drdbo * (dexp1*exp2 + exp1*dexp2);

              CEvd_vdW = Tap * de_vdW + dTap * e_vdW;
          }
          // lg energy
          if (flag_water) {
              r_ij5 = pow( r_ij, 5.0 );
              r_ij6 = pow( r_ij, 6.0 );
              re6 = pow( twbp->lgre, 6.0 );
              e_lg = -(twbp->lgcij/( r_ij6 + re6 ));
              CEvd_lg = drdbo * ( -6.0*e_lg *r_ij5/(r_ij6 + re6) );
          } else {
              e_lg = 0.0;
              CEvd_lg = 0.0;
          }

          e_b =  e_qm - ( e_vdW + e_lg );
          CEbo = CEvd_qm  - ( CEvd_vdW + CEvd_lg);

        }else{
          e_b = 0.0;
          CEbo = 0.0;
        }
 
      /*e_b = -twbp->De_s * exp(bos*bos) - bos*bos;
      CEbo = -twbp->De_s*2.0*bos*exp(bos*bos)-2.0*bos;*/

      /* calculate the Bond Energy */
      /*data->my_en.e_bond += ebond =
	-twbp->De_s * bo_ij->BO_s * exp_be12
	-twbp->De_p * bo_ij->BO_pi
	-twbp->De_pp * bo_ij->BO_pi2;*/

      data->my_en.e_bond += ebond = 
          e_b - twbp->De_p * bo_ij->BO_pi
          -twbp->De_pp * bo_ij->BO_pi2;


      /* tally into per-atom energy */
      if( system->pair_ptr->evflag)
	system->pair_ptr->ev_tally(i,j,natoms,1,ebond,0.0,0.0,0.0,0.0,0.0);

      /* calculate derivatives of Bond Orders */
      bo_ij->Cdbo += CEbo;
      bo_ij->Cdbopi -= (CEbo + twbp->De_p);
      bo_ij->Cdbopi2 -= (CEbo + twbp->De_pp);

      /* Stabilisation terminal triple bond */
      if( bo_ij->BO >= 1.00 ) {
	if( gp37 == 2 ||
	    (sbp_i->mass == 12.0000 && sbp_j->mass == 15.9990) ||
	    (sbp_j->mass == 12.0000 && sbp_i->mass == 15.9990) ) {
	  exphu = exp( -gp7 * SQR(bo_ij->BO - 2.50) );
	  exphua1 = exp(-gp3 * (workspace->total_bond_order[i]-bo_ij->BO));
	  exphub1 = exp(-gp3 * (workspace->total_bond_order[j]-bo_ij->BO));
	  exphuov = exp(gp4 * (workspace->Delta[i] + workspace->Delta[j]));
	  hulpov = 1.0 / (1.0 + 25.0 * exphuov);

	  estriph = gp10 * exphu * hulpov * (exphua1 + exphub1);
	  data->my_en.e_bond += estriph;

	  decobdbo = gp10 * exphu * hulpov * (exphua1 + exphub1) *
	    ( gp3 - 2.0 * gp7 * (bo_ij->BO-2.50) );
	  decobdboua = -gp10 * exphu * hulpov *
	    (gp3*exphua1 + 25.0*gp4*exphuov*hulpov*(exphua1+exphub1));
	  decobdboub = -gp10 * exphu * hulpov *
	    (gp3*exphub1 + 25.0*gp4*exphuov*hulpov*(exphua1+exphub1));

	  // tally into per-atom energy 
	  if( system->pair_ptr->evflag)
	    system->pair_ptr->ev_tally(i,j,natoms,1,estriph,0.0,0.0,0.0,0.0,0.0);

	  bo_ij->Cdbo += decobdbo;
	  workspace->CdDelta[i] += decobdboua;
	  workspace->CdDelta[j] += decobdboub;
	}
      }
    }
  }
}
