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
#include "rexpon_types.h"
#include "rexpon_nonbonded.h"
#include "rexpon_bond_orders.h"
#include "rexpon_list.h"
#include "rexpon_vector.h"
#include "iostream"

void vdW_Coulomb_Energy( rexpon_system *system, control_params *control,
                         simulation_data *data, storage *workspace,
                         rexpon_list **lists, output_controls *out_control )
{
  int i, j, pj, natoms;
  int start_i, end_i, flag;
  rc_tagint orig_i, orig_j;
  double p_vdW1, p_vdW1i;
  double powr_vdW1, powgi_vdW1;
  double tmp, r_ij, fn13, exp1, exp2;
  double Tap, dTap, dfn13, CEvd, CEclmb, de_core;
  double dr3gamij_1, dr3gamij_3;
  double e_ele, e_vdW, e_core, SMALL = 0.0001;
  double e_lg, de_lg, r_ij5, r_ij6, re6;
  double a0,r0,n0,scal0;
  double dexp1,dexp2,de_vdW;
  double d_vdW,alf_vdW,re_vdW,bet_vdW,gam_vdW,eta_vdW,del_vdW;
  double le,s0,s1,s2,s3,s4,s5,L2,L3,L4,L5;
  double dr,dr2,dr3,dr4,dr5;

  int type_i,type_j;
  int flag_water = 0;
  int flag_QMMM = 0;
  rvec temp, ext_press;
  two_body_parameters *twbp;
  far_neighbor_data *nbr_pj;
  rexpon_list *far_nbrs;

  // Tallying variables:
  double pe_vdw, f_tmp, delij[3];

  natoms = system->n;
  far_nbrs = (*lists) + FAR_NBRS;
  p_vdW1 = system->rexpon_param.gp.l[28];
  p_vdW1i = 1.0 / p_vdW1;
  e_core = 0;
  e_vdW = 0;
  e_lg = de_lg = 0.0;

  for( i = 0; i < natoms; ++i ) {
    if (system->my_atoms[i].type < 0) continue;
    start_i = Start_Index(i, far_nbrs);
    end_i   = End_Index(i, far_nbrs);
    orig_i  = system->my_atoms[i].orig_id;

    for( pj = start_i; pj < end_i; ++pj ) {
      nbr_pj = &(far_nbrs->select.far_nbr_list[pj]);
      j = nbr_pj->nbr;
      if (system->my_atoms[j].type < 0) continue;
      orig_j  = system->my_atoms[j].orig_id;

      flag = 0;
      if(nbr_pj->d <= control->nonb_cut) {
        if (j < natoms) flag = 1;
        else if (orig_i < orig_j) flag = 1;
        else if (orig_i == orig_j) {
          if (nbr_pj->dvec[2] > SMALL) flag = 1;
          else if (fabs(nbr_pj->dvec[2]) < SMALL) {
            if (nbr_pj->dvec[1] > SMALL) flag = 1;
            else if (fabs(nbr_pj->dvec[1]) < SMALL && nbr_pj->dvec[0] > SMALL)
              flag = 1;
          }
        }
      }

      if (flag) {

      r_ij = nbr_pj->d;
      twbp = &(system->rexpon_param.tbp[ system->my_atoms[i].type ]
                                       [ system->my_atoms[j].type ]);

      Tap = workspace->Tap[7] * r_ij + workspace->Tap[6];
      Tap = Tap * r_ij + workspace->Tap[5];
      Tap = Tap * r_ij + workspace->Tap[4];
      Tap = Tap * r_ij + workspace->Tap[3];
      Tap = Tap * r_ij + workspace->Tap[2];
      Tap = Tap * r_ij + workspace->Tap[1];
      Tap = Tap * r_ij + workspace->Tap[0];

      dTap = 7*workspace->Tap[7] * r_ij + 6*workspace->Tap[6];
      dTap = dTap * r_ij + 5*workspace->Tap[5];
      dTap = dTap * r_ij + 4*workspace->Tap[4];
      dTap = dTap * r_ij + 3*workspace->Tap[3];
      dTap = dTap * r_ij + 2*workspace->Tap[2];
      dTap += workspace->Tap[1]/r_ij;




      type_i  = system->my_atoms[i].type;
      type_j  = system->my_atoms[j].type;
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

     flag_QMMM = 0; 
     if ( ( !strcmp( system->rexpon_param.sbp[type_i].name, "HW" ) ||
            !strcmp( system->rexpon_param.sbp[type_i].name, "OW" ) ) ||
          ( !strcmp( system->rexpon_param.sbp[type_j].name, "HW" ) ||
            !strcmp( system->rexpon_param.sbp[type_j].name, "OW" ) ) )
         flag_QMMM = 1;  



      //std::cout<<" flag_water ============== "<<flag_water<<std::endl;

      /*vdWaals Calculations*/
      if(system->rexpon_param.gp.vdw_type==1 || system->rexpon_param.gp.vdw_type==3)
        { // shielding
          powr_vdW1 = pow(r_ij, p_vdW1);
          powgi_vdW1 = pow( 1.0 / twbp->gamma_w, p_vdW1);

          fn13 = pow( powr_vdW1 + powgi_vdW1, p_vdW1i );
          exp1 = exp( twbp->alpha * (1.0 - fn13 / twbp->r_vdW) );
          exp2 = exp( 0.5 * twbp->alpha * (1.0 - fn13 / twbp->r_vdW) );

          e_vdW = twbp->D * (exp1 - 2.0 * exp2);
          data->my_en.e_vdW += Tap * e_vdW;

          dfn13 = pow( powr_vdW1 + powgi_vdW1, p_vdW1i - 1.0) *
            pow(r_ij, p_vdW1 - 2.0);

          CEvd = dTap * e_vdW -
            Tap * twbp->D * (twbp->alpha / twbp->r_vdW) * (exp1 - exp2) * dfn13;
        }
      else{ // no shielding

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

              e_vdW = exp1*exp2;
              de_vdW = (dexp1*exp2 + exp1*dexp2)/r_ij;
    //          std::cout<<"==============e_vdW "<<e_vdW<<std::endl;

              data->my_en.e_vdW += Tap * e_vdW;
              CEvd = Tap * de_vdW + dTap * e_vdW;

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
            

              //std::cout<<"d_vdW re_vdW le alf_vdW s0 s1 s2 s3 s4 s5 "<<d_vdW<<" "<<re_vdW<<" "<< le<<" "<< alf_vdW<<" "<<
                //                                          s0<<" "<<s1<<" "<< s2<<" "<< s3<<" "<< s4<<" "<< s5<<std::endl;

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

              e_vdW = exp1 * exp2;
              de_vdW = (dexp1*exp2 + exp1*dexp2)/r_ij;

              // VERY DANGEROUS change
              // for REQM, we want to exclude i-j interactions that 
              // are all in QM region
              //e_vdW = flag_QMMM * e_vdW;
              //de_vdW = flag_QMMM * de_vdW;

              // temporary
              //e_vdW = 0.0;
              //de_vdW = 0.0;
              data->my_en.e_vdW += Tap * e_vdW;
              CEvd = Tap * de_vdW + dTap * e_vdW; 
              
          }
      }

      if(system->rexpon_param.gp.vdw_type==2 || system->rexpon_param.gp.vdw_type==3)
        { // innner wall
          e_core = twbp->ecore * exp(twbp->acore * (1.0-(r_ij/twbp->rcore)));
          data->my_en.e_vdW += Tap * e_core;

          de_core = -(twbp->acore/twbp->rcore) * e_core;
          CEvd += dTap * e_core + Tap * de_core / r_ij;

          //  lg correction, only if lgvdw is yes
          if (control->lgflag && flag_water) {
            r_ij5 = pow( r_ij, 5.0 );
            r_ij6 = pow( r_ij, 6.0 );
            re6 = pow( twbp->lgre, 6.0 );
            e_lg = -(twbp->lgcij/( r_ij6 + re6 ));
            data->my_en.e_vdW += Tap * e_lg;

            de_lg = -6.0 * e_lg *  r_ij5 / ( r_ij6 + re6 ) ;
            CEvd += dTap * e_lg + Tap * de_lg / r_ij;
          }else{
            e_lg = 0.0;
          }
 
        }

      /*Coulomb Calculations*/
      if (!control->coulomb_off_flag) {
        dr3gamij_1 = ( r_ij * r_ij * r_ij + twbp->gamma );
        dr3gamij_3 = pow( dr3gamij_1 , 0.33333333333333 );

        tmp = Tap / dr3gamij_3;
        data->my_en.e_ele += e_ele =
          C_ele * system->my_atoms[i].q * system->my_atoms[j].q * tmp;

        CEclmb = C_ele * system->my_atoms[i].q * system->my_atoms[j].q *
          ( dTap -  Tap * r_ij / dr3gamij_1 ) / dr3gamij_3;

      } else {
        CEclmb = 0.0;
        e_ele = 0.0;
      }

      /* tally into per-atom energy */
      if( system->pair_ptr->evflag || system->pair_ptr->vflag_atom) {
        pe_vdw = Tap * (e_vdW + e_core + e_lg);
        rvec_ScaledSum( delij, 1., system->my_atoms[i].x,
                              -1., system->my_atoms[j].x );
        f_tmp = -(CEvd + CEclmb);
        system->pair_ptr->ev_tally(i,j,natoms,1,pe_vdw,e_ele,
                        f_tmp,delij[0],delij[1],delij[2]);
      }

      if( control->virial == 0 ) {
        rvec_ScaledAdd( workspace->f[i], -(CEvd + CEclmb), nbr_pj->dvec );
        rvec_ScaledAdd( workspace->f[j], +(CEvd + CEclmb), nbr_pj->dvec );
      }
      else { /* NPT, iNPT or sNPT */
        rvec_Scale( temp, CEvd + CEclmb, nbr_pj->dvec );

        rvec_ScaledAdd( workspace->f[i], -1., temp );
        rvec_Add( workspace->f[j], temp );

        rvec_iMultiply( ext_press, nbr_pj->rel_box, temp );
        rvec_Add( data->my_ext_press, ext_press );
      }
      }
    }
  }

  if (!control->coulomb_off_flag)
    Compute_Polarization_Energy( system, data );
  else
    data->my_en.e_pol = 0.0;

}


void Tabulated_vdW_Coulomb_Energy( rexpon_system *system,control_params *control,
                                   simulation_data *data, storage *workspace,
                                   rexpon_list **lists,
                                   output_controls *out_control )
{
  int i, j, pj, r, natoms;
  int type_i, type_j, tmin, tmax;
  int start_i, end_i, flag;
  rc_tagint orig_i, orig_j;
  double r_ij, base, dif;
  double e_vdW, e_ele;
  double CEvd, CEclmb, SMALL = 0.0001;
  double f_tmp, delij[3];

  rvec temp, ext_press;
  far_neighbor_data *nbr_pj;
  rexpon_list *far_nbrs;
  LR_lookup_table *t;

  natoms = system->n;
  far_nbrs = (*lists) + FAR_NBRS;

  e_ele = e_vdW = 0;

  for( i = 0; i < natoms; ++i ) {
    type_i  = system->my_atoms[i].type;
    if (type_i < 0) continue;
    start_i = Start_Index(i,far_nbrs);
    end_i   = End_Index(i,far_nbrs);
    orig_i  = system->my_atoms[i].orig_id;

    for( pj = start_i; pj < end_i; ++pj ) {
      nbr_pj = &(far_nbrs->select.far_nbr_list[pj]);
      j = nbr_pj->nbr;
      type_j = system->my_atoms[j].type;
      if (type_j < 0) continue;
      orig_j  = system->my_atoms[j].orig_id;

      flag = 0;
      if(nbr_pj->d <= control->nonb_cut) {
        if (j < natoms) flag = 1;
        else if (orig_i < orig_j) flag = 1;
        else if (orig_i == orig_j) {
          if (nbr_pj->dvec[2] > SMALL) flag = 1;
          else if (fabs(nbr_pj->dvec[2]) < SMALL) {
            if (nbr_pj->dvec[1] > SMALL) flag = 1;
            else if (fabs(nbr_pj->dvec[1]) < SMALL && nbr_pj->dvec[0] > SMALL)
              flag = 1;
          }
        }
      }

      if (flag) {

      r_ij   = nbr_pj->d;
      tmin  = MIN( type_i, type_j );
      tmax  = MAX( type_i, type_j );
      t = &( LR[tmin][tmax] );

      /* Cubic Spline Interpolation */
      r = (int)(r_ij * t->inv_dx);
      if( r == 0 )  ++r;
      base = (double)(r+1) * t->dx;
      dif = r_ij - base;

      e_vdW = ((t->vdW[r].d*dif + t->vdW[r].c)*dif + t->vdW[r].b)*dif +
        t->vdW[r].a;

      if (!control->coulomb_off_flag) {
        e_ele = ((t->ele[r].d*dif + t->ele[r].c)*dif + t->ele[r].b)*dif +
          t->ele[r].a;
        e_ele *= system->my_atoms[i].q * system->my_atoms[j].q;
      } else {
        e_ele = 0.0;
      }

      data->my_en.e_vdW += e_vdW;
      data->my_en.e_ele += e_ele;

      CEvd = ((t->CEvd[r].d*dif + t->CEvd[r].c)*dif + t->CEvd[r].b)*dif +
        t->CEvd[r].a;

      if (!control->coulomb_off_flag) {
        CEclmb = ((t->CEclmb[r].d*dif+t->CEclmb[r].c)*dif+t->CEclmb[r].b)*dif +
          t->CEclmb[r].a;
        CEclmb *= system->my_atoms[i].q * system->my_atoms[j].q;
      } else {
        CEclmb = 0.0;
      }

      /* tally into per-atom energy */
      if( system->pair_ptr->evflag || system->pair_ptr->vflag_atom) {
        rvec_ScaledSum( delij, 1., system->my_atoms[i].x,
                              -1., system->my_atoms[j].x );
        f_tmp = -(CEvd + CEclmb);
        system->pair_ptr->ev_tally(i,j,natoms,1,e_vdW,e_ele,
                        f_tmp,delij[0],delij[1],delij[2]);
      }

      if( control->virial == 0 ) {
        rvec_ScaledAdd( workspace->f[i], -(CEvd + CEclmb), nbr_pj->dvec );
        rvec_ScaledAdd( workspace->f[j], +(CEvd + CEclmb), nbr_pj->dvec );
      }
      else { // NPT, iNPT or sNPT
        rvec_Scale( temp, CEvd + CEclmb, nbr_pj->dvec );

        rvec_ScaledAdd( workspace->f[i], -1., temp );
        rvec_Add( workspace->f[j], temp );

        rvec_iMultiply( ext_press, nbr_pj->rel_box, temp );
        rvec_Add( data->my_ext_press, ext_press );
      }
      }
    }
  }

  if (!control->coulomb_off_flag)
    Compute_Polarization_Energy( system, data );
  else
    data->my_en.e_pol = 0.0;
}



void Compute_Polarization_Energy( rexpon_system *system, simulation_data *data )
{
  int  i, type_i;
  double q, en_tmp;

  data->my_en.e_pol = 0.0;
  for( i = 0; i < system->n; i++ ) {
    type_i = system->my_atoms[i].type;
    if (type_i < 0) continue;
    q = system->my_atoms[i].q;

    en_tmp = KCALpMOL_to_EV * (system->rexpon_param.sbp[type_i].chi * q +
                (system->rexpon_param.sbp[type_i].eta / 2.) * SQR(q));
    data->my_en.e_pol += en_tmp;

    /* tally into per-atom energy */
    if( system->pair_ptr->evflag)
      system->pair_ptr->ev_tally(i,i,system->n,1,0.0,en_tmp,0.0,0.0,0.0,0.0);
  }
}

void LR_vdW_Coulomb( rexpon_system *system, storage *workspace,
        control_params *control, int i, int j, double r_ij, LR_data *lr )
{
  double p_vdW1 = system->rexpon_param.gp.l[28];
  double p_vdW1i = 1.0 / p_vdW1;
  double powr_vdW1, powgi_vdW1;
  double tmp, fn13, exp1, exp2;
  double Tap, dTap, dfn13;
  double dr3gamij_1, dr3gamij_3;
  double e_core, de_core;
  double e_lg, de_lg, r_ij5, r_ij6, re6;
  double a0,r0,n0,scal0;
  two_body_parameters *twbp;

  twbp = &(system->rexpon_param.tbp[i][j]);
  e_core = 0;
  de_core = 0;
  e_lg = de_lg = 0.0;

  /* calculate taper and its derivative */
  Tap = workspace->Tap[7] * r_ij + workspace->Tap[6];
  Tap = Tap * r_ij + workspace->Tap[5];
  Tap = Tap * r_ij + workspace->Tap[4];
  Tap = Tap * r_ij + workspace->Tap[3];
  Tap = Tap * r_ij + workspace->Tap[2];
  Tap = Tap * r_ij + workspace->Tap[1];
  Tap = Tap * r_ij + workspace->Tap[0];

  dTap = 7*workspace->Tap[7] * r_ij + 6*workspace->Tap[6];
  dTap = dTap * r_ij + 5*workspace->Tap[5];
  dTap = dTap * r_ij + 4*workspace->Tap[4];
  dTap = dTap * r_ij + 3*workspace->Tap[3];
  dTap = dTap * r_ij + 2*workspace->Tap[2];
  dTap += workspace->Tap[1]/r_ij;

  /*vdWaals Calculations*/
  if(system->rexpon_param.gp.vdw_type==1 || system->rexpon_param.gp.vdw_type==3)
    { // shielding
      powr_vdW1 = pow(r_ij, p_vdW1);
      powgi_vdW1 = pow( 1.0 / twbp->gamma_w, p_vdW1);

      fn13 = pow( powr_vdW1 + powgi_vdW1, p_vdW1i );
      exp1 = exp( twbp->alpha * (1.0 - fn13 / twbp->r_vdW) );
      exp2 = exp( 0.5 * twbp->alpha * (1.0 - fn13 / twbp->r_vdW) );

      lr->e_vdW = Tap * twbp->D * (exp1 - 2.0 * exp2);

      dfn13 = pow( powr_vdW1 + powgi_vdW1, p_vdW1i-1.0) * pow(r_ij, p_vdW1-2.0);

      lr->CEvd = dTap * twbp->D * (exp1 - 2.0 * exp2) -
        Tap * twbp->D * (twbp->alpha / twbp->r_vdW) * (exp1 - exp2) * dfn13;
    }
  else{ // no shielding
      r0 = twbp->repr0;
      a0 = twbp->repa0;
      n0 = twbp->repn;
      scal0 = twbp->repscal;
      exp1 = twbp->D*exp(-twbp->beta*r_ij+twbp->alpha);
      exp2 = exp(-a0*pow(r_ij,r0)+n0*r_ij+scal0);
      lr->e_vdW = Tap * exp1*exp2;
      lr->CEvd = dTap * exp1*exp2 - Tap*twbp->beta*exp1*exp2/r_ij - 
                 Tap*(a0*r0*pow(r_ij,(r0-1.0))-n0)*exp1*exp2/r_ij;

    /*lr->e_vdW = Tap * twbp->D * (exp1 - 2.0 * exp2);
    lr->CEvd = dTap * twbp->D * (exp1 - 2.0 * exp2) -
      Tap * twbp->D * (twbp->alpha / twbp->r_vdW) * (exp1 - exp2) / r_ij;*/
  }

  if(system->rexpon_param.gp.vdw_type==2 || system->rexpon_param.gp.vdw_type==3)
    { // innner wall
      e_core = twbp->ecore * exp(twbp->acore * (1.0-(r_ij/twbp->rcore)));
      lr->e_vdW += Tap * e_core;

      de_core = -(twbp->acore/twbp->rcore) * e_core;
      lr->CEvd += dTap * e_core + Tap * de_core / r_ij;

      //  lg correction, only if lgvdw is yes
      if (control->lgflag) {
        r_ij5 = pow( r_ij, 5.0 );
        r_ij6 = pow( r_ij, 6.0 );
        re6 = pow( twbp->lgre, 6.0 );
        e_lg = -(twbp->lgcij/( r_ij6 + re6 ));
        lr->e_vdW += Tap * e_lg;

        de_lg = -6.0 * e_lg *  r_ij5 / ( r_ij6 + re6 ) ;
        lr->CEvd += dTap * e_lg + Tap * de_lg/r_ij;
      }

    }


  /* Coulomb calculations */
  dr3gamij_1 = ( r_ij * r_ij * r_ij + twbp->gamma );
  dr3gamij_3 = pow( dr3gamij_1 , 0.33333333333333 );

  tmp = Tap / dr3gamij_3;
  lr->H = EV_to_KCALpMOL * tmp;
  lr->e_ele = C_ele * tmp;

  lr->CEclmb = C_ele * ( dTap -  Tap * r_ij / dr3gamij_1 ) / dr3gamij_3;
}
