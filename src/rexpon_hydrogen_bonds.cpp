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
#include "rexpon_hydrogen_bonds.h"
#include "rexpon_bond_orders.h"
#include "rexpon_list.h"
#include "rexpon_valence_angles.h"
#include "rexpon_vector.h"
#include "iostream"
void Hydrogen_Bonds( rexpon_system *system, control_params *control,
                     simulation_data *data, storage *workspace,
                     rexpon_list **lists, output_controls *out_control )
{
  int  i, j, k, pi, pk;
  int  type_i, type_j, type_k;
  int  start_j, end_j, hb_start_j, hb_end_j;
  int  hblist[MAX_BONDS];
  int  itr, top;
  int  num_hb_intrs = 0;
  ivec rel_jk;
  double r_jk, theta, cos_theta, sin_xhz2,sin_xhz4, sin_xhz6,cos_xhz1, sin_theta2;
  double e_hb, exp_hb1,exp_hb2, exp_hb3,exp_hb4, CEhb1, CEhb2, CEhb3,CEhb3_cut;
  rvec dcos_theta_di, dcos_theta_dj, dcos_theta_dk;
  rvec dvec_jk, force, ext_press;
  hbond_parameters *hbp;
  bond_order_data *bo_ij;
  bond_data *pbond_ij;
  far_neighbor_data *nbr_jk;
  rexpon_list *bonds, *hbonds;
  bond_data *bond_list;
  hbond_data *hbond_list;

  // tally variables
  double fi_tmp[3], fk_tmp[3], delij[3], delkj[3],delki[3];

  // HBLP
  int l, m, n, kiter, pk_b, start_k, end_k;
  int pi_b, start_i, end_i;

  double eng_tmp,e_lp, norm, d_km, d_kn, d_ir,d_is,d_alpha, d_beta, d_bis, d_per, r_ik;
  double r_alpha, r_beta;
  double sin_xhz8,sin_xhz10,sin_xhz12;
  double a1, a2, a3, a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,pow15;
  double ang_exp,hb_ang_exp, r1, c1, c2, a0,d0,r0,lp_scal,lp_exp,exp_lp,exp_lp1,exp_lp2;
  double exp_lp1_cut,exp_hb1_cut,exp5_cut,exp_lp_cut,exp_hb_cut;
  double term1,term2;
  double b_1,b_2,b_3,b_4,b_5,b_6,b_7,b_8,b_9,b_10,b_11,b_12,b_13,b_14,b_15,b_16,b_17,exp_hb,exp_angle;
  double exp_alpha,exp_alpha2,exp_beta;
  double CElp1, CElp2, CElp3,CElp4_cut, lp_cos,hb_exp,hb12,hb10,max_angle, min_angle, min_angle2 ;
  double alpha, beta, cos_alpha, cos_beta;
  double norm_bis_d,norm_bis_a,tetlp,cos_tetlp;
  rvec dcos_alpha_dj,dcos_alpha_dk, dcos_alpha_da ;
  rvec dcos_beta_dj,dcos_beta_dk, dcos_beta_da, dcos_beta_db ;
  rvec dvec_km, dvec_kn,dvec_alpha,dvec_beta,dvec_alpha1,dvec_beta1;
  rvec vec_A1A2,vec_A1B2,vec_B1A2,vec_B1B2;
  rvec vec_ab, vec_ab1, vec_ABm1, vec_ABm, vec_ab_inv,vec_ab1_inv;
  rvec dvec_is,dvec_ir;
  rvec dvec_bis, dvec_per;
  rvec dvec_bis1, dvec_per1;
  rvec dcos_alpha_dm, dcos_alpha_dn, dcos_beta_dm, dcos_beta_dn;
  rvec dvec_jk_inv, dvec_ik,dvec_ki;
  rvec dcos_alpha_di, dcos_beta_di;
  bond_order_data *bo_kl,*bo_il;
  bond_data *pbond_kl,*pbond_il;
  rvec dvec_kj,dvec_bis_d,dvec_bis_a;
  rvec dvec_ji;
  double alphaH1,betaH1,cos_alphaH1,cos_betaH1,sin_alpha6,sin_beta6,sin_alpha2,sin_beta2;
  double thetaLP1, thetaLP2, cos_thetaLP1, cos_thetaLP2;
  double fm_tmp[3], fn_tmp[3],delik[3], deljk[3], delnk[3], delmk[3];
  double fj_tmp[3] ;
  double CElp4,exp5,exp6,exp6lp,ehb,elp;
  double exp_rep1,exp_rep2;
  int loopi=0;
  double A1A2,A1B2,B1A2,B1B2;
  double norm_alpha,norm_beta,norm_alpha1,norm_beta1;
  double norm_ab,norm_ab1,norm_ABm, norm_ABm1;
  double r_cut = control->hbond_cut;
  double r_cut2;
  double TapO,dTapO,TapH,dTapH;
  int flag_water = 0;
  // END HBLP

  bonds = (*lists) + BONDS;
  bond_list = bonds->select.bond_list;
  hbonds = (*lists) + HBONDS;
  hbond_list = hbonds->select.hbond_list;

  for( j = 0; j < system->n; ++j )
    if( system->rexpon_param.sbp[system->my_atoms[j].type].p_hbond == 1 ) {
      type_j     = system->my_atoms[j].type;
      start_j    = Start_Index(j, bonds);
      end_j      = End_Index(j, bonds);
      hb_start_j = Start_Index( system->my_atoms[j].Hindex, hbonds );
      hb_end_j   = End_Index( system->my_atoms[j].Hindex, hbonds );
      if (type_j < 0) continue;

      top = 0;
      for( pi = start_j; pi < end_j; ++pi )  {
        pbond_ij = &( bond_list[pi] );
        i = pbond_ij->nbr;
        type_i = system->my_atoms[i].type;
	if (type_i < 0) continue;
        bo_ij = &(pbond_ij->bo_data);

        if( system->rexpon_param.sbp[type_i].p_hbond == 2 &&
            bo_ij->BO >= HB_THRESHOLD )
          hblist[top++] = pi;
      }

      for( pk = hb_start_j; pk < hb_end_j; ++pk ) {
        /* set k's varibles */
        k = hbond_list[pk].nbr;
        type_k = system->my_atoms[k].type;
	if (type_k < 0) continue;
        nbr_jk = hbond_list[pk].ptr;
        r_jk = nbr_jk->d;
        rvec_Scale( dvec_jk, hbond_list[pk].scl, nbr_jk->dvec );

        for( itr = 0; itr < top; ++itr ) {
          pi = hblist[itr];
          pbond_ij = &( bonds->select.bond_list[pi] );
          i = pbond_ij->nbr;

          if( system->my_atoms[i].orig_id != system->my_atoms[k].orig_id ) {
            bo_ij = &(pbond_ij->bo_data);
            type_i = system->my_atoms[i].type;
	    if (type_i < 0) continue;
            //hbp = &(system->rexpon_param.hbp[ type_i ][ type_j ][ type_k ]);
            hbp = &(system->rexpon_param.hbp[ type_j ][ type_i ][ type_k ]);

	    if (hbp->r0_hb <= 0.0) continue;
            ++num_hb_intrs;


            /* get Z-H vector*/
            rvec_Scale(dvec_ik, 1.0, dvec_jk);
            rvec_ScaledAdd(dvec_ik, -1.0,pbond_ij->dvec);
            r_ik = rvec_Norm( dvec_ik );

            rvec_Scale(dvec_ji, -1.0,pbond_ij->dvec);


            // taper
            TapO = workspace->TapOO[7] * r_jk + workspace->TapOO[6];
            TapO = TapO * r_jk + workspace->TapOO[5];
            TapO = TapO * r_jk + workspace->TapOO[4];
            TapO = TapO * r_jk + workspace->TapOO[3];
            TapO = TapO * r_jk + workspace->TapOO[2];
            TapO = TapO * r_jk + workspace->TapOO[1];
            TapO = TapO * r_jk + workspace->TapOO[0];

            TapH = workspace->TapHH[7] * r_ik + workspace->TapHH[6];
            TapH = TapH * r_ik + workspace->TapHH[5];
            TapH = TapH * r_ik + workspace->TapHH[4];
            TapH = TapH * r_ik + workspace->TapHH[3];
            TapH = TapH * r_ik + workspace->TapHH[2];
            TapH = TapH * r_ik + workspace->TapHH[1];
            TapH = TapH * r_ik + workspace->TapHH[0];

            dTapO = 7*workspace->TapOO[7] * r_jk + 6*workspace->TapOO[6];
            dTapO = dTapO * r_jk + 5*workspace->TapOO[5];
            dTapO = dTapO * r_jk + 4*workspace->TapOO[4];
            dTapO = dTapO * r_jk + 3*workspace->TapOO[3];
            dTapO = dTapO * r_jk + 2*workspace->TapOO[2];
            dTapO = dTapO * r_jk + 1*workspace->TapOO[1];

            dTapH = 7*workspace->TapHH[7] * r_ik + 6*workspace->TapHH[6];
            dTapH = dTapH * r_ik + 5*workspace->TapHH[5];
            dTapH = dTapH * r_ik + 4*workspace->TapHH[4];
            dTapH = dTapH * r_ik + 3*workspace->TapHH[3];
            dTapH = dTapH * r_ik + 2*workspace->TapHH[2];
            dTapH = dTapH * r_ik + 1*workspace->TapHH[1];


            // get HB params 
            b_1 = hbp->p_hb1*0.25;
            b_2 = hbp->p_hb2;
            b_3 = hbp->p_hb3;
            b_4 = hbp->p_hb4;
            b_5 = hbp->p_hb5;
            b_13 = hbp->p_hb6;
            b_14 = hbp->p_hb7;
            b_15 = hbp->p_hb8;
            b_16 = hbp->p_hb9;
            b_17 = hbp->p_hb10;

            a3 = hbp->p_hb16;
            a4 = hbp->p_hb17;
            a5 = hbp->p_hb18;
            a6 = hbp->p_hb19;
            a7 = hbp->p_hb20;
            a11 = hbp->p_hb21;
            a12 = hbp->p_hb22;
            a13 = hbp->p_hb23;



            /* get Z-H vector*/
            rvec_Scale(dvec_ki, -1.0, dvec_ik);


            flag_water = 0;
            /*if ( !strcmp( system->rexpon_param.sbp[type_k].name, "HW" ) &&
                 !strcmp( system->rexpon_param.sbp[type_j].name, "OW" )  )
                 flag_water = 1;
            if ( !strcmp( system->rexpon_param.sbp[type_k].name, "OW" ) &&
                 !strcmp( system->rexpon_param.sbp[type_j].name, "HW" ) )
                 flag_water = 1;
            if ( !strcmp( system->rexpon_param.sbp[type_k].name, "OW" ) &&
                 !strcmp( system->rexpon_param.sbp[type_j].name, "OW" ) )
                 flag_water = 1;
            if ( !strcmp( system->rexpon_param.sbp[type_k].name, "HW" ) &&
                 !strcmp( system->rexpon_param.sbp[type_j].name, "HW" ) )
                 flag_water = 1;*/


            //if ( flag_water) {

            start_k    = Start_Index(k, bonds);
            end_k      = End_Index(k, bonds);
            /* In water there are 2 H atoms bonded to each O atom*/
            // acceptor hydrogens 
            kiter = 0;
            for ( pk_b = start_k; pk_b < end_k; ++pk_b )  {
              pbond_kl = &( bonds->select.bond_list[pk_b] );
              l = pbond_kl->nbr;
              bo_kl = &(pbond_kl->bo_data);
              if ( bo_kl->BO > 0.0 ) {
                  if (kiter == 0 ) {
                     m = l;
                     rvec_Copy(dvec_km, pbond_kl->dvec);
                     d_km = pbond_kl->d;
                     kiter += 1;
                  } else if (kiter == 1) {
                     n = l;
                     rvec_Copy(dvec_kn, pbond_kl->dvec);
                     d_kn = pbond_kl->d;
                     kiter += 1;
                  }
              }
            }

            if (kiter == 2) {
                flag_water = 1;
            }

            if ( flag_water ) {
            /* Note: k--j-i
               pbond_ij->dvec = atom[i] - atom[j]
               dvec_jk = atom[k] - atom[j]
               dvec_km = atom[m] - atom[k]
               dvec_kn = atom[n] - atom[k]
               lone pair unit vectors */


            Calculate_Theta( dvec_ji, pbond_ij->d, dvec_ik, r_ik,
                             &theta, &cos_theta );

            /* the derivative of cos(theta) */
            Calculate_dCos_Theta( dvec_ji, pbond_ij->d, dvec_ik, r_ik,
                                  &dcos_theta_dj, &dcos_theta_di,
                                  &dcos_theta_dk );




            Calculate_LonePair_Vecs( dvec_kn, d_kn, dvec_km, d_km,
                                     &dvec_bis, &dvec_per,
                                     &dvec_alpha, &dvec_beta);

            /* get Z-H vector*/
            rvec_Scale(dvec_ki, -1.0, dvec_ik);

            r_alpha = rvec_Norm( dvec_alpha );
            r_beta = rvec_Norm( dvec_beta );
            /*calculate alpha and beta lone pair angles*/
            Calculate_Theta( dvec_alpha, r_alpha, dvec_ki, r_ik,
                             &alpha, &cos_alpha );
            Calculate_Theta( dvec_beta, r_beta, dvec_ki, r_ik,
                             &beta, &cos_beta );


            /* derivative of alpha and beta w.r.t k and j*/
            Calculate_dCos_Theta( dvec_ki, r_ik, dvec_alpha, r_alpha,
                                  &dcos_alpha_di, &dcos_alpha_dk,
                                  &dcos_alpha_da );
            Calculate_dCos_Theta( dvec_ki, r_ik, dvec_beta, r_beta,
                                  &dcos_beta_di, &dcos_beta_dk,
                                  &dcos_beta_db );

            /* since alpha and beta do not represent atoms the derivetive
                is translated to atoms m and n */
            Calculate_dCos_Alpha( dvec_bis, dvec_per, 1,
                                  dvec_km, d_km, dvec_kn, d_kn,
                                  dvec_alpha, dvec_ki, r_ik,
                                  &dcos_alpha_dm,&dcos_alpha_dn );
            Calculate_dCos_Alpha( dvec_bis, dvec_per, -1,
                                  dvec_km, d_km, dvec_kn, d_kn,
                                  dvec_beta, dvec_ki, r_ik,
                                  &dcos_beta_dm,&dcos_beta_dn );

            /* the force on atom k is equal to the sum of all forces
                on atom i, m and n but in opposite direction */
            rvec_Scale(dcos_alpha_dk, -1.0, dcos_alpha_di);
            rvec_ScaledAdd(dcos_alpha_dk, -1.0, dcos_alpha_dm);
            rvec_ScaledAdd(dcos_alpha_dk, -1.0, dcos_alpha_dn);
            rvec_Scale(dcos_beta_dk, -1.0, dcos_beta_di);
            rvec_ScaledAdd(dcos_beta_dk, -1.0, dcos_beta_dm);
            rvec_ScaledAdd(dcos_beta_dk, -1.0, dcos_beta_dn);


            /* hyrogen bond energy*/
            sin_theta2 = sin( theta/2.0 );
            sin_xhz2 = SQR(sin_theta2);
   
             

            // HBrep 

            exp_lp1 = b_1*exp(-b_2*r_jk+b_3)*exp(b_4/pow(r_jk,b_5));
            exp_lp2 = 1.0 + b_14*sin_xhz2;

            exp_lp = exp_lp1 * exp_lp2; 

            // HBatt
            exp5 = exp(-2*a5*(r_ik-a6));
            exp_hb1 = a3*(a4-pow(exp5-a7,2.0));
            sin_alpha2 = SQR(sin(alpha/2));
            sin_beta2 = SQR(sin(beta/2));
            sin_alpha6 = sin_alpha2*sin_alpha2*sin_alpha2;
            sin_beta6 = sin_beta2*sin_beta2*sin_beta2;
            exp6 = exp(-a12*(sin_alpha6+sin_beta6)+a13); 
            exp_hb2 = a11/(1+exp6); 

            exp_hb = exp_hb1 * exp_hb2;


            CEhb1 = 0.0;
            CEhb2 = -0.5*(b_14)* exp_lp1*TapO ;
            //CElp1 = -1.5*b_15/b_14*sin_alpha2*sin_alpha2*exp_lp2*exp_lp2*exp6lp* (exp_lp1*TapO );
            //CElp2 = -1.5*b_15/b_14*sin_beta2*sin_beta2*exp_lp2*exp_lp2*exp6lp*   (exp_lp1*TapO );
 
            CEhb3 = -exp_lp1*(b_2+b_4*b_5*pow(r_jk,-b_5-1))*exp_lp2*TapO + dTapO*exp_lp;

            //     std::cout<<"ehb10 "<< elp<<" "<<ehb<<" "<<r_ik<<std::endl;


           if ( r_ik > (r_cut-0.96)) {
              exp_hb = 0.0;
              exp_hb1 = 0.0;
              exp_hb2 = 0.0;
            }

            CElp1 = -1.5*a12/a11*sin_alpha2*sin_alpha2*exp_hb2*exp_hb2*exp6* (exp_hb1*TapH );
            CElp2 = -1.5*a12/a11*sin_beta2*sin_beta2*exp_hb2*exp_hb2*exp6*   (exp_hb1*TapH );
            CElp3 = 0.0;
            CElp4 = 4*a3*a5*exp5*(exp5-a7)*TapH*exp_hb2 + dTapH*exp_hb;


            // with Taper
            elp = exp_lp * TapO;
            ehb = exp_hb * TapH;

            //if (RAD2DEG(theta) < hbp->r0_hb) {
              //ehb = 0.0;
            //}


            //std::cout<<"ehb10 "<< elp<<" "<<ehb<<" "<<r_ik<<std::endl;


            data->my_en.e_hb += ehb;
            data->my_en.e_lp += elp;

            } else {

            Calculate_Theta( dvec_ji, pbond_ij->d, dvec_ik, r_ik,
                             &theta, &cos_theta );

            /* the derivative of cos(theta) */
            Calculate_dCos_Theta( dvec_ji, pbond_ij->d, dvec_ik, r_ik,
                                  &dcos_theta_dj, &dcos_theta_di,
                                  &dcos_theta_dk );


                 //std::cout<<"here1 "<<RAD2DEG(theta)<<" "<<i+1<<" "<<j+1<<" "<<k+1<<" "<<system->my_atoms[i].x[0]<<" "<<system->my_atoms[j].x[0]<<" "<<system->my_atoms[k].x[0]<<" "<<std::endl;
                 b_1 = 8.0; //hbp->p_hb1;
                 b_2 = 1.7241379; //hbp->p_hb2;
                 b_3 = 2.930; //hbp->p_hb3;
                 b_4 = 2.0; //hbp->p_hb4;

                 exp_lp1 = b_1 * ( exp(-2*b_2*(r_jk-b_3)) - 2*exp(-b_2*(r_jk-b_3)) );
                 exp_lp2 = pow(cos_theta,b_4);
                 exp_lp = exp_lp1 * exp_lp2;

                 elp = exp_lp * TapO;
                 ehb = 0.0;

                 //temporary
                 if (RAD2DEG(theta) < 90.0) {
                     exp_lp1 = 0.0;
                     exp_lp1 = 0.0;
                     exp_lp = 0.0;
                     elp = 0.0;
                 }

                 data->my_en.e_hb += elp;
                 data->my_en.e_lp += ehb;
                 //std::cout<<"ehb "<< elp<<std::endl;

                 CEhb1 = 0.0;
                 CEhb2 = b_4*pow(cos_theta,b_4-1)*exp_lp1*TapO;
                 CEhb3 = -2*b_1*b_2*(exp(-2*b_2*(r_jk-b_3)) - exp(-b_2*(r_jk-b_3)) )*TapO*exp_lp2 + dTapO*exp_lp;
            }



            /* hydrogen bond forces */
            bo_ij->Cdbo += CEhb1; // dbo term

            if( control->virial == 0 ) {
                 //std::cout<<"here2 "<<RAD2DEG(theta)<<std::endl;

              // dcos terms
              rvec_ScaledAdd( workspace->f[i], +CEhb2, dcos_theta_di );
              rvec_ScaledAdd( workspace->f[j], +CEhb2, dcos_theta_dj );
              rvec_ScaledAdd( workspace->f[k], +CEhb2, dcos_theta_dk );
              // dr terms
              rvec_ScaledAdd( workspace->f[j], -CEhb3/r_jk, dvec_jk );
              rvec_ScaledAdd( workspace->f[k], +CEhb3/r_jk, dvec_jk );
              if (flag_water) {
              // LP terms
              rvec_ScaledAdd( workspace->f[i], +CElp3, dcos_theta_di );
              rvec_ScaledAdd( workspace->f[j], +CElp3, dcos_theta_dj );
              rvec_ScaledAdd( workspace->f[k], +CElp3, dcos_theta_dk );
             
              rvec_ScaledAdd( workspace->f[i], +CElp1, dcos_alpha_di );
              rvec_ScaledAdd( workspace->f[k], +CElp1, dcos_alpha_dk );
              rvec_ScaledAdd( workspace->f[m], +CElp1, dcos_alpha_dm );
              rvec_ScaledAdd( workspace->f[n], +CElp1, dcos_alpha_dn );
              rvec_ScaledAdd( workspace->f[i], +CElp2, dcos_beta_di );
              rvec_ScaledAdd( workspace->f[k], +CElp2, dcos_beta_dk );
              rvec_ScaledAdd( workspace->f[m], +CElp2, dcos_beta_dm );
              rvec_ScaledAdd( workspace->f[n], +CElp2, dcos_beta_dn );
              rvec_ScaledAdd( workspace->f[i], -CElp4/r_ik, dvec_ik );
              rvec_ScaledAdd( workspace->f[k], +CElp4/r_ik, dvec_ik );
              }  

            }
            else {
                 //std::cout<<"here3 "<<RAD2DEG(theta)<<std::endl;

              rvec_Scale( force, +CEhb2, dcos_theta_di ); // dcos terms
              rvec_Add( workspace->f[i], force );
              rvec_iMultiply( ext_press, pbond_ij->rel_box, force );
              rvec_ScaledAdd( data->my_ext_press, 1.0, ext_press );

              rvec_ScaledAdd( workspace->f[j], +CEhb2, dcos_theta_dj );

              ivec_Scale( rel_jk, hbond_list[pk].scl, nbr_jk->rel_box );
              rvec_Scale( force, +CEhb2, dcos_theta_dk );
              rvec_Add( workspace->f[k], force );
              rvec_iMultiply( ext_press, rel_jk, force );
              rvec_ScaledAdd( data->my_ext_press, 1.0, ext_press );
              // dr terms
              rvec_ScaledAdd( workspace->f[j], -CEhb3/r_jk, dvec_jk );

              rvec_Scale( force, CEhb3/r_jk, dvec_jk );
              rvec_Add( workspace->f[k], force );
              rvec_iMultiply( ext_press, rel_jk, force );
              rvec_ScaledAdd( data->my_ext_press, 1.0, ext_press );
            }

            /* tally into per-atom virials */
            if (system->pair_ptr->vflag_atom || system->pair_ptr->evflag) {
              rvec_ScaledSum( delij,-1., system->my_atoms[i].x,
                                     1., system->my_atoms[j].x );
              rvec_ScaledSum( delki,-1., system->my_atoms[i].x,
                                     1., system->my_atoms[k].x );

              rvec_Scale(fj_tmp, CEhb2, dcos_theta_dj);
              rvec_Scale(fk_tmp, CEhb2, dcos_theta_dk);
              rvec_ScaledAdd(fk_tmp, CEhb3/r_jk, dvec_jk);
              rvec_ScaledAdd(fj_tmp, -CEhb3/r_jk, dvec_jk);

              system->pair_ptr->ev_tally3(j,i,k,elp,0.0,fj_tmp,fk_tmp,delij,delki);
                 //std::cout<<"here4 "<<RAD2DEG(theta)<<std::endl;


              if (flag_water) {
                 //std::cout<<"here5 "<<RAD2DEG(theta)<<std::endl;

              //rvec_ScaledSum( delik, 1., system->my_atoms[i].x,
                                    //-1., system->my_atoms[k].x );
              rvec_ScaledSum( delik, 1., system->my_atoms[i].x,
                                    -1., system->my_atoms[k].x );
              rvec_ScaledSum( delmk, 1., system->my_atoms[m].x,
                                    -1., system->my_atoms[k].x );
              rvec_ScaledSum( delnk, 1., system->my_atoms[n].x,
                                    -1., system->my_atoms[k].x );
              //rvec_Scale(fi_tmp, CElp3, dcos_theta_di);
              rvec_Scale(fi_tmp, -CElp4/r_ik, dvec_ik);
              rvec_ScaledAdd(fi_tmp, CElp3, dcos_theta_di);
              rvec_ScaledAdd(fi_tmp, CElp1, dcos_alpha_di);
              rvec_ScaledAdd(fi_tmp, CElp2, dcos_beta_di);
              rvec_Scale(fm_tmp, CElp1, dcos_alpha_dm);
              rvec_ScaledAdd(fm_tmp, CElp2, dcos_beta_dm);
              rvec_Scale(fn_tmp, CElp1, dcos_alpha_dn);
              rvec_ScaledAdd(fn_tmp, CElp2, dcos_beta_dn);

              //system->pair_ptr->ev_tally5(i,j,m,n,k,ehb,fi_tmp,fj_tmp,fm_tmp,fn_tmp,
                //                          delik,deljk,delmk,delnk);

              system->pair_ptr->ev_tally4(i,m,n,k,ehb,fi_tmp,fm_tmp,fn_tmp,
                                          delik,delmk,delnk);
              }
            }
          }
        }
      }
    }
}
/*------------------------------------------------------------------------------*/
void Calculate_LonePair_Vecs( rvec dvec_H1, double d_H1, rvec dvec_H2, double d_H2,
                              rvec* dvec_bis, rvec* dvec_per,
                              rvec* dvec_alpha, rvec* dvec_beta) 
{

  int t;
  double theta_0, theta_1, norm;
  rvec dvec_a,dvec_b,dvec_c, dvec_normal; 
  rvec dvec_bis_tmp, dvec_per_tmp;
  rvec dvec_alpha_tmp, dvec_beta_tmp ;
  double norm_b, norm_p;
  // unit vectors of a and b
  rvec_Scale( dvec_a, -1/d_H1, dvec_H1 );
  rvec_Scale( dvec_b, -1/d_H2, dvec_H2 );

  // bisector of a and b
  rvec_Sum( dvec_bis_tmp, dvec_a, dvec_b );
  norm_b = rvec_Norm( dvec_bis_tmp );
  //rvec_Scale( dvec_bis_tmp, 1/norm, dvec_bis_tmp );

  // cross product of a and b 
  rvec_Scale( dvec_a, 1/d_H1, dvec_H1 );
  rvec_Scale( dvec_b, 1/d_H2, dvec_H2 );
  rvec_Cross( dvec_per_tmp, dvec_a, dvec_b );
  norm_p = rvec_Norm( dvec_per_tmp );
  //rvec_Scale( dvec_per_tmp, 1/norm, dvec_per_tmp );

  // cross product of dvec_per qnd dvec_bis
  rvec_Cross( dvec_normal, dvec_per_tmp, dvec_bis_tmp );
  //norm = rvec_Norm( dvec_normal );
  //rvec_Scale( dvec_normal, 1/norm, dvec_normal );
  
  // Eq1 = angle btw lp and dvec_bis is 109.5/2   
  // Eq2 = angle btw lp and dvec_per is 90.0 - 109.5/2  
  // Eq3 = lp is in the plane of dvec_bis and dvec_per 
   
  theta_0 = DEG2RAD( 54.75 );
  theta_1 = DEG2RAD( 90.0 - 54.75 );

  // solve Eq1-3 with Guassian ellimination method AX = C
  dvec_c[0] = cos(theta_1) * norm_p;
  dvec_c[1] = cos(theta_0) * norm_b;
  dvec_c[2] = 0.0;

  // 1st lone pair vec
  Solve_Linear_Guassian(dvec_per_tmp, dvec_bis_tmp, 
                        dvec_normal, dvec_c, 
                        &dvec_alpha_tmp);
  // 2nd lone pair vec 
  rvec_Scale( dvec_per_tmp, -1.0, dvec_per_tmp );
  Solve_Linear_Guassian(dvec_per_tmp, dvec_bis_tmp,
                        dvec_normal, dvec_c,
                        &dvec_beta_tmp);

  // collect all vectors
  for( t = 0; t < 3; ++t ) {
    (*dvec_bis)[t] = dvec_bis_tmp[t];
    (*dvec_per)[t] = dvec_per_tmp[t];
    (*dvec_alpha)[t] = dvec_alpha_tmp[t];
    (*dvec_beta)[t] = dvec_beta_tmp[t];
  }

}
/*---------------------------------------------------------------------------*/
//void Solve_Linear_Guassian(rvec a1, rvec a2, rvec a3, rvec b, 
//                               rvec* x)
void Solve_Linear_Guassian(rvec a1, rvec a2, rvec a3, rvec a4,
                               rvec* x)
{
  //int n = 3;
  //int maxRow;
  //double A[n][n+1];
  //double maxEl, tmp, c;

  double D,a,b,c,d,l,m,n,k,p,q,r,s ;

  /* form 
  ax+by+cz+d=0
  lx+my+nz+k=0
  px+qy+rz+s=0
  */
 
  a = a1[0];
  b = a1[1];
  c = a1[2];
  l = a2[0];
  m = a2[1];
  n = a2[2];
  p = a3[0];
  q = a3[1];
  r = a3[2];
  d = -a4[0];
  k = -a4[1];
  s = -a4[2];


  D = (a*m*r+b*p*n+c*l*q)-(a*n*q+b*l*r+c*m*p);
  (*x)[0] = ((b*r*k+c*m*s+d*n*q)-(b*n*s+c*q*k+d*m*r))/D;
  (*x)[1]= ((a*n*s+c*p*k+d*l*r)-(a*r*k+c*l*s+d*n*p))/D;
  (*x)[2] = ((a*q*k+b*l*s+d*m*p)-(a*m*s+b*p*k+d*l*q))/D;


/*

  for (int i = 0; i < n; i++) {
    A[0][i] = a1[i];
    A[1][i] = a2[i];
    A[2][i] = a3[i];
  }

  for (int i = 0; i < n; i++) {
    A[i][3] = b[i];
  }

  for (int i = 0; i < n; i++) {
    // Search for maximum in this column
    maxEl = abs(A[i][i]);
    maxRow = i;
    for (int k = i+1; k < n; k++) {
      if (abs(A[k][i]) > maxEl) {
        maxEl = abs(A[k][i]);
        maxRow = k;
      }
    }

    // Swap maximum row with current row (column by column)
    for (int k = i; k < n+1; k++) {
      tmp = A[maxRow][k];
      A[maxRow][k] = A[i][k];
      A[i][k] = tmp;
    }

    // Make all rows below this one 0 in current column
    for (int k = i+1; k < n; k++) {
      c = -A[k][i]/A[i][i];
 std::cout<<"A "<<A[i][i]<<std::endl;
 
      for (int j = i; j < n+1; j++) {
        if (i == j) {
          A[k][j] = 0;
        } else {
          A[k][j] += c * A[i][j];
        }
      }
     }
  }
  // Solve equation Ax=b for an upper triangular matrix A
  for (int i = n-1; i >= 0; i--) {
   (*x)[i] = tmp = A[i][n]/A[i][i];
    for (int k = i-1; k>= 0; k--) {
      A[k][n] -= A[k][i] * tmp;
    }
  }
*/
}
/*------------------------------------------------------------------------------*/
void Calculate_dCos_Alpha( rvec v_b, rvec v_p, int dir, 
                           rvec v_m, double s_m, rvec v_n, double s_n,    
                           rvec v_l, rvec v_k, double s_k,
                           rvec* dcos_alpha_dm,
                           rvec* dcos_alpha_dn )
{
  int i, j;
  rvec v_c, Bm, Bn;
  rvec d_alpha_m, d_alpha_n ; 
  rvec v_cof1,v_cof2;
  rtensor d_b_m, d_b_n;
  rtensor d_c_m, d_c_n, d_p_m, d_p_n, A;
  double s2_m, s3_m, s2_n, s3_n,s_p, s2_p, s2_k, s3_l;
  double exp1m, exp2m, exp3m ;
  double exp1n, exp2n, exp3n ;
  double inprod;
  double s_l, s_b;
  double theta_0,theta_1,costheta1,costheta2;
  s_l = rvec_Norm(v_l);

  s2_m = SQR(s_m) ;
  s2_n = SQR(s_n) ;
  s3_m = s_m * s2_m ; 
  s3_n = s_m * s2_n ;
  s2_p = SQR(s_p) ;
  s3_l = SQR(s_l) * s_l ;
  s2_k = SQR(s_k);
  rvec v_m_unit, v_n_unit, v_sum; 

  theta_0 = DEG2RAD( 54.75 );
  theta_1 = DEG2RAD( 90.0 - 54.75 );

  costheta1 = cos(theta_0);
  costheta2 = cos(theta_1);
  

  // Step1: find the derevatives of x,y,z components of bisector vector
  // with respect to x,y,z components of the each H atom (m and n)
  // We store in a 3x3 materix  
  // v_b = - v_m/|v_m| - v_n/|v_n|
  // check later if this is diff than above:  v_b = - |v_n|* v_m - |v_m|*v_n   
  rvec_Scale( v_m_unit, -1/s_m, v_m ); // -v_m/|v_m|
  rvec_Scale( v_n_unit, -1/s_n, v_n ); // -v_n/|v_n|
  rvec_Sum( v_b, v_m_unit, v_n_unit ); // -v_m/|v_m| - v_n/|v_n| 
  s_b = rvec_Norm(v_b);
  // d_bis_i/d_m_j (i: x,y,z and j:x,y,z), m is one of the H
  // d_bis_i/d_n_j (i: x,y,z and j:x,y,z) n is the other H 
  for ( i = 0; i < 3; ++i ) {
   for ( j = 0; j < 3; ++j ) {
     if ( i == j )  {
       exp1m = -1/s_m + v_m[i]*v_m[i]/s3_m ;
       exp1n = -1/s_n + v_n[i]*v_n[i]/s3_n ;
     }
     else{
       exp1m = v_m[i]*v_m[j]/s3_m ;
       exp1n = v_n[i]*v_n[j]/s3_n ;
     }
     d_b_m[i][j] = exp1m;
     d_b_n[i][j] = exp1n;
   }
  }


  // Step 2. find the the derevatives of x,y,z components of perpendicular vector
  // with respect to x,y,z components of the each H atom (m and n)
  // We store in a 3x3 materix
  // v_p = v_m x v_n
  if (dir > 0 ) {
   rvec_Cross( v_p, v_n, v_m ) ;
   d_p_m[0][0] = d_p_m[1][1] = d_p_m[2][2] = 0.0 ;
   d_p_m[0][1] = -v_n[2];
   d_p_m[0][2] = v_n[1];
   d_p_m[1][2] = -v_n[0];
   d_p_m[1][0] = -d_p_m[0][1] ;
   d_p_m[2][0] = -d_p_m[0][2] ;
   d_p_m[2][1] = -d_p_m[1][2] ;
   d_p_n[0][0] = d_p_n[1][1] = d_p_n[2][2] = 0.0 ;
   d_p_n[0][1] = v_m[2];
   d_p_n[0][2] = -v_m[1];
   d_p_n[1][2] = v_m[0];
   d_p_n[1][0] = -d_p_n[0][1] ;
   d_p_n[2][0] = -d_p_n[0][2] ;
   d_p_n[2][1] = -d_p_n[1][2] ;
  }else {
   rvec_Cross( v_p, v_m, v_n ) ;
   d_p_m[0][0] = d_p_m[1][1] = d_p_m[2][2] = 0.0 ;
   d_p_m[0][1] = v_n[2];
   d_p_m[0][2] = -v_n[1];
   d_p_m[1][2] = v_n[0];
   d_p_m[1][0] = -d_p_m[0][1] ;
   d_p_m[2][0] = -d_p_m[0][2] ;
   d_p_m[2][1] = -d_p_m[1][2] ;
   d_p_n[0][0] = d_p_n[1][1] = d_p_n[2][2] = 0.0 ;
   d_p_n[0][1] = -v_m[2];
   d_p_n[0][2] = v_m[1];
   d_p_n[1][2] = -v_m[0];
   d_p_n[1][0] = -d_p_n[0][1] ;
   d_p_n[2][0] = -d_p_n[0][2] ;
   d_p_n[2][1] = -d_p_n[1][2] ;
  }
  s_p = rvec_Norm( v_p ) ;
  s2_p = SQR(s_p) ;


  // Step 3. find the the derevatives of x,y,z components of perpendicular vector
  // to v_p and v_b with respect to x,y,z components of m and n
  // We store in a 3x3 materix
  // v_c = v_b x v_p
  rvec_Cross( v_c, v_b, v_p ) ;
  for ( i = 0; i < 3; ++i ) {
    d_c_m[0][i] = -( d_p_m[1][i]*v_b[2] + v_p[1]*d_b_m[2][i] ) + ( d_p_m[2][i]*v_b[1] + v_p[2]*d_b_m[1][i] );
    d_c_m[1][i] = -( d_p_m[2][i]*v_b[0] + v_p[2]*d_b_m[0][i] ) + ( d_p_m[0][i]*v_b[2] + v_p[0]*d_b_m[2][i] );
    d_c_m[2][i] = -( d_p_m[0][i]*v_b[1] + v_p[0]*d_b_m[1][i] ) + ( d_p_m[1][i]*v_b[0] + v_p[1]*d_b_m[0][i] );
    d_c_n[0][i] = -( d_p_n[1][i]*v_b[2] + v_p[1]*d_b_n[2][i] ) + ( d_p_n[2][i]*v_b[1] + v_p[2]*d_b_n[1][i] );
    d_c_n[1][i] = -( d_p_n[2][i]*v_b[0] + v_p[2]*d_b_n[0][i] ) + ( d_p_n[0][i]*v_b[2] + v_p[0]*d_b_n[2][i] );
    d_c_n[2][i] = -( d_p_n[0][i]*v_b[1] + v_p[0]*d_b_n[1][i] ) + ( d_p_n[1][i]*v_b[0] + v_p[1]*d_b_n[0][i] );
  }
  

  // Step 4. Solve AX=B for x,y,z der of lone pair vector
  for ( i = 0; i < 3; ++i ) {
    for ( j = 0; j < 3; ++j ) {
      Bm[j] = Bn[j] = 0.0 ; 
    }
    for ( j = 0; j < 3; ++j ) {
      Bm[0] += -d_b_m[j][i] * v_l[j] + d_b_m[j][i] * v_b[j] * s_l/s_b * costheta1;
      Bm[1] += -d_p_m[j][i] * v_l[j] + d_p_m[j][i] * v_p[j] * s_l/s_p * costheta2;
      Bm[2] += -d_c_m[j][i] * v_l[j] ;
      Bn[0] += -d_b_n[j][i] * v_l[j] + d_b_n[j][i] * v_b[j] * s_l/s_b * costheta1;
      Bn[1] += -d_p_n[j][i] * v_l[j] + d_p_n[j][i] * v_p[j] * s_l/s_p * costheta2;
      Bn[2] += -d_c_n[j][i] * v_l[j] ;
      v_cof1[j] = v_b[j] - v_l[j] * s_b/s_l * costheta1;
      v_cof2[j] = v_p[j] - v_l[j] * s_p/s_l * costheta2;
    }

    //std::cout<<"Bm Bn "<<Bm[2]<<" "<<Bn[2]<<std::endl;
    Solve_Linear_Guassian(v_cof1, v_cof2, v_c, Bm, &d_alpha_m);
    Solve_Linear_Guassian(v_cof1, v_cof2, v_c, Bn, &d_alpha_n);
    //std::cout<<"d_alpha_m "<<d_alpha_m[0]<<" "<<d_alpha_m[1]<<" "<<d_alpha_m[2]<<" "<<std::endl;
    //std::cout<<"d_alpha_n "<<d_alpha_n[0]<<" "<<d_alpha_n[1]<<" "<<d_alpha_n[2]<<" "<<std::endl;
    exp1m = exp2m = 0.0;
    exp1n = exp2n = 0.0;
    inprod = 0.0;
    // v_k should not be normalized
    for ( j =0; j < 3; ++j ) {
      exp1m += d_alpha_m[j]*v_k[j];
      exp2m += d_alpha_m[j]*v_l[j];
      exp1n += d_alpha_n[j]*v_k[j];
      exp2n += d_alpha_n[j]*v_l[j];
      inprod += v_l[j]*v_k[j];
    }
    exp3m = exp1m/(s_k*s_l) - exp2m *inprod/(s_k * s3_l) ;
    exp3n = exp1n/(s_k*s_l) - exp2n *inprod/(s_k * s3_l) ;
    //std::cout<<"exp1m v_k "<<exp1m<<" "<<v_k[1]<<std::endl;
    //std::cout<<"exp3m "<<exp3m<<std::endl;
    //std::cout<<"exp3n "<<exp3n<<std::endl;
    //std::cout<<"----------------- "<<exp3n<<std::endl;
    (*dcos_alpha_dm)[i] = exp3m ;
    (*dcos_alpha_dn)[i] = exp3n ;
  }

}

