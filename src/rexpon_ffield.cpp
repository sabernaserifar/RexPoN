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
#include "error.h"
#include "rexpon_ffield.h"
#include "rexpon_tool_box.h"

char Read_Force_Field( FILE *fp, rexpon_interaction *rexpon,
                       control_params *control )
{
  char    *s;
  char   **tmp;
  char ****tor_flag;
  int      c, i, j, k, l, m, n, o, p, cnt;
  int lgflag = control->lgflag;
  int errorflag = 1;
  double     val;
  MPI_Comm comm;
  int me;

  comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &me);

  s = (char*) malloc(sizeof(char)*MAX_LINE);
  tmp = (char**) malloc(sizeof(char*)*MAX_TOKENS);
  for (i=0; i < MAX_TOKENS; i++)
    tmp[i] = (char*) malloc(sizeof(char)*MAX_TOKEN_LEN);

  /* reading first header comment */
  fgets( s, MAX_LINE, fp );

  /* line 2 is number of global parameters */
  fgets( s, MAX_LINE, fp );
  c = Tokenize( s, &tmp );

  /* reading the number of global parameters */
  n = atoi(tmp[0]);
  if (n < 1) {
    if (me == 0)
      fprintf( stderr, "WARNING: number of globals in ffield file is 0!\n" );
    fclose(fp);
    free(s);
    free(tmp);
    return 1;
  }

  rexpon->gp.n_global = n;
  rexpon->gp.l = (double*) malloc(sizeof(double)*n);

  /* see rexpon_types.h for mapping between l[i] and the lambdas used in ff */
  for (i=0; i < n; i++) {
    fgets(s,MAX_LINE,fp);
    c = Tokenize(s,&tmp);

    val = (double) atof(tmp[0]);
    rexpon->gp.l[i] = val;
  }

  control->bo_cut    = 0.01 * rexpon->gp.l[29];
  control->nonb_low  = rexpon->gp.l[11];
  control->nonb_cut  = rexpon->gp.l[12];

  /* next line is number of atom types and some comments */
  fgets( s, MAX_LINE, fp );
  c = Tokenize( s, &tmp );
  rexpon->num_atom_types = atoi(tmp[0]);

  /* 3 lines of comments */
  fgets(s,MAX_LINE,fp);
  fgets(s,MAX_LINE,fp);
  fgets(s,MAX_LINE,fp);

  /* Allocating structures in rexpon_interaction */
  rexpon->sbp = (single_body_parameters*)
    scalloc( rexpon->num_atom_types, sizeof(single_body_parameters), "sbp",
             comm );
  rexpon->tbp = (two_body_parameters**)
    scalloc( rexpon->num_atom_types, sizeof(two_body_parameters*), "tbp", comm );
  rexpon->thbp= (three_body_header***)
    scalloc( rexpon->num_atom_types, sizeof(three_body_header**), "thbp", comm );
  rexpon->hbp = (hbond_parameters***)
    scalloc( rexpon->num_atom_types, sizeof(hbond_parameters**), "hbp", comm );
  rexpon->fbp = (four_body_header****)
    scalloc( rexpon->num_atom_types, sizeof(four_body_header***), "fbp", comm );
  tor_flag  = (char****)
    scalloc( rexpon->num_atom_types, sizeof(char***), "tor_flag", comm );

  for( i = 0; i < rexpon->num_atom_types; i++ ) {
    rexpon->tbp[i] = (two_body_parameters*)
      scalloc( rexpon->num_atom_types, sizeof(two_body_parameters), "tbp[i]",
               comm );
    rexpon->thbp[i]= (three_body_header**)
      scalloc( rexpon->num_atom_types, sizeof(three_body_header*), "thbp[i]",
               comm );
    rexpon->hbp[i] = (hbond_parameters**)
      scalloc( rexpon->num_atom_types, sizeof(hbond_parameters*), "hbp[i]",
               comm );
    rexpon->fbp[i] = (four_body_header***)
      scalloc( rexpon->num_atom_types, sizeof(four_body_header**), "fbp[i]",
               comm );
    tor_flag[i]  = (char***)
      scalloc( rexpon->num_atom_types, sizeof(char**), "tor_flag[i]", comm );

    for( j = 0; j < rexpon->num_atom_types; j++ ) {
      rexpon->thbp[i][j]= (three_body_header*)
        scalloc( rexpon->num_atom_types, sizeof(three_body_header), "thbp[i,j]",
                 comm );
      rexpon->hbp[i][j] = (hbond_parameters*)
        scalloc( rexpon->num_atom_types, sizeof(hbond_parameters), "hbp[i,j]",
                 comm );
      rexpon->fbp[i][j] = (four_body_header**)
        scalloc( rexpon->num_atom_types, sizeof(four_body_header*), "fbp[i,j]",
                 comm );
      tor_flag[i][j]  = (char**)
        scalloc( rexpon->num_atom_types, sizeof(char*), "tor_flag[i,j]", comm );

      for (k=0; k < rexpon->num_atom_types; k++) {
        rexpon->fbp[i][j][k] = (four_body_header*)
          scalloc( rexpon->num_atom_types, sizeof(four_body_header), "fbp[i,j,k]",
                   comm );
        tor_flag[i][j][k]  = (char*)
          scalloc( rexpon->num_atom_types, sizeof(char), "tor_flag[i,j,k]",
                   comm );
      }
    }
  }

  rexpon->gp.vdw_type = 0;


  for( i = 0; i < rexpon->num_atom_types; i++ ) {
    /* line one */
    fgets( s, MAX_LINE, fp );
    c = Tokenize( s, &tmp );

    for( j = 0; j < (int)(strlen(tmp[0])); ++j )
      rexpon->sbp[i].name[j] = toupper( tmp[0][j] );

    val = atof(tmp[1]); rexpon->sbp[i].r_s        = val;
    val = atof(tmp[2]); rexpon->sbp[i].valency    = val;
    val = atof(tmp[3]); rexpon->sbp[i].mass       = val;
    val = atof(tmp[4]); rexpon->sbp[i].r_vdw      = val;
    val = atof(tmp[5]); rexpon->sbp[i].epsilon    = val;
    val = atof(tmp[6]); rexpon->sbp[i].gamma      = val;
    val = atof(tmp[7]); rexpon->sbp[i].r_pi       = val;
    val = atof(tmp[8]); rexpon->sbp[i].valency_e  = val;
    rexpon->sbp[i].nlp_opt = 0.5 * (rexpon->sbp[i].valency_e-rexpon->sbp[i].valency);

    /* line two */
    fgets( s, MAX_LINE, fp );
    c = Tokenize( s, &tmp );

    val = atof(tmp[0]); rexpon->sbp[i].alpha      = val;
    val = atof(tmp[1]); rexpon->sbp[i].gamma_w    = val;
    val = atof(tmp[2]); rexpon->sbp[i].valency_boc= val;
    val = atof(tmp[3]); rexpon->sbp[i].p_ovun5    = val;
    val = atof(tmp[4]);
    val = atof(tmp[5]); rexpon->sbp[i].chi        = val;
    val = atof(tmp[6]); rexpon->sbp[i].eta        = 2.0 * val;
    val = atof(tmp[7]); rexpon->sbp[i].p_hbond = (int) val;

    /* line 3 */
    fgets( s, MAX_LINE, fp );
    c = Tokenize( s, &tmp );

    val = atof(tmp[0]); rexpon->sbp[i].r_pi_pi    = val;
    val = atof(tmp[1]); rexpon->sbp[i].p_lp2      = val;
    val = atof(tmp[2]);
    val = atof(tmp[3]); rexpon->sbp[i].b_o_131    = val;
    val = atof(tmp[4]); rexpon->sbp[i].b_o_132    = val;
    val = atof(tmp[5]); rexpon->sbp[i].b_o_133    = val;
    val = atof(tmp[6]);
    val = atof(tmp[7]);

    /* line 4  */
    fgets( s, MAX_LINE, fp );
    c = Tokenize( s, &tmp );

    /* Sanity check */
    if (c < 3) {
      if (me == 0)
        fprintf(stderr, "Inconsistent ffield file (rexpon_ffield.cpp) \n");
      MPI_Abort( comm, FILE_NOT_FOUND );
    }

    val = atof(tmp[0]); rexpon->sbp[i].p_ovun2    = val;
    val = atof(tmp[1]); rexpon->sbp[i].p_val3     = val;
    val = atof(tmp[2]);
    val = atof(tmp[3]); rexpon->sbp[i].valency_val= val;
    val = atof(tmp[4]); rexpon->sbp[i].p_val5     = val;
    val = atof(tmp[5]); rexpon->sbp[i].rcore2     = val;
    val = atof(tmp[6]); rexpon->sbp[i].ecore2     = val;
    val = atof(tmp[7]); rexpon->sbp[i].acore2     = val;

    /* line 5, only if lgvdw is yes */
    if (lgflag) {
      fgets( s, MAX_LINE, fp );
      c = Tokenize( s, &tmp );

      /* Sanity check */
      if (c > 3) {
        if (me == 0)
          fprintf(stderr, "Inconsistent ffield file (rexpon_ffield.cpp) \n");
        MPI_Abort( comm, FILE_NOT_FOUND );
      }

      val = atof(tmp[0]); rexpon->sbp[i].lgcij           = val;
      val = atof(tmp[1]); rexpon->sbp[i].lgre           = val;
    }

    if( rexpon->sbp[i].rcore2>0.01 && rexpon->sbp[i].acore2>0.01 ){ // Inner-wall
      if( rexpon->sbp[i].gamma_w>0.5 ){ // Shielding vdWaals
        if( rexpon->gp.vdw_type != 0 && rexpon->gp.vdw_type != 3 ) {
          if (errorflag && (me == 0))
            fprintf( stderr, "Warning: inconsistent vdWaals-parameters\n"        \
                   "Force field parameters for element %s\n"                \
                   "indicate inner wall+shielding, but earlier\n"        \
                   "atoms indicate different vdWaals-method.\n"                \
                   "This may cause division-by-zero errors.\n"                \
                   "Keeping vdWaals-setting for earlier atoms.\n",
                   rexpon->sbp[i].name );
          errorflag = 0;
        } else {
          rexpon->gp.vdw_type = 3;
        }
      }
      else {  // No shielding vdWaals parameters present
        if( rexpon->gp.vdw_type != 0 && rexpon->gp.vdw_type != 2 ) {
          if (me == 0)
            fprintf( stderr, "Warning: inconsistent vdWaals-parameters\n" \
                   "Force field parameters for element %s\n"                \
                   "indicate inner wall without shielding, but earlier\n" \
                   "atoms indicate different vdWaals-method.\n"                \
                   "This may cause division-by-zero errors.\n"                \
                   "Keeping vdWaals-setting for earlier atoms.\n",
                   rexpon->sbp[i].name );
        } else {
          rexpon->gp.vdw_type = 2;
        }
      }
    }
    else{ // No Inner wall parameters present
      if( rexpon->sbp[i].gamma_w>0.5 ){ // Shielding vdWaals
        if( rexpon->gp.vdw_type != 0 && rexpon->gp.vdw_type != 1 ) {
          if (me == 0)
            fprintf( stderr, "Warning: inconsistent vdWaals-parameters\n"        \
                   "Force field parameters for element %s\n"                \
                   "indicate  shielding without inner wall, but earlier\n" \
                   "atoms indicate different vdWaals-method.\n"                \
                   "This may cause division-by-zero errors.\n"                \
                   "Keeping vdWaals-setting for earlier atoms.\n",
                   rexpon->sbp[i].name );
        } else {
          rexpon->gp.vdw_type = 1;
        }
      } else {
        if (me == 0)
          fprintf( stderr, "Error: inconsistent vdWaals-parameters\n"  \
                 "No shielding or inner-wall set for element %s\n",
                 rexpon->sbp[i].name );
        MPI_Abort( comm, INVALID_INPUT );
      }
    }
  }

  /* Equate vval3 to valf for first-row elements (25/10/2004) */
  for( i = 0; i < rexpon->num_atom_types; i++ )
    if( rexpon->sbp[i].mass < 21 &&
        rexpon->sbp[i].valency_val != rexpon->sbp[i].valency_boc ) {
      //if (me == 0)
        //fprintf(stderr,"Warning: changed valency_val to valency_boc for %s\n",
          //     rexpon->sbp[i].name );
      rexpon->sbp[i].valency_val = rexpon->sbp[i].valency_boc;
    }

  /* next line is number of two body combination and some comments */
  fgets(s,MAX_LINE,fp);
  c=Tokenize(s,&tmp);
  l = atoi(tmp[0]);

  /* a line of comments */
  fgets(s,MAX_LINE,fp);

  for (i=0; i < l; i++) {
    /* line 1 */
    fgets(s,MAX_LINE,fp);
    c=Tokenize(s,&tmp);

    j = atoi(tmp[0]) - 1;
    k = atoi(tmp[1]) - 1;

    if (j < rexpon->num_atom_types && k < rexpon->num_atom_types) {

      val = atof(tmp[2]); rexpon->tbp[j][k].De_s      = val;
      rexpon->tbp[k][j].De_s      = val;
      val = atof(tmp[3]); rexpon->tbp[j][k].De_p      = val;
      rexpon->tbp[k][j].De_p      = val;
      val = atof(tmp[4]); rexpon->tbp[j][k].De_pp     = val;
      rexpon->tbp[k][j].De_pp     = val;
      val = atof(tmp[5]); rexpon->tbp[j][k].p_be1     = val;
      rexpon->tbp[k][j].p_be1     = val;
      val = atof(tmp[6]); rexpon->tbp[j][k].p_bo5     = val;
      rexpon->tbp[k][j].p_bo5     = val;
      val = atof(tmp[7]); rexpon->tbp[j][k].v13cor    = val;
      rexpon->tbp[k][j].v13cor    = val;

      val = atof(tmp[8]); rexpon->tbp[j][k].p_bo6     = val;
      rexpon->tbp[k][j].p_bo6     = val;
      val = atof(tmp[9]); rexpon->tbp[j][k].p_ovun1 = val;
      rexpon->tbp[k][j].p_ovun1 = val;

      /* line 2 */
      fgets(s,MAX_LINE,fp);
      c=Tokenize(s,&tmp);

      val = atof(tmp[0]); rexpon->tbp[j][k].p_be2     = val;
      rexpon->tbp[k][j].p_be2     = val;
      val = atof(tmp[1]); rexpon->tbp[j][k].p_bo3     = val;
      rexpon->tbp[k][j].p_bo3     = val;
      val = atof(tmp[2]); rexpon->tbp[j][k].p_bo4     = val;
      rexpon->tbp[k][j].p_bo4     = val;
      val = atof(tmp[3]); rexpon->tbp[j][k].bom     = val;
      rexpon->tbp[k][j].bom     = val;
      val = atof(tmp[4]); rexpon->tbp[j][k].p_bo1     = val;
      rexpon->tbp[k][j].p_bo1     = val;
      val = atof(tmp[5]); rexpon->tbp[j][k].p_bo2     = val;
      rexpon->tbp[k][j].p_bo2     = val;
      val = atof(tmp[6]); rexpon->tbp[j][k].ovc       = val;
      rexpon->tbp[k][j].ovc       = val;
      val = atof(tmp[7]); rexpon->tbp[j][k].reqm       = val;
      rexpon->tbp[k][j].reqm       = val;
    }
  }

  for (i=0; i < rexpon->num_atom_types; i++)
    for (j=i; j < rexpon->num_atom_types; j++) {
      rexpon->tbp[i][j].r_s = 0.5 *
        (rexpon->sbp[i].r_s + rexpon->sbp[j].r_s);
      rexpon->tbp[j][i].r_s = 0.5 *
        (rexpon->sbp[j].r_s + rexpon->sbp[i].r_s);

      rexpon->tbp[i][j].r_p = 0.5 *
        (rexpon->sbp[i].r_pi + rexpon->sbp[j].r_pi);
      rexpon->tbp[j][i].r_p = 0.5 *
        (rexpon->sbp[j].r_pi + rexpon->sbp[i].r_pi);

      rexpon->tbp[i][j].r_pp = 0.5 *
        (rexpon->sbp[i].r_pi_pi + rexpon->sbp[j].r_pi_pi);
      rexpon->tbp[j][i].r_pp = 0.5 *
        (rexpon->sbp[j].r_pi_pi + rexpon->sbp[i].r_pi_pi);


      rexpon->tbp[i][j].p_boc3 =
        sqrt(rexpon->sbp[i].b_o_132 *
             rexpon->sbp[j].b_o_132);
      rexpon->tbp[j][i].p_boc3 =
        sqrt(rexpon->sbp[j].b_o_132 *
             rexpon->sbp[i].b_o_132);

      rexpon->tbp[i][j].p_boc4 =
        sqrt(rexpon->sbp[i].b_o_131 *
             rexpon->sbp[j].b_o_131);
      rexpon->tbp[j][i].p_boc4 =
        sqrt(rexpon->sbp[j].b_o_131 *
             rexpon->sbp[i].b_o_131);

      rexpon->tbp[i][j].p_boc5 =
        sqrt(rexpon->sbp[i].b_o_133 *
             rexpon->sbp[j].b_o_133);
      rexpon->tbp[j][i].p_boc5 =
        sqrt(rexpon->sbp[j].b_o_133 *
             rexpon->sbp[i].b_o_133);


      rexpon->tbp[i][j].D =
        sqrt(rexpon->sbp[i].epsilon *
             rexpon->sbp[j].epsilon);

      rexpon->tbp[j][i].D =
        sqrt(rexpon->sbp[j].epsilon *
             rexpon->sbp[i].epsilon);

      rexpon->tbp[i][j].alpha =
        sqrt(rexpon->sbp[i].alpha *
             rexpon->sbp[j].alpha);

      rexpon->tbp[j][i].alpha =
        sqrt(rexpon->sbp[j].alpha *
             rexpon->sbp[i].alpha);

      /*rexpon->tbp[i][j].r_vdW =
        2.0 * sqrt(rexpon->sbp[i].r_vdw * rexpon->sbp[j].r_vdw);

      rexpon->tbp[j][i].r_vdW =
        2.0 * sqrt(rexpon->sbp[j].r_vdw * rexpon->sbp[i].r_vdw);*/


      /* RexPoN */

      rexpon->tbp[i][j].r_vdW =
        1.0 * sqrt(rexpon->sbp[i].r_vdw * rexpon->sbp[j].r_vdw);

      rexpon->tbp[j][i].r_vdW =
        1.0 * sqrt(rexpon->sbp[j].r_vdw * rexpon->sbp[i].r_vdw);
      //

      rexpon->tbp[i][j].gamma_w =
        sqrt(rexpon->sbp[i].gamma_w *
             rexpon->sbp[j].gamma_w);

      rexpon->tbp[j][i].gamma_w =
        sqrt(rexpon->sbp[j].gamma_w *
             rexpon->sbp[i].gamma_w);

      rexpon->tbp[i][j].gamma =
        pow(rexpon->sbp[i].gamma *
            rexpon->sbp[j].gamma,-1.5);

      rexpon->tbp[j][i].gamma =
        pow(rexpon->sbp[j].gamma *
            rexpon->sbp[i].gamma,-1.5);


      /* RexPoN */
      rexpon->tbp[i][j].gamma =
        sqrt(rexpon->sbp[i].gamma *
             rexpon->sbp[j].gamma);

      rexpon->tbp[j][i].gamma =
        sqrt(rexpon->sbp[j].gamma *
             rexpon->sbp[i].gamma);
      //

      // additions for additional vdWaals interaction types - inner core

      rexpon->tbp[i][j].rcore = rexpon->tbp[j][i].rcore =
        sqrt( rexpon->sbp[i].rcore2 * rexpon->sbp[j].rcore2 );

      rexpon->tbp[i][j].ecore = rexpon->tbp[j][i].ecore =
        sqrt( rexpon->sbp[i].ecore2 * rexpon->sbp[j].ecore2 );

      rexpon->tbp[i][j].acore = rexpon->tbp[j][i].acore =
        sqrt( rexpon->sbp[i].acore2 * rexpon->sbp[j].acore2 );

      // additions for additional vdWalls interaction types lg correction

      rexpon->tbp[i][j].lgcij = rexpon->tbp[j][i].lgcij =
        sqrt( rexpon->sbp[i].lgcij * rexpon->sbp[j].lgcij );

      rexpon->tbp[i][j].lgre = rexpon->tbp[j][i].lgre = 2.0 * rexpon->gp.l[35] *
        sqrt( rexpon->sbp[i].lgre*rexpon->sbp[j].lgre );

    }

  fgets(s,MAX_LINE,fp);
  c=Tokenize(s,&tmp);
  l = atoi(tmp[0]);

  for (i=0; i < l; i++) {
    fgets(s,MAX_LINE,fp);
    c=Tokenize(s,&tmp);

    j = atoi(tmp[0]) - 1;
    k = atoi(tmp[1]) - 1;

    if (j < rexpon->num_atom_types && k < rexpon->num_atom_types)        {
      val = atof(tmp[2]);
      //if (val > 0.0) {
        rexpon->tbp[j][k].D = val;
        rexpon->tbp[k][j].D = val;
      //}

      val = atof(tmp[3]);
      //if (val > 0.0) {
        rexpon->tbp[j][k].r_vdW = 1 * val;
        rexpon->tbp[k][j].r_vdW = 1 * val;
      //}
        rexpon->tbp[j][k].beta = val;
        rexpon->tbp[k][j].beta = val;

      val = atof(tmp[4]);
      //if (val > 0.0) {
        rexpon->tbp[j][k].alpha = val;
        rexpon->tbp[k][j].alpha = val;
      //}

      val = atof(tmp[5]);
      //if (val > 0.0) {
        rexpon->tbp[j][k].r_s = val;
        rexpon->tbp[k][j].r_s = val;
      //}

      val = atof(tmp[6]);
      if (val > 0.0) {
        rexpon->tbp[j][k].r_p = val;
        rexpon->tbp[k][j].r_p = val;
      }

      val = atof(tmp[7]);
      if (val > 0.0) {
        rexpon->tbp[j][k].r_pp = val;
        rexpon->tbp[k][j].r_pp = val;
      }

      val = atof(tmp[8]);
        rexpon->tbp[j][k].repa0 = val;
        rexpon->tbp[k][j].repa0 = val;

      val = atof(tmp[9]);
        rexpon->tbp[j][k].repr0 = val;
        rexpon->tbp[k][j].repr0 = val;

      val = atof(tmp[10]);
        rexpon->tbp[j][k].repn = val;
        rexpon->tbp[k][j].repn = val;

      val = atof(tmp[11]);
        rexpon->tbp[j][k].repscal = val;
        rexpon->tbp[k][j].repscal = val;

      val = atof(tmp[12]);
        rexpon->tbp[j][k].reps = val;
        rexpon->tbp[k][j].reps = val;

      val = atof(tmp[13]);
        rexpon->tbp[j][k].lgre = 2*val;
        rexpon->tbp[k][j].lgre = 2*val;

      val = atof(tmp[14]);
      if (val >= 0.0) {
        rexpon->tbp[j][k].lgcij = val;
        rexpon->tbp[k][j].lgcij = val;
      }
    }
  }

  for( i = 0; i < rexpon->num_atom_types; ++i )
    for( j = 0; j < rexpon->num_atom_types; ++j )
      for( k = 0; k < rexpon->num_atom_types; ++k )
        rexpon->thbp[i][j][k].cnt = 0;

  fgets( s, MAX_LINE, fp );
  c = Tokenize( s, &tmp );
  l = atoi( tmp[0] );

  for( i = 0; i < l; i++ ) {
    fgets(s,MAX_LINE,fp);
    c=Tokenize(s,&tmp);

    j = atoi(tmp[0]) - 1;
    k = atoi(tmp[1]) - 1;
    m = atoi(tmp[2]) - 1;

    if (j < rexpon->num_atom_types && k < rexpon->num_atom_types &&
        m < rexpon->num_atom_types) {
      cnt = rexpon->thbp[j][k][m].cnt;
      rexpon->thbp[j][k][m].cnt++;
      rexpon->thbp[m][k][j].cnt++;

      val = atof(tmp[3]);
      rexpon->thbp[j][k][m].prm[cnt].theta_00 = val;
      rexpon->thbp[m][k][j].prm[cnt].theta_00 = val;

      val = atof(tmp[4]);
      rexpon->thbp[j][k][m].prm[cnt].p_val1 = val;
      rexpon->thbp[m][k][j].prm[cnt].p_val1 = val;

      val = atof(tmp[5]);
      rexpon->thbp[j][k][m].prm[cnt].p_val2 = val;
      rexpon->thbp[m][k][j].prm[cnt].p_val2 = val;

      val = atof(tmp[6]);
      rexpon->thbp[j][k][m].prm[cnt].p_coa1 = val;
      rexpon->thbp[m][k][j].prm[cnt].p_coa1 = val;

      val = atof(tmp[7]);
      rexpon->thbp[j][k][m].prm[cnt].p_val7 = val;
      rexpon->thbp[m][k][j].prm[cnt].p_val7 = val;

      val = atof(tmp[8]);
      rexpon->thbp[j][k][m].prm[cnt].p_pen1 = val;
      rexpon->thbp[m][k][j].prm[cnt].p_pen1 = val;

      val = atof(tmp[9]);
      rexpon->thbp[j][k][m].prm[cnt].p_val4 = val;
      rexpon->thbp[m][k][j].prm[cnt].p_val4 = val;
    }
  }

  /* clear all entries first */
  for( i = 0; i < rexpon->num_atom_types; ++i )
    for( j = 0; j < rexpon->num_atom_types; ++j )
      for( k = 0; k < rexpon->num_atom_types; ++k )
        for( m = 0; m < rexpon->num_atom_types; ++m ) {
          rexpon->fbp[i][j][k][m].cnt = 0;
          tor_flag[i][j][k][m] = 0;
        }

  /* next line is number of 4-body params and some comments */
  fgets( s, MAX_LINE, fp );
  c = Tokenize( s, &tmp );
  l = atoi( tmp[0] );

  for( i = 0; i < l; i++ ) {
    fgets( s, MAX_LINE, fp );
    c = Tokenize( s, &tmp );

    j = atoi(tmp[0]) - 1;
    k = atoi(tmp[1]) - 1;
    m = atoi(tmp[2]) - 1;
    n = atoi(tmp[3]) - 1;

    if (j >= 0 && n >= 0) { // this means the entry is not in compact form
      if (j < rexpon->num_atom_types && k < rexpon->num_atom_types &&
          m < rexpon->num_atom_types && n < rexpon->num_atom_types) {
        tor_flag[j][k][m][n] = 1;
        tor_flag[n][m][k][j] = 1;

        rexpon->fbp[j][k][m][n].cnt = 1;
        rexpon->fbp[n][m][k][j].cnt = 1;

        val = atof(tmp[4]);
        rexpon->fbp[j][k][m][n].prm[0].V1 = val;
        rexpon->fbp[n][m][k][j].prm[0].V1 = val;

        val = atof(tmp[5]);
        rexpon->fbp[j][k][m][n].prm[0].V2 = val;
        rexpon->fbp[n][m][k][j].prm[0].V2 = val;

        val = atof(tmp[6]);
        rexpon->fbp[j][k][m][n].prm[0].V3 = val;
        rexpon->fbp[n][m][k][j].prm[0].V3 = val;

        val = atof(tmp[7]);
        rexpon->fbp[j][k][m][n].prm[0].p_tor1 = val;
        rexpon->fbp[n][m][k][j].prm[0].p_tor1 = val;

        val = atof(tmp[8]);
        rexpon->fbp[j][k][m][n].prm[0].p_cot1 = val;
        rexpon->fbp[n][m][k][j].prm[0].p_cot1 = val;
      }
    }
    else { /* This means the entry is of the form 0-X-Y-0 */
      if( k < rexpon->num_atom_types && m < rexpon->num_atom_types )
        for( p = 0; p < rexpon->num_atom_types; p++ )
          for( o = 0; o < rexpon->num_atom_types; o++ ) {
            rexpon->fbp[p][k][m][o].cnt = 1;
            rexpon->fbp[o][m][k][p].cnt = 1;

            if (tor_flag[p][k][m][o] == 0) {
              rexpon->fbp[p][k][m][o].prm[0].V1 = atof(tmp[4]);
              rexpon->fbp[p][k][m][o].prm[0].V2 = atof(tmp[5]);
              rexpon->fbp[p][k][m][o].prm[0].V3 = atof(tmp[6]);
              rexpon->fbp[p][k][m][o].prm[0].p_tor1 = atof(tmp[7]);
              rexpon->fbp[p][k][m][o].prm[0].p_cot1 = atof(tmp[8]);
            }

            if (tor_flag[o][m][k][p] == 0) {
              rexpon->fbp[o][m][k][p].prm[0].V1 = atof(tmp[4]);
              rexpon->fbp[o][m][k][p].prm[0].V2 = atof(tmp[5]);
              rexpon->fbp[o][m][k][p].prm[0].V3 = atof(tmp[6]);
              rexpon->fbp[o][m][k][p].prm[0].p_tor1 = atof(tmp[7]);
              rexpon->fbp[o][m][k][p].prm[0].p_cot1 = atof(tmp[8]);
            }
          }
    }
  }



  /* next line is number of hydrogen bond params and some comments */
  fgets( s, MAX_LINE, fp );
  c = Tokenize( s, &tmp );
  l = atoi( tmp[0] );

  for( i = 0; i < rexpon->num_atom_types; ++i )
    for( j = 0; j < rexpon->num_atom_types; ++j )
      for( k = 0; k < rexpon->num_atom_types; ++k )
        rexpon->hbp[i][j][k].r0_hb = -1.0;

  for( i = 0; i < l; i++ ) {
    fgets( s, MAX_LINE, fp );
    c = Tokenize( s, &tmp );

    j = atoi(tmp[0]) - 1;
    k = atoi(tmp[1]) - 1;
    m = atoi(tmp[2]) - 1;


    if( j < rexpon->num_atom_types && m < rexpon->num_atom_types ) {
      val = atof(tmp[3]);
      rexpon->hbp[j][k][m].r0_hb = val;

      val = atof(tmp[4]);
      rexpon->hbp[j][k][m].p_hb1 = val;

      val = atof(tmp[5]);
      rexpon->hbp[j][k][m].p_hb2 = val;

      val = atof(tmp[6]);
      rexpon->hbp[j][k][m].p_hb3 = val;

      val = atof(tmp[7]);
      rexpon->hbp[j][k][m].p_hb4 = val;

      val = atof(tmp[8]);
      rexpon->hbp[j][k][m].p_hb5 = val;

      val = atof(tmp[9]);
      rexpon->hbp[j][k][m].p_hb6 = val;

      val = atof(tmp[10]);
      rexpon->hbp[j][k][m].p_hb7 = val;


      fgets(s,MAX_LINE,fp);
      c=Tokenize(s,&tmp);

      val = atof(tmp[0]);
      rexpon->hbp[j][k][m].p_hb8 = val;

      val = atof(tmp[1]);
      rexpon->hbp[j][k][m].p_hb9 = val;

      val = atof(tmp[2]);
      rexpon->hbp[j][k][m].p_hb10 = val;

      val = atof(tmp[3]);
      rexpon->hbp[j][k][m].p_hb11 = val;

      val = atof(tmp[4]);
      rexpon->hbp[j][k][m].p_hb12 = val;

      val = atof(tmp[5]);
      rexpon->hbp[j][k][m].p_hb13 = val;

      val = atof(tmp[6]);
      rexpon->hbp[j][k][m].p_hb14 = val;

      val = atof(tmp[7]);
      rexpon->hbp[j][k][m].p_hb15 = val;

      fgets(s,MAX_LINE,fp);
      c=Tokenize(s,&tmp);

      val = atof(tmp[0]);
      rexpon->hbp[j][k][m].p_hb16 = val;

      val = atof(tmp[1]);
      rexpon->hbp[j][k][m].p_hb17 = val;

      val = atof(tmp[2]);
      rexpon->hbp[j][k][m].p_hb18 = val;

      val = atof(tmp[3]);
      rexpon->hbp[j][k][m].p_hb19 = val;

      val = atof(tmp[4]);
      rexpon->hbp[j][k][m].p_hb20 = val;

      val = atof(tmp[5]);
      rexpon->hbp[j][k][m].p_hb21 = val;

      val = atof(tmp[6]);
      rexpon->hbp[j][k][m].p_hb22 = val;

      val = atof(tmp[7]);
      rexpon->hbp[j][k][m].p_hb23 = val;

    }
  }

  /* deallocate helper storage */
  for( i = 0; i < MAX_TOKENS; i++ )
    free( tmp[i] );
  free( tmp );
  free( s );


  /* deallocate tor_flag */
  for( i = 0; i < rexpon->num_atom_types; i++ ) {
    for( j = 0; j < rexpon->num_atom_types; j++ ) {
      for( k = 0; k < rexpon->num_atom_types; k++ ) {
        free( tor_flag[i][j][k] );
      }
      free( tor_flag[i][j] );
    }
    free( tor_flag[i] );
  }
  free( tor_flag );

  // close file

  fclose(fp);

  return SUCCESS;
}
