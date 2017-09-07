/*
 *                minimum free energy
 *                RNA secondary structure prediction
 *
 *                c Ivo Hofacker, Chrisoph Flamm
 *                original implementation by
 *                Walter Fontana
 *                g-quadruplex support and threadsafety
 *                by Ronny Lorenz
 *
 *                Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/utils.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/data_structures.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/loop_energies.h"
#include "ViennaRNA/mfe.h"

#define MAXSECTORS        500     /* dimension for a backtrack array */

/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE int
fill_arrays(vrna_fold_compound_t *vc);


PRIVATE int
postprocess_circular(vrna_fold_compound_t *vc,
                     sect                 bt_stack[],
                     int                  *bt);


PRIVATE void
backtrack(vrna_fold_compound_t  *vc,
          vrna_bp_stack_t       *bp_stack,
          sect                  bt_stack[],
          int                   s);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC float
vrna_mfe(vrna_fold_compound_t *vc,
         char                 *structure)
{
  char            *ss;
  int             length, energy, s;
  float           mfe;
  sect            bt_stack[MAXSECTORS]; /* stack of partial structures for backtracking */
  vrna_bp_stack_t *bp;

  s   = 0;
  mfe = (float)(INF / 100.);

  if (vc) {
    length = (int)vc->length;

    if (!vrna_fold_compound_prepare(vc, VRNA_OPTION_MFE)) {
      vrna_message_warning("vrna_mfe@mfe.c: Failed to prepare vrna_fold_compound");
      return mfe;
    }

    /* call user-defined recursion status callback function */
    if (vc->stat_cb)
      vc->stat_cb(VRNA_STATUS_MFE_PRE, vc->auxdata);

    energy = fill_arrays(vc);

    if (vc->params->model_details.circ)
      energy = postprocess_circular(vc, bt_stack, &s);

    /* call user-defined recursion status callback function */
    if (vc->stat_cb)
      vc->stat_cb(VRNA_STATUS_MFE_POST, vc->auxdata);

    if (structure && vc->params->model_details.backtrack) {
      /* add a guess of how many G's may be involved in a G quadruplex */
      bp = (vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (4 * (1 + length / 2)));

      backtrack(vc, bp, bt_stack, s);

      ss = vrna_db_from_bp_stack(bp, length);
      strncpy(structure, ss, length + 1);
      free(ss);
      free(bp);
    }

    switch (vc->params->model_details.backtrack_type) {
      case 'C':
        mfe = (float)vc->matrices->c[vc->jindx[length] + 1] / 100.;
        break;

      case 'M':
        mfe = (float)vc->matrices->fML[vc->jindx[length] + 1] / 100.;
        break;

      default:
        mfe = (float)energy / 100.;
        break;
    }

    if (vc->type == VRNA_FC_TYPE_COMPARATIVE)
      mfe /= (float)vc->n_seq;
  }

  return mfe;
}


/**
*** fill "c", "fML" and "f5" arrays and return  optimal energy
**/
PRIVATE int
fill_arrays(vrna_fold_compound_t *vc)
{
  unsigned char *hard_constraints;
  int           i, j, ij, length, energy, new_c, stackEnergy, turn,
                noLP, uniq_ML, dangle_model, *indx, *f5,
                *c, *fML, *fM1, hc_decompose, *cc, *cc1, *Fmi, *DMLi,
                *DMLi1, *DMLi2, *pscore;
  vrna_param_t  *P;
  vrna_mx_mfe_t *matrices;
  vrna_hc_t     *hc;
  vrna_ud_t     *domains_up;

  length            = (int)vc->length;
  pscore            = (vc->type == VRNA_FC_TYPE_COMPARATIVE) ? vc->pscore : NULL;
  indx              = vc->jindx;
  P                 = vc->params;
  noLP              = P->model_details.noLP;
  uniq_ML           = P->model_details.uniq_ML;
  dangle_model      = P->model_details.dangles;
  turn              = P->model_details.min_loop_size;
  hc                = vc->hc;
  hard_constraints  = hc->matrix;
  matrices          = vc->matrices;
  f5                = matrices->f5;
  c                 = matrices->c;
  fML               = matrices->fML;
  fM1               = matrices->fM1;
  domains_up        = vc->domains_up;

  /* allocate memory for all helper arrays */
  cc    = (int *)vrna_alloc(sizeof(int) * (length + 2));  /* auxilary arrays for canonical structures     */
  cc1   = (int *)vrna_alloc(sizeof(int) * (length + 2));  /* auxilary arrays for canonical structures     */
  Fmi   = (int *)vrna_alloc(sizeof(int) * (length + 1));  /* holds row i of fML (avoids jumps in memory)  */
  DMLi  = (int *)vrna_alloc(sizeof(int) * (length + 1));  /* DMLi[j] holds  MIN(fML[i,k]+fML[k+1,j])      */
  DMLi1 = (int *)vrna_alloc(sizeof(int) * (length + 1));  /*                MIN(fML[i+1,k]+fML[k+1,j])    */
  DMLi2 = (int *)vrna_alloc(sizeof(int) * (length + 1));  /*                MIN(fML[i+2,k]+fML[k+1,j])    */

  if ((turn < 0) || (turn > length))
    turn = length; /* does this make any sense? */

  /* pre-processing ligand binding production rule(s) */
  if (domains_up && domains_up->prod_cb)
    domains_up->prod_cb(vc, domains_up->data);

  /* prefill helper arrays */
  for (j = 0; j <= length; j++)
    Fmi[j] = DMLi[j] = DMLi1[j] = DMLi2[j] = INF;

  /* prefill matrices with init contributions */
  for (j = 1; j <= length; j++)
    for (i = (j > turn ? (j - turn) : 1); i <= j; i++) {
      c[indx[j] + i] = fML[indx[j] + i] = INF;
      if (uniq_ML)
        fM1[indx[j] + i] = INF;
    }

  /* start recursion */

  if (length <= turn) {
    /* clean up memory */
    free(cc);
    free(cc1);
    free(Fmi);
    free(DMLi);
    free(DMLi1);
    free(DMLi2);
    /* return free energy of unfolded chain */
    return 0;
  }

  for (i = length - turn - 1; i >= 1; i--) {
    /* i,j in [1..length] */

    for (j = i + turn + 1; j <= length; j++) {
      ij            = indx[j] + i;
      hc_decompose  = hard_constraints[ij];
      energy        = INF;

      /* do we evaluate this pair? */
      if (hc_decompose) {
        new_c = INF;

        /* check for hairpin loop */
        energy  = vrna_E_hp_loop(vc, i, j);
        new_c   = MIN2(new_c, energy);

        /* check for multibranch loops */
        energy  = vrna_E_mb_loop_fast(vc, i, j, DMLi1, DMLi2);
        new_c   = MIN2(new_c, energy);

        if (dangle_model == 3) {
          /* coaxial stacking */
          energy  = vrna_E_mb_loop_stack(vc, i, j);
          new_c   = MIN2(new_c, energy);
        }

        /* check for interior loops */
        energy  = vrna_E_int_loop(vc, i, j);
        new_c   = MIN2(new_c, energy);

        /* remember stack energy for --noLP option */
        if (noLP) {
          stackEnergy = vrna_E_stack(vc, i, j);
          new_c       = MIN2(new_c, cc1[j - 1] + stackEnergy);
          cc[j]       = new_c;
          c[ij]       = cc1[j - 1] + stackEnergy;
          if (vc->type == VRNA_FC_TYPE_COMPARATIVE && (cc[j] != INF))
            cc[j] -= pscore[indx[j] + i];
        } else {
          c[ij] = new_c;
        }

        if (vc->type == VRNA_FC_TYPE_COMPARATIVE && (c[ij] != INF))
          c[ij] -= pscore[indx[j] + i];
      } /* end >> if (pair) << */
      else {
        c[ij] = INF;
      }

      /* done with c[i,j], now compute fML[i,j] and fM1[i,j] */

      fML[ij] = vrna_E_ml_stems_fast(vc, i, j, Fmi, DMLi);

      if (uniq_ML)   /* compute fM1 for unique decomposition */
        fM1[ij] = E_ml_rightmost_stem(i, j, vc);
    } /* end of j-loop */

    {
      int *FF; /* rotate the auxilliary arrays */
      FF    = DMLi2;
      DMLi2 = DMLi1;
      DMLi1 = DMLi;
      DMLi  = FF;
      FF    = cc1;
      cc1   = cc;
      cc    = FF;
      for (j = 1; j <= length; j++)
        cc[j] = Fmi[j] = DMLi[j] = INF;
    }
  } /* end of i-loop */

  /* calculate energies of 5' fragments */
  (void)vrna_E_ext_loop_5(vc);

  /* clean up memory */
  free(cc);
  free(cc1);
  free(Fmi);
  free(DMLi);
  free(DMLi1);
  free(DMLi2);

  return f5[length];
}


/* post-processing step for circular RNAs */
PRIVATE int
postprocess_circular(vrna_fold_compound_t *vc,
                     sect                 bt_stack[],
                     int                  *bt)
{
  /*
   * auxiliarry arrays:
   * fM2 = multiloop region with exactly two stems, extending to 3' end
   * for stupid dangles=1 case we also need:
   * fM_d3 = multiloop region with >= 2 stems, starting at pos 2
   *         (a pair (k,n) will form 3' dangle with pos 1)
   * fM_d5 = multiloop region with >= 2 stems, extending to pos n-1
   *         (a pair (1,k) will form a 5' dangle with pos n)
   */
  unsigned char *hard_constraints, eval;
  char          *ptype;
  short         *S1;
  unsigned int  **a2s;
  int           Hi, Hj, Ii, Ij, Ip, Iq, ip, iq, Mi, *fM_d3, *fM_d5, Md3i,
                Md5i, FcMd3, FcMd5, FcH, FcI, FcM, Fc, *fM2, i, j, ij, u,
                length, new_c, fm, type, *my_c, *my_fML, *indx, FcO, tmp,
                dangle_model, turn, s;
  vrna_param_t  *P;
  vrna_hc_t     *hc;
  vrna_sc_t     *sc, **scs;

  length            = vc->length;
  P                 = vc->params;
  ptype             = (vc->type == VRNA_FC_TYPE_SINGLE) ? vc->ptype : NULL;
  indx              = vc->jindx;
  S1                = (vc->type == VRNA_FC_TYPE_SINGLE) ? vc->sequence_encoding : NULL;
  a2s               = (vc->type == VRNA_FC_TYPE_COMPARATIVE) ? vc->a2s : NULL;
  hc                = vc->hc;
  sc                = (vc->type == VRNA_FC_TYPE_SINGLE) ? vc->sc : NULL;
  scs               = (vc->type == VRNA_FC_TYPE_COMPARATIVE) ? vc->scs : NULL;
  dangle_model      = P->model_details.dangles;
  turn              = P->model_details.min_loop_size;
  hard_constraints  = hc->matrix;
  my_c              = vc->matrices->c;
  my_fML            = vc->matrices->fML;
  fM2               = vc->matrices->fM2;

  Fc = FcO = FcH = FcI = FcM = FcMd3 = FcMd5 = INF;

  /* unfolded state */
  eval = (hc->up_ext[1] >= length) ? 1 : 0;
  if (hc->f)
    eval = (hc->f(1, length, 1, length, VRNA_DECOMP_EXT_UP, hc->data)) ? eval : 0;

  if (eval) {
    Fc = 0; /* base line for unfolded state */

    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if(sc){
          if(sc->energy_up)
            Fc += sc->energy_up[1][length];

          if (sc->f)
            Fc += sc->f(1, length, 1, length, VRNA_DECOMP_EXT_UP, sc->data);
        }
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        if (scs) {
          for (s = 0; s < vc->n_seq; s++)
            if (scs[s]) {
              if (scs[s]->energy_up)
                Fc += scs[s]->energy_up[1][a2s[s][length]];
            }
        }

        break;
    }
    FcO = Fc;
  } else {
    Fc = INF;
  }

  for (i = 1; i < length; i++)
    for (j = i + turn + 1; j <= length; j++) {
      u = length - j + i - 1;
      if (u < turn)
        continue;

      ij = indx[j] + i;

      if (!hard_constraints[ij])
        continue;

      /* exterior hairpin case */
      new_c = vrna_E_hp_loop(vc, j, i);
      if (new_c != INF)
        new_c += my_c[ij];

      if (new_c < FcH) {
        FcH = new_c;
        Hi  = i;
        Hj  = j;
      }

      /* exterior interior loop case */
      ip = iq = 0;
      new_c = vrna_E_ext_int_loop(vc, i, j, &ip, &iq);
      if (new_c != INF)
        new_c += my_c[ij];

      if (ip != 0) {
        if (new_c < FcI) {
          FcI = new_c;
          Ii  = i;
          Ij  = j;
          Ip  = ip;
          Iq  = iq;
        }
      }
    } /* end of i,j loop */
  Fc = MIN2(Fc, FcH);
  Fc = MIN2(Fc, FcI);

  /*
   * compute the fM2 array (multi loops with exactly 2 helices)
   * to get a unique ML decomposition, just use fM1 instead of fML
   * below. However, that will not work with dangle_model==1
   */
  if ((sc) && (sc->f)) {
    if (hc->f) {
      for (i=1; i<length-turn; i++) {
        fM2[i] = INF;
        for (u=i+turn; u<length-turn; u++) {
          if (hc->f(i, length, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
            fm = sc->f(i, length, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
            if ((fm != INF) && (my_fML[indx[u]+i] != INF) && (my_fML[indx[length]+u+1] != INF)) {
              fm += my_fML[indx[u]+i] +
                    my_fML[indx[length]+u+1];
              fM2[i] = MIN2(fM2[i], fm);
            }
          }
        }
      }
    } else {
      for (i=1; i<length-turn; i++) {
        fM2[i] = INF;
        for (u=i+turn; u<length-turn; u++) {
          fm = sc->f(i, length, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
          if ((fm != INF) && (my_fML[indx[u]+i] != INF) && (my_fML[indx[length]+u+1] != INF)) {
            fm += my_fML[indx[u]+i] +
                  my_fML[indx[length]+u+1];
            fM2[i] = MIN2(fM2[i], fm);
          }
        }
      }
    }
  } else {
    if (hc->f) {
      for (i=1; i<length-turn; i++) {
        fM2[i] = INF;
        for (u=i+turn; u<length-turn; u++) {
          if (hc->f(i, length, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
            fm = my_fML[indx[u]+i];
            if ((fm != INF) && (my_fML[indx[length]+u+1] != INF)) {
              fm += my_fML[indx[length]+u+1];
              fM2[i] = MIN2(fM2[i], fm);
            }
          }
        }
      }
    } else {
      for (i=1; i<length-turn; i++) {
        fM2[i] = INF;
        for (u=i+turn; u<length-turn; u++) {
          fm = my_fML[indx[u]+i];
          if ((fm != INF) && (my_fML[indx[length]+u+1] != INF)) {
            fm += my_fML[indx[length]+u+1];
            fM2[i] = MIN2(fM2[i], fm);
          }
        }
      }
    }
  }

  if ((sc) && (sc->f)) {
    if (hc->f) {
      for (i=turn+1; i<length-2*turn; i++) {
        if (hc->f(1, length, i, i + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
          fm = sc->f(1, length, i, i + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
          if ((fm != INF) && (my_fML[indx[i]+1] != INF) && (fM2[i+1] != INF)) {
            fm += my_fML[indx[i]+1] +
                  fM2[i+1];

            if (fm<FcM) {
              FcM=fm; Mi=i;
            }
          }
        }
      }
    } else {
      for (i=turn+1; i<length-2*turn; i++) {
        fm = sc->f(1, length, i, i + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
        if ((fm != INF) && (my_fML[indx[i]+1] != INF) && (fM2[i+1] != INF)) {
          fm += my_fML[indx[i]+1] +
                fM2[i+1];

          if (fm<FcM) {
            FcM=fm; Mi=i;
          }
        }
      }
    }
  } else {
    if (hc->f) {
      for (i=turn+1; i<length-2*turn; i++) {
        if (hc->f(1, length, i, i + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
          fm = my_fML[indx[i]+1];
          if ((fm != INF) && (fM2[i+1] != INF)) {
            fm += fM2[i+1];

            if (fm<FcM) {
              FcM=fm; Mi=i;
            }
          }
        }
      }
    } else {
      for (i=turn+1; i<length-2*turn; i++) {
        fm = my_fML[indx[i]+1];
        if ((fm != INF) && (fM2[i+1] != INF)) {
          fm += fM2[i+1];

          if (fm<FcM) {
            FcM=fm; Mi=i;
          }
        }
      }
    }
  }

  if (FcM != INF) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        FcM += P->MLclosing;
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        FcM += vc->n_seq * P->MLclosing;
        break;
    }
  }
  Fc = MIN2(Fc, FcM);

  /*
      add multibranch loop configurations for odd dangle models
      not supported for comparative prediction (yet)
  */
  if ((vc->type == VRNA_FC_TYPE_SINGLE) && (dangle_model == 1) || (dangle_model == 3)) {
    fM_d3 =  (int *) vrna_alloc(sizeof(int)*(length+2));
    fM_d5 =  (int *) vrna_alloc(sizeof(int)*(length+2));

    for (i=turn+1; i<length-turn; i++)
      fM_d3[i] = INF;

    for (i=turn+1; i<length-turn; i++)
      fM_d5[i] = INF;

    if (hc->f) {
      if (hc->up_ml[1]) {
        /* 1st, construct fM_d3 array */
        for (i=turn+1; i<length-turn; i++) {
          if (hc->f(1, i, 2, i, VRNA_DECOMP_ML_ML, hc->data)) {
            for (u=2+turn; u<i-turn; u++) {
              if (hc->f(2, i, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
                fm = my_fML[indx[u]+2];

                /* skip configurations that violate (hard) constraints */
                if ((fm == INF) || (my_fML[indx[i]+u+1] == INF))
                  continue;

                fm += my_fML[indx[i]+u+1];

                if (sc) {
                  if (sc->energy_up)
                    fm += sc->energy_up[1][1];

                  if (sc->f) {
                    tmp = sc->f(1, i, 2, i, VRNA_DECOMP_ML_ML, sc->data);

                    if (tmp == INF)
                      continue;

                    fm += tmp;

                    tmp = sc->f(2, i, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

                    if (tmp == INF)
                      continue;

                    fm += tmp;
                  }
                }

                fM_d3[i] = MIN2(fM_d3[i], fm);
              }
            }
          }
        }

        /* 2nd, test ML loops with closing pair (i + 1, length), 3' dangle pos 1 */
        for (i=2*turn+1; i<length-turn; i++) {
          type = ptype[indx[length]+i+1];
          eval = (hard_constraints[indx[length] + i + 1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) ? 1 : 0;
          eval = hc->f(i + 1, length, 2, i, VRNA_DECOMP_PAIR_ML_EXT, hc->data) ? eval : 0;
          if (eval) {
            fm = fM_d3[i];
            if ((fm != INF) && (my_c[indx[length]+i+1] != INF)) {
              fm += my_c[indx[length]+i+1] +
                    E_MLstem(type, -1, S1[1], P) +
                    P->MLclosing;

              if (sc) {
                if (sc->f) {
                  tmp = sc->f(i+1, length, 2, i, VRNA_DECOMP_PAIR_ML_EXT, sc->data);

                  if (tmp == INF)
                    continue;

                  fm += tmp;
                }
              }

              if (fm<FcMd3) {
                FcMd3=fm; Md3i=i;
              }
            }
          }
        }

        /* 3rd, test ML loops with closing pair (i + 1, length), 5' dangle pos i, 3' dangle pos 1 */
        for (i=2*turn+1; i<length-turn; i++) {
          type = ptype[indx[length]+i+1];
          eval = (hard_constraints[indx[length] + i + 1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) ? 1 : 0;
          eval = hc->up_ml[i] ? eval : 0;
          eval = hc->f(i + 1, length, 2, i - 1, VRNA_DECOMP_PAIR_ML_EXT, hc->data) ? eval : 0;
          if (eval) {
            fm = fM_d3[i-1];
            if ((fm != INF) && (my_c[indx[length]+i+1] != INF)) {
              fm += my_c[indx[length]+i+1] +
                    E_MLstem(type, S1[i], S1[1], P) +
                    P->MLclosing;

              if (sc) {
                if (sc->energy_up)
                  fm += sc->energy_up[i][1];

                if (sc->f) {
                  tmp = sc->f(i + 1, length, 2, i - 1, VRNA_DECOMP_PAIR_ML_EXT, sc->data);

                  if (tmp == INF)
                    continue;

                  fm += tmp;
                }
              }

              if (fm<FcMd3) {
                FcMd3=fm; Md3i=-i;
              }
            }
          }
        }
      }
    } else {
      if (hc->up_ml[1]) {
        /* 1st, construct fM_d3 array */
        for (i=turn+1; i<length-turn; i++) {
          for (u=2+turn; u<i-turn; u++) {
            fm = my_fML[indx[u]+2];
            if ((fm != INF) && (my_fML[indx[i]+u+1] != INF)) {
              fm += my_fML[indx[i]+u+1];

              if (sc) {
                if (sc->energy_up)
                  fm += sc->energy_up[1][1];

                if (sc->f) {
                  tmp = sc->f(1, i, 2, i, VRNA_DECOMP_ML_ML, sc->data);

                  if (tmp == INF)
                    continue;

                  fm += tmp;

                  tmp = sc->f(2, i, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

                  if (tmp == INF)
                    continue;

                  fm += tmp;
                }
              }

              fM_d3[i] = MIN2(fM_d3[i], fm);
            }
          }
        }

        /* 2nd, test ML loops with closing pair (i + 1, length), 3' dangle pos 1 */
        for (i=2*turn+1; i<length-turn; i++) {
          if(hard_constraints[indx[length] + i + 1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP){
            fm = fM_d3[i];
            if ((fm != INF) && (my_c[indx[length]+i+1] != INF)) {
              type = ptype[indx[length]+i+1];
              fm += my_c[indx[length]+i+1]+
                    E_MLstem(type, -1, S1[1], P) +
                    P->MLclosing;

              if (sc) {
                if (sc->f) {
                  tmp = sc->f(i + 1, length, 2, i, VRNA_DECOMP_PAIR_ML_EXT, sc->data);

                  if (tmp == INF)
                    continue;

                  fm += tmp;
                }
              }

              if (fm<FcMd3) {
                FcMd3=fm; Md3i=i;
              }
            }
          }
        }

        /* 3rd, test ML loops with closing pair (i + 1, length), 5' dangle pos i, 3' dangle pos 1 */
        for (i=2*turn+1; i<length-turn; i++) {
          if(hard_constraints[indx[length] + i + 1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP){
            if(hc->up_ml[i]){
              fm = fM_d3[i-1];
              if ((fm != INF) && (my_c[indx[length]+i+1] != INF)) {
                type = ptype[indx[length]+i+1];
                fm += my_c[indx[length]+i+1]+
                      E_MLstem(type, S1[i], S1[1], P) +
                      P->MLclosing;

                if (sc) {
                  if (sc->energy_up)
                    fm += sc->energy_up[i][1];

                  if(sc->f) {
                    tmp = sc->f(i + 1, length, 2, i - 1, VRNA_DECOMP_PAIR_ML_EXT, sc->data);

                    if (tmp == INF)
                      continue;

                    fm += tmp;
                  }
                }

                if (fm<FcMd3) {
                  FcMd3=fm; Md3i=-i;
                }
              }
            }
          }
        }
      }
    }

    if (hc->f) {
      if (hc->up_ml[length]) {
        /* 1st, construct fM_d5 array */
        for (i=turn+1; i<length-turn; i++) {
          if (hc->f(i, length, i, length - 1, VRNA_DECOMP_ML_ML, hc->data)) {
            for (u=i+turn; u<length-turn; u++) {
              if (hc->f(i, length - 1, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
                fm = my_fML[indx[u]+i];
                if ((fm != INF) && (my_fML[indx[length-1]+u+1] != INF)) {
                  fm += my_fML[indx[length-1]+u+1];

                  if (sc) {
                    if (sc->energy_up)
                      fm += sc->energy_up[length][1];

                    if (sc->f) {
                      tmp = sc->f(i, length, i, length - 1, VRNA_DECOMP_ML_ML, sc->data);

                      if (tmp == INF)
                        continue;

                      fm += tmp;

                      tmp = sc->f(i, length - 1, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

                      if (tmp == INF)
                        continue;

                      fm += tmp;
                    }
                  }

                  fM_d5[i] = MIN2(fM_d5[i], fm);
                }
              }
            }
          }
        }

        /* 2nd, test ML loops with closing pair (1, i), 5' dangle pos n */
        for (i=turn+1; i<length-2*turn; i++) {
          eval = (hard_constraints[indx[i]+1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) ? 1 : 0;
          eval = hc->f(1, i, i + 1, length - 1, VRNA_DECOMP_PAIR_ML_EXT, hc->data) ? eval : 0;
          if(eval){
            fm = fM_d5[i+1];
            if ((fm != INF) && (my_c[indx[i]+1] != INF)) {
              type = ptype[indx[i]+1];
              fm += my_c[indx[i]+1] +
                    E_MLstem(type, S1[length], -1, P) +
                    P->MLclosing;

              if (sc) {
                if (sc->f) {
                  tmp = sc->f(1, i, i + 1, length - 1, VRNA_DECOMP_PAIR_ML_EXT, sc->data);

                  if (tmp == INF)
                    continue;

                  fm += tmp;
                }
              }

              if (fm<FcMd5) {
                FcMd5=fm; Md5i=i;
              }
            }
          }
        }

        /* 3rd, test ML loops with closing pair (1, i), 5' dangle pos n, 3' dangle pos i + 1 */
        for (i=turn+1; i<length-2*turn; i++) {
          eval = (hard_constraints[indx[i]+1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) ? 1 : 0;
          eval = hc->f(1, i, i + 1, length - 1, VRNA_DECOMP_PAIR_ML_EXT, hc->data) ? eval : 0;
          eval = hc->up_ml[i+1] ? eval : 0;
          if (eval) {
            fm = fM_d5[i+2];
            if ((fm != INF) && (my_c[indx[i]+1] != INF)) {
              type = ptype[indx[i]+1];
              fm += my_c[indx[i]+1] +
                    E_MLstem(type, S1[length], S1[i+1], P) +
                    P->MLclosing;

              if (sc) {
                if (sc->energy_up)
                  fm += sc->energy_up[i + 1][1];

                if (sc->f) {
                  tmp = sc->f(i, length - 1, i + 1, length - 1, VRNA_DECOMP_ML_ML, sc->data);

                  if (tmp == INF)
                    continue;

                  fm += tmp;

                  tmp = sc->f(1, i, i + 2, length - 1, VRNA_DECOMP_PAIR_ML_EXT, sc->data);

                  if (tmp == INF)
                    continue;

                  fm += tmp;
                }
              }

              if (fm<FcMd5) {
                FcMd5=fm; Md5i=-i;
              }
            }
          }
        }
      }
    } else {
      if (hc->up_ml[length]) {
        /* 1st, construct fM_d5 array */
        for (i=turn+1; i<length-turn; i++) {
          for (u=i+turn; u<length-turn; u++) {
            fm = my_fML[indx[u]+i];
            if ((fm != INF) && (my_fML[indx[length-1]+u+1] != INF)) {
              fm += my_fML[indx[length-1]+u+1];

              if (sc) {
                if (sc->energy_up)
                  fm += sc->energy_up[length][1];

                if (sc->f) {
                  tmp = sc->f(i, length, i, length - 1, VRNA_DECOMP_ML_ML, sc->data);

                  if (tmp == INF)
                    continue;

                  fm += tmp;

                  tmp = sc->f(i, length - 1, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

                  if (tmp == INF)
                    continue;

                  fm += tmp;
                }
              }

              fM_d5[i] = MIN2(fM_d5[i], fm);
            }
          }
        }

        /* 2nd, test ML loops with closing pair (1, i), 5' dangle pos n */
        for (i=turn+1; i<length-2*turn; i++) {
          if(hard_constraints[indx[i]+1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP){
            fm = fM_d5[i+1];
            if ((fm != INF) && (my_c[indx[i]+1] != INF)) {
              type = ptype[indx[i]+1];
              fm += my_c[indx[i]+1] +
                    E_MLstem(type, S1[length], -1, P) +
                    P->MLclosing;

              if (sc) {
                if (sc->f) {
                  tmp = sc->f(1, i, i + 1, length - 1, VRNA_DECOMP_PAIR_ML_EXT, sc->data);

                  if (tmp == INF)
                    continue;

                  fm += tmp;
                }
              }

              if (fm<FcMd5) {
                FcMd5=fm; Md5i=i;
              }
            }
          }
        }

        /* 3rd, test ML loops with closing pair (1, i), 5' dangle pos n, 3' dangle pos i + 1 */
        for (i=turn+1; i<length-2*turn; i++) {
          if(hard_constraints[indx[i]+1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP){
            if(hc->up_ml[i+1]){
              fm = fM_d5[i+2];
              if ((fm != INF) && (my_c[indx[i]+1] != INF)) {
                type = ptype[indx[i]+1];
                fm += my_c[indx[i]+1] + 
                      E_MLstem(type, S1[length], S1[i+1], P) +
                      P->MLclosing;

                if (sc) {
                  if (sc->energy_up)
                    fm += sc->energy_up[i + 1][1];

                  if (sc->f) {
                    tmp = sc->f(i, length - 1, i + 1, length - 1, VRNA_DECOMP_ML_ML, sc->data);

                    if (tmp == INF)
                      continue;

                    fm += tmp;

                    tmp = sc->f(1, i, i + 2, length - 1, VRNA_DECOMP_PAIR_ML_EXT, sc->data);

                    if (tmp == INF)
                      continue;

                    fm += tmp;
                  }
                }

                if (fm<FcMd5) {
                  FcMd5=fm; Md5i=-i;
                }
              }
            }
          }
        }
      }
    }

    if (FcMd5<MIN2(Fc,FcMd3)) {
      int real_i, sc_en = 0;

      /* looks like we have to do this ... */
      bt_stack[++(*bt)].i = 1;
      bt_stack[(*bt)].j = (Md5i>0)?Md5i:-Md5i;
      bt_stack[(*bt)].ml = 2;
      i = (Md5i>0)?Md5i+1 : -Md5i+2; /* let's backtrack fm_d5[Md5i+1] */
      real_i = (Md5i > 0) ? i : i - 1;

      if ((sc) && (sc->energy_up))
        sc_en += sc->energy_up[length][1];

      for (u=i+turn; u<length-turn; u++) {
        fm = my_fML[indx[u]+i] +
             my_fML[indx[length-1]+u+1] +
             sc_en;

        if (sc) {
          if (sc->energy_up)
            fm += sc->energy_up[real_i][i - real_i];

          if (sc->f)
            fm += sc->f(real_i, length, i, length - 1, VRNA_DECOMP_ML_ML, sc->data) +
                  sc->f(i, length - 1, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
        }

        if (fM_d5[i] == fm) {
          bt_stack[++(*bt)].i = i;
          bt_stack[(*bt)].j = u;
          bt_stack[(*bt)].ml = 1;
          bt_stack[++(*bt)].i =u+1;
          bt_stack[(*bt)].j = length-1;
          bt_stack[(*bt)].ml = 1;
          break;
        }
      }
      Fc = FcMd5;
    } else if (FcMd3<Fc) {
      int real_i, sc_en = 0;
      /* here we go again... */
      bt_stack[++(*bt)].i = (Md3i>0)?Md3i+1:-Md3i+1;
      bt_stack[(*bt)].j = length;
      bt_stack[(*bt)].ml = 2;
      i = (Md3i>0)? Md3i : -Md3i-1; /* let's backtrack fm_d3[Md3i] */
      real_i = (Md3i > 0) ? i : i + 1;

      if ((sc) && (sc->energy_up))
        sc_en += sc->energy_up[1][1];

      for (u=2+turn; u<i-turn; u++) {
        fm = my_fML[indx[u]+2] +
             my_fML[indx[i]+u+1] +
             sc_en;

        if (sc) {
          if (sc->energy_up)
            fm += sc->energy_up[real_i][real_i - i];

          if (sc->f)
            fm += sc->f(1, real_i, 2, i, VRNA_DECOMP_ML_ML, sc->data) +
                  sc->f(2, i, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
        }

        if (fM_d3[i] == fm) {
          bt_stack[++(*bt)].i = 2;
          bt_stack[(*bt)].j = u;
          bt_stack[(*bt)].ml = 1;
          bt_stack[++(*bt)].i =u+1;
          bt_stack[(*bt)].j = i;
          bt_stack[(*bt)].ml = 1;
          break;
        }
      }
      Fc = FcMd3;
    }
    free(fM_d3);
    free(fM_d5);
  }
  
  if(Fc < INF){
    if (FcH==Fc) {
      bt_stack[++(*bt)].i = Hi;
      bt_stack[(*bt)].j = Hj;
      bt_stack[(*bt)].ml = 2;
    }
    else if (FcI==Fc) {
      bt_stack[++(*bt)].i = Ii;
      bt_stack[(*bt)].j = Ij;
      bt_stack[(*bt)].ml = 2;
      bt_stack[++(*bt)].i = Ip;
      bt_stack[(*bt)].j = Iq;
      bt_stack[(*bt)].ml = 2;
    }
    else if (FcM==Fc) { /* grumpf we found a Multiloop */
      int eee;
      /* backtrack in fM2 */
      fm = fM2[Mi+1];
      for (u=Mi+turn+1; u<length-turn; u++) {
        eee = my_fML[indx[u]+Mi+1] +
              my_fML[indx[length]+u+1];

        if (sc) {
          if (sc->f)
            eee += sc->f(Mi + 1, length, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
        }

        if (fm == eee) {
                bt_stack[++(*bt)].i=Mi+1;
                bt_stack[(*bt)].j=u;
                bt_stack[(*bt)].ml = 1;
                bt_stack[++(*bt)].i=u+1;
                bt_stack[(*bt)].j=length;
                bt_stack[(*bt)].ml = 1;
                break;
        }
      }
      bt_stack[++(*bt)].i = 1;
      bt_stack[(*bt)].j = Mi;
      bt_stack[(*bt)].ml = 1;
    } else if (Fc == FcO) { /* unstructured */
      bt_stack[++(*bt)].i = 1;
      bt_stack[(*bt)].j = 1;
      bt_stack[(*bt)].ml = 0;
    }
  } else {
    /* forbidden, i.e. no configuration fulfills constraints */
  }
  vc->matrices->FcH = FcH;
  vc->matrices->FcI = FcI;
  vc->matrices->FcM = FcM;
  vc->matrices->Fc  = Fc;
  return Fc;
}


PUBLIC void
vrna_backtrack_from_intervals(vrna_fold_compound_t  *vc,
                              vrna_bp_stack_t       *bp_stack,
                              sect                  bt_stack[],
                              int                   s)
{
  if (vc)
    backtrack(vc, bp_stack, bt_stack, s);
}


/**
*** trace back through the "c", "f5" and "fML" arrays to get the
*** base pairing list. No search for equivalent structures is done.
*** This is fast, since only few structure elements are recalculated.
***
*** normally s=0.
*** If s>0 then s items have been already pushed onto the bt_stack
**/
PRIVATE void
backtrack(vrna_fold_compound_t  *vc,
          vrna_bp_stack_t       *bp_stack,
          sect                  bt_stack[],
          int                   s)
{
  char          backtrack_type;
  int           i, j, ij, k, length, b, *my_c, *indx, noLP, *pscore;
  vrna_param_t  *P;

  b               = 0;
  length          = vc->length;
  my_c            = vc->matrices->c;
  indx            = vc->jindx;
  P               = vc->params;
  noLP            = P->model_details.noLP;
  pscore          = vc->pscore;         /* covariance scores for comparative structure prediction */
  backtrack_type  = P->model_details.backtrack_type;

  if (s == 0) {
    bt_stack[++s].i = 1;
    bt_stack[s].j   = length;
    bt_stack[s].ml  = (backtrack_type == 'M') ? 1 : ((backtrack_type == 'C') ? 2 : 0);
  }

  while (s > 0) {
    int ml, cij;
    int canonical = 1;     /* (i,j) closes a canonical structure */

    /* pop one element from stack */
    i   = bt_stack[s].i;
    j   = bt_stack[s].j;
    ml  = bt_stack[s--].ml;

    switch (ml) {
      /* backtrack in f5 */
      case 0:
      {
        int p, q;
        if (vrna_BT_ext_loop_f5(vc, &j, &p, &q, bp_stack, &b)) {
          if (j > 0) {
            bt_stack[++s].i = 1;
            bt_stack[s].j   = j;
            bt_stack[s].ml  = 0;
          }

          if (p > 0) {
            i = p;
            j = q;
            goto repeat1;
          }

          continue;
        } else {
          vrna_message_error("backtracking failed in f5, segment [%d,%d]\n", i, j);
        }
      }
      break;

      /* trace back in fML array */
      case 1:
      {
        int p, q, comp1, comp2;
        if (vrna_BT_mb_loop_split(vc, &i, &j, &p, &q, &comp1, &comp2, bp_stack, &b)) {
          if (i > 0) {
            bt_stack[++s].i = i;
            bt_stack[s].j   = j;
            bt_stack[s].ml  = comp1;
          }

          if (p > 0) {
            bt_stack[++s].i = p;
            bt_stack[s].j   = q;
            bt_stack[s].ml  = comp2;
          }

          continue;
        } else {
          vrna_message_error("backtracking failed in fML, segment [%d,%d]\n", i, j);
        }
      }
      break;

      /* backtrack in c */
      case 2:
        bp_stack[++b].i = i;
        bp_stack[b].j   = j;
        goto repeat1;

        break;

      default:
        vrna_message_error("Backtracking failed due to unrecognized DP matrix!");
        break;
    }

repeat1:

    /*----- begin of "repeat:" -----*/
    ij = indx[j] + i;

    if (canonical)
      cij = my_c[ij];

    if (noLP) {
      if (vrna_BT_stack(vc, &i, &j, &cij, bp_stack, &b)) {
        canonical = 0;
        goto repeat1;
      }
    }

    canonical = 1;

    if (vc->type == VRNA_FC_TYPE_COMPARATIVE)
      cij += pscore[indx[j] + i];

    if (vrna_BT_hp_loop(vc, i, j, cij, bp_stack, &b))
      continue;

    if (vrna_BT_int_loop(vc, &i, &j, cij, bp_stack, &b)) {
      if (i < 0)
        continue;
      else
        goto repeat1;
    }

    /* (i.j) must close a multi-loop */
    int comp1, comp2;

    if (vrna_BT_mb_loop(vc, &i, &j, &k, cij, &comp1, &comp2)) {
      bt_stack[++s].i = i;
      bt_stack[s].j   = k;
      bt_stack[s].ml  = comp1;
      bt_stack[++s].i = k + 1;
      bt_stack[s].j   = j;
      bt_stack[s].ml  = comp2;
    } else {
      vrna_message_error("backtracking failed in repeat, segment [%d,%d]\n", i, j);
    }

    /* end of repeat: --------------------------------------------------*/
  } /* end of infinite while loop */

  bp_stack[0].i = b;    /* save the total number of base pairs */
}
