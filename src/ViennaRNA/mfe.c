/** \file **/

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

/* make this interface backward compatible with RNAlib < 2.2.0 */
#define VRNA_BACKWARD_COMPAT

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

PRIVATE int   fill_arrays(vrna_fold_compound_t *vc);


PRIVATE void  fill_arrays_circ(vrna_fold_compound_t *vc,
                               sect                 bt_stack[],
                               int                  *bt);


PRIVATE void  backtrack(vrna_fold_compound_t  *vc,
                        vrna_bp_stack_t       *bp_stack,
                        sect                  bt_stack[],
                        int                   s);


PRIVATE int   fill_arrays_comparative(vrna_fold_compound_t *vc);


PRIVATE void  fill_arrays_comparative_circ(vrna_fold_compound_t *vc,
                                           sect                 bt_stack[],
                                           int                  *bt);


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

    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        energy = fill_arrays(vc);
        if (vc->params->model_details.circ) {
          fill_arrays_circ(vc, bt_stack, &s);
          energy = vc->matrices->Fc;
        }

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        energy = fill_arrays_comparative(vc);
        if (vc->params->model_details.circ) {
          fill_arrays_comparative_circ(vc, bt_stack, &s);
          energy = vc->matrices->Fc;
        }

        break;

      default:
        vrna_message_warning("unrecognized fold compound type");
        return mfe;
        break;
    }


    /* call user-defined recursion status callback function */
    if (vc->stat_cb)
      vc->stat_cb(VRNA_STATUS_MFE_POST, vc->auxdata);

    if (structure && vc->params->model_details.backtrack) {
      bp = (vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (4 * (1 + length / 2))); /* add a guess of how many G's may be involved in a G quadruplex */

      backtrack(vc, bp, bt_stack, s);

      ss = vrna_db_from_bp_stack(bp, length);
      strncpy(structure, ss, length + 1);
      free(ss);
      free(bp);
    }

    if (vc->params->model_details.backtrack_type == 'C')
      mfe = (float)vc->matrices->c[vc->jindx[length] + 1] / 100.;
    else if (vc->params->model_details.backtrack_type == 'M')
      mfe = (float)vc->matrices->fML[vc->jindx[length] + 1] / 100.;
    else
      mfe = (float)energy / 100.;

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
  unsigned char type;
  char          *ptype;
  unsigned char *hard_constraints;
  int           i, j, ij, length, energy, new_c, stackEnergy, no_close, turn,
                noGUclosure, noLP, uniq_ML, dangle_model, *indx, *f5,
                *c, *fML, *fM1, hc_decompose, *cc, *cc1, *Fmi, *DMLi,
                *DMLi1, *DMLi2;
  vrna_param_t  *P;
  vrna_mx_mfe_t *matrices;
  vrna_hc_t     *hc;
  vrna_ud_t     *domains_up;

  length            = (int)vc->length;
  ptype             = vc->ptype;
  indx              = vc->jindx;
  P                 = vc->params;
  noGUclosure       = P->model_details.noGUclosure;
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
      type          = (unsigned char)ptype[ij];
      hc_decompose  = hard_constraints[ij];
      energy        = INF;

      no_close = (((type == 3) || (type == 4)) && noGUclosure);

      if (hc_decompose) {
        /* we evaluate this pair */
        new_c = INF;

        if (!no_close) {
          /* check for hairpin loop */
          energy  = vrna_E_hp_loop(vc, i, j);
          new_c   = MIN2(new_c, energy);

          /* check for multibranch loops */
          energy  = vrna_E_mb_loop_fast(vc, i, j, DMLi1, DMLi2);
          new_c   = MIN2(new_c, energy);
        }

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
        } else {
          c[ij] = new_c;
        }
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


#include "circfold.inc"


/**
*** the actual forward recursion to fill the energy arrays
**/
PRIVATE int
fill_arrays_comparative(vrna_fold_compound_t *vc)
{
  unsigned char *hard_constraints;
  short         **S;
  int           i, j, turn, energy, stackEnergy, new_c, s, *cc,
                *cc1, *Fmi, *DMLi, *DMLi1, *DMLi2, n_seq, length, *indx,
                *c, *f5, *fML, *pscore;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;

  n_seq             = vc->n_seq;
  length            = vc->length;
  S                 = vc->S;
  P                 = vc->params;
  md                = &(P->model_details);
  indx              = vc->jindx;          /* index for moving in the triangle matrices c[] and fMl[] */
  c                 = vc->matrices->c;    /* energy array, given that i-j pair */
  f5                = vc->matrices->f5;   /* energy of 5' end */
  fML               = vc->matrices->fML;  /* multi-loop auxiliary energy array */
  pscore            = vc->pscore;         /* precomputed array of pair types */
  turn              = md->min_loop_size;
  hc                = vc->hc;
  hard_constraints  = hc->matrix;

  /* allocate some memory for helper arrays */
  cc    = (int *)vrna_alloc(sizeof(int) * (length + 2));  /* linear array for calculating canonical structures */
  cc1   = (int *)vrna_alloc(sizeof(int) * (length + 2));  /*   "     "        */
  Fmi   = (int *)vrna_alloc(sizeof(int) * (length + 1));  /* holds row i of fML (avoids jumps in memory) */
  DMLi  = (int *)vrna_alloc(sizeof(int) * (length + 1));  /* DMLi[j] holds MIN(fML[i,k]+fML[k+1,j])  */
  DMLi1 = (int *)vrna_alloc(sizeof(int) * (length + 1));  /*             MIN(fML[i+1,k]+fML[k+1,j])  */
  DMLi2 = (int *)vrna_alloc(sizeof(int) * (length + 1));  /*             MIN(fML[i+2,k]+fML[k+1,j])  */


  if ((turn < 0) || (turn > length))
    turn = length;

  /* init energies */
  for (j = 1; j <= length; j++) {
    Fmi[j] = DMLi[j] = DMLi1[j] = DMLi2[j] = INF;
    for (i = (j > turn ? (j - turn) : 1); i < j; i++)
      c[indx[j] + i] = fML[indx[j] + i] = INF;
  }

  /* begin recursions */
  for (i = length - turn - 1; i >= 1; i--) {
    /* i,j in [1..length] */
    for (j = i + turn + 1; j <= length; j++) {
      int ij, psc;
      ij = indx[j] + i;

      psc = pscore[indx[j] + i];
      if (hard_constraints[ij]) {
        /* a pair to consider */
        new_c = INF;

        /* hairpin ----------------------------------------------*/
        energy  = vrna_E_hp_loop(vc, i, j);
        new_c   = MIN2(new_c, energy);

        /* check for multibranch loops */
        energy  = vrna_E_mb_loop_fast(vc, i, j, DMLi1, DMLi2);
        new_c   = MIN2(new_c, energy);

        /* check for interior loops */
        energy  = vrna_E_int_loop(vc, i, j);
        new_c   = MIN2(new_c, energy);

        /* remember stack energy for --noLP option */
        if (md->noLP) {
          stackEnergy = vrna_E_stack(vc, i, j);
          new_c       = MIN2(new_c, cc1[j - 1] + stackEnergy);
          cc[j]       = new_c - psc; /* add covariance bonnus/penalty */
          c[ij]       = cc1[j - 1] + stackEnergy - psc;
        } else {
          c[ij] = new_c - psc;       /* add covariance bonnus/penalty */
        }
      } /* end >> if (pair) << */
      else {
        c[ij] = INF;
      }

      /* done with c[i,j], now compute fML[i,j] */
      fML[ij] = vrna_E_ml_stems_fast(vc, i, j, Fmi, DMLi);
    } /* END for j */

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
  } /* END for i */

  /* calculate energies of 5' and 3' fragments */
  (void)vrna_E_ext_loop_5(vc);

  free(cc);
  free(cc1);
  free(Fmi);
  free(DMLi);
  free(DMLi1);
  free(DMLi2);

  return f5[length];
}


#include "ViennaRNA/alicircfold.inc"

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
  unsigned char type;
  char          *ptype, backtrack_type;
  int           i, j, ij, k, length, no_close, b, *my_c, *indx, noLP, noGUclosure, *pscore;
  vrna_param_t  *P;

  b               = 0;
  length          = vc->length;
  my_c            = vc->matrices->c;
  indx            = vc->jindx;
  P               = vc->params;
  noLP            = P->model_details.noLP;
  noGUclosure     = P->model_details.noGUclosure;
  ptype           = vc->ptype;
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

    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        type = (unsigned char)ptype[ij];
        if (type == 0)
          type = 7;

        no_close = (((type == 3) || (type == 4)) && noGUclosure);

        if (no_close) {
          if (cij == FORBIDDEN)
            continue;
        } else {
          if (vrna_BT_hp_loop(vc, i, j, cij, bp_stack, &b))
            continue;
        }

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        cij += pscore[indx[j] + i];
        if (vrna_BT_hp_loop(vc, i, j, cij, bp_stack, &b))
          continue;

        break;
    }

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
