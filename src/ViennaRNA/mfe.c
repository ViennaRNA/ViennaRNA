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

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/structures.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/mfe.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#define MAXSECTORS        500     /* dimension for a backtrack array */

struct aux_arrays {
  int *cc;    /* auxilary arrays for canonical structures     */
  int *cc1;   /* auxilary arrays for canonical structures     */
  int *Fmi;   /* holds row i of fML (avoids jumps in memory)  */
  int *DMLi;  /* DMLi[j] holds  MIN(fML[i,k]+fML[k+1,j])      */
  int *DMLi1; /*                MIN(fML[i+1,k]+fML[k+1,j])    */
  int *DMLi2; /*                MIN(fML[i+2,k]+fML[k+1,j])    */
};


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
fill_arrays(vrna_fold_compound_t *fc);


PRIVATE int
postprocess_circular(vrna_fold_compound_t *fc,
                     sect                 bt_stack[],
                     int                  *bt);


PRIVATE INLINE void
fill_fM_d5(vrna_fold_compound_t *fc,
           int                  *fM_d5);


PRIVATE INLINE void
fill_fM_d3(vrna_fold_compound_t *fc,
           int                  *fM_d3);


PRIVATE int
backtrack(vrna_fold_compound_t  *fc,
          vrna_bp_stack_t       *bp_stack,
          sect                  bt_stack[],
          int                   s);


PRIVATE INLINE int
decompose_pair(vrna_fold_compound_t *fc,
               int                  i,
               int                  j,
               struct aux_arrays    *aux);


PRIVATE INLINE struct aux_arrays *
get_aux_arrays(unsigned int length);


PRIVATE INLINE void
rotate_aux_arrays(struct aux_arrays *aux,
                  unsigned int      length);


PRIVATE INLINE void
free_aux_arrays(struct aux_arrays *aux);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC float
vrna_mfe(vrna_fold_compound_t *fc,
         char                 *structure)
{
  char            *ss;
  int             length, energy, s;
  float           mfe;
  sect            bt_stack[MAXSECTORS]; /* stack of partial structures for backtracking */
  vrna_bp_stack_t *bp;

  s   = 0;
  mfe = (float)(INF / 100.);

  if (fc) {
    length = (int)fc->length;

    if (!vrna_fold_compound_prepare(fc, VRNA_OPTION_MFE)) {
      vrna_message_warning("vrna_mfe@mfe.c: Failed to prepare vrna_fold_compound");
      return mfe;
    }

    /* call user-defined recursion status callback function */
    if (fc->stat_cb)
      fc->stat_cb(VRNA_STATUS_MFE_PRE, fc->auxdata);

    /* call user-defined grammar pre-condition callback function */
    if ((fc->aux_grammar) && (fc->aux_grammar->cb_proc))
      fc->aux_grammar->cb_proc(fc, VRNA_STATUS_MFE_PRE, fc->aux_grammar->data);

    energy = fill_arrays(fc);

    if (fc->params->model_details.circ)
      energy = postprocess_circular(fc, bt_stack, &s);

    if (structure && fc->params->model_details.backtrack) {
      /* add a guess of how many G's may be involved in a G quadruplex */
      bp = (vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (4 * (1 + length / 2)));

      if (backtrack(fc, bp, bt_stack, s) != 0) {
        ss = vrna_db_from_bp_stack(bp, length);
        strncpy(structure, ss, length + 1);
        free(ss);
      } else {
        memset(structure, '\0', sizeof(char) * (length + 1));
      }

      free(bp);
    }

    /* call user-defined recursion status callback function */
    if (fc->stat_cb)
      fc->stat_cb(VRNA_STATUS_MFE_POST, fc->auxdata);

    /* call user-defined grammar post-condition callback function */
    if ((fc->aux_grammar) && (fc->aux_grammar->cb_proc))
      fc->aux_grammar->cb_proc(fc, VRNA_STATUS_MFE_POST, fc->aux_grammar->data);

    switch (fc->params->model_details.backtrack_type) {
      case 'C':
        mfe = (float)fc->matrices->c[fc->jindx[length] + 1] / 100.;
        break;

      case 'M':
        mfe = (float)fc->matrices->fML[fc->jindx[length] + 1] / 100.;
        break;

      default:
        if (fc->type == VRNA_FC_TYPE_COMPARATIVE)
          mfe = (float)energy / (100. * (float)fc->n_seq);
        else
          mfe = (float)energy / 100.;

        break;
    }
  }

  return mfe;
}


PUBLIC int
vrna_backtrack_from_intervals(vrna_fold_compound_t  *fc,
                              vrna_bp_stack_t       *bp_stack,
                              sect                  bt_stack[],
                              int                   s)
{
  if (fc)
    return backtrack(fc, bp_stack, bt_stack, s);

  return 0;
}


PUBLIC float
vrna_backtrack5(vrna_fold_compound_t  *fc,
                unsigned int          length,
                char                  *structure)
{
  char            *ss;
  int             s;
  float           mfe;
  sect            bt_stack[MAXSECTORS]; /* stack of partial structures for backtracking */
  vrna_bp_stack_t *bp;

  s   = 0;
  mfe = (float)(INF / 100.);

  if ((fc) && (structure) && (fc->matrices) && (fc->matrices->f5) &&
      (!fc->params->model_details.circ)) {
    memset(structure, '\0', sizeof(char) * (length + 1));

    if (length > fc->length)
      return mfe;

    /* add a guess of how many G's may be involved in a G quadruplex */
    bp = (vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (4 * (1 + length / 2)));

    bt_stack[++s].i = 1;
    bt_stack[s].j   = length;
    bt_stack[s].ml  = 0;


    if (backtrack(fc, bp, bt_stack, s) != 0) {
      ss = vrna_db_from_bp_stack(bp, length);
      strncpy(structure, ss, length + 1);
      free(ss);

      if (fc->type == VRNA_FC_TYPE_COMPARATIVE)
        mfe = (float)fc->matrices->f5[length] / (100. * (float)fc->n_seq);
      else
        mfe = (float)fc->matrices->f5[length] / 100.;
    }

    free(bp);
  }

  return mfe;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */

/* fill DP matrices */
PRIVATE int
fill_arrays(vrna_fold_compound_t *fc)
{
  int               i, j, ij, length, turn, uniq_ML, *indx, *f5, *c, *fML, *fM1;
  vrna_param_t      *P;
  vrna_mx_mfe_t     *matrices;
  vrna_ud_t         *domains_up;
  struct aux_arrays *helper_arrays;

  length      = (int)fc->length;
  indx        = fc->jindx;
  P           = fc->params;
  uniq_ML     = P->model_details.uniq_ML;
  turn        = P->model_details.min_loop_size;
  matrices    = fc->matrices;
  f5          = matrices->f5;
  c           = matrices->c;
  fML         = matrices->fML;
  fM1         = matrices->fM1;
  domains_up  = fc->domains_up;

  /* allocate memory for all helper arrays */
  helper_arrays = get_aux_arrays(length);

  if ((turn < 0) || (turn > length))
    turn = length; /* does this make any sense? */

  /* pre-processing ligand binding production rule(s) */
  if (domains_up && domains_up->prod_cb)
    domains_up->prod_cb(fc, domains_up->data);

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
    free_aux_arrays(helper_arrays);

    /* return free energy of unfolded chain */
    return 0;
  }

  for (i = length - turn - 1; i >= 1; i--) {
    for (j = i + turn + 1; j <= length; j++) {
      ij = indx[j] + i;

      /* decompose subsegment [i, j] with pair (i, j) */
      c[ij] = decompose_pair(fc, i, j, helper_arrays);

      /* decompose subsegment [i, j] that is multibranch loop part with at least one branch */
      fML[ij] = vrna_E_ml_stems_fast(fc, i, j, helper_arrays->Fmi, helper_arrays->DMLi);

      /* decompose subsegment [i, j] that is multibranch loop part with exactly one branch */
      if (uniq_ML)
        fM1[ij] = E_ml_rightmost_stem(i, j, fc);

      if ((fc->aux_grammar) && (fc->aux_grammar->cb_aux))
        fc->aux_grammar->cb_aux(fc, i, j, fc->aux_grammar->data);
    } /* end of j-loop */

    rotate_aux_arrays(helper_arrays, length);
  } /*
     * end of i-loop
     * calculate energies of 5' fragments
     */
  (void)vrna_E_ext_loop_5(fc);

  /* clean up memory */
  free_aux_arrays(helper_arrays);

  return f5[length];
}


/* post-processing step for circular RNAs */
PRIVATE int
postprocess_circular(vrna_fold_compound_t *fc,
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
  short         *S1, **SS, **S5, **S3;
  unsigned int  **a2s;
  int           Hi, Hj, Ii, Ij, Ip, Iq, ip, iq, Mi, *fM_d3, *fM_d5, Md3i,
                Md5i, FcMd3, FcMd5, FcH, FcI, FcM, Fc, *fM2, i, j, ij, u,
                length, new_c, fm, type, *my_c, *my_fML, *indx, FcO, tmp,
                dangle_model, turn, s, n_seq;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;
  vrna_sc_t     *sc, **scs;

  length            = fc->length;
  n_seq             = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  P                 = fc->params;
  md                = &(P->model_details);
  ptype             = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->ptype : NULL;
  indx              = fc->jindx;
  S1                = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  SS                = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5                = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3                = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  a2s               = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->a2s;
  hc                = fc->hc;
  sc                = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sc : NULL;
  scs               = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->scs;
  dangle_model      = md->dangles;
  turn              = md->min_loop_size;
  hard_constraints  = hc->mx;
  my_c              = fc->matrices->c;
  my_fML            = fc->matrices->fML;
  fM2               = fc->matrices->fM2;

  Fc  = FcO = FcH = FcI = FcM = FcMd3 = FcMd5 = INF;
  Mi  = Md5i = Md3i = Iq = Ip = Ij = Ii = Hj = Hi = 0;

  /* unfolded state */
  eval = (hc->up_ext[1] >= length) ? 1 : 0;
  if (hc->f)
    eval = (hc->f(1, length, 1, length, VRNA_DECOMP_EXT_UP, hc->data)) ? eval : 0;

  if (eval) {
    Fc = 0; /* base line for unfolded state */

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (sc) {
          if (sc->energy_up)
            Fc += sc->energy_up[1][length];

          if (sc->f)
            Fc += sc->f(1, length, 1, length, VRNA_DECOMP_EXT_UP, sc->data);
        }

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        if (scs) {
          for (s = 0; s < fc->n_seq; s++)
            if (scs[s])
              if (scs[s]->energy_up)
                Fc += scs[s]->energy_up[1][a2s[s][length]];
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

      if (!hard_constraints[length * i + j])
        continue;

      /* exterior hairpin case */
      new_c = vrna_E_hp_loop(fc, j, i);
      if (new_c != INF)
        new_c += my_c[ij];

      if (new_c < FcH) {
        FcH = new_c;
        Hi  = i;
        Hj  = j;
      }

      /* exterior interior loop case */
      ip    = iq = 0;
      new_c = vrna_E_ext_int_loop(fc, i, j, &ip, &iq);
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
  Fc  = MIN2(Fc, FcH);
  Fc  = MIN2(Fc, FcI);

  /*
   * compute the fM2 array (multi loops with exactly 2 helices)
   * to get a unique ML decomposition, just use fM1 instead of fML
   * below. However, that will not work with dangle_model==1
   */
  int *fml_tmp = my_fML;

  /* some pre-processing to reduce redundant code */
  if ((hc->f) ||
      ((fc->type == VRNA_FC_TYPE_SINGLE) && (sc)) ||
      ((fc->type == VRNA_FC_TYPE_COMPARATIVE) && (scs)))
    fml_tmp = vrna_alloc(sizeof(int) * (length + 2));
  else
    fml_tmp += indx[length];

  for (i = 1; i < length - turn; i++) {
    fM2[i] = INF;
    /* some preparations in case we have to comply/add certain constraints */
    if ((fml_tmp - indx[length]) != my_fML) {
      /* copy original data */
      for (u = i + turn; u < length - turn; u++)
        fml_tmp[u + 1] = my_fML[indx[length] + u + 1];

      /* apply hard constraints */
      if (hc->f) {
        for (u = i + turn; u < length - turn; u++)
          if (!hc->f(i, length, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data))
            fml_tmp[u + 1] = INF;
      }

      switch (fc->type) {
        case VRNA_FC_TYPE_SINGLE:
          if ((sc) && (sc->f)) {
            for (u = i + turn; u < length - turn; u++) {
              fm = sc->f(i, length, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
              if ((fm != INF) && (fml_tmp[u + 1] != INF))
                fml_tmp[u + 1] += fm;
            }
          }

          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          if (scs) {
            for (u = i + turn; u < length - turn; u++) {
              fm = 0;
              for (s = 0; s < n_seq; s++)
                if ((scs[s]) && (scs[s]->f))
                  fm += scs[s]->f(i, length, u, u + 1, VRNA_DECOMP_ML_ML_ML, scs[s]->data);

              if ((fm != INF) && (fml_tmp[u + 1] != INF))
                fml_tmp[u + 1] += fm;
            }
          }

          break;
      }
    }

    /* actual decomposition */
    for (u = i + turn; u < length - turn; u++) {
      fm = my_fML[indx[u] + i];
      if ((fm != INF) && (fml_tmp[u + 1] != INF)) {
        fm      += fml_tmp[u + 1];
        fM2[i]  = MIN2(fM2[i], fm);
      }
    }
  }

  if ((fml_tmp - indx[length]) != my_fML)
    free(fml_tmp);

  /*
   *  Now, process all exterior multibranch loop configurations
   */
  int *fm2_tmp = fM2;

  if (hc->f) {
    fm2_tmp = vrna_alloc(sizeof(int) * (length + 1));

    /* copy data */
    for (i = turn + 1; i < length - 2 * turn; i++)
      fm2_tmp[i + 1] = fM2[i + 1];

    /* apply hard constraints */
    for (i = turn + 1; i < length - 2 * turn; i++)
      if (!hc->f(1, length, i, i + 1, VRNA_DECOMP_ML_ML_ML, hc->data))
        fm2_tmp[i + 1] = INF;
  }

  if ((fc->type == VRNA_FC_TYPE_SINGLE) && (sc) && (sc->f)) {
    if (fm2_tmp == fM2) {
      fm2_tmp = vrna_alloc(sizeof(int) * (length + 1));

      /* copy data */
      for (i = turn + 1; i < length - 2 * turn; i++)
        fm2_tmp[i + 1] = fM2[i + 1];
    }

    /* add soft constraints */
    for (i = turn + 1; i < length - 2 * turn; i++) {
      fm = sc->f(1, length, i, i + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
      if ((fm != INF) && (fm2_tmp[i + 1] != INF))
        fm2_tmp[i + 1] += fm;
    }
  }

  if ((fc->type == VRNA_FC_TYPE_COMPARATIVE) && (scs)) {
    if (fm2_tmp == fM2) {
      fm2_tmp = vrna_alloc(sizeof(int) * (length + 1));

      /* copy data */
      for (i = turn + 1; i < length - 2 * turn; i++)
        fm2_tmp[i + 1] = fM2[i + 1];
    }

    /* add soft constraints */
    for (i = turn + 1; i < length - 2 * turn; i++) {
      fm = 0;
      for (s = 0; s < n_seq; s++)
        if ((scs[s]) && (scs[s]->f))
          fm += scs[s]->f(1, length, i, i + 1, VRNA_DECOMP_ML_ML_ML, scs[s]->data);

      if ((fm != INF) && (fm2_tmp[i + 1] != INF))
        fm2_tmp[i + 1] += fm;
    }
  }

  /* actual decomposition */
  for (i = turn + 1; i < length - 2 * turn; i++) {
    if ((my_fML[indx[i] + 1] != INF) && (fm2_tmp[i + 1] != INF)) {
      fm = my_fML[indx[i] + 1] +
           fm2_tmp[i + 1];

      if (fm < FcM) {
        FcM = fm;
        Mi  = i;
      }
    }
  }

  if (fm2_tmp != fM2)
    free(fm2_tmp);

  if (FcM != INF) {
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        FcM += P->MLclosing;
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        FcM += n_seq * P->MLclosing;
        break;
    }
  }

  Fc = MIN2(Fc, FcM);

  /* add multibranch loop configurations for odd dangle models */
  if ((dangle_model == 1) || (dangle_model == 3)) {
    fM_d5 = (int *)vrna_alloc(sizeof(int) * (length + 2));
    fM_d3 = (int *)vrna_alloc(sizeof(int) * (length + 2));

    for (i = turn + 1; i < length - turn; i++)
      fM_d5[i] = INF;

    for (i = turn + 1; i < length - turn; i++)
      fM_d3[i] = INF;

    if (hc->up_ml[1]) {
      fill_fM_d3(fc, fM_d3);

      /*
       **********************************************************************
       *  1. test ML loops with closing pair (i + 1, length), 3' dangle pos 1
       **********************************************************************
       */
      int *c_tmp, *c_tmp2;

      c_tmp   = vrna_alloc(sizeof(int) * (length + 2));
      c_tmp2  = my_c + indx[length];

      /* copy contributions for enclosed structure part */
      for (i = 2 * turn + 1; i < length - turn; i++)
        c_tmp[i + 1] = c_tmp2[i + 1];

      /* obey hard constraints */
      if (hc->f) {
        for (i = 2 * turn + 1; i < length - turn; i++)
          if (!hc->f(i + 1, length, 2, i, VRNA_DECOMP_PAIR_ML_EXT, hc->data))
            c_tmp[i + 1] = INF;
      }

      /* add soft constraints */
      if ((fc->type == VRNA_FC_TYPE_SINGLE) && (sc) && (sc->f)) {
        for (i = 2 * turn + 1; i < length - turn; i++) {
          tmp = sc->f(i + 1, length, 2, i,
                      VRNA_DECOMP_PAIR_ML_EXT,
                      sc->data);
          if (tmp != INF) {
            if (c_tmp[i + 1] != INF)
              c_tmp[i + 1] += tmp;
          } else {
            c_tmp[i + 1] = INF;
          }
        }
      }

      if ((fc->type == VRNA_FC_TYPE_COMPARATIVE) && (scs)) {
        for (i = 2 * turn + 1; i < length - turn; i++) {
          tmp = 0;
          for (s = 0; s < n_seq; s++)
            if ((scs[s]) && (scs[s]->f))
              tmp += scs[s]->f(i + 1, length, 2, i,
                               VRNA_DECOMP_PAIR_ML_EXT,
                               scs[s]->data);

          if (tmp != INF) {
            if (c_tmp[i + 1] != INF)
              c_tmp[i + 1] += tmp;
          } else {
            c_tmp[i + 1] = INF;
          }
        }
      }

      /* add contributions for enclosing pair */
      for (i = 2 * turn + 1; i < length - turn; i++) {
        if (c_tmp[i + 1] != INF) {
          /* obey internal hard constraints */
          if (hard_constraints[length * length + i + 1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
            tmp = 0;
            switch (fc->type) {
              case VRNA_FC_TYPE_SINGLE:
                type  = vrna_get_ptype(indx[length] + i + 1, ptype);
                tmp   = E_MLstem(type, -1, S1[1], P) +
                        P->MLclosing;
                break;

              case VRNA_FC_TYPE_COMPARATIVE:
                tmp = P->MLclosing * n_seq;
                for (s = 0; s < n_seq; s++) {
                  type  = vrna_get_ptype_md(SS[s][i + 1], SS[s][length], md);
                  tmp   += E_MLstem(type, -1, S3[s][length], P);
                }
                break;
            }

            c_tmp[i + 1] += tmp;
          } else {
            c_tmp[i + 1] = INF;
          }
        }
      }

      /* actual decomposition */
      for (i = 2 * turn + 1; i < length - turn; i++) {
        fm = fM_d3[i];
        if ((fm != INF) && (c_tmp[i + 1] != INF)) {
          fm += c_tmp[i + 1];
          if (fm < FcMd3) {
            FcMd3 = fm;
            Md3i  = i;
          }
        }
      }

      /*
       **********************************************************************
       *  2. test ML loops with closing pair (i + 1, length), 5' dangle
       *  pos i, 3' dangle pos 1
       **********************************************************************
       */

      /* copy contributions for enclosed structure part */
      for (i = 2 * turn + 1; i < length - turn; i++)
        c_tmp[i + 1] = c_tmp2[i + 1];


      /* obey hard constraints */
      if (hc->f) {
        for (i = 2 * turn + 1; i < length - turn; i++)
          if (!hc->f(i + 1, length, 2, i - 1, VRNA_DECOMP_PAIR_ML_EXT, hc->data))
            c_tmp[i + 1] = INF;
      }

      /* add soft constraints (static and user-defined) */
      if ((fc->type == VRNA_FC_TYPE_SINGLE) && (sc)) {
        if (sc->energy_up) {
          for (i = 2 * turn + 1; i < length - turn; i++)
            if (c_tmp[i + 1] != INF)
              c_tmp[i + 1] += sc->energy_up[i][1];
        }

        if (sc->f) {
          for (i = 2 * turn + 1; i < length - turn; i++) {
            tmp = sc->f(i + 1, length, 2, i - 1,
                        VRNA_DECOMP_PAIR_ML_EXT,
                        sc->data);
            if (tmp == INF)
              c_tmp[i + 1] = INF;
            else if (c_tmp[i + 1] != INF)
              c_tmp[i + 1] += tmp;
          }
        }
      }

      if ((fc->type == VRNA_FC_TYPE_COMPARATIVE) && (scs)) {
        for (i = 2 * turn + 1; i < length - turn; i++) {
          if (c_tmp[i + 1] != INF) {
            for (s = 0; s < n_seq; s++) {
              if (scs[s]) {
                if (scs[s]->energy_up)
                  c_tmp[i + 1] += scs[s]->energy_up[a2s[s][i]][1];

                if (scs[s]->f) {
                  c_tmp[i + 1] += scs[s]->f(i + 1, length, 2, i - 1,
                                            VRNA_DECOMP_PAIR_ML_EXT,
                                            scs[s]->data);
                }
              }
            }
          }
        }
      }

      /* add contributions for enclosing pair */
      for (i = 2 * turn + 1; i < length - turn; i++) {
        if (c_tmp[i + 1] != INF) {
          /* obey internal hard constraints */
          if ((hard_constraints[length * length + i + 1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) &&
              (hc->up_ml[i])) {
            tmp = 0;
            switch (fc->type) {
              case VRNA_FC_TYPE_SINGLE:
                type  = vrna_get_ptype(indx[length] + i + 1, ptype);
                tmp   = E_MLstem(type, S1[i], S1[1], P) +
                        P->MLclosing;
                break;

              case VRNA_FC_TYPE_COMPARATIVE:
                tmp = P->MLclosing * n_seq;
                for (s = 0; s < n_seq; s++) {
                  type  = vrna_get_ptype_md(SS[s][i + 1], SS[s][length], md);
                  tmp   += E_MLstem(type, S5[s][i + 1], S3[s][length], P);
                }
                break;
            }

            c_tmp[i + 1] += tmp;
          } else {
            c_tmp[i + 1] = INF;
          }
        }
      }

      /* actual decomposition */
      for (i = 2 * turn + 1; i < length - turn; i++) {
        fm = fM_d3[i - 1];
        if ((fm != INF) && (c_tmp[i + 1] != INF)) {
          fm += c_tmp[i + 1];
          if (fm < FcMd3) {
            FcMd3 = fm;
            Md3i  = -i;
          }
        }
      }

      free(c_tmp);
    }

    if (hc->up_ml[length]) {
      fill_fM_d5(fc, fM_d5);

      /*
       **********************************************************************
       * 1. test ML loops with closing pair (1, i), 5' dangle pos n
       **********************************************************************
       */
      int *fmd5_tmp = vrna_alloc(sizeof(int) * (length + 2));

      /* copy contributions */
      for (i = turn + 1; i < length - turn; i++)
        fmd5_tmp[i + 1] = fM_d5[i + 1];

      /* obey hard constraints */
      if (hc->f) {
        for (i = turn + 1; i < length - turn; i++)
          if (!hc->f(1, i, i + 1, length - 1, VRNA_DECOMP_PAIR_ML_EXT, hc->data))
            fmd5_tmp[i + 1] = INF;
      }

      /* add soft constraints */
      if ((fc->type == VRNA_FC_TYPE_SINGLE) && (sc) && (sc->f)) {
        for (i = turn + 1; i < length - turn; i++) {
          fm = sc->f(1, i, i + 1, length - 1,
                     VRNA_DECOMP_PAIR_ML_EXT,
                     sc->data);
          if ((fm != INF) && (fmd5_tmp[i + 1] != INF))
            fmd5_tmp += fm;
        }
      }

      if ((fc->type == VRNA_FC_TYPE_COMPARATIVE) && (scs)) {
        for (i = turn + 1; i < length - turn; i++) {
          fm = 0;
          for (s = 0; s < n_seq; s++)
            if ((scs[s]) && (scs[s]->f))
              fm += scs[s]->f(1, i, i + 1, length - 1,
                              VRNA_DECOMP_PAIR_ML_EXT,
                              scs[s]->data);

          if ((fm != INF) && (fmd5_tmp[i + 1] != INF))
            fmd5_tmp += fm;
        }
      }

      /* add contributions for enclosing pair */
      for (i = turn + 1; i < length - turn; i++) {
        if (fmd5_tmp[i + 1] != INF) {
          if (hard_constraints[length + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
            tmp = 0;
            switch (fc->type) {
              case VRNA_FC_TYPE_SINGLE:
                type  = vrna_get_ptype(indx[i] + 1, ptype);
                tmp   = E_MLstem(type, S1[length], -1, P) +
                        P->MLclosing;
                break;

              case VRNA_FC_TYPE_COMPARATIVE:
                tmp = P->MLclosing * n_seq;
                for (s = 0; s < n_seq; s++) {
                  type  = vrna_get_ptype_md(SS[s][1], SS[s][i], md);
                  tmp   += E_MLstem(type, S5[s][1], -1, P);
                }
                break;
            }

            fmd5_tmp[i + 1] += tmp;
          } else {
            fmd5_tmp[i + 1] = INF;
          }
        }
      }
      /* actual decomposition */
      for (i = turn + 1; i < length - turn; i++) {
        fm = fmd5_tmp[i + 1];
        if ((fm != INF) && (my_c[indx[i] + 1] != INF)) {
          fm += my_c[indx[i] + 1];
          if (fm < FcMd5) {
            FcMd5 = fm;
            Md5i  = i;
          }
        }
      }

      /*
       **********************************************************************
       * 2. test ML loops with closing pair (1, i), 5' dangle pos n, 3' dangle
       * pos i + 1
       **********************************************************************
       */

      /* copy contributions */
      for (i = turn + 1; i < length - turn; i++)
        fmd5_tmp[i + 2] = fM_d5[i + 2];

      /* obey hard constraints */
      if (hc->f) {
        for (i = turn + 1; i < length - turn; i++)
          if (!hc->f(1, i, i + 2, length - 1, VRNA_DECOMP_PAIR_ML_EXT, hc->data))
            fmd5_tmp[i + 2] = INF;
      }

      /* add soft constraints (static and user-defined) */
      if ((fc->type == VRNA_FC_TYPE_SINGLE) && (sc)) {
        if (sc->energy_up) {
          for (i = turn + 1; i < length - turn; i++)
            if (fmd5_tmp[i + 2] != INF)
              fmd5_tmp[i + 2] += sc->energy_up[i + 1][1];
        }

        if (sc->f) {
          for (i = turn + 1; i < length - turn; i++) {
            tmp = sc->f(1, i, i + 2, length - 1,
                        VRNA_DECOMP_PAIR_ML_EXT,
                        sc->data);
            if (tmp == INF)
              fmd5_tmp[i + 2] = INF;
            else if (fmd5_tmp[i + 2] != INF)
              fmd5_tmp[i + 2] += tmp;
          }
        }
      }

      if ((fc->type == VRNA_FC_TYPE_COMPARATIVE) && (scs)) {
        for (i = turn + 1; i < length - turn; i++) {
          if (fmd5_tmp[i + 2] != INF) {
            for (s = 0; s < n_seq; s++) {
              if (scs[s]) {
                if (scs[s]->energy_up)
                  fmd5_tmp[i + 2] += scs[s]->energy_up[a2s[s][i + 1]][1];

                if (scs[s]->f) {
                  fmd5_tmp[i + 2] += scs[s]->f(1, i, i + 2, length - 1,
                                               VRNA_DECOMP_PAIR_ML_EXT,
                                               scs[s]->data);
                }
              }
            }
          }
        }
      }

      /* add contributions for enclosing pair */
      for (i = turn + 1; i < length - turn; i++) {
        if (fmd5_tmp[i + 2] != INF) {
          /* obey internal hard constraints */
          if ((hard_constraints[length + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) &&
              (hc->up_ml[i + 1])) {
            tmp = 0;
            switch (fc->type) {
              case VRNA_FC_TYPE_SINGLE:
                type  = vrna_get_ptype(indx[i] + 1, ptype);
                tmp   = E_MLstem(type, S1[length], S1[i + 1], P) +
                        P->MLclosing;
                break;

              case VRNA_FC_TYPE_COMPARATIVE:
                tmp = P->MLclosing * n_seq;
                for (s = 0; s < n_seq; s++) {
                  type  = vrna_get_ptype_md(SS[s][1], SS[s][i], md);
                  tmp   += E_MLstem(type, S5[s][1], S3[s][i], P);
                }
                break;
            }

            fmd5_tmp[i + 2] += tmp;
          } else {
            fmd5_tmp[i + 2] = INF;
          }
        }
      }

      /* actual decomposition */
      for (i = turn + 1; i < length - turn; i++) {
        fm = fmd5_tmp[i + 2];
        if ((fm != INF) && (my_c[indx[i] + 1] != INF)) {
          fm += my_c[indx[i] + 1];

          if (fm < FcMd5) {
            FcMd5 = fm;
            Md5i  = -i;
          }
        }
      }

      free(fmd5_tmp);
    }

    if (FcMd5 < MIN2(Fc, FcMd3)) {
      int real_i, sc_en = 0;

      /* looks like we have to do this ... */
      bt_stack[++(*bt)].i = 1;
      bt_stack[(*bt)].j   = (Md5i > 0) ? Md5i : -Md5i;
      bt_stack[(*bt)].ml  = 2;
      i                   = (Md5i > 0) ? Md5i + 1 : -Md5i + 2; /* let's backtrack fm_d5[Md5i+1] */
      real_i              = (Md5i > 0) ? i : i - 1;

      if ((fc->type == VRNA_FC_TYPE_SINGLE) && (sc) && (sc->energy_up)) {
        sc_en += sc->energy_up[length][1];
      } else if ((fc->type == VRNA_FC_TYPE_COMPARATIVE) && (scs)) {
        for (s = 0; s < n_seq; s++)
          if ((scs[s]) && (scs[s]->energy_up))
            sc_en += scs[s]->energy_up[a2s[s][length]][1];
      }

      for (u = i + turn; u < length - turn; u++) {
        fm = my_fML[indx[u] + i] +
             my_fML[indx[length - 1] + u + 1] +
             sc_en;

        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            if (sc) {
              if (sc->energy_up)
                fm += sc->energy_up[real_i][i - real_i];

              if (sc->f) {
                fm += sc->f(real_i, length, i, length - 1,
                            VRNA_DECOMP_ML_ML,
                            sc->data) +
                      sc->f(i, length - 1, u, u + 1,
                            VRNA_DECOMP_ML_ML_ML,
                            sc->data);
              }
            }

            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            if (scs) {
              for (s = 0; s < n_seq; s++) {
                if (scs[s]) {
                  if (scs[s]->energy_up)
                    fm += scs[s]->energy_up[a2s[s][real_i]][i - real_i];

                  if (scs[s]->f) {
                    fm += scs[s]->f(real_i, length, i, length - 1,
                                    VRNA_DECOMP_ML_ML,
                                    scs[s]->data) +
                          scs[s]->f(i, length - 1, u, u + 1,
                                    VRNA_DECOMP_ML_ML_ML,
                                    scs[s]->data);
                  }
                }
              }
            }

            break;
        }

        if (fM_d5[i] == fm) {
          bt_stack[++(*bt)].i = i;
          bt_stack[(*bt)].j   = u;
          bt_stack[(*bt)].ml  = 1;
          bt_stack[++(*bt)].i = u + 1;
          bt_stack[(*bt)].j   = length - 1;
          bt_stack[(*bt)].ml  = 1;
          break;
        }
      }
      Fc = FcMd5;
    } else if (FcMd3 < Fc) {
      int real_i, sc_en = 0;
      /* here we go again... */
      bt_stack[++(*bt)].i = (Md3i > 0) ? Md3i + 1 : -Md3i + 1;
      bt_stack[(*bt)].j   = length;
      bt_stack[(*bt)].ml  = 2;
      i                   = (Md3i > 0) ? Md3i : -Md3i - 1; /* let's backtrack fm_d3[Md3i] */
      real_i              = (Md3i > 0) ? i : i + 1;

      if ((fc->type == VRNA_FC_TYPE_SINGLE) && (sc) && (sc->energy_up)) {
        sc_en += sc->energy_up[1][1];
      } else if ((fc->type == VRNA_FC_TYPE_COMPARATIVE) && (scs)) {
        for (s = 0; s < n_seq; s++)
          if ((scs[s]) && (scs[s]->energy_up))
            sc_en += scs[s]->energy_up[a2s[s][1]][1];
      }

      for (u = 2 + turn; u < i - turn; u++) {
        fm = my_fML[indx[u] + 2] +
             my_fML[indx[i] + u + 1] +
             sc_en;

        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            if (sc) {
              if (sc->energy_up)
                fm += sc->energy_up[real_i][real_i - i];

              if (sc->f) {
                fm += sc->f(1, real_i, 2, i,
                            VRNA_DECOMP_ML_ML,
                            sc->data) +
                      sc->f(2, i, u, u + 1,
                            VRNA_DECOMP_ML_ML_ML,
                            sc->data);
              }
            }

            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            if (scs) {
              for (s = 0; s < n_seq; s++) {
                if (scs[s]) {
                  if (scs[s]->energy_up)
                    fm += scs[s]->energy_up[a2s[s][real_i]][real_i - i];

                  if (scs[s]->f) {
                    fm += scs[s]->f(1, real_i, 2, i,
                                    VRNA_DECOMP_ML_ML,
                                    scs[s]->data) +
                          scs[s]->f(2, i, u, u + 1,
                                    VRNA_DECOMP_ML_ML_ML,
                                    scs[s]->data);
                  }
                }
              }
            }

            break;
        }

        if (fM_d3[i] == fm) {
          bt_stack[++(*bt)].i = 2;
          bt_stack[(*bt)].j   = u;
          bt_stack[(*bt)].ml  = 1;
          bt_stack[++(*bt)].i = u + 1;
          bt_stack[(*bt)].j   = i;
          bt_stack[(*bt)].ml  = 1;
          break;
        }
      }
      Fc = FcMd3;
    }

    free(fM_d3);
    free(fM_d5);
  }

  if (Fc < INF) {
    if (FcH == Fc) {
      bt_stack[++(*bt)].i = Hi;
      bt_stack[(*bt)].j   = Hj;
      bt_stack[(*bt)].ml  = 2;
    } else if (FcI == Fc) {
      bt_stack[++(*bt)].i = Ii;
      bt_stack[(*bt)].j   = Ij;
      bt_stack[(*bt)].ml  = 2;
      bt_stack[++(*bt)].i = Ip;
      bt_stack[(*bt)].j   = Iq;
      bt_stack[(*bt)].ml  = 2;
    } else if (FcM == Fc) {
      /* grumpf we found a Multiloop */
      int eee;
      /* backtrack in fM2 */
      fm = fM2[Mi + 1];
      for (u = Mi + turn + 1; u < length - turn; u++) {
        eee = my_fML[indx[u] + Mi + 1] +
              my_fML[indx[length] + u + 1];

        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            if (sc) {
              if (sc->f)
                eee += sc->f(Mi + 1, length, u, u + 1,
                             VRNA_DECOMP_ML_ML_ML,
                             sc->data);
            }

            break;

          case VRNA_FC_TYPE_COMPARATIVE:
            if (scs) {
              for (s = 0; s < n_seq; s++)
                if ((scs[s]) && (scs[s]->f))
                  eee += scs[s]->f(Mi + 1, length, u, u + 1,
                                   VRNA_DECOMP_ML_ML_ML,
                                   scs[s]->data);
            }

            break;
        }

        if (fm == eee) {
          bt_stack[++(*bt)].i = Mi + 1;
          bt_stack[(*bt)].j   = u;
          bt_stack[(*bt)].ml  = 1;
          bt_stack[++(*bt)].i = u + 1;
          bt_stack[(*bt)].j   = length;
          bt_stack[(*bt)].ml  = 1;
          break;
        }
      }
      bt_stack[++(*bt)].i = 1;
      bt_stack[(*bt)].j   = Mi;
      bt_stack[(*bt)].ml  = 1;
    } else if (Fc == FcO) {
      /* unstructured */
      bt_stack[++(*bt)].i = 1;
      bt_stack[(*bt)].j   = 1;
      bt_stack[(*bt)].ml  = 0;
    }
  } else {
    /* forbidden, i.e. no configuration fulfills constraints */
  }

  fc->matrices->FcH = FcH;
  fc->matrices->FcI = FcI;
  fc->matrices->FcM = FcM;
  fc->matrices->Fc  = Fc;
  return Fc;
}


/*
 **************************
 *  Construct fM_d5 array
 **************************
 */
PRIVATE INLINE void
fill_fM_d5(vrna_fold_compound_t *fc,
           int                  *fM_d5)
{
  unsigned int  s, n_seq, **a2s;
  int           *fm_tmp, *fm_tmp2, sc_base, fm, tmp, i, u,
                length, turn, *my_fML, *indx;
  vrna_md_t     *md;
  vrna_param_t  *P;
  vrna_hc_t     *hc;
  vrna_sc_t     *sc, **scs;

  n_seq   = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  length  = fc->length;
  a2s     = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->a2s;
  P       = fc->params;
  md      = &(P->model_details);
  my_fML  = fc->matrices->fML;
  hc      = fc->hc;
  sc      = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sc : NULL;
  scs     = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->scs;
  indx    = fc->jindx;
  turn    = md->min_loop_size;

  fm_tmp2 = vrna_alloc(sizeof(int) * (length + 2));
  sc_base = 0;

  if ((fc->type == VRNA_FC_TYPE_SINGLE) && (sc) && (sc->energy_up))
    sc_base += sc->energy_up[length][1];
  else if ((fc->type == VRNA_FC_TYPE_COMPARATIVE) && (scs))
    for (s = 0; s < n_seq; s++)
      if ((scs[s]) && (scs[s]->energy_up))
        sc_base += scs[s]->energy_up[a2s[s][length]][1];

  for (i = turn + 1; i < length - turn; i++) {
    fm_tmp = my_fML + indx[length - 1];

    if (sc_base != 0) {
      fm_tmp = fm_tmp2;

      for (u = 2 + turn; u < i - turn; u++)
        fm_tmp[u + 1] = my_fML[indx[length - 1] + u + 1] +
                        sc_base;
    }

    if (hc->f) {
      if (!hc->f(i, length, i, length - 1, VRNA_DECOMP_ML_ML, hc->data))
        continue;

      if (fm_tmp != fm_tmp2) {
        fm_tmp = fm_tmp2;

        for (u = 2 + turn; u < i - turn; u++)
          fm_tmp[u + 1] = my_fML[indx[length - 1] + u + 1];
      }

      for (u = 2 + turn; u < i - turn; u++)
        if (!hc->f(i, length - 1, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data))
          fm_tmp[u + 1] = INF;
    }

    if ((fc->type == VRNA_FC_TYPE_SINGLE) && (sc) && (sc->f)) {
      if (fm_tmp != fm_tmp2) {
        fm_tmp = fm_tmp2;

        for (u = 2 + turn; u < i - turn; u++)
          fm_tmp[u + 1] = my_fML[indx[length - 1] + u + 1];
      }

      fm = sc->f(i, length, i, length - 1,
                 VRNA_DECOMP_ML_ML,
                 sc->data);

      if (fm != INF) {
        for (u = 2 + turn; u < i - turn; u++) {
          if (fm_tmp[u + 1] != INF) {
            tmp = sc->f(i, length - 1, u, u + 1,
                        VRNA_DECOMP_ML_ML_ML,
                        sc->data);
            if (tmp != INF)
              tmp += fm;

            fm_tmp[u + 1] += tmp;
          }
        }
      } else {
        for (u = 2 + turn; u < i - turn; u++)
          fm_tmp[u + 1] = INF;
      }
    }

    if ((fc->type == VRNA_FC_TYPE_COMPARATIVE) && (scs)) {
      if (fm_tmp != fm_tmp2) {
        fm_tmp = fm_tmp2;

        for (u = 2 + turn; u < i - turn; u++)
          fm_tmp[u + 1] = my_fML[indx[length - 1] + u + 1];
      }

      fm = 0;
      for (s = 0; s < n_seq; s++)
        if ((scs[s]) && (scs[s]->f))
          fm += scs[s]->f(i, length, i, length - 1,
                          VRNA_DECOMP_ML_ML,
                          scs[s]->data);

      if (fm != INF) {
        for (u = 2 + turn; u < i - turn; u++) {
          if (fm_tmp[u + 1] != INF) {
            tmp = 0;
            for (s = 0; s < n_seq; s++)
              if ((scs[s]) && (scs[s]->f))
                tmp += scs[s]->f(i, length - 1, u, u + 1,
                                 VRNA_DECOMP_ML_ML_ML,
                                 scs[s]->data);

            tmp           += fm;
            fm_tmp[u + 1] += tmp;
          }
        }
      } else {
        for (u = 2 + turn; u < i - turn; u++)
          fm_tmp[u + 1] = INF;
      }
    }

    /* actually compute the entries for fM_d3 array */
    for (u = i + turn; u < length - turn; u++) {
      fm = my_fML[indx[u] + i];

      /* skip configurations that violate (hard) constraints */
      if ((fm == INF) || (fm_tmp[u + 1] == INF))
        continue;

      fm += fm_tmp[u + 1];

      fM_d5[i] = MIN2(fM_d5[i], fm);
    }
  }

  free(fm_tmp2);
}


/*
 **************************
 *  Construct fM_d3 array
 **************************
 */
PRIVATE INLINE void
fill_fM_d3(vrna_fold_compound_t *fc,
           int                  *fM_d3)
{
  unsigned int  s, n_seq, **a2s;
  int           *fm_tmp, *fm_tmp2, sc_base, fm, tmp, i, u,
                length, turn, *my_fML, *indx;
  vrna_md_t     *md;
  vrna_param_t  *P;
  vrna_hc_t     *hc;
  vrna_sc_t     *sc, **scs;

  n_seq   = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  length  = fc->length;
  a2s     = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->a2s;
  P       = fc->params;
  md      = &(P->model_details);
  my_fML  = fc->matrices->fML;
  hc      = fc->hc;
  sc      = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sc : NULL;
  scs     = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->scs;
  indx    = fc->jindx;
  turn    = md->min_loop_size;

  fm_tmp2 = vrna_alloc(sizeof(int) * (length + 2));
  sc_base = 0;

  if ((fc->type == VRNA_FC_TYPE_SINGLE) && (sc) && (sc->energy_up))
    sc_base += sc->energy_up[1][1];
  else if ((fc->type == VRNA_FC_TYPE_COMPARATIVE) && (scs))
    for (s = 0; s < n_seq; s++)
      if ((scs[s]) && (scs[s]->energy_up))
        sc_base += scs[s]->energy_up[a2s[s][1]][1];

  for (i = turn + 1; i < length - turn; i++) {
    fm_tmp = my_fML + indx[i];

    if (sc_base != 0) {
      fm_tmp = fm_tmp2;

      for (u = 2 + turn; u < i - turn; u++)
        fm_tmp[u + 1] = my_fML[indx[i] + u + 1] +
                        sc_base;
    }

    if (hc->f) {
      if (!hc->f(1, i, 2, i, VRNA_DECOMP_ML_ML, hc->data))
        continue;

      if (fm_tmp != fm_tmp2) {
        fm_tmp = fm_tmp2;

        for (u = 2 + turn; u < i - turn; u++)
          fm_tmp[u + 1] = my_fML[indx[i] + u + 1];
      }

      for (u = 2 + turn; u < i - turn; u++)
        if (!hc->f(2, i, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data))
          fm_tmp[u + 1] = INF;
    }

    if ((fc->type == VRNA_FC_TYPE_SINGLE) && (sc) && (sc->f)) {
      if (fm_tmp != fm_tmp2) {
        fm_tmp = fm_tmp2;

        for (u = 2 + turn; u < i - turn; u++)
          fm_tmp[u + 1] = my_fML[indx[i] + u + 1];
      }

      fm = sc->f(1, i, 2, i, VRNA_DECOMP_ML_ML, sc->data);

      if (fm != INF) {
        for (u = 2 + turn; u < i - turn; u++) {
          if (fm_tmp[u + 1] != INF) {
            tmp = sc->f(2, i, u, u + 1,
                        VRNA_DECOMP_ML_ML_ML,
                        sc->data);
            if (tmp != INF) {
              tmp           += fm;
              fm_tmp[u + 1] += tmp;
            } else {
              fm_tmp[u + 1] = INF;
            }
          }
        }
      } else {
        for (u = 2 + turn; u < i - turn; u++)
          fm_tmp[u + 1] = INF;
      }
    }

    if ((fc->type == VRNA_FC_TYPE_COMPARATIVE) && (scs)) {
      if (fm_tmp != fm_tmp2) {
        fm_tmp = fm_tmp2;

        for (u = 2 + turn; u < i - turn; u++)
          fm_tmp[u + 1] = my_fML[indx[i] + u + 1];
      }

      fm = 0;

      for (s = 0; s < n_seq; s++)
        if ((scs[s]) && (scs[s]->f))
          fm += scs[s]->f(1, i, 2, i,
                          VRNA_DECOMP_ML_ML,
                          scs[s]->data);

      for (u = 2 + turn; u < i - turn; u++) {
        if (fm_tmp[u + 1] != INF) {
          tmp = fm;
          for (s = 0; s < n_seq; s++)
            if ((scs[s]) && (scs[s]->f))
              tmp += scs[s]->f(2, i, u, u + 1,
                               VRNA_DECOMP_ML_ML_ML,
                               scs[s]->data);

          fm_tmp[u + 1] += tmp;
        }
      }
    }

    /* actually compute the entries for fM_d3 array */
    for (u = 2 + turn; u < i - turn; u++) {
      fm = my_fML[indx[u] + 2];

      /* skip configurations that violate (hard) constraints */
      if ((fm == INF) || (fm_tmp[u + 1] == INF))
        continue;

      fm += fm_tmp[u + 1];

      fM_d3[i] = MIN2(fM_d3[i], fm);
    }
  }

  free(fm_tmp2);
}


/**
*** trace back through the "c", "f5" and "fML" arrays to get the
*** base pairing list. No search for equivalent structures is done.
*** This is fast, since only few structure elements are recalculated.
***
*** normally s=0.
*** If s>0 then s items have been already pushed onto the bt_stack
**/
PRIVATE int
backtrack(vrna_fold_compound_t  *fc,
          vrna_bp_stack_t       *bp_stack,
          sect                  bt_stack[],
          int                   s)
{
  char          backtrack_type;
  int           i, j, ij, k, length, b, *my_c, *indx, noLP, *pscore, ret;
  vrna_param_t  *P;

  ret             = 1;
  b               = 0;
  length          = fc->length;
  my_c            = fc->matrices->c;
  indx            = fc->jindx;
  P               = fc->params;
  noLP            = P->model_details.noLP;
  pscore          = fc->pscore;         /* covariance scores for comparative structure prediction */
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
        if (vrna_BT_ext_loop_f5(fc, &j, &p, &q, bp_stack, &b)) {
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
          vrna_message_warning("backtracking failed in f5, segment [%d,%d]\n", i, j);
          ret = 0;
          goto backtrack_exit;
        }
      }
      break;

      /* trace back in fML array */
      case 1:
      {
        int p, q, comp1, comp2;
        if (vrna_BT_mb_loop_split(fc, &i, &j, &p, &q, &comp1, &comp2, bp_stack, &b)) {
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
          ret = 0;
          goto backtrack_exit;
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
        ret = 0;
        goto backtrack_exit;
    }

repeat1:

    /*----- begin of "repeat:" -----*/
    ij = indx[j] + i;

    if (canonical)
      cij = my_c[ij];

    if (noLP) {
      if (vrna_BT_stack(fc, &i, &j, &cij, bp_stack, &b)) {
        canonical = 0;
        goto repeat1;
      }
    }

    canonical = 1;

    if (fc->type == VRNA_FC_TYPE_COMPARATIVE)
      cij += pscore[indx[j] + i];

    if (vrna_BT_hp_loop(fc, i, j, cij, bp_stack, &b))
      continue;

    if (vrna_BT_int_loop(fc, &i, &j, cij, bp_stack, &b)) {
      if (i < 0)
        continue;
      else
        goto repeat1;
    }

    /* (i.j) must close a multi-loop */
    int comp1, comp2;

    if (vrna_BT_mb_loop(fc, &i, &j, &k, cij, &comp1, &comp2)) {
      bt_stack[++s].i = i;
      bt_stack[s].j   = k;
      bt_stack[s].ml  = comp1;
      bt_stack[++s].i = k + 1;
      bt_stack[s].j   = j;
      bt_stack[s].ml  = comp2;
    } else {
      vrna_message_warning("backtracking failed in repeat, segment [%d,%d]\n", i, j);
      ret = 0;
      goto backtrack_exit;
    }

    /* end of repeat: --------------------------------------------------*/
  } /* end of infinite while loop */

backtrack_exit:

  bp_stack[0].i = b;    /* save the total number of base pairs */

  return ret;
}


PRIVATE INLINE int
decompose_pair(vrna_fold_compound_t *fc,
               int                  i,
               int                  j,
               struct aux_arrays    *aux)
{
  unsigned char hc_decompose;
  unsigned int  n;
  int           e, new_c, energy, stackEnergy, ij, dangle_model, noLP,
                *DMLi1, *DMLi2, *cc, *cc1;

  n             = fc->length;
  ij            = fc->jindx[j] + i;
  dangle_model  = fc->params->model_details.dangles;
  noLP          = fc->params->model_details.noLP;
  hc_decompose  = fc->hc->mx[n * i + j];
  DMLi1         = aux->DMLi1;
  DMLi2         = aux->DMLi2;
  cc            = aux->cc;
  cc1           = aux->cc1;
  e             = INF;

  /* do we evaluate this pair? */
  if (hc_decompose) {
    new_c = INF;

    /* check for hairpin loop */
    energy  = vrna_E_hp_loop(fc, i, j);
    new_c   = MIN2(new_c, energy);

    /* check for multibranch loops */
    energy  = vrna_E_mb_loop_fast(fc, i, j, DMLi1, DMLi2);
    new_c   = MIN2(new_c, energy);

    if (dangle_model == 3) {
      /* coaxial stacking */
      energy  = vrna_E_mb_loop_stack(fc, i, j);
      new_c   = MIN2(new_c, energy);
    }

    /* check for interior loops */
    energy  = vrna_E_int_loop(fc, i, j);
    new_c   = MIN2(new_c, energy);

    /* remember stack energy for --noLP option */
    if (noLP) {
      stackEnergy = vrna_E_stack(fc, i, j);
      new_c       = MIN2(new_c, cc1[j - 1] + stackEnergy);
      cc[j]       = new_c;
      if ((fc->type == VRNA_FC_TYPE_COMPARATIVE) &&
          (cc[j] != INF))
        cc[j] -= fc->pscore[ij];

      e = cc1[j - 1] + stackEnergy;
    } else {
      e = new_c;
    }

    /* finally, check for auxiliary grammar rule(s) */
    if ((fc->aux_grammar) && (fc->aux_grammar->cb_aux_c)) {
      energy  = fc->aux_grammar->cb_aux_c(fc, i, j, fc->aux_grammar->data);
      new_c   = MIN2(new_c, energy);
    }

    if ((fc->type == VRNA_FC_TYPE_COMPARATIVE) &&
        (e != INF))
      e -= fc->pscore[ij];
  } /* end >> if (pair) << */

  return e;
}


PRIVATE INLINE struct aux_arrays *
get_aux_arrays(unsigned int length)
{
  unsigned int      j;
  struct aux_arrays *aux = (struct aux_arrays *)vrna_alloc(sizeof(struct aux_arrays));

  aux->cc     = (int *)vrna_alloc(sizeof(int) * (length + 2));  /* auxilary arrays for canonical structures     */
  aux->cc1    = (int *)vrna_alloc(sizeof(int) * (length + 2));  /* auxilary arrays for canonical structures     */
  aux->Fmi    = (int *)vrna_alloc(sizeof(int) * (length + 1));  /* holds row i of fML (avoids jumps in memory)  */
  aux->DMLi   = (int *)vrna_alloc(sizeof(int) * (length + 1));  /* DMLi[j] holds  MIN(fML[i,k]+fML[k+1,j])      */
  aux->DMLi1  = (int *)vrna_alloc(sizeof(int) * (length + 1));  /*                MIN(fML[i+1,k]+fML[k+1,j])    */
  aux->DMLi2  = (int *)vrna_alloc(sizeof(int) * (length + 1));  /*                MIN(fML[i+2,k]+fML[k+1,j])    */

  /* prefill helper arrays */
  for (j = 0; j <= length; j++)
    aux->Fmi[j] = aux->DMLi[j] = aux->DMLi1[j] = aux->DMLi2[j] = INF;

  return aux;
}


PRIVATE INLINE void
rotate_aux_arrays(struct aux_arrays *aux,
                  unsigned int      length)
{
  unsigned int  j;
  int           *FF;

  FF          = aux->DMLi2;
  aux->DMLi2  = aux->DMLi1;
  aux->DMLi1  = aux->DMLi;
  aux->DMLi   = FF;
  FF          = aux->cc1;
  aux->cc1    = aux->cc;
  aux->cc     = FF;
  for (j = 1; j <= length; j++)
    aux->cc[j] = aux->Fmi[j] = aux->DMLi[j] = INF;
}


PRIVATE INLINE void
free_aux_arrays(struct aux_arrays *aux)
{
  free(aux->cc);
  free(aux->cc1);
  free(aux->Fmi);
  free(aux->DMLi);
  free(aux->DMLi1);
  free(aux->DMLi2);
  free(aux);
}
