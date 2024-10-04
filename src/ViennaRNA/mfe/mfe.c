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
#include "ViennaRNA/structures/pairtable.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/mfe/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/eval/hairpin.h"
#include "ViennaRNA/eval/internal.h"
#include "ViennaRNA/eval/multibranch.h"
#include "ViennaRNA/mfe/exterior.h"
#include "ViennaRNA/mfe/internal.h"
#include "ViennaRNA/mfe/multibranch.h"
#include "ViennaRNA/backtrack/global.h"
#include "ViennaRNA/backtrack/exterior.h"
#include "ViennaRNA/backtrack/hairpin.h"
#include "ViennaRNA/backtrack/internal.h"
#include "ViennaRNA/backtrack/multibranch.h"
#include "ViennaRNA/backtrack/gquad.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/mfe/global.h"

#include "ViennaRNA/intern/grammar_dat.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#define MAXSECTORS        500     /* dimension for a backtrack array */

#define   ADD_OR_INF(a, b)     (((a) != INF) && ((b) != INF) ?  (a) + (b) : INF)

#define M2_FORWARD

struct aux_arrays {
  int                   *cc;  /* auxilary arrays for canonical structures     */
  int                   *cc1; /* auxilary arrays for canonical structures     */
  vrna_mx_mfe_aux_ml_t  ml_helpers;
};

#include "ViennaRNA/constraints/exterior_hc.inc"
#include "ViennaRNA/constraints/exterior_sc.inc"
#include "ViennaRNA/constraints/multibranch_hc.inc"
#include "ViennaRNA/constraints/multibranch_sc.inc"

struct ms_helpers {
  vrna_hc_eval_f        evaluate;
  struct hc_ext_def_dat hc_dat_local;
  struct sc_f5_dat      sc_wrapper;
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
fill_arrays(vrna_fold_compound_t  *fc,
            struct ms_helpers     *ms_dat);


PRIVATE int
postprocess_circular(vrna_fold_compound_t *fc,
                     vrna_bts_t           bt_stack);


PRIVATE INLINE void
fill_fM_d5(vrna_fold_compound_t *fc,
           int                  *fM_d5);


PRIVATE INLINE void
fill_fM_d3(vrna_fold_compound_t *fc,
           int                  *fM_d3);


PRIVATE int
backtrack(vrna_fold_compound_t  *fc,
          vrna_bps_t            bp_stack,
          vrna_bts_t            bt_stack,
          struct ms_helpers     *ms_dat);


PRIVATE INLINE int
decompose_pair(vrna_fold_compound_t *fc,
               int                  i,
               int                  j,
               struct aux_arrays    *aux,
               struct ms_helpers    *ms_dat);


PRIVATE INLINE struct aux_arrays *
get_aux_arrays(unsigned int length);


PRIVATE INLINE void
rotate_aux_arrays(struct aux_arrays *aux,
                  unsigned int      length);


PRIVATE INLINE void
free_aux_arrays(struct aux_arrays *aux);


PRIVATE struct ms_helpers *
get_ms_helpers(vrna_fold_compound_t *fc);


PRIVATE void
free_ms_helpers(struct ms_helpers *ms_dat,
                size_t            strands);


PRIVATE void
update_fms5_arrays(vrna_fold_compound_t *fc,
                   int                  i,
                   struct ms_helpers    *ms_dat);


PRIVATE void
update_fms3_arrays(vrna_fold_compound_t *fc,
                   unsigned int         s,
                   struct ms_helpers    *ms_dat);


PRIVATE int
pair_multi_strand(vrna_fold_compound_t  *fc,
                  int                   i,
                  int                   j,
                  struct ms_helpers     *ms_dat);


PRIVATE int
BT_multi_strand(vrna_fold_compound_t  *fc,
                int                   *i,
                int                   *j,
                unsigned int          *sn1,
                unsigned int          *sn2,
                int                   en,
                struct ms_helpers     *ms_dat);


PRIVATE int
BT_fms5_split(vrna_fold_compound_t  *fc,
              unsigned int          strand,
              int                   *i,
              int                   *k,
              int                   *l,
              struct ms_helpers     *ms_dat);


PRIVATE int
BT_fms3_split(vrna_fold_compound_t  *fc,
              unsigned int          strand,
              int                   *j,
              int                   *k,
              int                   *l,
              struct ms_helpers     *ms_dat);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC float
vrna_mfe(vrna_fold_compound_t *fc,
         char                 *structure)
{
  char              *ss;
  unsigned int      length;
  int               energy;
  float             mfe;
  vrna_bts_t        bt_stack; /* stack of partial structures for backtracking */
  vrna_bps_t        bp;
  struct ms_helpers *ms_dat;

  bt_stack  = vrna_bts_init(MAXSECTORS);
  bp        = NULL;

  mfe = (float)(INF / 100.);

  if (fc) {
    length  = fc->length;
    ms_dat  = NULL;

    if (!vrna_fold_compound_prepare(fc, VRNA_OPTION_MFE)) {
      vrna_log_warning("vrna_mfe@mfe.c: Failed to prepare vrna_fold_compound");
      vrna_bts_free(bt_stack);
      return mfe;
    }

    /* call user-defined recursion status callback function */
    if (fc->stat_cb)
      fc->stat_cb(fc, VRNA_STATUS_MFE_PRE, fc->auxdata);

    /* call user-defined grammar pre-condition callback function */
    if (fc->aux_grammar) {
      for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->cbs_status); i++)
        if (fc->aux_grammar->cbs_status[i])
          fc->aux_grammar->cbs_status[i](fc, VRNA_STATUS_MFE_PRE, fc->aux_grammar->datas[i]);
    }

    if (fc->strands > 1)
      ms_dat = get_ms_helpers(fc);

    energy = fill_arrays(fc, ms_dat);

    if (fc->params->model_details.circ)
      energy = postprocess_circular(fc, bt_stack);

    if (structure && fc->params->model_details.backtrack) {
      /* add a guess of how many G's may be involved in a G quadruplex */
      bp = vrna_bps_init(4 * (1 + length / 2));
      if (backtrack(fc, bp, bt_stack, ms_dat) != 0) {
        if ((fc->aux_grammar) &&
            (fc->aux_grammar->serialize_bp)) {
          ss = fc->aux_grammar->serialize_bp(fc, bp, fc->aux_grammar->serialize_bp_data);
        } else {
          ss = vrna_db_from_bps(bp, length);
        }

        strncpy(structure, ss, length + 1);
        free(ss);
      } else {
        memset(structure, '\0', sizeof(char) * (length + 1));
      }

      vrna_bps_free(bp);
    }

    /* call user-defined recursion status callback function */
    if (fc->stat_cb)
      fc->stat_cb(fc, VRNA_STATUS_MFE_POST, fc->auxdata);

    /* call user-defined grammar post-condition callback function */
    if (fc->aux_grammar) {
      for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->cbs_status); i++)
        if (fc->aux_grammar->cbs_status[i])
          fc->aux_grammar->cbs_status[i](fc, VRNA_STATUS_MFE_POST, fc->aux_grammar->datas[i]);
    }

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

    free_ms_helpers(ms_dat, fc->strands);
  }

  vrna_bts_free(bt_stack);

  return mfe;
}


PUBLIC int
vrna_backtrack_from_intervals(vrna_fold_compound_t  *fc,
                              vrna_bp_stack_t       *bp_stack,
                              sect                  bt_stack[],
                              int                   s)
{
  int i, ret = 0;

  if (fc) {
    vrna_bts_t bts;
    if (s > 0) {
      bts = vrna_bts_init((unsigned int)s);
      for (i = 0; i < s; i++)
        vrna_bts_push(bts,
                      (vrna_sect_t){
                        .i = bt_stack[i].i,
                        .j = bt_stack[i].j,
                        .ml = bt_stack[i].ml
                      });
    } else {
      bts = vrna_bts_init(0);
    }

    vrna_bps_t bps = vrna_bps_init(0);
    ret = backtrack(fc, bps, bts, NULL);

    /* copy bps elements to bp_stack?! */
    if (bp_stack) {
      unsigned int j = bp_stack[0].i;
      while (vrna_bps_size(bps) > 0) {
        vrna_bp_t bp = vrna_bps_pop(bps);
        bp_stack[++j].i = bp.i;
        bp_stack[j].j   = bp.j;
      }
      bp_stack[0].i = j;
    }

    vrna_bts_free(bts);
    vrna_bps_free(bps);
  }

  return ret;
}


PUBLIC float
vrna_backtrack5(vrna_fold_compound_t  *fc,
                unsigned int          length,
                char                  *structure)
{
  char        *ss;
  float       mfe;
  vrna_bts_t  bt_stack;     /* stack of partial structures for backtracking */
  vrna_bps_t  bp_stack;

  mfe = (float)(INF / 100.);

  if ((fc) && (structure) && (fc->matrices) && (fc->matrices->f5) &&
      (!fc->params->model_details.circ)) {
    memset(structure, '\0', sizeof(char) * (length + 1));

    if (length > fc->length)
      return mfe;

    /* add a guess of how many G's may be involved in a G quadruplex */
    bp_stack  = vrna_bps_init(4 * (1 + length / 2));
    bt_stack  = vrna_bts_init(MAXSECTORS);
    vrna_bts_push(bt_stack,
                  ((vrna_sect_t){
                    .i = 1,
                    .j = length,
                    .ml = VRNA_MX_FLAG_F5
                  }));

    if (backtrack(fc, bp_stack, bt_stack, NULL) != 0) {
      if ((fc->aux_grammar) &&
          (fc->aux_grammar->serialize_bp)) {
        ss = fc->aux_grammar->serialize_bp(fc,
                                           bp_stack,
                                           fc->aux_grammar->serialize_bp_data);
      } else {
        ss = vrna_db_from_bps(bp_stack, length);
      }

      strncpy(structure, ss, length + 1);
      free(ss);

      if (fc->type == VRNA_FC_TYPE_COMPARATIVE)
        mfe = (float)fc->matrices->f5[length] / (100. * (float)fc->n_seq);
      else
        mfe = (float)fc->matrices->f5[length] / 100.;
    }

    vrna_bts_free(bt_stack);
    vrna_bps_free(bp_stack);
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
fill_arrays(vrna_fold_compound_t  *fc,
            struct ms_helpers     *ms_dat)
{
  unsigned int      *sn;
  unsigned int      i, j, length, uniq_ML;
  int               ij, *indx, *f5, *c, *fML, *fM1, *fM2_real;
  vrna_param_t      *P;
  vrna_md_t         *md;
  vrna_mx_mfe_t     *matrices;
  vrna_ud_t         *domains_up;
  struct aux_arrays *helper_arrays;

  length      = (int)fc->length;
  indx        = fc->jindx;
  P           = fc->params;
  md          = &(P->model_details);
  uniq_ML     = md->uniq_ML;
  matrices    = fc->matrices;
  f5          = matrices->f5;
  c           = matrices->c;
  fML         = matrices->fML;
  fM1         = matrices->fM1;
  fM2_real    = matrices->fM2_real;
  domains_up  = fc->domains_up;
  sn          = fc->strand_number;

  /* allocate memory for all helper arrays */
  helper_arrays = get_aux_arrays(length);

  /* pre-processing ligand binding production rule(s) */
  if (domains_up && domains_up->prod_cb)
    domains_up->prod_cb(fc, domains_up->data);

  /* prefill matrices with init contributions */
  for (i = 1; i <= length; i++) {
    c[indx[i] + i] = fML[indx[i] + i] = INF;

    if (fM1)
      fM1[indx[i] + i] = INF;

    if (fM2_real)
      for (j = i; j <= length; j++)
        fM2_real[indx[j] + i] = INF;
  }

  /* start recursion */
  if (length <= ((fc->strands > 1) ? fc->strands : (unsigned int)md->min_loop_size)) {
    /* clean up memory */
    free_aux_arrays(helper_arrays);

    /* return free energy of unfolded chain */
    return 0;
  }

  for (i = length - 1; i >= 1; i--) {
    if ((fc->strands > 1) &&
        (sn[i] != sn[i + 1]))
      update_fms3_arrays(fc, sn[i + 1], ms_dat);

    for (j = i + 1; j <= length; j++) {
      ij = indx[j] + i;

      /* decompose subsegment [i, j] with pair (i, j) */
      c[ij] = decompose_pair(fc, i, j, helper_arrays, ms_dat);

      /* decompose subsegment [i, j] that is multibranch loop part with at least two branches */
      if (fM2_real)
        fM2_real[ij] = vrna_mfe_multibranch_m2_fast(fc,
                                                    i,
                                                    j,
                                                    helper_arrays->ml_helpers);

      /* decompose subsegment [i, j] that is multibranch loop part with at least one branch */
      fML[ij] = vrna_mfe_multibranch_stems_fast(fc, i, j, helper_arrays->ml_helpers);


      /* decompose subsegment [i, j] that is multibranch loop part with exactly one branch */
      if (fM1)
        fM1[ij] = vrna_mfe_multibranch_m1(fc, i, j);

      if (fc->aux_grammar)
        /* call auxiliary grammar rules */
        for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->aux); i++)
          if (fc->aux_grammar->aux[i].cb)
            (void)fc->aux_grammar->aux[i].cb(fc, i, j, fc->aux_grammar->aux[i].data);
    } /* end of j-loop */

    rotate_aux_arrays(helper_arrays, length);

    if (fc->strands > 1)
      update_fms5_arrays(fc, i, ms_dat);
  } /*
     * end of i-loop
     * calculate energies of 5' fragments
     */
  (void)vrna_mfe_exterior_f5(fc);

  /* clean up memory */
  free_aux_arrays(helper_arrays);

  return f5[length];
}


/* post-processing step for circular RNAs */
PRIVATE int
postprocess_circular(vrna_fold_compound_t *fc,
                     vrna_bts_t           bt_stack)
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
  unsigned char     *hard_constraints, eval;
  char              *ptype;
  short             *S, *S1, **SS, **S5, **S3;
  unsigned int      **a2s, u1, u2, us, us1, us2, s1, s2, p, q, si, sj,
                    Hgi, Hgj, Igi, Igj, Igp, Igq, Igg, Mgi, Mgj, Hi, Hj,
                    Ii, Ij, Ip, Iq, Mi, i, j, u, length, s, n_seq, turn,
                    dangle_model, with_gquad;
  int               *fM_d3, *fM_d5, Md3i, Md5i, FcMd3, FcMd5, FcH, FcI, FcM, Fc, ij,
                    new_c, fm, type, *my_c, *my_fML, *indx, FcO, tmp, FgH, FgI,
                    FgM, e, *fM2_real, *fM1_new;
  vrna_param_t      *P;
  vrna_md_t         *md;
  vrna_hc_t         *hc;
  vrna_sc_t         *sc, **scs;
  vrna_smx_csr(int) *c_gq;
  struct sc_mb_dat  sc_mb_wrapper;

  length            = fc->length;
  n_seq             = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  P                 = fc->params;
  md                = &(P->model_details);
  ptype             = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->ptype : NULL;
  indx              = fc->jindx;
  S                 = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding2 : NULL;
  S1                = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  SS                = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5                = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3                = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  a2s               = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->a2s;
  hc                = fc->hc;
  sc                = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sc : NULL;
  scs               = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->scs;
  dangle_model      = md->dangles;
  with_gquad        = md->gquad;
  turn              = md->min_loop_size;
  hard_constraints  = hc->mx;
  my_c              = fc->matrices->c;
  my_fML            = fc->matrices->fML;
  fM1_new           = fc->matrices->fM1_new;
  fM2_real          = fc->matrices->fM2_real;
  c_gq              = fc->matrices->c_gq;

  init_sc_mb(fc, &sc_mb_wrapper);

  Fc  = FcO = FcH = FcI = FcM = FcMd3 = FcMd5 = INF;
  Mi  = Md5i = Md3i = Iq = Ip = Ij = Ii = Hj = Hi = 0;

  /* explicit gquadruplex cases */
  FgH = FgI = FgM = INF;
  Hgi = Hgj = Igi = Igj = Igp = Igq = Igg = Mgi = Mgj = 0;

  /* unfolded state */
  eval = (hc->up_ext[1] >= length) ? 1 : 0;
  if (hc->f)
    eval = (hc->f(1, length, 1, length, VRNA_DECOMP_EXT_UP, hc->data)) ? eval : 0;

  if (eval) {
    Fc = 0;

    switch (fc->type) {
      case VRNA_FC_TYPE_COMPARATIVE:
        for (s = 0; s < n_seq; s++)
          Fc += vrna_E_exterior_loop(a2s[s][length], md);
        break;

      default:
        Fc += vrna_E_exterior_loop(length, md);
        break;
    }

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
            if (scs[s]) {
              if (scs[s]->energy_up)
                Fc += scs[s]->energy_up[1][a2s[s][length]];

              if (scs[s]->f)
                Fc += scs[s]->f(1, length, 1, length, VRNA_DECOMP_EXT_UP, scs[s]->data);
            }
        }

        break;
    }
    FcO = Fc;
  } else {
    Fc = INF;
  }

  /* fill fM1 */
  for (j = length; j > 0; j--)
    fM1_new[j] = INF;

  for (j = MIN2(turn + 2, VRNA_GQUAD_MIN_BOX_SIZE); j <= length; j++) {
    /* regular base pairs */
    for (u = j - turn - 1; u >= 1; u--) {
      eval = (hc->up_ml[1] >= (u - 1)) ?
             (hc->mx[length * j + u] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) :
             0;
      if ((hc->f) && (!hc->f(1, j, u, j, VRNA_DECOMP_ML_ML, hc->data)))
        eval = 0;

      if (eval) {
        e = my_c[indx[j] + u] +
            (n_seq) * (u - 1) * P->MLbase;

        switch (fc->type) {
          case VRNA_FC_TYPE_COMPARATIVE:
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][u], SS[s][j], md);
              e     += vrna_E_multibranch_stem(type, S5[s][u], S3[s][j], P);
            }
            break;

          default:
            type  = vrna_get_ptype_md(S[u], S[j], md);
            e     += vrna_E_multibranch_stem(type, S1[u - 1], S1[j + 1], P);
            break;
        }

        if (sc_mb_wrapper.red_ml)
          e += sc_mb_wrapper.red_ml(1, j, u, j, &sc_mb_wrapper);

        fM1_new[j] = MIN2(fM1_new[j], e);
      }
    }

    if (with_gquad) {
      /* g-quads */
      if (j >= VRNA_GQUAD_MIN_BOX_SIZE) {
        for (u = j - VRNA_GQUAD_MIN_BOX_SIZE + 1; u >= 1; u--) {
          eval = (hc->up_ml[1] >= (u - 1)) ? 1 : 0;
          if ((hc->f) && (!hc->f(1, j, u, j, VRNA_DECOMP_ML_ML, hc->data)))
            eval = 0;

          if (eval) {
#ifndef VRNA_DISABLE_C11_FEATURES
            e = vrna_smx_csr_get(c_gq, u, j, INF);
#else
            e = vrna_smx_csr_int_get(c_gq, u, j, INF);
#endif
            if (e != INF) {
              e += (vrna_E_multibranch_stem(0, -1, -1, P) + (u - 1) * P->MLbase) *
                   (int)n_seq;

              if (sc_mb_wrapper.red_ml)
                e += sc_mb_wrapper.red_ml(1, j, u, j, &sc_mb_wrapper);

              fM1_new[j] = MIN2(fM1_new[j], e);
            }
          }
        }
      }
    }
  }


  if (with_gquad) {
    /* consider all configurations where a G-quadruplex spans over the artificial cutpoint */
    unsigned int start = 1;
    if (length > VRNA_GQUAD_MAX_BOX_SIZE)
      start = length - VRNA_GQUAD_MAX_BOX_SIZE;

    /* loop over each possible start of a cutpoint-spanning gquad */
    for (i = start; i <= length; i++) {
      unsigned int  start_j = 1;
      unsigned int  stop_j  = length - 1;
      if ((length - i + 1) < VRNA_GQUAD_MIN_BOX_SIZE)
        start_j = VRNA_GQUAD_MIN_BOX_SIZE + i - length - 1;

      if (length > VRNA_GQUAD_MAX_BOX_SIZE)
        stop_j = VRNA_GQUAD_MAX_BOX_SIZE + i - length - 1;

      if (stop_j >= i)
        stop_j = i - 1;

      for (j = start_j; j <= stop_j; j++) {
#ifndef VRNA_DISABLE_C11_FEATURES
        new_c = vrna_smx_csr_get(c_gq, i, j, INF);
#else
        new_c = vrna_smx_csr_int_get(c_gq, i, j, INF);
#endif

        if (new_c != INF) {
          /* case 1: gquad is the only structure, rest is unpaired */
          if (i - j > 3) {
            /* keep at least 3 unpaired bases between start and end of gquad */
            unsigned int u;
            u = i - j - 1;
            /* 1st, obey hard constraints */
            if (hc->up_ext[j + 1] >= u)
              eval = 1;

            if (hc->f)
              eval = (hc->f(j + 1, length, i, length, VRNA_DECOMP_EXT_EXT, hc->data)) ? eval : 0;

            if (eval) {
              int e = new_c;

              if (md->circ_penalty)
                e += vrna_E_hairpin(u, 0, -1, -1, NULL, P) * (int)n_seq;

              switch (fc->type) {
                case VRNA_FC_TYPE_SINGLE:
                  if (sc) {
                    if (sc->energy_up)
                      e += sc->energy_up[j + 1][u];

                    if (sc->f)
                      e += sc->f(j + 1, length, i, length, VRNA_DECOMP_EXT_EXT, sc->data);
                  }

                  break;
                case VRNA_FC_TYPE_COMPARATIVE:
                  if (scs) {
                    for (s = 0; s < n_seq; s++) {
                      if (scs[s]) {
                        if (scs[s]->energy_up) {
                          s   = a2s[s][j] + 1;
                          us  = a2s[s][i] - s;
                          if (us > 0)
                            e += scs[s]->energy_up[s][us];
                        }

                        if (scs[s]->f)
                          e +=
                            scs[s]->f(j + 1, length, i, length, VRNA_DECOMP_EXT_EXT, scs[s]->data);
                      }
                    }
                  }

                  break;
              }

              if (e < FgH) {
                FgH = e;
                Hgi = i;
                Hgj = j;
              }
            }
          } /* end case 1 */

          /* !!! TODO: Make sure that we leave at least 3 unpaired bases between the two components! */

          /* case 2.1: gquad forms an internal-loop like structure with another gquadruplex */
          for (u1 = 0, p = j + 1; p + VRNA_GQUAD_MIN_BOX_SIZE - 1 < i; p++, u1++) {
            /* obey hard constraints */
            if (hc->up_int[j + 1] < u1)
              break;

            if (u1 > MAXLOOP)
              break;

            for (u2 = 0, q = i - 1; q >= p + VRNA_GQUAD_MIN_BOX_SIZE - 1; q--, u2++) {
              /* obey hard constraints */
              if (hc->up_int[q + 1] < u2)
                break;

              if (u1 + u2 > MAXLOOP)
                break;

              /* obey user-defined hard constraints */
              if (hc->f) {
                vrna_log_debug(
                  "user-defined hard constraints not implemented for int-loop type gquads yet!");
              }

#ifndef VRNA_DISABLE_C11_FEATURES
              e = vrna_smx_csr_get(c_gq, p, q, INF);
#else
              e = vrna_smx_csr_int_get(c_gq, p, q, INF);
#endif
              if (e != INF) {
                switch (fc->type) {
                  case VRNA_FC_TYPE_SINGLE:
                    e += P->internal_loop[u1 + u2];
                    if (sc) {
                      if (sc->energy_up) {
                        if (u1 > 0)
                          e += sc->energy_up[j + 1][u1];

                        if (u2 > 0)
                          e += sc->energy_up[q + 1][u2];
                      }

                      if (sc->f) {
                        vrna_log_debug(
                          "user-defined soft constraints not fully implemented for int-loop type gquads yet");
                      }
                    }

                    break;
                  case VRNA_FC_TYPE_COMPARATIVE:
                    for (s = 0; s < n_seq; s++) {
                      s1  = a2s[s][j] + 1;
                      us1 = a2s[s][p] - s1;
                      s2  = a2s[s][q] + 1;
                      us2 = a2s[s][i] - s2;
                      e   += P->internal_loop[us1 + us2];
                    }

                    if (scs) {
                      for (s = 0; s < n_seq; s++) {
                        if (scs[s]) {
                          if (scs[s]->energy_up) {
                            s1  = a2s[s][j] + 1;
                            us1 = a2s[s][p] - s1;
                            s2  = a2s[s][q] + 1;
                            us2 = a2s[s][i] - s2;
                            if (us1 > 0)
                              e += sc->energy_up[s1][us1];

                            if (us2 > 0)
                              e += sc->energy_up[s2][us2];
                          }

                          if (scs[s]->f) {
                            vrna_log_debug(
                              "user-defined soft constraints not fully implemented for int-loop type gquads yet");
                          }
                        }
                      }
                    }

                    break;
                }

                if (new_c + e < FgI) {
                  FgI = new_c + e;
                  Igi = i;
                  Igj = j;
                  Igp = p;
                  Igq = q;
                  Igg = 1;
                }
              }
            }
          } /* end case 2.1 */

          /* case 2.2: gquad forms an internal-loop like structure with another base pair */
          for (u1 = 0, p = j + 1; p + turn + 1 < i; p++, u1++) {
            /* obey hard constraints */
            if (hc->up_int[j + 1] < u1)
              break;

            if (u1 > MAXLOOP)
              break;

            for (u2 = 0, q = i - 1; q > p + turn; q--, u2++) {
              /* obey hard constraints */
              if (hc->up_int[q + 1] < u2)
                break;

              if (!(hc->mx[length * p + q] &
                    (VRNA_CONSTRAINT_CONTEXT_INT_LOOP | VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC)))
                continue;

              if (u1 + u2 > MAXLOOP)
                break;

              /* obey user-defined hard constraints */
              if (hc->f) {
                vrna_log_debug(
                  "user-defined hard constraints not implemented for int-loop type gquads yet!");
              }

              e = my_c[indx[q] + p];

              if (e != INF) {
                switch (fc->type) {
                  case VRNA_FC_TYPE_SINGLE:
                    type = vrna_get_ptype_md(S[q], S[p], md);
                    if (md->dangles == 2)
                      e += P->mismatchI[type][S1[q + 1]][S1[p - 1]];

                    if (type > 2)
                      e += P->TerminalAU;

                    e += P->internal_loop[u1 + u2];

                    if (sc) {
                      if (sc->energy_up) {
                        if (u1 > 0)
                          e += sc->energy_up[j + 1][u1];

                        if (u2 > 0)
                          e += sc->energy_up[q + 1][u2];
                      }

                      if (sc->f) {
                        vrna_log_debug(
                          "user-defined soft constraints not fully implemented for int-loop type gquads yet");
                      }
                    }

                    break;
                  case VRNA_FC_TYPE_COMPARATIVE:
                    for (s = 0; s < n_seq; s++) {
                      type = vrna_get_ptype_md(SS[s][q], SS[s][p], md);
                      if (md->dangles == 2)
                        e += P->mismatchI[type][S3[s][q]][S5[s][p]];

                      if (type > 2)
                        e += P->TerminalAU;

                      s1  = a2s[s][j] + 1;
                      us1 = a2s[s][p] - s1;
                      s2  = a2s[s][q] + 1;
                      us2 = a2s[s][i] - s2;
                      e   += P->internal_loop[us1 + us2];
                    }

                    if (scs) {
                      for (s = 0; s < n_seq; s++) {
                        if (scs[s]) {
                          if (scs[s]->energy_up) {
                            s1  = a2s[s][j] + 1;
                            us1 = a2s[s][p] - s1;
                            s2  = a2s[s][q] + 1;
                            us2 = a2s[s][i] - s2;
                            if (us1 > 0)
                              e += sc->energy_up[s1][us1];

                            if (us2 > 0)
                              e += sc->energy_up[s2][us2];
                          }

                          if (scs[s]->f) {
                            vrna_log_debug(
                              "user-defined soft constraints not fully implemented for int-loop type gquads yet");
                          }
                        }
                      }
                    }

                    break;
                }

                if (new_c + e < FgI) {
                  FgI = new_c + e;
                  Igi = i;
                  Igj = j;
                  Igp = p;
                  Igq = q;
                  Igg = 0;
                }
              }
            }
          } /* end case 2.2 */

          /* case 3, gquad forms a multi-branch loop like structure with other base pairs or gquadruplexes */
          if (fM2_real[indx[i - 1] + j + 1] != INF) {
            e = new_c +
                fM2_real[indx[i - 1] + j + 1] +
                n_seq *
                (P->MLclosing + vrna_E_multibranch_stem(0, -1, -1, P));
            if (e < FgM) {
              FgM = e;
              Mgi = i;
              Mgj = j;
            }
          }
        }
      }
    }


    /* last case explicitely handled here: everything unpaired except for one gquad somewhere not spanning artifical cutpoint */
    for (i = 1; i + VRNA_GQUAD_MIN_BOX_SIZE - 1 <= length; i++)
      for (j = i + VRNA_GQUAD_MIN_BOX_SIZE - 1;
           (j <= length) && j <= (i + VRNA_GQUAD_MAX_BOX_SIZE - 1);
           j++) {
#ifndef VRNA_DISABLE_C11_FEATURES
        e = vrna_smx_csr_get(c_gq, i, j, INF);
#else
        e = vrna_smx_csr_int_get(c_gq, i, j, INF);
#endif
        if (e != INF) {
          /* obey constraints */
          u1  = i - 1;
          u2  = length - j;
          if (u1 + u2 < 3)
            continue;

          eval = (hc->up_ext[1] >= u1) ? 1 : 0;
          if (u2 > 0)
            eval = (hc->up_ext[j + 1] >= u2) ? eval : 0;

          if (hc->f) {
            if (u1 > 0)
              eval = (hc->f(1, i - 1, 1, i - 1, VRNA_DECOMP_EXT_UP, hc->data)) ? eval : 0;

            if (u2 > 0)
              eval = (hc->f(j + 1, length, j + 1, length, VRNA_DECOMP_EXT_UP, hc->data)) ? eval : 0;
          }

          if (eval) {
            if (md->circ_penalty)
              e += vrna_E_hairpin(u1 + u2, 0, -1, -1, NULL, P) * (int)n_seq;

            /* apply soft constraints, if any */
            switch (fc->type) {
              case VRNA_FC_TYPE_SINGLE:
                if (sc) {
                  if (sc->energy_up) {
                    if (u1 > 0)
                      e += sc->energy_up[1][u1];

                    if (u2 > 0)
                      e += sc->energy_up[j + 1][u2];
                  }

                  if (sc->f) {
                    if (u1 > 0)
                      e += sc->f(1, i - 1, 1, i - 1, VRNA_DECOMP_EXT_UP, sc->data);

                    if (u2 > 0)
                      e += sc->f(j + 1, length, j + 1, length, VRNA_DECOMP_EXT_UP, sc->data);
                  }
                }

                break;

              case VRNA_FC_TYPE_COMPARATIVE:
                if (scs) {
                  for (s = 0; s < n_seq; s++) {
                    s1  = a2s[s][1];
                    us1 = a2s[s][i - 1];
                    s2  = a2s[s][j] + 1;
                    us2 = a2s[s][length] - a2s[s][j];
                    if (scs[s]->energy_up) {
                      if (us1 > 0)
                        e += scs[s]->energy_up[s1][us1];

                      if (us2 > 0)
                        e += scs[s]->energy_up[s2][us2];
                    }

                    if (scs[s]->f) {
                      if (i > 1)
                        e += scs[s]->f(1, i - 1, 1, i - 1, VRNA_DECOMP_EXT_UP, scs[s]->data);

                      if (j < length)
                        e += scs[s]->f(j + 1,
                                       length,
                                       j + 1,
                                       length,
                                       VRNA_DECOMP_EXT_UP,
                                       scs[s]->data);
                    }
                  }
                }

                break;
            }

            if (e < FgH) {
              FgH = e;
              Hgi = i;
              Hgj = j;
            }
          }
        }
      }

    /* internal loop cases with at least one gquad below */
    for (i = 1; i < MAXLOOP + 1; i++) {
      u1 = i - 1;
      if (hc->up_int[1] + 1 >= i) {
        /* [gquad] + (basepair) */
        for (j = i + VRNA_GQUAD_MIN_BOX_SIZE - 1; j + turn + 2 <= length; j++) {
#ifndef VRNA_DISABLE_C11_FEATURES
          e = vrna_smx_csr_get(c_gq, i, j, INF);
#else
          e = vrna_smx_csr_int_get(c_gq, i, j, INF);
#endif
          if (e != INF) {
            for (p = j + 1; p + u1 <= j + MAXLOOP + 1; p++) {
              u2 = p - j - 1;
              unsigned int  u3, us3;
              unsigned int  stop = p + turn + 1;

              if (stop + MAXLOOP < length + u1 + u2)
                stop = length + u1 + u2 - MAXLOOP;

              if (hc->up_int[j + 1] >= u2) {
                for (u3 = 0, q = length; q >= stop; q--, u3++) {
                  if (((u2 == 0) && (u1 + u3 < 3)) ||
                      ((u1 + u3 == 0) && (u2 < 3)))
                    continue;

                  unsigned int sq, sp;
                  sq    = S1[q + 1];
                  sp    = S1[p - 1];
                  eval  = (hc->up_int[q] >= u3) ? 1 : 0;
                  eval  =
                    (hc->mx[length * p + q] &
                     (VRNA_CONSTRAINT_CONTEXT_INT_LOOP |
                      VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC)) ? eval : 0;
                  if (eval) {
                    int pq      = indx[q] + p;
                    int energy  = my_c[pq];
                    if (energy != INF) {
                      switch (fc->type) {
                        case VRNA_FC_TYPE_SINGLE:
                          type = vrna_get_ptype_md(S[q], S[p], md);
                          if (dangles == 2)
                            energy += P->mismatchI[type][sq][sp];

                          if (type > 2)
                            energy += P->TerminalAU;

                          energy += P->internal_loop[u1 + u2 + u3];
                          break;

                        case VRNA_FC_TYPE_COMPARATIVE:
                          for (s = 0; s < n_seq; s++) {
                            type = vrna_get_ptype_md(SS[s][q], SS[s][p], md);
                            if (md->dangles == 2)
                              energy += P->mismatchI[type][S3[s][q]][S5[s][p]];

                            if (type > 2)
                              energy += P->TerminalAU;

                            us1     = (i > 1) ? a2s[s][i - 1] - a2s[s][1] : 0;
                            us2     = a2s[s][p - 1] - a2s[s][j];
                            us3     = a2s[s][length] - a2s[s][q];
                            energy  += P->internal_loop[us1 + us2 + us3];
                          }
                          break;
                      }

                      if (e + energy < FgI) {
                        FgI = e + energy;
                        Igi = i;
                        Igj = j;
                        Igp = p;
                        Igq = q;
                        Igg = 0;
                      }
                    }
                  }
                }
              }
            }
          }
        }

        /* [gquad] + [gquad] */
        for (j = i + VRNA_GQUAD_MIN_BOX_SIZE - 1; j + VRNA_GQUAD_MIN_BOX_SIZE + 1 <= length; j++) {
#ifndef VRNA_DISABLE_C11_FEATURES
          e = vrna_smx_csr_get(c_gq, i, j, INF);
#else
          e = vrna_smx_csr_int_get(c_gq, i, j, INF);
#endif
          if (e != INF) {
            for (p = j + 1; p + VRNA_GQUAD_MIN_BOX_SIZE - 1 <= length; p++) {
              u2 = p - j - 1;
              if (hc->up_int[j + 1] >= u2) {
                unsigned int  stop = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;

                if (stop + MAXLOOP < length + u1 + u2)
                  stop = length + u1 + u2 - MAXLOOP;

                unsigned int  u3, us3;
                for (u3 = 0, q = length; q >= stop; q--, u3++) {
                  if (((u2 == 0) && (u1 + u3 < 3)) ||
                      ((u1 + u3 == 0) && (u2 < 3)))
                    continue;

#ifndef VRNA_DISABLE_C11_FEATURES
                  int energy = vrna_smx_csr_get(c_gq, p, q, INF);
#else
                  int energy = vrna_smx_csr_int_get(c_gq, p, q, INF);
#endif
                  if (energy != INF) {
                    switch (fc->type) {
                      case VRNA_FC_TYPE_SINGLE:
                        energy += P->internal_loop[u1 + u2 + u3];
                        break;

                      case VRNA_FC_TYPE_COMPARATIVE:
                        for (s = 0; s < n_seq; s++) {
                          us1     = (i > 1) ? a2s[s][i - 1] - a2s[s][1] : 0;
                          us2     = a2s[s][p - 1] - a2s[s][j];
                          us3     = a2s[s][length] - a2s[s][q];
                          energy  += P->internal_loop[us1 + us2 + us3];
                        }
                        break;
                    }

                    if (e + energy < FgI) {
                      FgI = e + energy;
                      Igi = i;
                      Igj = j;
                      Igp = p;
                      Igq = q;
                      Igg = 1;
                    }
                  }
                }
              }
            }
          }
        }

        /* (basepair) + [gquad] */
        for (j = i + turn + 1; j + VRNA_GQUAD_MIN_BOX_SIZE <= length; j++) {
          eval  = (hc->mx[length * i + j] & (VRNA_CONSTRAINT_CONTEXT_INT_LOOP | VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC)) ? 1 : 0;
          ij    = indx[j] + i;
          e     = my_c[ij];
          if ((eval) &&
              (e != INF)) {
            switch (fc->type) {
              case VRNA_FC_TYPE_SINGLE:
                si    = S1[i - 1];
                sj    = S1[j + 1];
                type  = vrna_get_ptype_md(S[j], S[i], md);
                if (dangles == 2)
                  e += P->mismatchI[type][sj][si];

                if (type > 2)
                  e += P->TerminalAU;

                break;

              case VRNA_FC_TYPE_COMPARATIVE:
                for (s = 0; s < n_seq; s++) {
                  type = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
                  if (md->dangles == 2)
                    e += P->mismatchI[type][S3[s][j]][S5[s][i]];

                  if (type > 2)
                    e += P->TerminalAU;
                }
                break;
            }

            for (p = j + 1; p + VRNA_GQUAD_MIN_BOX_SIZE - 1 <= length; p++) {
              u2 = p - j - 1;
              if (hc->up_int[j + 1] < u2)
                break;

              unsigned int  stop = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
              if (stop + MAXLOOP < length + u1 + u2)
                stop = length + u1 + u2 - MAXLOOP;

              unsigned int  u3, us3;
              for (u3 = 0, q = length; q >= stop; q--, u3++) {
                if (((u2 == 0) && (u1 + u3 < 3)) ||
                    ((u1 + u3 == 0) && (u2 < 3)))
                  continue;

                eval = (hc->up_int[q] >= u3) ? 1 : 0;
                if (eval) {
#ifndef VRNA_DISABLE_C11_FEATURES
                  int energy = vrna_smx_csr_get(c_gq, p, q, INF);
#else
                  int energy = vrna_smx_csr_int_get(c_gq, p, q, INF);
#endif
                  if (energy != INF) {
                    switch (fc->type) {
                      case VRNA_FC_TYPE_SINGLE:
                        energy += P->internal_loop[u1 + u2 + u3];
                        break;

                      case VRNA_FC_TYPE_COMPARATIVE:
                        for (s = 0; s < n_seq; s++) {
                          us1     = (i > 1) ? a2s[s][i - 1] - a2s[s][1] : 0;
                          us2     = a2s[s][p - 1] - a2s[s][j];
                          us3     = a2s[s][length] - a2s[s][q];
                          energy  += P->internal_loop[us1 + us2 + us3];
                        }
                        break;
                    }

                    if (e + energy < FgI) {
                      FgI = e + energy;
                      Igi = p;
                      Igj = q;
                      Igp = i;
                      Igq = j;
                      Igg = 0;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    vrna_log_debug("FgH = %d, FgI = %d, FgM = %d", FgH, FgI, FgM);
    Fc  = MIN2(Fc, FgH);
    Fc  = MIN2(Fc, FgI);
    Fc  = MIN2(Fc, FgM);
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
      new_c = vrna_eval_hairpin(fc, j, i, VRNA_EVAL_LOOP_DEFAULT);
      if (new_c != INF)
        new_c += my_c[ij];

      if (new_c < FcH) {
        FcH = new_c;
        Hi  = i;
        Hj  = j;
      }

      /* exterior internal loop case */
      new_c = vrna_mfe_internal_ext(fc, i, j, &p, &q);
      if (new_c != INF)
        new_c += my_c[ij];

      if (p != 0) {
        if (new_c < FcI) {
          FcI = new_c;
          Ii  = i;
          Ij  = j;
          Ip  = p;
          Iq  = q;
        }
      }
    } /* end of i,j loop */
  Fc  = MIN2(Fc, FcH);
  Fc  = MIN2(Fc, FcI);

  /* use fM1_new and fM2 to construct segments with at least 3 branches */
  unsigned int  space3 = 2 * MIN2(turn + 2, VRNA_GQUAD_MIN_BOX_SIZE);

  int           *fM1_tmp = fM1_new;

  /* apply hard constraints if necessary */
  if (hc->f) {
    fM1_tmp = (int *)vrna_alloc(sizeof(int) * (length + 2));
    fM1_tmp = memcpy(fM1_tmp, fM1_new, sizeof(int) * (length + 2));

    for (u = turn + 2; u <= length; u++)
      if (!hc->f(1, length, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data))
        fM1_tmp[u] = INF;
  }

  /* apply soft constraints if necessary */
  if (sc_mb_wrapper.decomp_ml) {
    if (fM1_tmp == fM1_new) {
      fM1_tmp = (int *)vrna_alloc(sizeof(int) * (length + 2));
      fM1_tmp = memcpy(fM1_tmp, fM1_new, sizeof(int) * (length + 2));
    }

    for (u = turn + 2; u <= length; u++)
      fM1_tmp[u] = ADD_OR_INF(fM1_tmp[u], sc_mb_wrapper.decomp_ml(1, length, u, u + 1, &sc_mb_wrapper));
  }

  for (u = MIN2(turn + 2, VRNA_GQUAD_MIN_BOX_SIZE); u + space3 <= length; u++) {
    if ((fM1_tmp[u] != INF) &&
        (fM2_real[indx[length] + u + 1] != INF)) {
      new_c = fM1_tmp[u] +
              fM2_real[indx[length] + u + 1];

      if (new_c < FcM) {
        FcM = new_c;
        Mi  = u;
      }
    }
  }

  if (FcM != INF)
    FcM += n_seq * P->MLclosing;

  if (fM1_tmp != fM1_new)
    free(fM1_tmp);

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
                tmp   = vrna_E_multibranch_stem(type, -1, S1[1], P) +
                        P->MLclosing;
                break;

              case VRNA_FC_TYPE_COMPARATIVE:
                tmp = P->MLclosing * n_seq;
                for (s = 0; s < n_seq; s++) {
                  type  = vrna_get_ptype_md(SS[s][i + 1], SS[s][length], md);
                  tmp   += vrna_E_multibranch_stem(type, -1, S3[s][length], P);
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
                tmp   = vrna_E_multibranch_stem(type, S1[i], S1[1], P) +
                        P->MLclosing;
                break;

              case VRNA_FC_TYPE_COMPARATIVE:
                tmp = P->MLclosing * n_seq;
                for (s = 0; s < n_seq; s++) {
                  type  = vrna_get_ptype_md(SS[s][i + 1], SS[s][length], md);
                  tmp   += vrna_E_multibranch_stem(type, S5[s][i + 1], S3[s][length], P);
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
            Md3i  = -(int)i;
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
                tmp   = vrna_E_multibranch_stem(type, S1[length], -1, P) +
                        P->MLclosing;
                break;

              case VRNA_FC_TYPE_COMPARATIVE:
                tmp = P->MLclosing * n_seq;
                for (s = 0; s < n_seq; s++) {
                  type  = vrna_get_ptype_md(SS[s][1], SS[s][i], md);
                  tmp   += vrna_E_multibranch_stem(type, S5[s][1], -1, P);
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
                tmp   = vrna_E_multibranch_stem(type, S1[length], S1[i + 1], P) +
                        P->MLclosing;
                break;

              case VRNA_FC_TYPE_COMPARATIVE:
                tmp = P->MLclosing * n_seq;
                for (s = 0; s < n_seq; s++) {
                  type  = vrna_get_ptype_md(SS[s][1], SS[s][i], md);
                  tmp   += vrna_E_multibranch_stem(type, S5[s][1], S3[s][i], P);
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
            Md5i  = -(int)i;
          }
        }
      }

      free(fmd5_tmp);
    }

    /* Start backtracing in exterior loops and prepare for everything else */

    if (FcMd5 < MIN2(Fc, FcMd3)) {
      int real_i, sc_en = 0;

      /* looks like we have to do this ... */
      vrna_bts_push(bt_stack,
                    ((vrna_sect_t){
                      .i = 1,
                      .j = (Md5i > 0) ? Md5i : -Md5i,
                      .ml = VRNA_MX_FLAG_C
                    }));

      i       = (Md5i > 0) ? Md5i + 1 : -Md5i + 2;             /* let's backtrack fm_d5[Md5i+1] */
      real_i  = (Md5i > 0) ? i : i - 1;

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
          vrna_bts_push(bt_stack,
                        ((vrna_sect_t){
            .i = i,
            .j = u,
            .ml = VRNA_MX_FLAG_M}));
          vrna_bts_push(bt_stack,
                        ((vrna_sect_t){
                        .i = u + 1,
                        .j   = length - 1,
                        .ml  = VRNA_MX_FLAG_M}));
          break;
        }
      }
      Fc = FcMd5;
    } else if (FcMd3 < Fc) {
      int real_i, sc_en = 0;
      /* here we go again... */
      vrna_bts_push(bt_stack,
                    ((vrna_sect_t){
                    .i = (Md3i > 0) ? Md3i + 1 : -Md3i + 1,
                    .j   = length,
                    .ml  = VRNA_MX_FLAG_C}));

      i       = (Md3i > 0) ? Md3i : -Md3i - 1;             /* let's backtrack fm_d3[Md3i] */
      real_i  = (Md3i > 0) ? i : i + 1;

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
          vrna_bts_push(bt_stack,
                        ((vrna_sect_t){
                        .i = 2,
                        .j = u,
                        .ml = VRNA_MX_FLAG_M}));
          vrna_bts_push(bt_stack,
                        ((vrna_sect_t){
                        .i = u + 1,
                        .j   = i,
                        .ml  = VRNA_MX_FLAG_M}));
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
      vrna_bts_push(bt_stack, ((vrna_sect_t){
        .i = Hi,
        .j = Hj,
        .ml = VRNA_MX_FLAG_C}));
    } else if (FcI == Fc) {
      vrna_bts_push(bt_stack, ((vrna_sect_t){
        .i = Ii,
        .j = Ij,
        .ml = VRNA_MX_FLAG_C}));
      vrna_bts_push(bt_stack, ((vrna_sect_t){
        .i = Ip,
        .j = Iq,
        .ml = VRNA_MX_FLAG_C}));
    } else if (FcM == Fc) {
      /* 1. find component in fM1_new */
      for (u = 1; u + MIN2(turn + 1, VRNA_GQUAD_MIN_BOX_SIZE - 1) <= Mi; u++) {
        if (hc->up_ml[1] < (u - 1))
          break;

        eval = ((hc->f) && (!hc->f(1, Mi, u, Mi, VRNA_DECOMP_ML_ML, hc->data))) ? 0 : 1;

        if ((eval) &&
            (hc->mx[length * Mi + u] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP)) {
          e = my_c[indx[Mi] + u] +
              (n_seq) * (u - 1) * P->MLbase;

          switch (fc->type) {
            case VRNA_FC_TYPE_COMPARATIVE:
              for (s = 0; s < n_seq; s++) {
                type  = vrna_get_ptype_md(SS[s][u], SS[s][Mi], md);
                e     += vrna_E_multibranch_stem(type, S5[s][u], S3[s][Mi], P);
              }
              break;

            default:
              type  = vrna_get_ptype_md(S[u], S[Mi], md);
              e     += vrna_E_multibranch_stem(type, S1[u - 1], S1[Mi + 1], P);
              break;
          }

          if (sc_mb_wrapper.red_ml)
            e += sc_mb_wrapper.red_ml(1, Mi, u, Mi, &sc_mb_wrapper);

          if (e == fM1_new[Mi]) {
            vrna_bts_push(bt_stack, ((vrna_sect_t){
                    .i = u,
                    .j = Mi,
                    .ml = VRNA_MX_FLAG_C}));
            break;
          }
        }

        if ((with_gquad) &&
            (eval) &&
            (u + VRNA_GQUAD_MIN_BOX_SIZE - 1 <= Mi)) {
#ifndef VRNA_DISABLE_C11_FEATURES
          e = vrna_smx_csr_get(c_gq, u, Mi, INF);
#else
          e = vrna_smx_csr_int_get(c_gq, u, Mi, INF);
#endif
          if (e != INF) {
            e += (vrna_E_multibranch_stem(0, -1, -1, P) + (u - 1) * P->MLbase) *
                 (int)n_seq;

            if (sc_mb_wrapper.red_ml)
              e += sc_mb_wrapper.red_ml(1, j, u, j, &sc_mb_wrapper);

            if (e == fM1_new[Mi]) {
              vrna_bts_push(bt_stack, ((vrna_sect_t){
                    .i = u,
                    .j = Mi,
                    .ml = VRNA_MX_FLAG_G}));
              break;
            }
          }
        }
      }

      if (u > Mi)
        vrna_log_error("Backtrack failed in fM1_new[%d] = %d", Mi, fM1_new[Mi]);

      /* 2. find split-point in fM2_real */
      for (u = Mi + 1 + MIN2(turn + 1, VRNA_GQUAD_MIN_BOX_SIZE - 1);
           u + MIN2(turn + 2, VRNA_GQUAD_MIN_BOX_SIZE) <= length;
           u++) {
        if ((my_fML[indx[u] + Mi + 1] == INF) ||
            (my_fML[indx[length] + u + 1] == INF) ||
            ((hc->f) && (!hc->f(Mi + 1, length, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data))))
          continue;

        new_c = my_fML[indx[u] + Mi + 1] +
                my_fML[indx[length] + u + 1];

        if (sc_mb_wrapper.decomp_ml)
          new_c += sc_mb_wrapper.decomp_ml(Mi + 1, length, u, u + 1, &sc_mb_wrapper);

        if (new_c == fM2_real[indx[length] + Mi + 1]) {
            vrna_bts_push(bt_stack, ((vrna_sect_t){
                  .i = Mi + 1,
                  .j = u,
                  .ml = VRNA_MX_FLAG_M}));
            vrna_bts_push(bt_stack, ((vrna_sect_t){
                  .i = u + 1,
                  .j = length,
                  .ml = VRNA_MX_FLAG_M}));
            break;
        }
      }

      if (u + MIN2(turn + 2, VRNA_GQUAD_MIN_BOX_SIZE) > length)
        vrna_log_error("Backtrack failed in fM2_real[%d][%d] = %d", Mi + 1, length,
                       fM2_real[indx[length] + Mi + 1]);
    } else if (Fc == FcO) {
      /* unstructured */
      vrna_bts_push(bt_stack, ((vrna_sect_t){
        .i = 1,
        .j = 1,
        .ml = VRNA_MX_FLAG_F5}));
    } else if (with_gquad) {
      if (Fc == FgH) {
        vrna_bts_push(bt_stack, ((vrna_sect_t){
          .i = Hgi,
          .j = Hgj,
          .ml = VRNA_MX_FLAG_G}));
      } else if (Fc == FgI) {
        vrna_bts_push(bt_stack, ((vrna_sect_t){
          .i = Igi,
          .j = Igj,
          .ml = VRNA_MX_FLAG_G}));
        if (Igg) {
          vrna_bts_push(bt_stack, ((vrna_sect_t){
            .i = Igp,
            .j = Igq,
            .ml = VRNA_MX_FLAG_G}));
        } else {
          vrna_bts_push(bt_stack, ((vrna_sect_t){
            .i = Igp,
            .j = Igq,
            .ml = VRNA_MX_FLAG_C}));
        }
      } else if (Fc == FgM) {
        vrna_bts_push(bt_stack, ((vrna_sect_t){
          .i = Mgi,
          .j = Mgj,
          .ml = VRNA_MX_FLAG_G}));

        /* determine split index for fM2_real */
        unsigned int  ii, jj;
        int           ee;
        ii  = Mgj + 1;
        jj  = Mgi - 1;
        ee  = fM2_real[indx[jj] + ii];

        for (u = ii + turn + 1; u + 1 + turn + 1 <= jj; u++) {
          if ((my_fML[indx[u] + ii] == INF) ||
              (my_fML[indx[jj] + u + 1] == INF) ||
              ((hc->f) && (!hc->f(ii, jj, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data))))
            continue;

          e = my_fML[indx[u] + ii] + my_fML[indx[jj] + u + 1];

          if (sc_mb_wrapper.decomp_ml)
            e += sc_mb_wrapper.decomp_ml(ii, jj, u, u + 1, &sc_mb_wrapper);

          if (ee == e) {
            vrna_bts_push(bt_stack, ((vrna_sect_t){
              .i = ii,
              .j = u,
              .ml = VRNA_MX_FLAG_M}));
            vrna_bts_push(bt_stack, ((vrna_sect_t){
              .i = u + 1,
              .j = jj,
              .ml = VRNA_MX_FLAG_M}));
            break;
          }
        }
        if (u + 1 + turn + 1 > jj)
          vrna_log_error("Backtracking failed for gquad");
      }
    }
  } else {
    /* forbidden, i.e. no configuration fulfills constraints */
  }

  fc->matrices->FcH = FcH;
  fc->matrices->FcI = FcI;
  fc->matrices->FcM = FcM;
  fc->matrices->Fc  = Fc;

  free_sc_mb(&sc_mb_wrapper);

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


#if 1

PRIVATE struct ms_helpers *
get_ms_helpers(vrna_fold_compound_t *fc)
{
  struct ms_helpers *dat = (struct ms_helpers *)vrna_alloc(sizeof(struct ms_helpers));

  dat->evaluate = prepare_hc_ext_def(fc, &(dat->hc_dat_local));

  init_sc_f5(fc, &(dat->sc_wrapper));

  return dat;
}


PRIVATE void
free_ms_helpers(struct ms_helpers *ms_dat,
                size_t            strands)
{
  if (strands > 1)
    free_sc_f5(&(ms_dat->sc_wrapper));

  free(ms_dat);
}


PRIVATE void
update_fms5_arrays(vrna_fold_compound_t *fc,
                   int                  i,
                   struct ms_helpers    *ms_dat)
{
  short                 *S1, *S2, s5, s3;
  unsigned int          *sn, *se, type;
  int                   e, tmp, **fms5, *c, *idx, n, end, dangle_model, base;
  vrna_param_t          *params;
  vrna_md_t             *md;
  vrna_hc_eval_f        evaluate;
  struct hc_ext_def_dat *hc_dat_local;
  struct sc_f5_dat      *sc_wrapper;
  sc_ext_red_cb         sc_spl;
  sc_ext_red_cb         sc_red_stem;
  sc_ext_red_cb         sc_red_ext;

  n             = fc->length;
  S1            = fc->sequence_encoding;
  S2            = fc->sequence_encoding2;
  idx           = fc->jindx;
  c             = fc->matrices->c;
  fms5          = fc->matrices->fms5;
  sn            = fc->strand_number;
  se            = fc->strand_end;
  params        = fc->params;
  md            = &(params->model_details);
  dangle_model  = md->dangles;
  evaluate      = ms_dat->evaluate;
  hc_dat_local  = &(ms_dat->hc_dat_local);
  sc_wrapper    = &(ms_dat->sc_wrapper);
  sc_spl        = sc_wrapper->decomp;
  sc_red_stem   = sc_wrapper->red_stem;
  sc_red_ext    = sc_wrapper->red_ext;


  for (size_t s = 0; s < fc->strands; s++) {
    e   = tmp = INF;
    end = se[s];

    if (i < end) {
      if (evaluate(i, end, i + 1, end, VRNA_DECOMP_EXT_EXT, hc_dat_local)) {
        tmp = fms5[s][i + 1];
        if ((sc_red_ext) && (tmp != INF))
          tmp += sc_red_ext(i, end, i + 1, end, sc_wrapper);
      }
    } else {
      tmp = 0;
    }

    e = MIN2(e, tmp);

    for (int k = i + 1; k < end; k++) {
      if ((evaluate(i, end, k, k + 1, VRNA_DECOMP_EXT_STEM_EXT, hc_dat_local)) &&
          (fms5[s][k + 1] != INF) &&
          (c[idx[k] + i] != INF)) {
        type = vrna_get_ptype_md(S2[i], S2[k], md);

        switch (dangle_model) {
          case 2:
            s5  = ((i > 1) && (sn[i - 1] == sn[i])) ? S1[i - 1] : -1;
            s3  = ((k < n) && (sn[k] == sn[k + 1])) ? S1[k + 1] : -1;
            break;

          default:
            s5 = s3 = -1;
            break;
        }

        base = vrna_E_exterior_stem(type, s5, s3, params);

        tmp = base +
              fms5[s][k + 1] +
              c[idx[k] + i];

        if (sc_red_stem)
          tmp += sc_red_stem(i, k, i, k, sc_wrapper);

        if (sc_spl)
          tmp += sc_spl(i, end, k, k + 1, sc_wrapper);

        e = MIN2(e, tmp);
      }
    }

    if ((evaluate(i, end, i, end, VRNA_DECOMP_EXT_STEM, hc_dat_local)) &&
        (c[idx[end] + i] != INF)) {
      type = vrna_get_ptype_md(S2[i], S2[end], md);
      switch (dangle_model) {
        case 2:
          s5  = ((i > 1) && (sn[i - 1] == sn[i])) ? S1[i - 1] : -1;
          s3  = -1;
          break;

        default:
          s5 = s3 = -1;
          break;
      }

      base = vrna_E_exterior_stem(type, s5, s3, params);

      tmp = base +
            c[idx[end] + i];

      if (sc_red_stem)
        tmp += sc_red_stem(i, end, i, end, sc_wrapper);

      e = MIN2(e, tmp);
    }

    /* damn odd dangles below */
    if (dangle_model % 2) {
      s5 = S1[i];

      for (int k = i + 2; k + 1 < end; k++) {
        if ((evaluate(i, end, k, k + 1, VRNA_DECOMP_EXT_STEM_EXT1, hc_dat_local)) &&
            (fms5[s][k + 1] != INF) &&
            (c[idx[k] + i + 1] != INF)) {
          type  = vrna_get_ptype_md(S2[i + 1], S2[k], md);
          base  = vrna_E_exterior_stem(type, s5, -1, params);
          tmp   = base +
                  fms5[s][k + 1] +
                  c[idx[k] + i + 1];

          if (sc_red_stem)
            tmp += sc_red_stem(i, k, i + 1, k, sc_wrapper);

          if (sc_spl)
            tmp += sc_spl(i, end, k, k + 1, sc_wrapper);

          e = MIN2(e, tmp);
        }
      }

      for (int k = i + 1; k + 1 < end; k++) {
        s3 = S1[k + 1];

        if ((evaluate(i, end, k, k + 2, VRNA_DECOMP_EXT_STEM_EXT, hc_dat_local)) &&
            (fms5[s][k + 2] != INF) &&
            (c[idx[k] + i] != INF)) {
          type  = vrna_get_ptype_md(S2[i], S2[k], md);
          base  = vrna_E_exterior_stem(type, -1, s3, params);
          tmp   = base +
                  fms5[s][k + 2] +
                  c[idx[k] + i];

          if (sc_red_stem)
            tmp += sc_red_stem(i, k + 1, i, k, sc_wrapper);

          if (sc_spl)
            tmp += sc_spl(i, end, k + 1, k + 2, sc_wrapper);

          e = MIN2(e, tmp);
        }

        if ((evaluate(i, end, k, k + 2, VRNA_DECOMP_EXT_STEM_EXT1, hc_dat_local)) &&
            (fms5[s][k + 2] != INF) &&
            (c[idx[k] + i + 1] != INF) &&
            (i + 1 < k)) {
          type  = vrna_get_ptype_md(S2[i + 1], S2[k], md);
          base  = vrna_E_exterior_stem(type, s5, s3, params);
          tmp   = base +
                  fms5[s][k + 2] +
                  c[idx[k] + i + 1];

          if (sc_red_stem)
            tmp += sc_red_stem(i, k + 1, i + 1, k, sc_wrapper);

          if (sc_spl)
            tmp += sc_spl(i, end, k + 1, k + 2, sc_wrapper);

          e = MIN2(e, tmp);
        }
      } /* end for (k... */

      s5  = S1[i];
      s3  = S1[end];

      if ((evaluate(i, end, i, end - 1, VRNA_DECOMP_EXT_STEM, hc_dat_local)) &&
          (c[idx[end - 1] + i] != INF) &&
          (i + 1 < end)) {
        type  = vrna_get_ptype_md(S2[i], S2[end - 1], md);
        base  = vrna_E_exterior_stem(type, -1, s3, params);
        tmp   = base +
                c[idx[end - 1] + i];

        if (sc_red_stem)
          tmp += sc_red_stem(i, end, i, end - 1, sc_wrapper);

        e = MIN2(e, tmp);
      }

      if ((evaluate(i, end, i + 1, end, VRNA_DECOMP_EXT_STEM, hc_dat_local)) &&
          (c[idx[end] + i + 1] != INF) &&
          (i + 1 < end)) {
        type  = vrna_get_ptype_md(S2[i + 1], S2[end], md);
        base  = vrna_E_exterior_stem(type, s5, -1, params);
        tmp   = base +
                c[idx[end] + i + 1];

        if (sc_red_stem)
          tmp += sc_red_stem(i, end, i + 1, end, sc_wrapper);

        e = MIN2(e, tmp);
      }

      if ((evaluate(i, end, i + 1, end - 1, VRNA_DECOMP_EXT_STEM, hc_dat_local)) &&
          (c[idx[end - 1] + i + 1] != INF) &&
          (i + 2 < end)) {
        type  = vrna_get_ptype_md(S2[i + 1], S2[end - 1], md);
        base  = vrna_E_exterior_stem(type, s5, s3, params);
        tmp   = base +
                c[idx[end - 1] + i + 1];

        if (sc_red_stem)
          tmp += sc_red_stem(i, end, i + 1, end - 1, sc_wrapper);

        e = MIN2(e, tmp);
      }
    }

    fms5[s][i] = e;
  }
}


PRIVATE void
update_fms3_arrays(vrna_fold_compound_t *fc,
                   unsigned int         s,
                   struct ms_helpers    *ms_dat)
{
  short                 *S1, *S2, s5, s3;
  unsigned int          *sn, *ss, type;
  int                   *c, **fms3, base, e, tmp, j, k, start, n, *idx, dangle_model;
  vrna_param_t          *params;
  vrna_md_t             *md;
  vrna_hc_eval_f        evaluate;
  struct hc_ext_def_dat *hc_dat_local;
  struct sc_f5_dat      *sc_wrapper;
  sc_ext_red_cb         sc_spl;
  sc_ext_red_cb         sc_red_stem;
  sc_ext_red_cb         sc_red_ext;

  n             = fc->length;
  S1            = fc->sequence_encoding;
  S2            = fc->sequence_encoding2;
  idx           = fc->jindx;
  c             = fc->matrices->c;
  fms3          = fc->matrices->fms3;
  sn            = fc->strand_number;
  ss            = fc->strand_start;
  params        = fc->params;
  md            = &(params->model_details);
  dangle_model  = md->dangles;
  evaluate      = ms_dat->evaluate;
  hc_dat_local  = &(ms_dat->hc_dat_local);
  sc_wrapper    = &(ms_dat->sc_wrapper);
  sc_spl        = sc_wrapper->decomp;
  sc_red_stem   = sc_wrapper->red_stem;
  sc_red_ext    = sc_wrapper->red_ext;

  start = ss[s];

  for (j = start; j <= n; j++) {
    e = tmp = INF;

    if (start < j) {
      if (evaluate(start, j, start, j - 1, VRNA_DECOMP_EXT_EXT, hc_dat_local)) {
        tmp = fms3[s][j - 1];

        if (sc_red_ext)
          tmp += sc_red_ext(start, j, start, j - 1, sc_wrapper);

        e = MIN2(e, tmp);
      }
    } else {
      e = 0;
    }

    if (evaluate(start, j, start, j, VRNA_DECOMP_EXT_STEM, hc_dat_local) &&
        (c[idx[j] + start] != INF)) {
      type = vrna_get_ptype_md(S2[start], S2[j], md);

      switch (dangle_model) {
        case 2:
          s5  = -1;
          s3  = (sn[j] == sn[j + 1]) ? S1[j + 1] : -1;
          break;

        default:
          s5 = s3 = -1;
          break;
      }

      base = vrna_E_exterior_stem(type, s5, s3, params);

      tmp = base +
            c[idx[j] + start];

      if (sc_red_stem)
        tmp += sc_red_stem(start, j, start, j, sc_wrapper);

      e = MIN2(e, tmp);
    }

    for (k = start; k < j; k++) {
      if (evaluate(start, j, k, k + 1, VRNA_DECOMP_EXT_EXT_STEM, hc_dat_local) &&
          (fms3[s][k] != INF) &&
          (c[idx[j] + k + 1] != INF)) {
        type = vrna_get_ptype_md(S2[k + 1], S2[j], md);

        switch (dangle_model) {
          case 2:
            s5  = (sn[k] == sn[k + 1]) ? S1[k] : -1;
            s3  = (sn[j] == sn[j + 1]) ? S1[j + 1] : -1;
            break;

          default:
            s5 = s3 = -1;
            break;
        }

        base = vrna_E_exterior_stem(type, s5, s3, params);

        tmp = base +
              fms3[s][k] +
              c[idx[j] + k + 1];

        if (sc_red_stem)
          tmp += sc_red_stem(k + 1, j, k + 1, j, sc_wrapper);

        if (sc_spl)
          tmp += sc_spl(start, j, k, k + 1, sc_wrapper);

        e = MIN2(e, tmp);
      }
    }

    if (dangle_model % 2) {
      s5  = S1[start];
      s3  = S1[j];

      if (evaluate(start, j, start, j - 1, VRNA_DECOMP_EXT_STEM, hc_dat_local) &&
          (c[idx[j - 1] + start] != INF) &&
          (start + 1 < j)) {
        type  = vrna_get_ptype_md(S2[start], S2[j - 1], md);
        base  = vrna_E_exterior_stem(type, -1, s3, params);
        tmp   = base +
                c[idx[j - 1] + start];

        if (sc_red_stem)
          tmp += sc_red_stem(start, j, start, j - 1, sc_wrapper);

        e = MIN2(e, tmp);
      }

      if (evaluate(start, j, start + 1, j, VRNA_DECOMP_EXT_STEM, hc_dat_local) &&
          (c[idx[j] + start + 1] != INF) &&
          (start + 1 < j)) {
        type  = vrna_get_ptype_md(S2[start + 1], S2[j], md);
        base  = vrna_E_exterior_stem(type, s5, -1, params);
        tmp   = base +
                c[idx[j] + start + 1];

        if (sc_red_stem)
          tmp += sc_red_stem(start, j, start + 1, j, sc_wrapper);

        e = MIN2(e, tmp);
      }

      if (evaluate(start, j, start + 1, j - 1, VRNA_DECOMP_EXT_STEM, hc_dat_local) &&
          (c[idx[j - 1] + start + 1] != INF) &&
          (start + 2 < j)) {
        type  = vrna_get_ptype_md(S2[start + 1], S2[j - 1], md);
        base  = vrna_E_exterior_stem(type, s5, s3, params);
        tmp   = base +
                c[idx[j - 1] + start + 1];

        if (sc_red_stem)
          tmp += sc_red_stem(start, j, start + 1, j - 1, sc_wrapper);

        e = MIN2(e, tmp);
      }

      for (k = start; k < j; k++) {
        s5 = S1[k + 1];

        if (evaluate(start, j, k, k + 2, VRNA_DECOMP_EXT_EXT_STEM, hc_dat_local) &&
            (fms3[s][k] != INF) &&
            (c[idx[j] + k + 2] != INF) &&
            k + 2 < j) {
          type  = vrna_get_ptype_md(S2[k + 2], S2[j], md);
          base  = vrna_E_exterior_stem(type, s5, -1, params);
          tmp   = base +
                  fms3[s][k] +
                  c[idx[j] + k + 2];

          if (sc_red_stem)
            tmp += sc_red_stem(k + 1, j, k + 2, j, sc_wrapper);

          if (sc_spl)
            tmp += sc_spl(start, j, k, k + 1, sc_wrapper);

          e = MIN2(e, tmp);
        }

        if (evaluate(start, j, k, k + 1, VRNA_DECOMP_EXT_EXT_STEM1, hc_dat_local) &&
            (fms3[s][k] != INF) &&
            (c[idx[j - 1] + k + 1] != INF) &&
            (k + 2 < j)) {
          type  = vrna_get_ptype_md(S2[k + 1], S2[j - 1], md);
          base  = vrna_E_exterior_stem(type, -1, s3, params);
          tmp   = base +
                  fms3[s][k] +
                  c[idx[j - 1] + k + 1];

          if (sc_red_stem)
            tmp += sc_red_stem(k + 1, j, k + 1, j - 1, sc_wrapper);

          if (sc_spl)
            tmp += sc_spl(start, j, k, k + 1, sc_wrapper);

          e = MIN2(e, tmp);
        }

        if (evaluate(start, j, k, k + 2, VRNA_DECOMP_EXT_EXT_STEM1, hc_dat_local) &&
            (fms3[s][k] != INF) &&
            (c[idx[j - 1] + k + 2] != INF) &&
            (k + 3 < j)) {
          type  = vrna_get_ptype_md(S2[k + 2], S2[j - 1], md);
          base  = vrna_E_exterior_stem(type, s5, s3, params);
          tmp   = base +
                  fms3[s][k] +
                  c[idx[j - 1] + k + 2];

          if (sc_red_stem)
            tmp += sc_red_stem(k + 1, j, k + 2, j - 1, sc_wrapper);

          if (sc_spl)
            tmp += sc_spl(start, j, k, k + 1, sc_wrapper);

          e = MIN2(e, tmp);
        }
      }
    }

    fms3[s][j] = e;
  }
}


PRIVATE int
pair_multi_strand(vrna_fold_compound_t  *fc,
                  int                   i,
                  int                   j,
                  struct ms_helpers     *ms_dat)
{
  short                 *S1, *S2, s5, s3;
  unsigned int          *sn, *ends, type, nick;
  int                   start, contribution, **fms5, **fms3, base, tmp, tmp2, dangle_model;
  vrna_param_t          *params;
  vrna_md_t             *md;
  vrna_hc_eval_f        evaluate;
  struct hc_ext_def_dat *hc_dat_local;
  struct sc_f5_dat      *sc_wrapper;
  sc_ext_red_cb         sc_red_stem;

  contribution  = INF;
  S1            = fc->sequence_encoding;
  S2            = fc->sequence_encoding2;
  params        = fc->params;
  md            = &(params->model_details);
  dangle_model  = md->dangles;
  sn            = fc->strand_number;
  ends          = fc->strand_end;
  fms5          = fc->matrices->fms5;
  fms3          = fc->matrices->fms3;
  evaluate      = ms_dat->evaluate;
  hc_dat_local  = &(ms_dat->hc_dat_local);
  sc_wrapper    = &(ms_dat->sc_wrapper);
  sc_red_stem   = sc_wrapper->red_stem;

  if ((sn[i] != sn[j]) &&
      (evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
    /* most obious strand nick is at end of sn[i] and start of sn[j] */
    type = vrna_get_ptype_md(S2[j], S2[i], md);

    switch (dangle_model) {
      case 2:
        s5  = (sn[j] == sn[j - 1]) ? S1[j - 1] : -1;
        s3  = (sn[i] == sn[i + 1]) ? S1[i + 1] : -1;
        break;

      default:
        s5 = s3 = -1;
        break;
    }

    base = vrna_E_exterior_stem(type, s5, s3, params) +
           params->DuplexInit;

    if (sc_red_stem)
      base += sc_red_stem(j, i, j, i, sc_wrapper);

    tmp = INF;

    /*
     *  if (evaluate(i + 1,
     *               j - 1,
     *               ends[sn[i]],
     *               ends[sn[i]] + 1,
     *               VRNA_DECOMP_EXT_EXT_EXT,
     *               &hc_dat_local))
     */
    if (sn[i] != sn[i + 1]) {
      if ((sn[j - 1] != sn[j]) &&
          (i + 1 == j)) {
        tmp2  = 0;
        tmp   = MIN2(tmp, tmp2);
      } else if (sn[j - 1] == sn[j]) {
        tmp2  = fms3[sn[i + 1]][j - 1];
        tmp   = MIN2(tmp, tmp2);
      }
    } else if (sn[j - 1] != sn[j]) {
      tmp2  = fms5[sn[j - 1]][i + 1];
      tmp   = MIN2(tmp, tmp2);
    } else {
      if ((fms5[sn[i]][i + 1] != INF) &&
          (fms3[sn[ends[sn[i]] + 1]][j - 1] != INF)) {
        tmp2 = 0;

        if (ends[sn[i]] > (unsigned int)i)
          tmp2 += fms5[sn[i]][i + 1];

        if ((unsigned int)(j - 1) > ends[sn[i]])
          tmp2 += fms3[sn[ends[sn[i]] + 1]][j - 1];

        tmp = MIN2(tmp, tmp2);
      }

      /* check whether we find more strand nicks between i and j */
      nick = ends[sn[i]] + 1;
      while (sn[nick] != sn[j]) {
        /*
         *    if (evaluate(i + 1,
         *                 j - 1,
         *                 ends[sn[nick]],
         *                 ends[sn[nick]] + 1,
         *                 VRNA_DECOMP_EXT_EXT_EXT,
         *                 &hc_dat_local))
         */
        if ((fms5[sn[nick]][i + 1] != INF) &&
            (fms3[sn[ends[sn[nick]] + 1]][j - 1] != INF)) {
          tmp2 = 0;
          if ((unsigned int)(i + 1) <= ends[sn[nick]])
            tmp2 += fms5[sn[nick]][i + 1];

          if (ends[sn[nick]] + 1 <= (unsigned int)(j - 1))
            tmp2 += fms3[sn[ends[sn[nick]] + 1]][j - 1];

          tmp = MIN2(tmp, tmp2);
        }

        nick = ends[sn[nick]] + 1;
      }
    }

    if (tmp != INF)
      contribution = base + tmp;

    /* odd dangles below */
    if (dangle_model % 2) {
      s5  = (sn[j] == sn[j - 1]) ? S1[j - 1] : -1;
      s3  = (sn[i] == sn[i + 1]) ? S1[i + 1] : -1;

      if ((j > i + 1) &&
          (sn[i] != sn[i + 1]) &&
          (sn[j - 1] == sn[j])) {
        tmp = vrna_E_exterior_stem(type, s5, -1, params) +
              params->DuplexInit;

        if (sc_red_stem)
          tmp += sc_red_stem(j - 1, i, j, i, sc_wrapper);

        if (sn[j - 2] == sn[j]) {
          if (fms3[sn[i + 1]][j - 2] != INF) {
            tmp           += fms3[sn[i + 1]][j - 2];
            contribution  = MIN2(contribution, tmp);
          }
        } else {
          contribution = MIN2(contribution, tmp);
        }
      } else if ((i + 1 < j) &&
                 (sn[j - 1] != sn[j]) &&
                 (sn[i] == sn[i + 1])) {
        tmp = vrna_E_exterior_stem(type, -1, s3, params) +
              params->DuplexInit;

        if (sc_red_stem)
          tmp += sc_red_stem(j, i + 1, j, i, sc_wrapper);

        if (sn[i] == sn[i + 2]) {
          if (fms5[sn[j - 1]][i + 2] != INF) {
            tmp           += fms5[sn[j - 1]][i + 2];
            contribution  = MIN2(contribution, tmp);
          }
        } else {
          contribution = MIN2(contribution, tmp);
        }
      } else if ((sn[i] == sn[i + 1]) &&
                 (sn[j - 1] == sn[j])) {
        /* 1st case, i + 1 and j - 1 dangle on enclosing pair */
        tmp   = INF;
        base  = vrna_E_exterior_stem(type, s5, s3, params) +
                params->DuplexInit;

        if (sc_red_stem)
          base += sc_red_stem(j - 1, i + 1, j, i, sc_wrapper);

        start = i;
        nick  = ends[sn[start]] + 1;
        do {
          if ((fms5[sn[start]][i + 2] != INF) &&
              (fms3[sn[nick]][j - 2] != INF)) {
            tmp2 = 0;

            if (nick > (unsigned int)(i + 2))
              tmp2 += fms5[sn[start]][i + 2];

            if ((unsigned int)j > nick + 1)
              tmp2 += fms3[sn[nick]][j - 2];

            tmp = MIN2(tmp, tmp2);
          }

          start = nick;
          nick  = ends[sn[nick]] + 1;
        } while (sn[nick] != sn[j]);

        if (tmp != INF) {
          tmp           += base;
          contribution  = MIN2(contribution, tmp);
        }

        /* 2nd case, only i + 1 may dangle on enclosing pair */
        tmp   = INF;
        base  = vrna_E_exterior_stem(type, -1, s3, params) +
                params->DuplexInit;

        if (sc_red_stem)
          base += sc_red_stem(j, i + 1, j, i, sc_wrapper);

        start = i;
        nick  = ends[sn[start]] + 1;
        do {
          if ((fms5[sn[start]][i + 2] != INF) &&
              (fms3[sn[nick]][j - 1] != INF)) {
            tmp2 = 0;

            if (nick > (unsigned int)(i + 2))
              tmp2 += fms5[sn[start]][i + 2];

            if ((unsigned int)j > nick)
              tmp2 += fms3[sn[nick]][j - 1];

            tmp = MIN2(tmp, tmp2);
          }

          start = nick;
          nick  = ends[sn[nick]] + 1;
        } while (sn[nick] != sn[j]);

        if (tmp != INF) {
          tmp           += base;
          contribution  = MIN2(contribution, tmp);
        }

        /* 3rd case, only j - 1 may dangle on enclosing pair */
        tmp   = INF;
        base  = vrna_E_exterior_stem(type, s5, -1, params) +
                params->DuplexInit;

        if (sc_red_stem)
          base += sc_red_stem(j - 1, i, j, i, sc_wrapper);

        start = i;
        nick  = ends[sn[start]] + 1;
        do {
          if ((fms5[sn[start]][i + 1] != INF) &&
              (fms3[sn[nick]][j - 2] != INF)) {
            tmp2 = 0;

            if (nick > (unsigned int)(i + 1))
              tmp2 += fms5[sn[start]][i + 1];

            if ((unsigned int)j > nick + 1)
              tmp2 += fms3[sn[nick]][j - 2];

            tmp = MIN2(tmp, tmp2);
          }

          start = nick;
          nick  = ends[sn[nick]] + 1;
        } while (sn[nick] != sn[j]);

        if (tmp != INF) {
          tmp           += base;
          contribution  = MIN2(contribution, tmp);
        }
      }
    }
  }

  return contribution;
}


PRIVATE int
BT_multi_strand(vrna_fold_compound_t  *fc,
                int                   *i,
                int                   *j,
                unsigned int          *sn1,
                unsigned int          *sn2,
                int                   en,
                struct ms_helpers     *ms_dat)
{
  short                 *S1, *S2, s5, s3;
  unsigned int          *sn, *ends, type, nick;
  int                   **fms5, **fms3;
  int                   start, base, tmp, dangle_model;
  vrna_param_t          *params;
  vrna_md_t             *md;
  vrna_hc_eval_f        evaluate;
  struct hc_ext_def_dat *hc_dat_local;
  struct sc_f5_dat      *sc_wrapper;
  sc_ext_red_cb         sc_red_stem;

  if (fc) {
    S1            = fc->sequence_encoding;
    S2            = fc->sequence_encoding2;
    params        = fc->params;
    md            = &(params->model_details);
    dangle_model  = md->dangles;
    sn            = fc->strand_number;
    ends          = fc->strand_end;
    fms5          = fc->matrices->fms5;
    fms3          = fc->matrices->fms3;
    evaluate      = ms_dat->evaluate;
    hc_dat_local  = &(ms_dat->hc_dat_local);
    sc_wrapper    = &(ms_dat->sc_wrapper);
    sc_red_stem   = sc_wrapper->red_stem;

    if ((sn[*i] != sn[*j]) &&
        (evaluate(*i, *j, *i, *j, VRNA_DECOMP_EXT_STEM, hc_dat_local))) {
      /* most obious strand nick is at end of sn[i] and start of sn[j] */
      type = vrna_get_ptype_md(S2[*j], S2[*i], md);

      switch (dangle_model) {
        case 2:
          s5  = (sn[*j] == sn[*j - 1]) ? S1[*j - 1] : -1;
          s3  = (sn[*i] == sn[*i + 1]) ? S1[*i + 1] : -1;
          break;

        default:
          s5 = s3 = -1;
          break;
      }

      base = vrna_E_exterior_stem(type, s5, s3, params) +
             params->DuplexInit;

      if (sc_red_stem)
        base += sc_red_stem(*j, *i, *j, *i, sc_wrapper);

      if (sn[*i] != sn[*i + 1]) {
        if ((sn[*j - 1] != sn[*j]) &&
            (*i + 1 == *j)) {
          tmp = 0;
          if (tmp + base == en) {
            *sn1  = 0;
            *sn2  = 0;
            *i    = 0;
            *j    = 0;
            return 1;
          }
        } else if (sn[*j - 1] == sn[*j]) {
          tmp = fms3[sn[*i + 1]][*j - 1];
          if (tmp + base == en) {
            *sn1  = 0;
            *sn2  = sn[*i + 1];
            *i    = 0;
            *j    = *j - 1;
            return 1;
          }
        }
      } else if (sn[*j - 1] != sn[*j]) {
        tmp = fms5[sn[*j - 1]][*i + 1];
        if (tmp + base == en) {
          *sn1  = sn[*j - 1];
          *sn2  = 0;
          *i    = *i + 1;
          *j    = 0;
          return 1;
        }
      } else {
        tmp = 0;

        if (ends[sn[*i]] > (unsigned int)(*i))
          tmp += fms5[sn[*i]][*i + 1];

        if ((unsigned int)(*j - 1) > ends[sn[*i]])
          tmp += fms3[sn[ends[sn[*i]] + 1]][*j - 1];

        if (tmp + base == en) {
          *sn1  = sn[*i];
          *sn2  = sn[ends[sn[*i]] + 1];
          *i    = (ends[sn[*i]] > (unsigned int)(*i)) ? *i + 1 : 0;
          *j    = ((unsigned int)(*j - 1) > ends[sn[*i]]) ? *j - 1 : 0;
          return 1;
        }

        /* check whether we find more strand nicks between i and j */
        nick = ends[sn[*i]] + 1;
        while (sn[nick] != sn[*j]) {
          tmp = 0;
          if ((unsigned int)(*i + 1) <= ends[sn[nick]])
            tmp += fms5[sn[nick]][*i + 1];

          if (ends[sn[nick]] + 1 <= (unsigned int)(*j - 1))
            tmp += fms3[sn[ends[sn[nick]] + 1]][*j - 1];

          if (tmp + base == en) {
            *sn1  = sn[nick];
            *sn2  = sn[ends[sn[nick]] + 1];
            *i    = ((unsigned int)(*i + 1) <= ends[sn[nick]]) ? *i + 1 : 0;
            *j    = (ends[sn[nick]] + 1 <= (unsigned int)(*j - 1)) ? *j - 1 : 0;
            return 1;
          }

          nick = ends[sn[nick]] + 1;
        }
      }

      /* odd dangle cases below */
      if (dangle_model % 2) {
        s5  = (sn[*j] == sn[*j - 1]) ? S1[*j - 1] : -1;
        s3  = (sn[*i] == sn[*i + 1]) ? S1[*i + 1] : -1;

        if ((*j > *i + 1) &&
            (sn[*i] != sn[*i + 1]) &&
            (sn[*j - 1] == sn[*j])) {
          tmp = vrna_E_exterior_stem(type, s5, -1, params) +
                params->DuplexInit;

          if (sc_red_stem)
            tmp += sc_red_stem(*j - 1, *i, *j, *i, sc_wrapper);

          if (sn[*j - 2] == sn[*j]) {
            if (fms3[sn[*i + 1]][*j - 2] != INF) {
              tmp += fms3[sn[*i + 1]][*j - 2];

              if (tmp == en) {
                *sn1  = 0;
                *sn2  = sn[*i + 1];
                *i    = 0;
                *j    = *j - 2;
                return 1;
              }
            }
          } else if (tmp == en) {
            *sn1  = 0;
            *sn2  = 0;
            *i    = 0;
            *j    = 0;
            return 1;
          }
        } else if ((*i + 1 < *j) &&
                   (sn[*j - 1] != sn[*j]) &&
                   (sn[*i] == sn[*i + 1])) {
          tmp = vrna_E_exterior_stem(type, -1, s3, params) +
                params->DuplexInit;

          if (sc_red_stem)
            tmp += sc_red_stem(*j, *i + 1, *j, *i, sc_wrapper);

          if (sn[*i] == sn[*i + 2]) {
            if (fms5[sn[*j - 1]][*i + 2] != INF) {
              tmp += fms5[sn[*j - 1]][*i + 2];

              if (tmp == en) {
                *sn1  = sn[*j - 1];
                *sn2  = 0;
                *i    = *i + 2;
                *j    = 0;
                return 1;
              }
            }
          } else if (tmp == en) {
            *sn1  = 0;
            *sn2  = 0;
            *i    = 0;
            *j    = 0;
            return 1;
          }
        } else if ((sn[*i] == sn[*i + 1]) &&
                   (sn[*j - 1] == sn[*j])) {
          /* 1st case, i + 1 and j - 1 dangle on enclosing pair */
          base = vrna_E_exterior_stem(type, s5, s3, params) +
                 params->DuplexInit;

          if (sc_red_stem)
            base += sc_red_stem(*j - 1, *i + 1, *j, *i, sc_wrapper);

          start = *i;
          nick  = ends[sn[start]] + 1;
          do {
            if ((fms5[sn[start]][*i + 2] != INF) &&
                (fms3[sn[nick]][*j - 2] != INF)) {
              tmp = 0;

              if (nick > (unsigned int)(*i + 2))
                tmp += fms5[sn[start]][*i + 2];

              if ((unsigned int)(*j) > nick + 1)
                tmp += fms3[sn[nick]][*j - 2];

              if (tmp + base == en) {
                *sn1  = sn[start];
                *sn2  = sn[nick];
                *i    = (nick > (unsigned int)(*i + 2)) ? *i + 2 : 0;
                *j    = ((unsigned int)(*j) > nick + 1) ? *j - 2 : 0;
                return 1;
              }
            }

            start = nick;
            nick  = ends[sn[nick]] + 1;
          } while (sn[nick] != sn[*j]);

          /* 2nd case, i + 1 may dangle on enclosing pair */
          base = vrna_E_exterior_stem(type, -1, s3, params) +
                 params->DuplexInit;

          if (sc_red_stem)
            tmp += sc_red_stem(*j, *i + 1, *j, *i, sc_wrapper);

          start = *i;
          nick  = ends[sn[start]] + 1;
          do {
            if ((fms5[sn[start]][*i + 2] != INF) &&
                (fms3[sn[nick]][*j - 1] != INF)) {
              tmp = 0;

              if (nick > (unsigned int)(*i + 2))
                tmp += fms5[sn[start]][*i + 2];

              if ((unsigned int)(*j) > nick)
                tmp += fms3[sn[nick]][*j - 1];

              if (tmp + base == en) {
                *sn1  = sn[start];
                *sn2  = sn[nick];
                *i    = (nick > (unsigned int)(*i + 2)) ? *i + 2 : 0;
                *j    = ((unsigned int)(*j) > nick) ? *j - 1 : 0;
                return 1;
              }
            }

            start = nick;
            nick  = ends[sn[nick]] + 1;
          } while (sn[nick] != sn[*j]);

          /* 3rd case, j - 1 may dangle on enclosing pair */
          base = vrna_E_exterior_stem(type, s5, -1, params) +
                 params->DuplexInit;

          if (sc_red_stem)
            tmp += sc_red_stem(*j - 1, *i, *j, *i, sc_wrapper);

          start = *i;
          nick  = ends[sn[start]] + 1;
          do {
            if ((fms5[sn[start]][*i + 1] != INF) &&
                (fms3[sn[nick]][*j - 2] != INF)) {
              tmp = 0;

              if (nick > (unsigned int)(*i + 1))
                tmp += fms5[sn[start]][*i + 1];

              if ((unsigned int)(*j) > nick + 1)
                tmp += fms3[sn[nick]][*j - 2];

              if (tmp + base == en) {
                *sn1  = sn[start];
                *sn2  = sn[nick];
                *i    = (nick > (unsigned int)(*i + 1)) ? *i + 1 : 0;
                *j    = ((unsigned int)(*j) > nick + 1) ? *j - 2 : 0;
                return 1;
              }
            }

            start = nick;
            nick  = ends[sn[nick]] + 1;
          } while (sn[nick] != sn[*j]);
        }
      }
    }
  }

  return 0;
}


PRIVATE int
BT_fms5_split(vrna_fold_compound_t  *fc,
              unsigned int          strand,
              int                   *i,
              int                   *k,
              int                   *l,
              struct ms_helpers     *ms_dat)
{
  short                 *S1, *S2, s5, s3;
  unsigned int          *sn, *se, type;
  int                   u, *idx, end, *c, **fms5, base, tmp, dangle_model;
  vrna_param_t          *params;
  vrna_md_t             *md;
  vrna_hc_eval_f        evaluate;
  struct hc_ext_def_dat *hc_dat_local;
  struct sc_f5_dat      *sc_wrapper;
  sc_ext_red_cb         sc_spl;
  sc_ext_red_cb         sc_red_stem;
  sc_ext_red_cb         sc_red_ext;

  if (!ms_dat)
    return 0;

  S1            = fc->sequence_encoding;
  S2            = fc->sequence_encoding2;
  sn            = fc->strand_number;
  se            = fc->strand_end;
  end           = (int)se[strand];
  idx           = fc->jindx;
  params        = fc->params;
  md            = &(params->model_details);
  dangle_model  = md->dangles;
  c             = fc->matrices->c;
  fms5          = fc->matrices->fms5;
  evaluate      = ms_dat->evaluate;
  hc_dat_local  = &(ms_dat->hc_dat_local);
  sc_wrapper    = &(ms_dat->sc_wrapper);
  sc_spl        = sc_wrapper->decomp;
  sc_red_stem   = sc_wrapper->red_stem;
  sc_red_ext    = sc_wrapper->red_ext;

  if (*i == end) {
    *i  = 0;
    *k  = 0;
    *l  = 0;
    return 1;
  }

  if (evaluate(*i, end, *i + 1, end, VRNA_DECOMP_EXT_EXT, hc_dat_local)) {
    tmp = fms5[strand][*i + 1];

    if (sc_red_ext)
      tmp += sc_red_ext(*i, end, *i + 1, end, sc_wrapper);

    if (fms5[strand][*i] == tmp) {
      *i  = *i + 1;
      *k  = 0;
      *l  = 0;
      return 1;
    }
  }

  if (evaluate(*i, end, *i, end, VRNA_DECOMP_EXT_STEM, hc_dat_local)) {
    type = vrna_get_ptype_md(S2[*i], S2[end], md);

    switch (dangle_model) {
      case 2:
        s5  = ((*i > 1) && (sn[*i - 1] == sn[*i])) ? S1[*i - 1] : -1;
        s3  = -1;
        break;

      default:
        s5 = s3 = -1;
        break;
    }

    base = vrna_E_exterior_stem(type, s5, s3, params);

    if (sc_red_stem)
      base += sc_red_stem(*i, end, *i, end, sc_wrapper);

    if (fms5[strand][*i] == c[idx[end] + *i] + base) {
      *k  = end;
      *l  = end;
      return 1;
    }
  }

  for (u = *i + 1; u < end; u++) {
    if (evaluate(*i, end, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT, hc_dat_local)) {
      type = vrna_get_ptype_md(S2[*i], S2[u], md);

      switch (dangle_model) {
        case 2:
          s5  = ((*i > 1) && (sn[*i - 1] == sn[*i])) ? S1[*i - 1] : -1;
          s3  = (sn[u] == sn[u + 1]) ? S1[u + 1] : -1;
          break;

        default:
          s5 = s3 = -1;
          break;
      }

      tmp = fms5[strand][u + 1] +
            c[idx[u] + *i] +
            vrna_E_exterior_stem(type, s5, s3, params);

      if (sc_red_stem)
        tmp += sc_red_stem(*i, u, *i, u, sc_wrapper);

      if (sc_spl)
        tmp += sc_spl(*i, end, u, u + 1, sc_wrapper);

      if (tmp == fms5[strand][*i]) {
        *k  = u;
        *l  = u + 1;
        return 1;
      }
    }
  }

  if (dangle_model % 2) {
    s5  = S1[*i];
    s3  = S1[end];

    if (evaluate(*i, end, *i + 1, end, VRNA_DECOMP_EXT_STEM, hc_dat_local)) {
      type  = vrna_get_ptype_md(S2[*i + 1], S2[end], md);
      base  = vrna_E_exterior_stem(type, s5, -1, params);

      if (sc_red_stem)
        base += sc_red_stem(*i, end, *i + 1, end, sc_wrapper);

      if (fms5[strand][*i] == c[idx[end] + *i + 1] + base) {
        *i  = *i + 1;
        *k  = end;
        *l  = end;
        return 1;
      }
    }

    if (evaluate(*i, end, *i, end - 1, VRNA_DECOMP_EXT_STEM, hc_dat_local)) {
      type  = vrna_get_ptype_md(S2[*i], S2[end - 1], md);
      base  = vrna_E_exterior_stem(type, -1, s3, params);

      if (sc_red_stem)
        base += sc_red_stem(*i, end, *i, end - 1, sc_wrapper);

      if (fms5[strand][*i] == c[idx[end - 1] + *i] + base) {
        *k  = end - 1;
        *l  = end;
        return 1;
      }
    }

    if (evaluate(*i, end, *i + 1, end - 1, VRNA_DECOMP_EXT_STEM, hc_dat_local)) {
      type  = vrna_get_ptype_md(S2[*i + 1], S2[end - 1], md);
      base  = vrna_E_exterior_stem(type, s5, s3, params);

      if (sc_red_stem)
        base += sc_red_stem(*i, end, *i + 1, end - 1, sc_wrapper);

      if (fms5[strand][*i] == c[idx[end - 1] + *i + 1] + base) {
        *i  = *i + 1;
        *k  = end - 1;
        *l  = end;
        return 1;
      }
    }

    for (u = *i + 1; u < end; u++) {
      if (evaluate(*i, end, u, u + 1, VRNA_DECOMP_EXT_STEM_EXT1, hc_dat_local)) {
        type  = vrna_get_ptype_md(S2[*i + 1], S2[u], md);
        tmp   = fms5[strand][u + 1] +
                c[idx[u] + *i + 1] +
                vrna_E_exterior_stem(type, s5, -1, params);

        if (sc_red_stem)
          tmp += sc_red_stem(*i, u, *i + 1, u, sc_wrapper);

        if (sc_spl)
          tmp += sc_spl(*i, end, u, u + 1, sc_wrapper);

        if (tmp == fms5[strand][*i]) {
          *i  = *i + 1;
          *k  = u;
          *l  = u + 1;
          return 1;
        }
      }
    }

    for (u = *i + 1; u + 1 < end; u++) {
      s3 = (sn[u] == sn[u + 1]) ? S1[u + 1] : -1;

      if (evaluate(*i, end, u, u + 2, VRNA_DECOMP_EXT_STEM_EXT, hc_dat_local)) {
        type  = vrna_get_ptype_md(S2[*i], S2[u], md);
        tmp   = fms5[strand][u + 2] +
                c[idx[u] + *i] +
                vrna_E_exterior_stem(type, -1, s3, params);

        if (sc_red_stem)
          tmp += sc_red_stem(*i, u + 1, *i, u, sc_wrapper);

        if (sc_spl)
          tmp += sc_spl(*i, end, u + 1, u + 2, sc_wrapper);

        if (tmp == fms5[strand][*i]) {
          *k  = u;
          *l  = u + 2;
          return 1;
        }
      }

      if (evaluate(*i, end, u, u + 2, VRNA_DECOMP_EXT_STEM_EXT1, hc_dat_local)) {
        type  = vrna_get_ptype_md(S2[*i + 1], S2[u], md);
        tmp   = fms5[strand][u + 2] +
                c[idx[u] + *i + 1] +
                vrna_E_exterior_stem(type, s5, s3, params);

        if (sc_red_stem)
          tmp += sc_red_stem(*i, u + 1, *i + 1, u, sc_wrapper);

        if (sc_spl)
          tmp += sc_spl(*i, end, u + 1, u + 2, sc_wrapper);

        if (tmp == fms5[strand][*i]) {
          *i  = *i + 1;
          *k  = u;
          *l  = u + 2;
          return 1;
        }
      }
    }
  }

  return 0;
}


PRIVATE int
BT_fms3_split(vrna_fold_compound_t  *fc,
              unsigned int          strand,
              int                   *j,
              int                   *k,
              int                   *l,
              struct ms_helpers     *ms_dat)
{
  short                 *S1, *S2, s5, s3;
  unsigned int          *sn, *ss, type;
  int                   u, *idx, start, n, *c, **fms3, base, dangle_model;
  vrna_param_t          *params;
  vrna_md_t             *md;
  vrna_hc_eval_f        evaluate;
  struct hc_ext_def_dat *hc_dat_local;
  struct sc_f5_dat      *sc_wrapper;
  sc_ext_red_cb         sc_spl;
  sc_ext_red_cb         sc_red_stem;
  sc_ext_red_cb         sc_red_ext;

  n             = fc->length;
  S1            = fc->sequence_encoding;
  S2            = fc->sequence_encoding2;
  sn            = fc->strand_number;
  ss            = fc->strand_start;
  start         = (int)ss[strand];
  idx           = fc->jindx;
  params        = fc->params;
  md            = &(params->model_details);
  dangle_model  = md->dangles;
  c             = fc->matrices->c;
  fms3          = fc->matrices->fms3;
  evaluate      = ms_dat->evaluate;
  hc_dat_local  = &(ms_dat->hc_dat_local);
  sc_wrapper    = &(ms_dat->sc_wrapper);
  sc_spl        = sc_wrapper->decomp;
  sc_red_stem   = sc_wrapper->red_stem;
  sc_red_ext    = sc_wrapper->red_ext;

  if (*j == start) {
    *j  = 0;
    *k  = 0;
    *l  = 0;
    return 1;
  }

  if (evaluate(start, *j, start, *j - 1, VRNA_DECOMP_EXT_EXT, hc_dat_local)) {
    base = fms3[strand][*j - 1];

    if (sc_red_ext)
      base += sc_red_ext(start, *j, start, *j - 1, sc_wrapper);

    if (fms3[strand][*j] == base) {
      *j  = *j - 1;
      *k  = 0;
      *l  = 0;
      return 1;
    }
  }

  if (evaluate(start, *j, start, *j, VRNA_DECOMP_EXT_STEM, hc_dat_local)) {
    type = vrna_get_ptype_md(S2[start], S2[*j], md);

    switch (dangle_model) {
      case 2:
        s5  = -1;
        s3  = ((*j < n) && (sn[*j] == sn[*j + 1])) ? S1[*j + 1] : -1;
        break;

      default:
        s5 = s3 = -1;
        break;
    }

    base = vrna_E_exterior_stem(type, s5, s3, params);

    if (sc_red_stem)
      base += sc_red_stem(start, *j, start, *j, sc_wrapper);

    if (fms3[strand][*j] == c[idx[*j] + start] + base) {
      *k  = start;
      *l  = start;
      return 1;
    }
  }

  for (u = start; u < *j; u++) {
    if (evaluate(start, *j, u, u + 1, VRNA_DECOMP_EXT_EXT_STEM, hc_dat_local)) {
      type = vrna_get_ptype_md(S2[u + 1], S2[*j], md);

      switch (dangle_model) {
        case 2:
          s5  = (sn[u] == sn[u + 1]) ? S1[u] : -1;
          s3  = ((*j < n) && (sn[*j] == sn[*j + 1])) ? S1[*j + 1] : -1;
          break;

        default:
          s5 = s3 = -1;
          break;
      }

      base = vrna_E_exterior_stem(type, s5, s3, params);

      if (sc_red_stem)
        base += sc_red_stem(u + 1, *j, u + 1, *j, sc_wrapper);

      if (sc_spl)
        base += sc_spl(start, *j, u, u + 1, sc_wrapper);

      if (fms3[strand][*j] == fms3[strand][u] + c[idx[*j] + u + 1] + base) {
        *k  = u + 1;
        *l  = u;
        return 1;
      }
    }
  }

  if (dangle_model % 2) {
    s5  = S1[start];
    s3  = S1[*j];

    if (evaluate(start, *j, start + 1, *j, VRNA_DECOMP_EXT_STEM, hc_dat_local)) {
      type = vrna_get_ptype_md(S2[start + 1], S2[*j], md);

      base = vrna_E_exterior_stem(type, s5, -1, params);

      if (sc_red_stem)
        base += sc_red_stem(start, *j, start + 1, *j, sc_wrapper);

      if (fms3[strand][*j] == c[idx[*j] + start + 1] + base) {
        *k  = start + 1;
        *l  = start;
        return 1;
      }
    }

    if (evaluate(start, *j, start, *j - 1, VRNA_DECOMP_EXT_STEM, hc_dat_local)) {
      type  = vrna_get_ptype_md(S2[start], S2[*j - 1], md);
      base  = vrna_E_exterior_stem(type, -1, s3, params);

      if (sc_red_stem)
        base += sc_red_stem(start, *j, start, *j - 1, sc_wrapper);

      if (fms3[strand][*j] == c[idx[*j - 1] + start] + base) {
        *j  = *j - 1;
        *k  = start;
        *l  = start;
        return 1;
      }
    }

    if (evaluate(start, *j, start + 1, *j - 1, VRNA_DECOMP_EXT_STEM, hc_dat_local)) {
      type  = vrna_get_ptype_md(S2[start + 1], S2[*j - 1], md);
      base  = vrna_E_exterior_stem(type, s5, s3, params);

      if (sc_red_stem)
        base += sc_red_stem(start, *j, start + 1, *j - 1, sc_wrapper);

      if (fms3[strand][*j] == c[idx[*j - 1] + start + 1] + base) {
        *j  = *j - 1;
        *k  = start + 1;
        *l  = start;
        return 1;
      }
    }

    for (u = start; u < *j; u++) {
      s5 = S1[u + 1];

      if (evaluate(start, *j, u, u + 1, VRNA_DECOMP_EXT_EXT_STEM1, hc_dat_local)) {
        type  = vrna_get_ptype_md(S2[u + 1], S2[*j - 1], md);
        base  = vrna_E_exterior_stem(type, -1, s3, params);

        if (sc_red_stem)
          base += sc_red_stem(u + 1, *j, u + 1, *j - 1, sc_wrapper);

        if (sc_spl)
          base += sc_spl(start, *j, u, u + 1, sc_wrapper);

        if (fms3[strand][*j] == fms3[strand][u] + c[idx[*j - 1] + u + 1] + base) {
          *j  = *j - 1;
          *k  = u + 1;
          *l  = u;
          return 1;
        }
      }

      if (evaluate(start, *j, u, u + 2, VRNA_DECOMP_EXT_EXT_STEM, hc_dat_local)) {
        type  = vrna_get_ptype_md(S2[u + 2], S2[*j], md);
        base  = vrna_E_exterior_stem(type, s5, -1, params);

        if (sc_red_stem)
          base += sc_red_stem(u + 1, *j, u + 2, *j, sc_wrapper);

        if (sc_spl)
          base += sc_spl(start, *j, u, u + 1, sc_wrapper);

        if (fms3[strand][*j] == fms3[strand][u] + c[idx[*j] + u + 2] + base) {
          *k  = u + 2;
          *l  = u;
          return 1;
        }
      }

      if (evaluate(start, *j, u, u + 2, VRNA_DECOMP_EXT_EXT_STEM1, hc_dat_local)) {
        type  = vrna_get_ptype_md(S2[u + 2], S2[*j - 1], md);
        base  = vrna_E_exterior_stem(type, s5, s3, params);

        if (sc_red_stem)
          base += sc_red_stem(u + 1, *j, u + 2, *j - 1, sc_wrapper);

        if (sc_spl)
          base += sc_spl(start, *j, u, u + 1, sc_wrapper);

        if (fms3[strand][*j] == fms3[strand][u] + c[idx[*j - 1] + u + 2] + base) {
          *j  = *j - 1;
          *k  = u + 2;
          *l  = u;
          return 1;
        }
      }
    }
  }

  return 0;
}


#endif


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
          vrna_bps_t            bp_stack,
          vrna_bts_t            bt_stack,
          struct ms_helpers     *ms_dat)
{
  char          backtrack_type;
  unsigned int  L, ll[3];
  int           i, j, ij, k, l, length, *my_c, *indx, noLP, *pscore, ret;
  vrna_param_t  *P;
  vrna_gr_aux_t aux_grammar;

  ret             = 1;
  length          = fc->length;
  my_c            = fc->matrices->c;
  indx            = fc->jindx;
  P               = fc->params;
  noLP            = P->model_details.noLP;
  pscore          = fc->pscore;         /* covariance scores for comparative structure prediction */
  backtrack_type  = P->model_details.backtrack_type;
  aux_grammar     = fc->aux_grammar;


  if (vrna_bts_size(bt_stack) == 0) {
    vrna_bts_push(bt_stack, ((vrna_sect_t){
      .i = 1,
      .j = length,
      .ml = (backtrack_type == 'M') ? VRNA_MX_FLAG_M : ((backtrack_type == 'C') ? VRNA_MX_FLAG_C : VRNA_MX_FLAG_F5)}));
  }

  while (vrna_bts_size(bt_stack) > 0) {
    int         ml, cij, canonical;
    vrna_sect_t e;

    /* pop one element from stack */
    canonical = 1;
    e         = vrna_bts_pop(bt_stack);
    i         = e.i;
    j         = e.j;
    ml        = e.ml;

    switch (ml) {
      /* backtrack in f5 */
      case VRNA_MX_FLAG_F5:
      {
        if (vrna_bt_f(fc, 1, j, bp_stack, bt_stack)) {
          continue;
        } else {
          vrna_log_warning("backtracking failed in f, segment [%d,%d], e = %d\n",
                           i,
                           j,
                           fc->matrices->f5[j]);
          ret = 0;
          goto backtrack_exit;
        }
      }
      break;

      /* trace back in fML array */
      case VRNA_MX_FLAG_M:
      {
        if (vrna_bt_m(fc, i, j, bp_stack, bt_stack)) {
          continue;
        } else {
          vrna_log_warning("backtracking failed in fML, segment [%d,%d]",
                           i,
                           j);
          ret = 0;
          goto backtrack_exit;
        }
      }
      break;

      /* backtrack in c */
      case VRNA_MX_FLAG_C:
        vrna_bps_push(bp_stack,
                      (vrna_bp_t){
                        .i = i,
                        .j = j
                      });
        goto repeat1;

        break;

      /* backtrack in fms5 */
      case VRNA_MX_FLAG_MS5:
      {
        unsigned int strand = j;

        if (BT_fms5_split(fc, strand, &i, &k, &l, ms_dat)) {
          if (k > 0) {
            vrna_bts_push(bt_stack, ((vrna_sect_t){
              .i = i,
              .j = k,
              .ml = VRNA_MX_FLAG_C}));

            if ((unsigned int)k < fc->strand_end[strand]) {
              vrna_bts_push(bt_stack, ((vrna_sect_t){
                .i = l,
                .j = strand,
                .ml = VRNA_MX_FLAG_MS5}));
            }
          } else if (i > 0) {
            vrna_bts_push(bt_stack, ((vrna_sect_t){
              .i = i,
              .j = strand,
              .ml = VRNA_MX_FLAG_MS5}));
          }

          continue;
        } else {
          vrna_log_warning("backtracking failed in fsm5[%d][%d] (%d:%d)\n",
                           strand,
                           i,
                           fc->strand_start[strand],
                           fc->strand_end[strand]);
          ret = 0;
          goto backtrack_exit;
        }
      }
      break;

      /* backtrack in fms3 */
      case VRNA_MX_FLAG_MS3:
      {
        unsigned int strand = i;

        if (BT_fms3_split(fc, strand, &j, &k, &l, ms_dat)) {
          if (k > 0) {
            vrna_bts_push(bt_stack, ((vrna_sect_t){
              .i = k,
              .j = j,
              .ml = VRNA_MX_FLAG_C}));

            if ((unsigned int)k > fc->strand_start[strand]) {
              vrna_bts_push(bt_stack, ((vrna_sect_t){
                .i = strand,
                .j = l,
                .ml = VRNA_MX_FLAG_MS3}));
            }
          } else if (j > 0) {
            vrna_bts_push(bt_stack, ((vrna_sect_t){
              .i = strand,
              .j = j,
              .ml = VRNA_MX_FLAG_MS3}));
          }

          continue;
        } else {
          vrna_log_warning("backtracking failed in fsm3[%d][%d] (%d:%d)\n",
                           strand,
                           j,
                           fc->strand_start[strand],
                           fc->strand_end[strand]);
          ret = 0;
          goto backtrack_exit;
        }
      }
      break;

      case VRNA_MX_FLAG_G:
        if (vrna_bt_gquad(fc, i, j, &L, ll)) {
          vrna_bps_push(bp_stack,
                        (vrna_bp_t){
                          .i = i,
                          .j = i,
                          .L = L,
                          .l = { ll[0], ll[1], ll[2] }
                        });
          continue;
        } else {
          vrna_log_warning("backtracking failed in G, segment [%d,%d]\n", i, j);
        }

        break;

      default:
        /* catch auxiliary grammar backtracking here */
        if ((aux_grammar) &&
            ((unsigned int)ml > VRNA_MX_FLAG_MAX)) {
          unsigned int flag = ml -
                              VRNA_MX_FLAG_MAX -
                              1;
          if ((flag < vrna_array_size(aux_grammar->aux)) &&
              (aux_grammar->aux[flag].cb_bt)) {
            if (aux_grammar->aux[flag].cb_bt(fc, i, j, INF, bp_stack, bt_stack,
                                             aux_grammar->aux[flag].data)) {
              continue;
            } else {
              vrna_log_warning(
                "backtracking failed in auxiliary grammar backtrack %u, segment [%d, %d]\n",
                ml,
                i,
                j);
            }
          } else {
            vrna_log_error(
              "backtracking requested but unavailable for auxiliary grammar %u, segment [%d, %d]\n",
              ml,
              i,
              j);
          }
        }

        ret = 0;
        goto backtrack_exit;
    }

repeat1:

    /*----- begin of "repeat:" -----*/
    ij = indx[j] + i;

    if (canonical)
      cij = my_c[ij];

    if (noLP) {
      if (vrna_bt_stacked_pairs(fc, i, j, &cij, bp_stack, bt_stack)) {
        canonical = 0;
        /* remove enclosed element from backtrack stack to
         * allow for immediately going back to repeat1
         */
        vrna_sect_t trash = vrna_bts_pop(bt_stack);
        i = trash.i;
        j = trash.j;

        goto repeat1;
      }
    }

    canonical = 1;

    if (fc->type == VRNA_FC_TYPE_COMPARATIVE)
      cij += pscore[indx[j] + i];

    if (vrna_bt_hairpin(fc, i, j, cij, bp_stack, bt_stack))
      continue;

    if (vrna_bt_internal_loop(fc, i, j, cij, bp_stack, bt_stack))
      continue;

    if (fc->strands > 1) {
      unsigned int sn1, sn2;

      if ((ms_dat) &&
          (BT_multi_strand(fc, &i, &j, &sn1, &sn2, cij, ms_dat))) {
        if (i > 0) {
          vrna_bts_push(bt_stack, ((vrna_sect_t){
            .i = i,
            .j = sn1,
            .ml = VRNA_MX_FLAG_MS5}));
        }

        if (j > 0) {
          vrna_bts_push(bt_stack, ((vrna_sect_t){
            .i = sn2,
            .j = j,
            .ml = VRNA_MX_FLAG_MS3}));
        }

        continue;
      }
    }

    /* (i.j) must close a multi-loop */
    if (vrna_bt_multibranch_loop(fc, i, j, cij, bp_stack, bt_stack)) {
      continue;
    } else if (aux_grammar) {
      ret = 0;
      /* go through each user-provided backtrack callback and try finding the solution */
      for (size_t c = 0; c < vrna_array_size(aux_grammar->c); c++)
        if ((aux_grammar->c[c].cb_bt) &&
            (ret = aux_grammar->c[c].cb_bt(fc, i, j, cij, bp_stack, bt_stack, aux_grammar->c[c].data)))
          break;

      if (ret)
        continue;
    } else {
      vrna_log_warning("backtracking failed in repeat, segment [%d,%d]\n", i, j);
      ret = 0;
      goto backtrack_exit;
    }

    /* end of repeat: --------------------------------------------------*/
  } /* end of infinite while loop */

backtrack_exit:

  return ret;
}


PRIVATE INLINE int
decompose_pair(vrna_fold_compound_t *fc,
               int                  i,
               int                  j,
               struct aux_arrays    *aux,
               struct ms_helpers    *ms_dat)
{
  unsigned char hc_decompose;
  unsigned int  n;
  int           e, new_c, energy, stackEnergy, ij, dangle_model, noLP, *cc, *cc1;

  n             = fc->length;
  ij            = fc->jindx[j] + i;
  dangle_model  = fc->params->model_details.dangles;
  noLP          = fc->params->model_details.noLP;
  hc_decompose  = fc->hc->mx[n * i + j];
  cc            = aux->cc;
  cc1           = aux->cc1;
  e             = INF;

  /* do we evaluate this pair? */
  if (hc_decompose) {
    new_c = INF;

    /* check for hairpin loop */
    energy  = vrna_eval_hairpin(fc, i, j, VRNA_EVAL_LOOP_DEFAULT);
    new_c   = MIN2(new_c, energy);

    /* check for multibranch loops */
    energy  = vrna_mfe_multibranch_loop_fast(fc, i, j, aux->ml_helpers);
    new_c   = MIN2(new_c, energy);

    if (dangle_model == 3) {
      /* coaxial stacking */
      energy  = vrna_mfe_multibranch_loop_stack(fc, i, j);
      new_c   = MIN2(new_c, energy);
    }

    /* check for internal loops */
    energy  = vrna_mfe_internal(fc, i, j);
    new_c   = MIN2(new_c, energy);

    /* multi-strand decomposition */
    if (fc->strands > 1) {
      energy  = pair_multi_strand(fc, i, j, ms_dat);
      new_c   = MIN2(new_c, energy);
    }

    /* remember stack energy for --noLP option */
    if (noLP) {
      stackEnergy = vrna_eval_stack(fc, i, j, VRNA_EVAL_LOOP_DEFAULT);
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
    if (fc->aux_grammar) {
      for (size_t c = 0; c < vrna_array_size(fc->aux_grammar->c); c++) {
        if (fc->aux_grammar->c[c].cb) {
          energy  = fc->aux_grammar->c[c].cb(fc, i, j, fc->aux_grammar->c[c].data);
          e       = MIN2(e, energy);
        }
      }
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
  struct aux_arrays *aux = (struct aux_arrays *)vrna_alloc(sizeof(struct aux_arrays));

  aux->cc   = (int *)vrna_alloc(sizeof(int) * (length + 2));    /* auxilary arrays for canonical structures     */
  aux->cc1  = (int *)vrna_alloc(sizeof(int) * (length + 2));    /* auxilary arrays for canonical structures     */

  /* get fast multiloop decomposition helpers */
  aux->ml_helpers = vrna_mfe_multibranch_fast_init(length);

  return aux;
}


PRIVATE INLINE void
rotate_aux_arrays(struct aux_arrays *aux,
                  unsigned int      length)
{
  unsigned int  j;
  int           *FF;

  FF        = aux->cc1;
  aux->cc1  = aux->cc;
  aux->cc   = FF;
  for (j = 1; j <= length; j++)
    aux->cc[j] = INF;

  vrna_mfe_multibranch_fast_rotate(aux->ml_helpers);
}


PRIVATE INLINE void
free_aux_arrays(struct aux_arrays *aux)
{
  free(aux->cc);
  free(aux->cc1);
  vrna_mfe_multibranch_fast_free(aux->ml_helpers);

  free(aux);
}
