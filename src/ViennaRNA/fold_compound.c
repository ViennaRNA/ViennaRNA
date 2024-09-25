/** \file fold_compound.c **/

/*
 *                Fold Compound API
 *
 *                This file contains everything necessary to create
 *                and destroy data structures of type vrna_fold_compound_s,
 *                the most fundamental data structure throughout RNAlib
 *
 *                c Ronny Lorenz
 *
 *                ViennaRNA package
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
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/sequences/alignments.h"
#include "ViennaRNA/params/ribosum.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/partfunc/global.h"
#include "ViennaRNA/cofold.h"
#include "ViennaRNA/mm.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/structures/utils.h"
#include "ViennaRNA/fold_compound.h"


#ifndef INLINE
#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif
#endif

/*
 #################################
 # PRIVATE MACROS                #
 #################################
 */

#define WITH_PTYPE          1L    /* passed to set_fold_compound() to indicate that we need to set fc->ptype */
#define WITH_PTYPE_COMPAT   2L    /* passed to set_fold_compound() to indicate that we need to set fc->ptype_compat */

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
PRIVATE void
set_fold_compound(vrna_fold_compound_t  *fc,
                  unsigned int          options,
                  unsigned int          aux);


PRIVATE void
make_pscores(vrna_fold_compound_t *fc);


PRIVATE void
sanitize_bp_span(vrna_fold_compound_t *fc,
                 unsigned int         options);


PRIVATE void
add_params(vrna_fold_compound_t *fc,
           vrna_md_t            *md_p,
           unsigned int         options);


PRIVATE vrna_fold_compound_t *
init_fc_single(void);


PRIVATE vrna_fold_compound_t *
init_fc_comparative(void);


PRIVATE INLINE void
nullify(vrna_fold_compound_t *fc);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC void
vrna_fold_compound_free(vrna_fold_compound_t *fc)
{
  unsigned int s;

  if (fc) {
    /* first destroy common attributes */
    vrna_mx_mfe_free(fc);
    vrna_mx_pf_free(fc);
    free(fc->iindx);
    free(fc->jindx);
    free(fc->params);
    free(fc->exp_params);

    vrna_hc_free(fc->hc);
    vrna_ud_remove(fc);
    vrna_gr_reset(fc);
    vrna_sequence_remove_all(fc);

    /* now distinguish the fc type */
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        free(fc->sequence);
        free(fc->sequence_encoding);
        free(fc->sequence_encoding2);
        free(fc->ptype);
        free(fc->ptype_pf_compat);
        vrna_sc_free(fc->sc);
        break;
      case VRNA_FC_TYPE_COMPARATIVE:
        for (s = 0; s < fc->n_seq; s++) {
          free(fc->sequences[s]);
          free(fc->S[s]);
          free(fc->S5[s]);
          free(fc->S3[s]);
          free(fc->Ss[s]);
          free(fc->a2s[s]);
        }
        free(fc->sequences);
        free(fc->cons_seq);
        free(fc->S_cons);
        free(fc->S);
        free(fc->S5);
        free(fc->S3);
        free(fc->Ss);
        free(fc->a2s);
        free(fc->pscore);
        free(fc->pscore_pf_compat);
        if (fc->scs) {
          for (s = 0; s < fc->n_seq; s++)
            vrna_sc_free(fc->scs[s]);
          free(fc->scs);
        }

        break;
      default:                      /* do nothing */
        break;
    }

    /* free Distance Class Partitioning stuff (should be NULL if not used) */
    free(fc->reference_pt1);
    free(fc->reference_pt2);
    free(fc->referenceBPs1);
    free(fc->referenceBPs2);
    free(fc->bpdist);
    free(fc->mm1);
    free(fc->mm2);

    /* free local folding related stuff (should be NULL if not used) */
    free(fc->ptype_local);
    free(fc->pscore_local);
#ifdef VRNA_WITH_SVM
    vrna_zsc_filter_free(fc);
#endif

    if (fc->free_auxdata)
      fc->free_auxdata(fc->auxdata);

    free(fc);
  }
}


PUBLIC vrna_fold_compound_t *
vrna_fold_compound(const char       *sequence,
                   const vrna_md_t  *md_p,
                   unsigned int     options)
{
  unsigned int          length, aux_options;
  vrna_fold_compound_t  *fc;
  vrna_md_t             md;

  if (sequence == NULL)
    return NULL;

  /* sanity check */
  length = strlen(sequence);
  if (length == 0) {
    vrna_log_warning("vrna_fold_compound@data_structures.c: "
                         "sequence length must be greater 0");
    return NULL;
  }

  if (length > vrna_sequence_length_max(options)) {
    vrna_log_warning("vrna_fold_compound@data_structures.c: "
                         "sequence length of %d exceeds addressable range",
                         length);
    return NULL;
  }

  fc = init_fc_single();

  fc->length    = length;
  fc->sequence  = strdup(sequence);

  aux_options = 0L;


  /* get a copy of the model details */
  if (md_p)
    md = *md_p;
  else
    vrna_md_set_default(&md);

  /* now for the energy parameters */
  add_params(fc, &md, options);

  sanitize_bp_span(fc, options);

  if (options & VRNA_OPTION_WINDOW) {
    set_fold_compound(fc, options, aux_options);

    if (!(options & VRNA_OPTION_EVAL_ONLY)) {
      /* add minimal hard constraint data structure */
      vrna_hc_init_window(fc);

      /* add DP matrices */
      vrna_mx_add(fc, VRNA_MX_WINDOW, options);
    }
  } else {
    /* regular global structure prediction */
    aux_options |= WITH_PTYPE;

    if (options & VRNA_OPTION_PF)
      aux_options |= WITH_PTYPE_COMPAT;

    set_fold_compound(fc, options, aux_options);

    if (!(options & VRNA_OPTION_EVAL_ONLY)) {
      /* add default hard constraints */
      vrna_hc_init(fc);

      /* add DP matrices (if required) */
      vrna_mx_add(fc, VRNA_MX_DEFAULT, options);
    }
  }

  return fc;
}


PUBLIC vrna_fold_compound_t *
vrna_fold_compound_comparative(const char   **sequences,
                               vrna_md_t    *md_p,
                               unsigned int options)
{
  return vrna_fold_compound_comparative2(sequences,
                                         NULL,
                                         NULL,
                                         NULL,
                                         NULL,
                                         md_p,
                                         options);
}


PUBLIC vrna_fold_compound_t *
vrna_fold_compound_comparative2(const char                **sequences,
                                const char                **names,
                                const unsigned char       *orientation,
                                const unsigned long long  *start,
                                const unsigned long long  *genome_size,
                                vrna_md_t                 *md_p,
                                unsigned int              options)
{
  unsigned int          s, n_seq, length;
  vrna_fold_compound_t  *fc;
  vrna_md_t             md;
  unsigned int          aux_options;

  aux_options = 0U;

  if (sequences == NULL)
    return NULL;

  for (s = 0; sequences[s]; s++);  /* count the sequences */

  n_seq = s;

  length = strlen(sequences[0]);
  /* sanity check */
  if (length == 0) {
    vrna_log_warning("vrna_fold_compound_comparative: "
                         "sequence length must be greater 0");
  } else if (length > vrna_sequence_length_max(options)) {
    vrna_log_warning("vrna_fold_compound_comparative: "
                         "sequence length of %d exceeds addressable range",
                         length);
  }

  for (s = 0; s < n_seq; s++)
    if (strlen(sequences[s]) != length) {
      vrna_log_warning("vrna_fold_compound_comparative: "
                           "uneqal sequence lengths in alignment");
      return NULL;
    }

  fc = init_fc_comparative();

  if (fc) {
    fc->n_seq   = n_seq;
    fc->length  = length;

    /* get a copy of the model details */
    if (md_p)
      md = *md_p;
    else /* this fallback relies on global parameters and thus is not threadsafe */
      vrna_md_set_default(&md);

    /* now for the energy parameters */
    add_params(fc, &md, options);

    sanitize_bp_span(fc, options);

    vrna_msa_add(fc,
                 sequences,
                 names,
                 orientation,
                 start,
                 genome_size,
                 VRNA_SEQUENCE_RNA);

    fc->sequences = vrna_alloc(sizeof(char *) * (fc->n_seq + 1));
    for (s = 0; sequences[s]; s++)
      fc->sequences[s] = strdup(sequences[s]);

    if (options & VRNA_OPTION_WINDOW) {
      set_fold_compound(fc, options, aux_options);

      fc->pscore_local = vrna_alloc(sizeof(int *) * (fc->length + 1));

#if 0
      for (i = (int)fc->length; (i > (int)fc->length - fc->window_size - 5) && (i >= 0); i--)
        fc->pscore_local[i] = vrna_alloc(sizeof(int) * (fc->window_size + 5));
#endif

      if (!(options & VRNA_OPTION_EVAL_ONLY)) {
        /* add minimal hard constraint data structure */
        vrna_hc_init_window(fc);

        /* add DP matrices */
        vrna_mx_add(fc, VRNA_MX_WINDOW, options);
      }
    } else {
      /* regular global structure prediction */

      aux_options |= WITH_PTYPE;

      if (options & VRNA_OPTION_PF)
        aux_options |= WITH_PTYPE_COMPAT;

      set_fold_compound(fc, options, aux_options);

      make_pscores(fc);

      if (!(options & VRNA_OPTION_EVAL_ONLY)) {
        /* add default hard constraints */
        vrna_hc_init(fc);

        /* add DP matrices (if required) */
        vrna_mx_add(fc, VRNA_MX_DEFAULT, options);
      }
    }
  }

  return fc;
}


PUBLIC vrna_fold_compound_t *
vrna_fold_compound_TwoD(const char    *sequence,
                        const char    *s1,
                        const char    *s2,
                        vrna_md_t     *md_p,
                        unsigned int  options)
{
  unsigned int          length, l, turn;
  vrna_fold_compound_t  *fc;
  vrna_md_t             md;


  if (sequence == NULL)
    return NULL;

  /* sanity check */
  length = strlen(sequence);
  if (length == 0) {
    vrna_log_warning("vrna_fold_compound_TwoD: "
                         "sequence length must be greater 0");
    return NULL;
  } else if (length > vrna_sequence_length_max(options)) {
    vrna_log_warning("vrna_fold_compound_TwoD: "
                         "sequence length of %d exceeds addressable range",
                         length);
    return NULL;
  }

  l = strlen(s1);
  if (l != length) {
    vrna_log_warning("vrna_fold_compound_TwoD: "
                         "sequence and s1 differ in length");
    return NULL;
  }

  l = strlen(s2);
  if (l != length) {
    vrna_log_warning("vrna_fold_compound_TwoD: "
                         "sequence and s2 differ in length");
    return NULL;
  }

  fc = init_fc_single();
  if (fc) {
    fc->length    = length;
    fc->sequence  = strdup(sequence);

    /* get a copy of the model details */
    if (md_p)
      md = *md_p;
    else /* this fallback relies on global parameters and thus is not threadsafe */
      vrna_md_set_default(&md);

    /* always make uniq ML decomposition ! */
    md.uniq_ML      = 1;
    md.compute_bpp  = 0;

    /* now for the energy parameters */
    add_params(fc, &md, options);

    set_fold_compound(fc, options, WITH_PTYPE | WITH_PTYPE_COMPAT);

    if (!(options & VRNA_OPTION_EVAL_ONLY)) {
      vrna_hc_init(fc); /* add default hard constraints */

      /* add DP matrices */
      vrna_mx_add(fc, VRNA_MX_2DFOLD, options);
    }

    /* set all fields that are unique to Distance class partitioning... */
    turn              = fc->params->model_details.min_loop_size;
    fc->reference_pt1 = vrna_ptable(s1);
    fc->reference_pt2 = vrna_ptable(s2);
    fc->referenceBPs1 = vrna_refBPcnt_matrix(fc->reference_pt1, turn);
    fc->referenceBPs2 = vrna_refBPcnt_matrix(fc->reference_pt2, turn);
    fc->bpdist        = vrna_refBPdist_matrix(fc->reference_pt1, fc->reference_pt2, turn);
    /* compute maximum matching with reference structure 1 disallowed */
    fc->mm1 = maximumMatchingConstraint(fc->sequence, fc->reference_pt1);
    /* compute maximum matching with reference structure 2 disallowed */
    fc->mm2 = maximumMatchingConstraint(fc->sequence, fc->reference_pt2);

    fc->maxD1 = fc->mm1[fc->iindx[1] - length] + fc->referenceBPs1[fc->iindx[1] - length];
    fc->maxD2 = fc->mm2[fc->iindx[1] - length] + fc->referenceBPs2[fc->iindx[1] - length];
  }

  return fc;
}


PUBLIC void
vrna_fold_compound_add_auxdata(vrna_fold_compound_t       *fc,
                               void                       *data,
                               vrna_auxdata_free_f f)
{
  if (fc && data) {
    if (fc->free_auxdata) /* free pre-existing auxdata */
      fc->free_auxdata(fc->auxdata);

    fc->auxdata       = data;
    fc->free_auxdata  = f;
  }
}


PUBLIC void
vrna_fold_compound_add_callback(vrna_fold_compound_t            *fc,
                                vrna_recursion_status_f  f)
{
  if (fc && f)
    fc->stat_cb = f;
}


PUBLIC int
vrna_fold_compound_prepare(vrna_fold_compound_t *fc,
                           unsigned int         options)
{
  int ret = 1; /* success */

  /* check maximum sequence length restrictions */
  if (fc->length > vrna_sequence_length_max(options)) {
    vrna_log_warning(
      "vrna_fold_compound_prepare@data_structures.c: sequence length of %d exceeds addressable range",
      fc->length);
    return 0;
  }

  /* make sure to always provide sane bp-span settings */
  sanitize_bp_span(fc, options);

  /* prepare Boltzmann factors if required */
  vrna_params_prepare(fc, options);

  /* prepare ptype array(s) */
  vrna_ptypes_prepare(fc, options);

  if (options & VRNA_OPTION_PF) {
    /* prepare for partition function computation */

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:     /* get pre-computed Boltzmann factors if not present*/
        if (fc->domains_up)         /* turn on unique ML decomposition with qm1 array */
          fc->exp_params->model_details.uniq_ML = 1;

        break;

      default:
        /* not doing anything here... */
        break;
    }
  }

  /* prepare hard constraints */
  vrna_hc_prepare(fc, options);

  /* prepare soft constraints data structure, if required */
  vrna_sc_prepare(fc, options);

  /* Add DP matrices, if not they are not present or do not fit current settings */
  vrna_mx_prepare(fc, options);

  /* prepare auxiliary grammar rules data structure, if required */
  ret &= vrna_gr_prepare(fc, options);

  return ret;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE void
sanitize_bp_span(vrna_fold_compound_t *fc,
                 unsigned int         options)
{
  vrna_md_t *md;

  md = &(fc->params->model_details);

  /* make sure that min_loop_size, max_bp_span, and window_size are sane */
  if (options & VRNA_OPTION_WINDOW) {
    if (md->window_size <= 0)
      md->window_size = (int)fc->length;
    else if (md->window_size > (int)fc->length)
      md->window_size = (int)fc->length;

    fc->window_size = md->window_size;
  } else {
    /* non-local fold mode */
    md->window_size = (int)fc->length;
  }

  if ((md->max_bp_span <= 0) || (md->max_bp_span > md->window_size))
    md->max_bp_span = md->window_size;
}


PRIVATE void
add_params(vrna_fold_compound_t *fc,
           vrna_md_t            *md_p,
           unsigned int         options)
{
  /*
   * ALWAYS provide regular energy parameters
   * remove previous parameters if present and they differ from current model
   */
  if (fc->params) {
    if (memcmp(md_p, &(fc->params->model_details), sizeof(vrna_md_t)) != 0) {
      free(fc->params);
      fc->params = NULL;
    }
  }

  if (!fc->params)
    fc->params = vrna_params(md_p);

  vrna_params_prepare(fc, options);
}


PRIVATE void
set_fold_compound(vrna_fold_compound_t  *fc,
                  unsigned int          options,
                  unsigned int          aux)
{
  char          *sequence, **sequences, **ptr;
  unsigned int  length, s;
  vrna_md_t     *md_p;

  sequence  = NULL;
  sequences = NULL;

  md_p = &(fc->params->model_details);

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      sequence = fc->sequence;

      fc->sequence  = NULL;
      fc->length    = 0;

      /* split input sequences at default delimiter '&' */
      sequences = vrna_strsplit(sequence, NULL);

      /* add individual sequences to fold compound */
      for (ptr = sequences; *ptr; ptr++) {
        vrna_sequence_add(fc, *ptr, VRNA_SEQUENCE_RNA);
        free(*ptr);
      }

      free(sequences);
      free(sequence);

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
      if (fc->strands > 1) {
        fc->cutpoint = fc->nucleotides[0].length + 1;

#if 0
        if (md_p->min_loop_size == TURN)
          md_p->min_loop_size = 0;                              /* is it safe to set this here? */
#endif
      }

#endif

      if (!(options & VRNA_OPTION_EVAL_ONLY)) {
        /* temporary hack for multi-strand case */
        if (fc->strands > 1) {
          int min_loop_size = md_p->min_loop_size;
          md_p->min_loop_size = 0;
          fc->ptype = (aux & WITH_PTYPE) ? vrna_ptypes(fc->sequence_encoding2, md_p) : NULL;
          md_p->min_loop_size = min_loop_size;
        } else {
          fc->ptype = (aux & WITH_PTYPE) ? vrna_ptypes(fc->sequence_encoding2, md_p) : NULL;
        }
#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
        /* backward compatibility ptypes */
        fc->ptype_pf_compat =
          (aux & WITH_PTYPE_COMPAT) ? get_ptypes(fc->sequence_encoding2, md_p, 1) : NULL;
#endif
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      sequences = fc->sequences;

      fc->length = length = fc->length;

      fc->cons_seq  = vrna_aln_consensus_sequence((const char **)sequences, md_p);
      fc->S_cons    = vrna_seq_encode_simple(fc->cons_seq, md_p);

      fc->pscore = vrna_alloc(sizeof(int) * ((length * (length + 1)) / 2 + 2));
      /* backward compatibility ptypes */
      fc->pscore_pf_compat =
        (aux & WITH_PTYPE_COMPAT) ? vrna_alloc(sizeof(int) *
                                               ((length * (length + 1)) / 2 + 2)) : NULL;

      oldAliEn = fc->oldAliEn = md_p->oldAliEn;

      fc->S   = vrna_alloc((fc->n_seq + 1) * sizeof(short *));
      fc->S5  = vrna_alloc((fc->n_seq + 1) * sizeof(short *));
      fc->S3  = vrna_alloc((fc->n_seq + 1) * sizeof(short *));
      fc->a2s = vrna_alloc((fc->n_seq + 1) * sizeof(unsigned int *));
      fc->Ss  = vrna_alloc((fc->n_seq + 1) * sizeof(char *));

      for (s = 0; s < fc->n_seq; s++) {
        vrna_aln_encode(fc->sequences[s],
                        &(fc->S[s]),
                        &(fc->S5[s]),
                        &(fc->S3[s]),
                        &(fc->Ss[s]),
                        &(fc->a2s[s]),
                        md_p);
      }
      fc->S5[fc->n_seq]   = NULL;
      fc->S3[fc->n_seq]   = NULL;
      fc->a2s[fc->n_seq]  = NULL;
      fc->Ss[fc->n_seq]   = NULL;
      fc->S[fc->n_seq]    = NULL;

      break;

    default:                      /* do nothing ? */
      break;
  }

  vrna_sequence_prepare(fc);

  if (!(options & VRNA_OPTION_WINDOW) && (fc->length <= vrna_sequence_length_max(options))) {
    fc->iindx = vrna_idx_row_wise(fc->length);
    fc->jindx = vrna_idx_col_wise(fc->length);
  }
}


PRIVATE void
make_pscores(vrna_fold_compound_t *fc)
{
  /*
   * calculate co-variance bonus for each pair depending on
   * compensatory/consistent mutations and incompatible seqs
   * should be 0 for conserved pairs, >0 for good pairs
   */

#define NONE -10000 /* score for forbidden pairs */

  int       i, j, k, l, s, max_span, turn;
  short     **S   = fc->S;
  char      **AS  = fc->sequences;
  int       n_seq = fc->n_seq;
  vrna_md_t *md   =
    (fc->params) ? &(fc->params->model_details) : &(fc->exp_params->model_details);
  int       *pscore   = fc->pscore;             /* precomputed array of pair types */
  int       *indx     = fc->jindx;
  int       *my_iindx = fc->iindx;
  int       n         = fc->length;

  turn = md->min_loop_size;

  if (md->ribo) {
    float **dm = NULL;

    if (RibosumFile != NULL)
      dm = readribosum(RibosumFile);
    else
      dm = get_ribosum((const char **)AS, n_seq, n);

    if (dm) {
      for (i = 0; i < 7; i++) {
        for (j = 0; j < 7; j++)
          md->pair_dist[i][j] = dm[i][j];

        free(dm[i]);
      }
      free(dm);
    }
  }

  max_span = md->max_bp_span;
  if ((max_span < turn + 2) || (max_span > n))
    max_span = n;

  for (i = 1; i < n; i++) {
    for (j = i + 1; (j < i + turn + 1) && (j <= n); j++)
      pscore[indx[j] + i] = NONE;
    for (j = i + turn + 1; j <= n; j++) {
      if ((j - i + 1) > max_span) {
        pscore[indx[j] + i] = NONE;
        continue;
      }

      unsigned int pfreq[8] = {
        0, 0, 0, 0, 0, 0, 0, 0
      };

      for (s = 0; s < n_seq; s++) {
        int type;
        if (S[s][i] == 0 && S[s][j] == 0) {
          type = 7;                             /* gap-gap  */
        } else {
          if ((AS[s][i] == '~') || (AS[s][j] == '~')) {
            type = 7;
          } else {
            type = md->pair[S[s][i]][S[s][j]];
            if ((md->noGU) && ((type == 3) || (type == 4)))
              type = 0;
          }
        }

        pfreq[type]++;
      }

      pscore[indx[j] + i] = vrna_pscore_freq(fc, &pfreq[0], 6);
    }
  }

  if (md->noLP) {
    /* remove unwanted pairs */
    for (k = 1; k < n - turn - 1; k++)
      for (l = 1; l <= 2; l++) {
        int type, ntype = 0, otype = 0;
        i     = k;
        j     = i + turn + l;
        type  = pscore[indx[j] + i];
        while ((i >= 1) && (j <= n)) {
          if ((i > 1) && (j < n))
            ntype = pscore[indx[j + 1] + i - 1];

          if ((otype < md->cv_fact * MINPSCORE) &&
              (ntype < md->cv_fact * MINPSCORE))  /* too many counterexamples */
            pscore[indx[j] + i] = NONE;           /* i.j can only form isolated pairs */

          otype = type;
          type  = ntype;
          i--;
          j++;
        }
      }
  }

  /* copy over pscores for backward compatibility */
  if (fc->pscore_pf_compat) {
    for (i = 1; i < n; i++)
      for (j = i; j <= n; j++)
        fc->pscore_pf_compat[my_iindx[i] - j] = (short)pscore[indx[j] + i];
  }
}


PRIVATE vrna_fold_compound_t *
init_fc_single(void)
{
  vrna_fold_compound_t  init = {
    .type = VRNA_FC_TYPE_SINGLE
  };
  vrna_fold_compound_t  *fc = vrna_alloc(sizeof(vrna_fold_compound_t));

  if (fc) {
    memcpy(fc, &init, sizeof(vrna_fold_compound_t));
    nullify(fc);
  }

  return fc;
}


PRIVATE vrna_fold_compound_t *
init_fc_comparative(void)
{
  vrna_fold_compound_t  init = {
    .type = VRNA_FC_TYPE_COMPARATIVE
  };
  vrna_fold_compound_t  *fc = vrna_alloc(sizeof(vrna_fold_compound_t));

  if (fc) {
    memcpy(fc, &init, sizeof(vrna_fold_compound_t));
    nullify(fc);
  }

  return fc;
}


PRIVATE INLINE void
nullify(vrna_fold_compound_t *fc)
{
  if (fc) {
    fc->length            = 0;
    fc->strands           = 0;
    fc->cutpoint          = -1;
    fc->strand_number     = NULL;
    fc->strand_order      = NULL;
    fc->strand_order_uniq = NULL;
    fc->strand_start      = NULL;
    fc->strand_end        = NULL;
    fc->nucleotides       = NULL;
    fc->alignment         = NULL;

    fc->hc            = NULL;
    fc->matrices      = NULL;
    fc->exp_matrices  = NULL;
    fc->params        = NULL;
    fc->exp_params    = NULL;
    fc->iindx         = NULL;
    fc->jindx         = NULL;

    fc->stat_cb       = NULL;
    fc->auxdata       = NULL;
    fc->free_auxdata  = NULL;

    fc->domains_struc = NULL;
    fc->domains_up    = NULL;
    fc->aux_grammar   = NULL;

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        fc->sequence            = NULL;
        fc->sequence_encoding   = NULL;
        fc->encoding5           = NULL;
        fc->encoding3           = NULL;
        fc->sequence_encoding2  = NULL;
        fc->ptype               = NULL;
        fc->ptype_pf_compat     = NULL;
        fc->sc                  = NULL;

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        fc->sequences         = NULL;
        fc->n_seq             = 0;
        fc->cons_seq          = NULL;
        fc->S_cons            = NULL;
        fc->S                 = NULL;
        fc->S5                = NULL;
        fc->S3                = NULL;
        fc->Ss                = NULL;
        fc->a2s               = NULL;
        fc->pscore            = NULL;
        fc->pscore_local      = NULL;
        fc->pscore_pf_compat  = NULL;
        fc->scs               = NULL;
        fc->oldAliEn          = 0;

        break;
    }

    fc->maxD1         = 0;
    fc->maxD2         = 0;
    fc->reference_pt1 = NULL;
    fc->reference_pt2 = NULL;
    fc->referenceBPs1 = NULL;
    fc->referenceBPs2 = NULL;
    fc->bpdist        = NULL;
    fc->mm1           = NULL;
    fc->mm2           = NULL;

    fc->window_size = -1;
    fc->ptype_local = NULL;
#ifdef VRNA_WITH_SVM
    fc->zscore_data = NULL;
#endif
  }
}
