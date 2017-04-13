/** \file data_structures.c **/

/*
 *                Data structure creation/destruction
 *
 *                This file contains everything which is necessary to
 *                obtain and destroy datastructures used in the folding
 *                recurrences throughout the VienneRNA paclage
 *
 *                c Ronny Lorenx
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
#include "ViennaRNA/structure_utils.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/data_structures.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/aln_util.h"
#include "ViennaRNA/ribo.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/cofold.h"
#include "ViennaRNA/mm.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/*
 #################################
 # PRIVATE MACROS                #
 #################################
 */

#define WITH_PTYPE          1L    /* passed to set_fold_compound() to indicate that we need to set vc->ptype */
#define WITH_PTYPE_COMPAT   2L    /* passed to set_fold_compound() to indicate that we need to set vc->ptype_compat */

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
PRIVATE void  set_fold_compound(vrna_fold_compound_t  *vc,
                                vrna_md_t             *md_p,
                                unsigned int          options,
                                unsigned int          aux);


PRIVATE void  make_pscores(vrna_fold_compound_t *vc);


PRIVATE void  add_params(vrna_fold_compound_t *vc,
                         vrna_md_t            *md_p,
                         unsigned int         options);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC void
vrna_fold_compound_free(vrna_fold_compound_t *vc)
{
  int s;

  if (vc) {
    /* first destroy common attributes */
    vrna_mx_mfe_free(vc);
    vrna_mx_pf_free(vc);
    free(vc->iindx);
    free(vc->jindx);
    free(vc->params);
    free(vc->exp_params);
    free(vc->strand_number);
    vrna_hc_free(vc->hc);
    vrna_ud_remove(vc);

    /* now distinguish the vc type */
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        free(vc->sequence);
        free(vc->sequence_encoding);
        free(vc->sequence_encoding2);
        free(vc->ptype);
        free(vc->ptype_pf_compat);
        vrna_sc_free(vc->sc);
        break;
      case VRNA_FC_TYPE_COMPARATIVE:
        for (s = 0; s < vc->n_seq; s++) {
          free(vc->sequences[s]);
          free(vc->S[s]);
          free(vc->S5[s]);
          free(vc->S3[s]);
          free(vc->Ss[s]);
          free(vc->a2s[s]);
        }
        free(vc->sequences);
        free(vc->cons_seq);
        free(vc->S_cons);
        free(vc->S);
        free(vc->S5);
        free(vc->S3);
        free(vc->Ss);
        free(vc->a2s);
        free(vc->pscore);
        free(vc->pscore_pf_compat);
        if (vc->scs) {
          for (s = 0; s < vc->n_seq; s++)
            vrna_sc_free(vc->scs[s]);
          free(vc->scs);
        }

        break;
      default:                      /* do nothing */
        break;
    }

    /* free Distance Class Partitioning stuff (should be NULL if not used) */
    free(vc->reference_pt1);
    free(vc->reference_pt2);
    free(vc->referenceBPs1);
    free(vc->referenceBPs2);
    free(vc->bpdist);
    free(vc->mm1);
    free(vc->mm2);

    /* free local folding related stuff (should be NULL if not used) */
    free(vc->ptype_local);
    free(vc->pscore_local);

    if (vc->free_auxdata)
      vc->free_auxdata(vc->auxdata);

    free(vc);
  }
}


PUBLIC vrna_fold_compound_t *
vrna_fold_compound(const char   *sequence,
                   vrna_md_t    *md_p,
                   unsigned int options)
{
  int                   i;
  unsigned int          length, aux_options;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  if (sequence == NULL)
    return NULL;

  /* sanity check */
  length = strlen(sequence);
  if (length == 0)
    vrna_message_error("vrna_fold_compound@data_structures.c: sequence length must be greater 0");

  if (length > vrna_sequence_length_max(options))
    vrna_message_error("vrna_fold_compound@data_structures.c: sequence length of %d exceeds addressable range", length);

  vc            = vrna_alloc(sizeof(vrna_fold_compound_t));
  vc->type      = VRNA_FC_TYPE_SINGLE;
  vc->length    = length;
  vc->sequence  = strdup(sequence);
  aux_options   = 0L;


  /* get a copy of the model details */
  if (md_p)
    md = *md_p;
  else
    vrna_md_set_default(&md);

  if (options & VRNA_OPTION_WINDOW) {
    /* sliding window structure prediction */
    if (md.window_size <= 0)
      md.window_size = (int)vc->length;
    else if (md.window_size > (int)vc->length)
      md.window_size = (int)vc->length;

    vc->window_size = md.window_size;

    if ((md.max_bp_span <= 0) || (md.max_bp_span > md.window_size))
      md.max_bp_span = md.window_size;

    set_fold_compound(vc, &md, options, aux_options);

    vc->ptype_local = vrna_alloc(sizeof(char *) * (vc->length + 1));

    if (!(options & VRNA_OPTION_EVAL_ONLY)) {
      /* add default hard constraints */
      /* vrna_hc_init(vc); */ /* no hard constraints in Lfold, yet! */

      /* add DP matrices */
      vrna_mx_add(vc, VRNA_MX_WINDOW, options);
    }
  } else {
    /* regular global structure prediction */

    /* set window size to entire sequence */
    md.window_size = (int)vc->length;

    aux_options |= WITH_PTYPE;

    if (options & VRNA_OPTION_PF)
      aux_options |= WITH_PTYPE_COMPAT;

    set_fold_compound(vc, &md, options, aux_options);

    if (!(options & VRNA_OPTION_EVAL_ONLY)) {
      /* add default hard constraints */
      vrna_hc_init(vc);

      /* add DP matrices (if required) */
      vrna_mx_add(vc, VRNA_MX_DEFAULT, options);
    }
  }

  return vc;
}


PUBLIC vrna_fold_compound_t *
vrna_fold_compound_comparative(const char   **sequences,
                               vrna_md_t    *md_p,
                               unsigned int options)
{
  int                   i, s, n_seq, length;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;
  unsigned int          aux_options;

  aux_options = 0U;

  if (sequences == NULL)
    return NULL;

  for (s = 0; sequences[s]; s++) ; /* count the sequences */

  n_seq = s;

  length = strlen(sequences[0]);
  /* sanity check */
  if (length == 0)
    vrna_message_error("vrna_fold_compound_comparative@data_structures.c: sequence length must be greater 0");
  else if (length > vrna_sequence_length_max(options))
    vrna_message_error("vrna_fold_compound_comparative@data_structures.c: sequence length of %d exceeds addressable range", length);

  for (s = 0; s < n_seq; s++)
    if (strlen(sequences[s]) != length)
      vrna_message_error("vrna_fold_compound_comparative@data_structures.c: uneqal sequence lengths in alignment");

  vc        = vrna_alloc(sizeof(vrna_fold_compound_t));
  vc->type  = VRNA_FC_TYPE_COMPARATIVE;

  vc->n_seq     = n_seq;
  vc->length    = length;
  vc->sequences = vrna_alloc(sizeof(char *) * (vc->n_seq + 1));
  for (s = 0; sequences[s]; s++)
    vc->sequences[s] = strdup(sequences[s]);

  /* get a copy of the model details */
  if (md_p)
    md = *md_p;
  else /* this fallback relies on global parameters and thus is not threadsafe */
    vrna_md_set_default(&md);

  if (options & VRNA_OPTION_WINDOW) {
    /* sliding window structure prediction */
    if (md.window_size <= 0)
      md.window_size = (int)vc->length;
    else if (md.window_size > (int)vc->length)
      md.window_size = (int)vc->length;

    vc->window_size = md.window_size;

    if ((md.max_bp_span <= 0) || (md.max_bp_span > md.window_size))
      md.max_bp_span = md.window_size;

    set_fold_compound(vc, &md, options, aux_options);

    vc->pscore_local = vrna_alloc(sizeof(int *) * (vc->length + 1));

#if 0
    for (i = (int)vc->length; (i > (int)vc->length - vc->window_size - 5) && (i >= 0); i--)
      vc->pscore_local[i] = vrna_alloc(sizeof(int) * (vc->window_size + 5));
#endif

    if (!(options & VRNA_OPTION_EVAL_ONLY)) {
      /* add default hard constraints */
      /* vrna_hc_init(vc); */ /* no hard constraints in aliLfold, yet! */

      /* add DP matrices */
      vrna_mx_add(vc, VRNA_MX_WINDOW, options);
    }
  } else {
    /* regular global structure prediction */

    aux_options |= WITH_PTYPE;

    if (options & VRNA_OPTION_PF)
      aux_options |= WITH_PTYPE_COMPAT;

    set_fold_compound(vc, &md, options, aux_options);

    make_pscores(vc);

    if (!(options & VRNA_OPTION_EVAL_ONLY)) {
      /* add default hard constraints */
      vrna_hc_init(vc);

      /* add DP matrices (if required) */
      vrna_mx_add(vc, VRNA_MX_DEFAULT, options);
    }
  }

  return vc;
}


PUBLIC vrna_fold_compound_t *
vrna_fold_compound_TwoD(const char    *sequence,
                        const char    *s1,
                        const char    *s2,
                        vrna_md_t     *md_p,
                        unsigned int  options)
{
  int                   length, l, turn;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;


  if (sequence == NULL)
    return NULL;

  /* sanity check */
  length = strlen(sequence);
  if (length == 0)
    vrna_message_error("vrna_fold_compound_TwoD: sequence length must be greater 0");
  else if (length > vrna_sequence_length_max(options))
    vrna_message_error("vrna_fold_compound_TwoD@data_structures.c: sequence length of %d exceeds addressable range", length);

  l = strlen(s1);
  if (l != length)
    vrna_message_error("vrna_fold_compound_TwoD: sequence and s1 differ in length");

  l = strlen(s2);
  if (l != length)
    vrna_message_error("vrna_fold_compound_TwoD: sequence and s2 differ in length");

  vc            = vrna_alloc(sizeof(vrna_fold_compound_t));
  vc->type      = VRNA_FC_TYPE_SINGLE;
  vc->length    = length;
  vc->sequence  = strdup(sequence);

  /* get a copy of the model details */
  if (md_p)
    md = *md_p;
  else /* this fallback relies on global parameters and thus is not threadsafe */
    vrna_md_set_default(&md);

  /* always make uniq ML decomposition ! */
  md.uniq_ML      = 1;
  md.compute_bpp  = 0;

  set_fold_compound(vc, &md, options, WITH_PTYPE | WITH_PTYPE_COMPAT);

  if (!(options & VRNA_OPTION_EVAL_ONLY)) {
    vrna_hc_init(vc); /* add default hard constraints */

    /* add DP matrices */
    vrna_mx_add(vc, VRNA_MX_2DFOLD, options);
  }

  /* set all fields that are unique to Distance class partitioning... */
  turn              = vc->params->model_details.min_loop_size;
  vc->reference_pt1 = vrna_ptable(s1);
  vc->reference_pt2 = vrna_ptable(s2);
  vc->referenceBPs1 = vrna_refBPcnt_matrix(vc->reference_pt1, turn);
  vc->referenceBPs2 = vrna_refBPcnt_matrix(vc->reference_pt2, turn);
  vc->bpdist        = vrna_refBPdist_matrix(vc->reference_pt1, vc->reference_pt2, turn);
  /* compute maximum matching with reference structure 1 disallowed */
  vc->mm1 = maximumMatchingConstraint(vc->sequence, vc->reference_pt1);
  /* compute maximum matching with reference structure 2 disallowed */
  vc->mm2 = maximumMatchingConstraint(vc->sequence, vc->reference_pt2);

  vc->maxD1 = vc->mm1[vc->iindx[1] - length] + vc->referenceBPs1[vc->iindx[1] - length];
  vc->maxD2 = vc->mm2[vc->iindx[1] - length] + vc->referenceBPs2[vc->iindx[1] - length];

  return vc;
}


PUBLIC void
vrna_fold_compound_add_auxdata(vrna_fold_compound_t       *vc,
                               void                       *data,
                               vrna_callback_free_auxdata *f)
{
  if (vc && data) {
    if (vc->free_auxdata) /* free pre-existing auxdata */
      vc->free_auxdata(vc->auxdata);

    vc->auxdata       = data;
    vc->free_auxdata  = f;
  }
}


PUBLIC void
vrna_fold_compound_add_callback(vrna_fold_compound_t            *vc,
                                vrna_callback_recursion_status  *f)
{
  if (vc && f)
    vc->stat_cb = f;
}


PUBLIC int
vrna_fold_compound_prepare(vrna_fold_compound_t *vc,
                           unsigned int         options)
{
  int ret = 1; /* success */

  /* check maximum sequence length restrictions */
  if (vc->length > vrna_sequence_length_max(options)) {
    vrna_message_warning("vrna_fold_compound_prepare@data_structures.c: sequence length of %d exceeds addressable range", vc->length);
    return 0;
  }

  if (options & VRNA_OPTION_MFE) {
    /* prepare for MFE computation */
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (!vc->ptype)
          if (!(options & VRNA_OPTION_WINDOW))
            vc->ptype = vrna_ptypes(vc->sequence_encoding2,
                                    &(vc->params->model_details));

        break;
      case VRNA_FC_TYPE_COMPARATIVE:
        break;
      default:
        break;
    }
  }

  if (options & VRNA_OPTION_PF) {
    /* prepare for partition function computation */

    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:     /* get pre-computed Boltzmann factors if not present*/
        if (!vc->exp_params)
          vc->exp_params = vrna_exp_params(&(vc->params->model_details));

        if (!vc->ptype)
          vc->ptype = vrna_ptypes(vc->sequence_encoding2, &(vc->exp_params->model_details));

#ifdef VRNA_BACKWARD_COMPAT
        /* backward compatibility ptypes */
        if (!vc->ptype_pf_compat)
          vc->ptype_pf_compat = get_ptypes(vc->sequence_encoding2, &(vc->exp_params->model_details), 1);

#endif
        /* get precomputed Boltzmann factors for soft-constraints (if any) */
        if (vc->sc) {
          if (!vc->sc->exp_energy_up)
            vrna_sc_set_up(vc, NULL, VRNA_OPTION_PF);

          if (!vc->sc->exp_energy_bp)
            vrna_sc_set_bp(vc, NULL, VRNA_OPTION_PF);

          if (!vc->sc->exp_energy_stack)
            vrna_sc_add_SHAPE_deigan(vc, NULL, 0, 0, VRNA_OPTION_PF);
        }

        if (vc->domains_up)                            /* turn on unique ML decomposition with qm1 array */
          vc->exp_params->model_details.uniq_ML = 1;

        break;

      case VRNA_FC_TYPE_COMPARATIVE:  /* get pre-computed Boltzmann factors if not present*/
        if (!vc->exp_params)
          vc->exp_params = vrna_exp_params_comparative(vc->n_seq, &(vc->params->model_details));

        break;

      default:
        break;
    }
  }

  /* Add DP matrices, if not they are not present or do not fit current settings */
  vrna_mx_prepare(vc, options);

  return ret;
}


#ifndef VRNA_DISABLE_C11_FEATURES
PUBLIC void
vrna_C11_features(void)
{
  __asm("nop");
}


#endif

/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE void
add_params(vrna_fold_compound_t *vc,
           vrna_md_t            *md_p,
           unsigned int         options)
{
  /* ALWAYS add regular energy parameters */
  vc->params = vrna_params(md_p);

  if (options & VRNA_OPTION_PF) {
    vc->exp_params = (vc->type == VRNA_FC_TYPE_SINGLE) ? \
                     vrna_exp_params(md_p) : \
                     vrna_exp_params_comparative(vc->n_seq, md_p);
  }
}


PRIVATE void
set_fold_compound(vrna_fold_compound_t  *vc,
                  vrna_md_t             *md_p,
                  unsigned int          options,
                  unsigned int          aux)
{
  char          *sequence, **sequences;
  unsigned int  length, s, i;
  int           cp;                           /* cut point for cofold */
  char          *seq, *seq2;

  sequence  = NULL;
  sequences = NULL;
  cp        = -1;

  /* some default init values */
  vc->params        = NULL;
  vc->exp_params    = NULL;
  vc->matrices      = NULL;
  vc->exp_matrices  = NULL;
  vc->hc            = NULL;
  vc->auxdata       = NULL;
  vc->free_auxdata  = NULL;

  vc->strand_number = NULL;
  vc->domains_struc = NULL;
  vc->domains_up    = NULL;
  vc->aux_grammar   = NULL;

  switch (vc->type) {
    case VRNA_FC_TYPE_SINGLE:
      sequence = vc->sequence;

      seq2          = strdup(sequence);
      seq           = vrna_cut_point_remove(seq2, &cp);                   /*  splice out the '&' if concatenated sequences and
                                                                           * reset cp... this should also be safe for
                                                                           * single sequences */
      vc->cutpoint  = cp;

      if ((cp > 0) && (md_p->min_loop_size == TURN))
        md_p->min_loop_size = 0;                              /* is it safe to set this here? */

      free(vc->sequence);
      vc->sequence            = seq;
      vc->length              = length = strlen(seq);
      vc->sequence_encoding   = vrna_seq_encode(seq, md_p);
      vc->sequence_encoding2  = vrna_seq_encode_simple(seq, md_p);

      vc->strand_number = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (vc->length + 2));
      if (cp > 0) {
        for (s = i = 0; i <= vc->length + 1; i++) {
          /* this sets pos. 0 and n+1 as well */
          if (i == vc->cutpoint)
            s++;

          vc->strand_number[i] = s;
        }
      }

      if (!(options & VRNA_OPTION_EVAL_ONLY)) {
        vc->ptype = (aux & WITH_PTYPE) ? vrna_ptypes(vc->sequence_encoding2, md_p) : NULL;
        /* backward compatibility ptypes */
        vc->ptype_pf_compat = (aux & WITH_PTYPE_COMPAT) ? get_ptypes(vc->sequence_encoding2, md_p, 1) : NULL;
      } else {
        vc->ptype           = NULL;
        vc->ptype_pf_compat = NULL;
      }

      vc->sc = NULL;
      free(seq2);
      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      sequences = vc->sequences;

      vc->length = length = vc->length;

      vc->strand_number = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (vc->length + 1));

      vc->cons_seq  = consensus((const char **)sequences);
      vc->S_cons    = vrna_seq_encode_simple(vc->cons_seq, md_p);

      vc->pscore = vrna_alloc(sizeof(int) * ((length * (length + 1)) / 2 + 2));
      /* backward compatibility ptypes */
      vc->pscore_pf_compat = (aux & WITH_PTYPE_COMPAT) ? vrna_alloc(sizeof(int) * ((length * (length + 1)) / 2 + 2)) : NULL;

      oldAliEn = vc->oldAliEn = md_p->oldAliEn;

      vc->S   = vrna_alloc((vc->n_seq + 1) * sizeof(short *));
      vc->S5  = vrna_alloc((vc->n_seq + 1) * sizeof(short *));
      vc->S3  = vrna_alloc((vc->n_seq + 1) * sizeof(short *));
      vc->a2s = vrna_alloc((vc->n_seq + 1) * sizeof(unsigned short *));
      vc->Ss  = vrna_alloc((vc->n_seq + 1) * sizeof(char *));

      for (s = 0; s < vc->n_seq; s++) {
        vrna_aln_encode(vc->sequences[s],
                        &(vc->S[s]),
                        &(vc->S5[s]),
                        &(vc->S3[s]),
                        &(vc->Ss[s]),
                        &(vc->a2s[s]),
                        md_p);
      }
      vc->S5[vc->n_seq]   = NULL;
      vc->S3[vc->n_seq]   = NULL;
      vc->a2s[vc->n_seq]  = NULL;
      vc->Ss[vc->n_seq]   = NULL;
      vc->S[vc->n_seq]    = NULL;

      vc->scs = NULL;
      break;

    default:                      /* do nothing ? */
      break;
  }

  if (vc->length <= vrna_sequence_length_max(options)) {
    vc->iindx = vrna_idx_row_wise(vc->length);
    vc->jindx = vrna_idx_col_wise(vc->length);
  } else {
    vc->iindx = NULL;
    vc->jindx = NULL;
  }

  /* now come the energy parameters */
  add_params(vc, md_p, options);
}


PRIVATE void
make_pscores(vrna_fold_compound_t *vc)
{
  /* calculate co-variance bonus for each pair depending on  */
  /* compensatory/consistent mutations and incompatible seqs */
  /* should be 0 for conserved pairs, >0 for good pairs      */

#define NONE -10000 /* score for forbidden pairs */

  char      *structure = NULL;
  int       i, j, k, l, s, max_span, turn;
  float     **dm;
  int       olddm[7][7] = { { 0, 0, 0, 0, 0, 0, 0 }, /* hamming distance between pairs */
                            { 0, 0, 2, 2, 1, 2, 2 } /* CG */,
                            { 0, 2, 0, 1, 2, 2, 2 } /* GC */,
                            { 0, 2, 1, 0, 2, 1, 2 } /* GU */,
                            { 0, 1, 2, 2, 0, 2, 1 } /* UG */,
                            { 0, 2, 2, 1, 2, 0, 2 } /* AU */,
                            { 0, 2, 2, 2, 1, 2, 0 } /* UA */ };

  short     **S       = vc->S;
  char      **AS      = vc->sequences;
  int       n_seq     = vc->n_seq;
  vrna_md_t *md       = (vc->params) ? &(vc->params->model_details) : &(vc->exp_params->model_details);
  int       *pscore   = vc->pscore;             /* precomputed array of pair types */
  int       *indx     = vc->jindx;
  int       *my_iindx = vc->iindx;
  int       n         = vc->length;

  turn = md->min_loop_size;

  if (md->ribo) {
    if (RibosumFile != NULL)
      dm = readribosum(RibosumFile);
    else
      dm = get_ribosum((const char **)AS, n_seq, n);
  } else {
    /*use usual matrix*/
    dm = vrna_alloc(7 * sizeof(float *));
    for (i = 0; i < 7; i++) {
      dm[i] = vrna_alloc(7 * sizeof(float));
      for (j = 0; j < 7; j++)
        dm[i][j] = (float)olddm[i][j];
    }
  }

  max_span = md->max_bp_span;
  if ((max_span < turn + 2) || (max_span > n))
    max_span = n;

  for (i = 1; i < n; i++) {
    for (j = i + 1; (j < i + turn + 1) && (j <= n); j++)
      pscore[indx[j] + i] = NONE;
    for (j = i + turn + 1; j <= n; j++) {
      int     pfreq[8] = {
        0, 0, 0, 0, 0, 0, 0, 0
      };
      double  score;
      for (s = 0; s < n_seq; s++) {
        int type;
        if (S[s][i] == 0 && S[s][j] == 0) {
          type = 7;                             /* gap-gap  */
        } else {
          if ((AS[s][i] == '~') || (AS[s][j] == '~'))
            type = 7;
          else
            type = md->pair[S[s][i]][S[s][j]];
        }

        pfreq[type]++;
      }
      if (pfreq[0] * 2 + pfreq[7] > n_seq) {
        pscore[indx[j] + i] = NONE;
        continue;
      }

      for (k = 1, score = 0; k <= 6; k++) /* ignore pairtype 7 (gap-gap) */
        for (l = k; l <= 6; l++)
          score += pfreq[k] * pfreq[l] * dm[k][l];
      /* counter examples score -1, gap-gap scores -0.25   */
      pscore[indx[j] + i] = md->cv_fact *
                            ((UNIT * score) / n_seq - md->nc_fact * UNIT * (pfreq[0] + pfreq[7] * 0.25));

      if ((j - i + 1) > max_span)
        pscore[indx[j] + i] = NONE;
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

          if ((otype < md->cv_fact * MINPSCORE) && (ntype < md->cv_fact * MINPSCORE)) /* too many counterexamples */
            pscore[indx[j] + i] = NONE;                                               /* i.j can only form isolated pairs */

          otype = type;
          type  = ntype;
          i--;
          j++;
        }
      }
  }

  if (fold_constrained && (structure != NULL)) {
    int psij, hx, hx2, *stack, *stack2;
    stack   = vrna_alloc(sizeof(int) * (n + 1));
    stack2  = vrna_alloc(sizeof(int) * (n + 1));

    for (hx = hx2 = 0, j = 1; j <= n; j++) {
      switch (structure[j - 1]) {
        case 'x': /* can't pair */
          for (l = 1; l < j - turn; l++)
            pscore[indx[j] + l] = NONE;
          for (l = j + turn + 1; l <= n; l++)
            pscore[indx[l] + j] = NONE;
          break;
        case '(':
          stack[hx++] = j;
        /* fallthrough */
        case '[':
          stack2[hx2++] = j;
        /* fallthrough */
        case '<': /* pairs upstream */
          for (l = 1; l < j - turn; l++)
            pscore[indx[j] + l] = NONE;
          break;
        case ']':
          if (hx2 <= 0)
            vrna_message_error("unbalanced brackets in constraints\n%s", structure);

          i                   = stack2[--hx2];
          pscore[indx[j] + i] = NONE;
          break;
        case ')':
          if (hx <= 0)
            vrna_message_error("unbalanced brackets in constraints\n%s", structure);

          i     = stack[--hx];
          psij  = pscore[indx[j] + i]; /* store for later */
          for (k = j; k <= n; k++)
            for (l = i; l <= j; l++)
              pscore[indx[k] + l] = NONE;
          for (l = i; l <= j; l++)
            for (k = 1; k <= i; k++)
              pscore[indx[l] + k] = NONE;
          for (k = i + 1; k < j; k++)
            pscore[indx[k] + i] = pscore[indx[j] + k] = NONE;
          pscore[indx[j] + i] = (psij > 0) ? psij : 0;
        /* fallthrough */
        case '>': /* pairs downstream */
          for (l = j + turn + 1; l <= n; l++)
            pscore[indx[l] + j] = NONE;
          break;
      }
    }
    if (hx != 0)
      vrna_message_error("unbalanced brackets in constraint string\n%s", structure);

    free(stack);
    free(stack2);
  }

  /*free dm */
  for (i = 0; i < 7; i++)
    free(dm[i]);
  free(dm);

  /* copy over pscores for backward compatibility */
  if (vc->pscore_pf_compat) {
    for (i = 1; i < n; i++)
      for (j = i; j <= n; j++)
        vc->pscore_pf_compat[my_iindx[i] - j] = (short)pscore[indx[j] + i];
  }
}
