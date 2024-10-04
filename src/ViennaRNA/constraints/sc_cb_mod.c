#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/fold_compound.h"
#include "ViennaRNA/datastructures/string.h"
#include "ViennaRNA/datastructures/array.h"
#include "ViennaRNA/io/utils.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/params/default.h"
#include "json/json.h"

#include "ViennaRNA/constraints/soft_special.h"
#include "ViennaRNA/constraints/sc_cb_intern.h"

#ifndef INLINE
# ifdef __GNUC__
#   define INLINE inline
# else
#   define INLINE
# endif
#endif

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
PRIVATE INLINE void
init_stacks(struct vrna_sc_mod_param_s  *params,
            energy_corrections          *diffs,
            vrna_param_t                *P);


PRIVATE INLINE void
init_mismatches(struct vrna_sc_mod_param_s  *params,
                energy_corrections          *diffs,
                vrna_param_t                *P);


PRIVATE INLINE void
init_dangles(struct vrna_sc_mod_param_s *params,
             energy_corrections         *diffs,
             vrna_param_t               *P);


PRIVATE INLINE void
init_terminal(struct vrna_sc_mod_param_s  *params,
              energy_corrections          *diffs,
              vrna_param_t                *P);


PRIVATE INLINE int
mismatch(vrna_fold_compound_t *fc,
         unsigned int         i,
         unsigned int         j,
         energy_corrections   *data);


PRIVATE INLINE int
terminal(unsigned int       i,
         unsigned int       j,
         energy_corrections *data);


PRIVATE INLINE int
sc_PAIR_HP_terminal(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   k,
                    int                   l,
                    void                  *d);


PRIVATE INLINE int
sc_PAIR_HP_mismatch(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   k,
                    int                   l,
                    void                  *d);


PRIVATE int
sc_PAIR_HP(vrna_fold_compound_t *fc,
           int                  i,
           int                  j,
           int                  k,
           int                  l,
           void                 *d);


PRIVATE INLINE int
sc_PAIR_IL_stack(vrna_fold_compound_t *fc,
                 int                  i,
                 int                  j,
                 int                  k,
                 int                  l,
                 void                 *d);


PRIVATE INLINE int
sc_PAIR_IL_terminal(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   k,
                    int                   l,
                    void                  *d);


PRIVATE INLINE int
sc_PAIR_IL_mismatch(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   k,
                    int                   l,
                    void                  *d);


PRIVATE INLINE int
sc_PAIR_IL_mismatch_terminal(vrna_fold_compound_t *fc,
                             int                  i,
                             int                  j,
                             int                  k,
                             int                  l,
                             void                 *d);


PRIVATE INLINE int
sc_PAIR_IL_stack_terminal(vrna_fold_compound_t  *fc,
                          int                   i,
                          int                   j,
                          int                   k,
                          int                   l,
                          void                  *d);


PRIVATE INLINE int
sc_PAIR_IL_stack_mismatch(vrna_fold_compound_t  *fc,
                          int                   i,
                          int                   j,
                          int                   k,
                          int                   l,
                          void                  *d);


PRIVATE INLINE int
sc_PAIR_IL(vrna_fold_compound_t *fc,
           int                  i,
           int                  j,
           int                  k,
           int                  l,
           void                 *d);


PRIVATE INLINE int
sc_PAIR_ML_terminal(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   k,
                    int                   l,
                    void                  *d);


PRIVATE INLINE int
sc_PAIR_ML_mismatch(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   k,
                    int                   l,
                    void                  *d);


PRIVATE INLINE int
sc_PAIR_ML(vrna_fold_compound_t *fc,
           int                  i,
           int                  j,
           int                  k,
           int                  l,
           void                 *d);


PRIVATE INLINE int
sc_STEM_terminal(vrna_fold_compound_t *fc,
                 int                  i,
                 int                  j,
                 int                  k,
                 int                  l,
                 void                 *d);


PRIVATE INLINE int
sc_STEM_mismatch(vrna_fold_compound_t *fc,
                 int                  i,
                 int                  j,
                 int                  k,
                 int                  l,
                 void                 *d);


PRIVATE INLINE int
sc_STEM(vrna_fold_compound_t  *fc,
        int                   i,
        int                   j,
        int                   k,
        int                   l,
        void                  *d);


PRIVATE INLINE int
sc_EXT_STEM_EXT_terminal(vrna_fold_compound_t *fc,
                         int                  i,
                         int                  j,
                         int                  k,
                         int                  l,
                         void                 *d);


PRIVATE INLINE int
sc_EXT_STEM_EXT_mismatch(vrna_fold_compound_t *fc,
                         int                  i,
                         int                  j,
                         int                  k,
                         int                  l,
                         void                 *d);


PRIVATE INLINE int
sc_EXT_STEM_EXT(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j,
                int                   k,
                int                   l,
                void                  *d);


PRIVATE INLINE int
sc_EXT_EXT_STEM_terminal(vrna_fold_compound_t *fc,
                         int                  i,
                         int                  j,
                         int                  k,
                         int                  l,
                         void                 *d);


PRIVATE INLINE int
sc_EXT_EXT_STEM_mismatch(vrna_fold_compound_t *fc,
                         int                  i,
                         int                  j,
                         int                  k,
                         int                  l,
                         void                 *d);


PRIVATE INLINE int
sc_EXT_EXT_STEM(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j,
                int                   k,
                int                   l,
                void                  *d);


PRIVATE INLINE int
sc_EXT_STEM_OUTSIDE_terminal(vrna_fold_compound_t *fc,
                             int                  i,
                             int                  j,
                             int                  k,
                             int                  l,
                             void                 *d);


PRIVATE INLINE int
sc_EXT_STEM_OUTSIDE_mismatch(vrna_fold_compound_t *fc,
                             int                  i,
                             int                  j,
                             int                  k,
                             int                  l,
                             void                 *d);


PRIVATE INLINE int
sc_EXT_STEM_OUTSIDE(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   k,
                    int                   l,
                    void                  *d);


PRIVATE INLINE int
sc_ML_ML_STEM_terminal(vrna_fold_compound_t *fc,
                       int                  i,
                       int                  j,
                       int                  k,
                       int                  l,
                       void                 *d);


PRIVATE INLINE int
sc_ML_ML_STEM_mismatch(vrna_fold_compound_t *fc,
                       int                  i,
                       int                  j,
                       int                  k,
                       int                  l,
                       void                 *d);


PRIVATE INLINE int
sc_ML_ML_STEM(vrna_fold_compound_t  *fc,
              int                   i,
              int                   j,
              int                   k,
              int                   l,
              void                  *d);


PRIVATE int
prepare_mod_data(vrna_fold_compound_t *fc,
                 void                 *data,
                 unsigned int         event,
                 void                 *event_data);


PRIVATE void
free_energy_corrections(void *d);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_sc_mod_json(vrna_fold_compound_t *fc,
                 const char           *json,
                 const unsigned int   *modification_sites,
                 unsigned int         options)
{
  int                 ret;
  vrna_sc_mod_param_t params;

  ret = 0;

  if ((fc) &&
      (json) &&
      (modification_sites)) {
    params = vrna_sc_mod_read_from_json(json,
                                        &(fc->params->model_details));

    ret = vrna_sc_mod(fc, params, modification_sites, options);

    vrna_sc_mod_parameters_free(params);
  }

  return ret;
}


PUBLIC int
vrna_sc_mod_jsonfile(vrna_fold_compound_t *fc,
                     const char           *jsonfile,
                     const unsigned int   *modification_sites,
                     unsigned int         options)
{
  int                 ret;
  vrna_sc_mod_param_t params;

  ret = 0;

  if ((fc) &&
      (jsonfile) &&
      (modification_sites)) {
    params = vrna_sc_mod_read_from_jsonfile(jsonfile,
                                            &(fc->params->model_details));

    ret = vrna_sc_mod(fc, params, modification_sites, options);

    vrna_sc_mod_parameters_free(params);
  }

  return ret;
}


PUBLIC int
vrna_sc_mod(vrna_fold_compound_t      *fc,
            const vrna_sc_mod_param_t params,
            const unsigned int        *modification_sites,
            unsigned int              options)
{
  int ret = 0;

  if ((fc) &&
      (params) &&
      (modification_sites)) {
    unsigned int        *sn       = fc->strand_number;
    unsigned int        *ss       = fc->strand_start;
    vrna_md_t           *md       = &(fc->params->model_details);
    char                bases[8]  = "_ACGUTM";
    energy_corrections  *diffs;

    bases[6]  = params->one_letter_code;
    diffs     = (energy_corrections *)vrna_alloc(sizeof(energy_corrections));

    /* copy ptypes */
    memcpy(&(diffs->ptypes[0][0]), &(params->ptypes[0][0]), sizeof(params->ptypes));

    diffs->enc = NULL;

     /* store (current) number of strands and prepare per-strand modification site lists */
    diffs->strands = fc->strands;
    vrna_array_init_size(diffs->modification_sites, diffs->strands);
    for (size_t i = 0; i < diffs->strands; i++) {
      vrna_array_make(unsigned int, msites);
      vrna_array_append(diffs->modification_sites, msites);
    }

    /* store modification sites on a per-strand basis */
    for (size_t i = 0; modification_sites[i]; i++) {
      unsigned int msite        = modification_sites[i];
      unsigned int strand       = sn[msite];
      unsigned int actual_msite = msite - ss[strand] + 1;
      unsigned int enc          = fc->sequence_encoding[msite];
      unsigned int unmod        = params->unmodified_encoding;
      unsigned int fallback     = params->fallback_encoding;
      unsigned int pass         = 1;

      if (msite > fc->length) {
        if (!(options & VRNA_SC_MOD_SILENT))
          vrna_log_warning("modification site %u after sequence length (%u)",
                               msite,
                               fc->length);
        continue;
      }

      if (options & (VRNA_SC_MOD_CHECK_UNMOD | VRNA_SC_MOD_CHECK_FALLBACK))
        pass = 0;

      if ((options & VRNA_SC_MOD_CHECK_UNMOD) && (enc == unmod))
        pass = 1;
      else if ((options & VRNA_SC_MOD_CHECK_FALLBACK) && (enc == fallback))
        pass = 1;

      if (!pass) {
        if (!(options & VRNA_SC_MOD_SILENT))
          vrna_log_warning(
            "modification site %u lists wrong unmodified base %c (should be %c)",
            msite,
            bases[fc->sequence_encoding[msite]],
            params->fallback);

        continue;
      }

      vrna_array_append(diffs->modification_sites[strand], actual_msite);

      /* increase return value per modified base we found */
      ret++;

      /* allow for all pairing partners specified in the input */
      for (unsigned int j = 1; j < msite; j++) {
        if ((sn[msite] != sn[j]) ||
            ((msite - j - 1) >= (unsigned int)md->min_loop_size)) {
          for (unsigned int cnt = 0; cnt < params->num_ptypes / 2; cnt++) {
            unsigned int pp_enc = params->pairing_partners_encoding[cnt];
            if ((unsigned int)fc->sequence_encoding[j] == pp_enc) {
              vrna_hc_add_bp(fc,
                             j,
                             msite,
                             VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS | VRNA_CONSTRAINT_CONTEXT_NO_REMOVE);
            }
          }
        }
      }
      for (unsigned int j = msite + 1; j <= fc->length; j++) {
        if ((sn[msite] != sn[j]) ||
            ((j - msite - 1) >= (unsigned int)md->min_loop_size)) {
          for (unsigned int cnt = 0; cnt < params->num_ptypes / 2; cnt++) {
            unsigned int pp_enc = params->pairing_partners_encoding[cnt];
            if ((unsigned int)fc->sequence_encoding[j] == pp_enc) {
              vrna_hc_add_bp(fc,
                             msite,
                             j,
                             VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS | VRNA_CONSTRAINT_CONTEXT_NO_REMOVE);
            }
          }
        }
      }
    }

    init_stacks(params, diffs, fc->params);
    init_terminal(params, diffs, fc->params);
    init_mismatches(params, diffs, fc->params);
    init_dangles(params, diffs, fc->params);

    unsigned int available = params->available;

    /* bind callbacks depending on what data is provided for this modified base */
    if (available & MOD_PARAMS_TERMINAL_dG) {
      if (available & MOD_PARAMS_MISMATCH_dG) {
        vrna_sc_multi_cb_add(fc,
                             &sc_PAIR_HP,
                             NULL,
                             (void *)diffs,
                             &prepare_mod_data,
                             &free_energy_corrections,
                             VRNA_DECOMP_PAIR_HP);

        vrna_sc_multi_cb_add(fc,
                             (available & MOD_PARAMS_STACK_dG) ? &sc_PAIR_IL : sc_PAIR_IL_mismatch_terminal,
                             NULL,
                             (void *)diffs,
                             NULL,
                             NULL,
                             VRNA_DECOMP_PAIR_IL);

        vrna_sc_multi_cb_add(fc,
                             &sc_PAIR_ML,
                             NULL,
                             (void *)diffs,
                             NULL,
                             NULL,
                             VRNA_DECOMP_PAIR_ML);

        vrna_sc_multi_cb_add(fc,
                             &sc_STEM,
                             NULL,
                             (void *)diffs,
                             NULL,
                             NULL,
                             VRNA_DECOMP_EXT_STEM);

        vrna_sc_multi_cb_add(fc,
                             &sc_EXT_STEM_EXT,
                             NULL,
                             (void *)diffs,
                             NULL,
                             NULL,
                             VRNA_DECOMP_EXT_STEM_EXT);

        vrna_sc_multi_cb_add(fc,
                             &sc_EXT_EXT_STEM,
                             NULL,
                             (void *)diffs,
                             NULL,
                             NULL,
                             VRNA_DECOMP_EXT_EXT_STEM);

        vrna_sc_multi_cb_add(fc,
                             &sc_EXT_STEM_OUTSIDE,
                             NULL,
                             (void *)diffs,
                             NULL,
                             NULL,
                             VRNA_DECOMP_EXT_STEM_OUTSIDE);

        vrna_sc_multi_cb_add(fc,
                             &sc_STEM,
                             NULL,
                             (void *)diffs,
                             NULL,
                             NULL,
                             VRNA_DECOMP_ML_STEM);

        vrna_sc_multi_cb_add(fc,
                             &sc_ML_ML_STEM,
                             NULL,
                             (void *)diffs,
                             NULL,
                             NULL,
                             VRNA_DECOMP_ML_ML_STEM);
      } else {
        /* no mismatch energies available */
        vrna_sc_multi_cb_add(fc,
                             &sc_PAIR_HP_terminal,
                             NULL,
                             (void *)diffs,
                             &prepare_mod_data,
                             &free_energy_corrections,
                             VRNA_DECOMP_PAIR_HP);

        vrna_sc_multi_cb_add(fc,
                             (available & MOD_PARAMS_STACK_dG) ? &sc_PAIR_IL_stack_terminal : &sc_PAIR_IL_terminal,
                             NULL,
                             (void *)diffs,
                             NULL,
                             NULL,
                             VRNA_DECOMP_PAIR_IL);

        vrna_sc_multi_cb_add(fc,
                             &sc_PAIR_ML_terminal,
                             NULL,
                             (void *)diffs,
                             NULL,
                             NULL,
                             VRNA_DECOMP_PAIR_ML);

        vrna_sc_multi_cb_add(fc,
                             &sc_STEM_terminal,
                             NULL,
                             (void *)diffs,
                             NULL,
                             NULL,
                             VRNA_DECOMP_EXT_STEM);

        vrna_sc_multi_cb_add(fc,
                             &sc_EXT_STEM_EXT_terminal,
                             NULL,
                             (void *)diffs,
                             NULL,
                             NULL,
                             VRNA_DECOMP_EXT_STEM_EXT);

        vrna_sc_multi_cb_add(fc,
                             &sc_EXT_EXT_STEM_terminal,
                             NULL,
                             (void *)diffs,
                             NULL,
                             NULL,
                             VRNA_DECOMP_EXT_EXT_STEM);

        vrna_sc_multi_cb_add(fc,
                             &sc_EXT_STEM_OUTSIDE_terminal,
                             NULL,
                             (void *)diffs,
                             NULL,
                             NULL,
                             VRNA_DECOMP_EXT_STEM_OUTSIDE);

        vrna_sc_multi_cb_add(fc,
                             &sc_STEM_terminal,
                             NULL,
                             (void *)diffs,
                             NULL,
                             NULL,
                             VRNA_DECOMP_ML_STEM);

        vrna_sc_multi_cb_add(fc,
                             &sc_ML_ML_STEM_terminal,
                             NULL,
                             (void *)diffs,
                             NULL,
                             NULL,
                             VRNA_DECOMP_ML_ML_STEM);
      }
    } else if (available & MOD_PARAMS_MISMATCH_dG) {
      /* no terminalAU-like parameters available */
      vrna_sc_multi_cb_add(fc,
                           &sc_PAIR_HP_mismatch,
                           NULL,
                           (void *)diffs,
                           &prepare_mod_data,
                           &free_energy_corrections,
                           VRNA_DECOMP_PAIR_HP);

      vrna_sc_multi_cb_add(fc,
                           (available & MOD_PARAMS_STACK_dG) ? &sc_PAIR_IL_stack_mismatch : &sc_PAIR_IL_mismatch,
                           NULL,
                           (void *)diffs,
                           NULL,
                           NULL,
                           VRNA_DECOMP_PAIR_IL);

      vrna_sc_multi_cb_add(fc,
                           &sc_PAIR_ML_mismatch,
                           NULL,
                           (void *)diffs,
                           NULL,
                           NULL,
                           VRNA_DECOMP_PAIR_ML);

      vrna_sc_multi_cb_add(fc,
                           &sc_STEM_mismatch,
                           NULL,
                           (void *)diffs,
                           NULL,
                           NULL,
                           VRNA_DECOMP_EXT_STEM);

      vrna_sc_multi_cb_add(fc,
                           &sc_EXT_STEM_EXT_mismatch,
                           NULL,
                           (void *)diffs,
                           NULL,
                           NULL,
                           VRNA_DECOMP_EXT_STEM_EXT);

      vrna_sc_multi_cb_add(fc,
                           &sc_EXT_EXT_STEM_mismatch,
                           NULL,
                           (void *)diffs,
                           NULL,
                           NULL,
                           VRNA_DECOMP_EXT_EXT_STEM);

      vrna_sc_multi_cb_add(fc,
                           &sc_EXT_STEM_OUTSIDE_mismatch,
                           NULL,
                           (void *)diffs,
                           NULL,
                           NULL,
                           VRNA_DECOMP_EXT_STEM_OUTSIDE);

      vrna_sc_multi_cb_add(fc,
                           &sc_STEM_mismatch,
                           NULL,
                           (void *)diffs,
                           NULL,
                           NULL,
                           VRNA_DECOMP_ML_STEM);

      vrna_sc_multi_cb_add(fc,
                           &sc_ML_ML_STEM_mismatch,
                           NULL,
                           (void *)diffs,
                           NULL,
                           NULL,
                           VRNA_DECOMP_ML_ML_STEM);
    } else if (available & MOD_PARAMS_STACK_dG) {
      /* just stacking parameters available */
      vrna_sc_multi_cb_add(fc,
                           &sc_PAIR_IL_stack,
                           NULL,
                           (void *)diffs,
                           &prepare_mod_data,
                           &free_energy_corrections,
                           VRNA_DECOMP_PAIR_IL);
    }
  }

  return ret;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE INLINE void
init_stacks(struct vrna_sc_mod_param_s  *params,
            energy_corrections          *diffs,
            vrna_param_t                *P)
{
  char          nt[MAX_ALPHABET] = {
    '\0', 'A', 'C', 'G', 'U', 'M'
  };
  unsigned int  i, si, sj, enc_unmod, enc_pp, tt, pair_MP, pair_PM;
  int           e,
                (*dG)[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET],
                (*dH)[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET];
  double        tempf;
  vrna_md_t     *md;

  md        = &(P->model_details);
  tempf     = (md->temperature + K0) / (37. + K0);
  enc_unmod = params->fallback_encoding;
  dG        = &(params->stack_dG);
  dH        = &(params->stack_dH);
  nt[5]     = params->one_letter_code;

  if (params->available & MOD_PARAMS_STACK_dG) {
    for (i = 1; i <= params->num_ptypes; i += 2) {
      enc_pp = params->pairing_partners_encoding[(i - 1) / 2];
      /* pair type of unmodified version as encoded in RNAlib */
      pair_MP = md->pair[enc_unmod][enc_pp];
      pair_PM = md->pair[enc_pp][enc_unmod];

      if (pair_MP == 0)
        pair_MP = 7;

      if (pair_PM == 0)
        pair_PM = 7;

      for (si = 1; si < MAX_ALPHABET; si++) {
        for (sj = 1; sj < MAX_ALPHABET; sj++) {
          if (si == 5) {
            if (sj == 5)
              tt = md->pair[enc_unmod][enc_unmod];
            else
              tt = md->pair[sj][enc_unmod];
          } else if (sj == 5) {
            tt = md->pair[enc_unmod][si];
          } else {
            tt = md->pair[sj][si];
          }

          if (tt == 0)
            tt = 7;

          if ((*dG)[i][sj][si] != INF) {
            if (params->available & MOD_PARAMS_STACK_dH)
              e = (*dH)[i][sj][si] - ((*dH)[i][sj][si] - (*dG)[i][sj][si]) * tempf;
            else
              e = (*dG)[i][sj][si];

            diffs->stack_diff[i][sj][si] = e - P->stack[pair_MP][tt];

            vrna_log_debug("d_stack(%c%c, %c%c) = %d = %d - %d",
                   nt[5],
                   nt[enc_pp],
                   nt[sj],
                   nt[si],
                   diffs->stack_diff[i][sj][si],
                   e,
                   P->stack[pair_MP][tt]);
          }

          if ((*dG)[i + 1][sj][si] != INF) {
            if (params->available & MOD_PARAMS_STACK_dH)
              e = (*dH)[i + 1][sj][si] - ((*dH)[i + 1][sj][si] - (*dG)[i + 1][sj][si]) * tempf;
            else
              e = (*dG)[i + 1][sj][si];

            diffs->stack_diff[i + 1][sj][si] = e - P->stack[pair_PM][tt];

            vrna_log_debug("d_stack(%c%c, %c%c) = %d = %d - %d",
                   nt[enc_pp],
                   nt[5],
                   nt[sj],
                   nt[si],
                   diffs->stack_diff[i + 1][sj][si],
                   e,
                   P->stack[pair_PM][tt]);
          }
        }
      }
    }
  }
}


PRIVATE INLINE void
init_mismatches(struct vrna_sc_mod_param_s  *params,
                energy_corrections          *diffs,
                vrna_param_t                *P)
{
  char          nt[MAX_ALPHABET] = {
    '\0', 'A', 'C', 'G', 'U', 'M'
  };

  char          bp[3] = {
    0
  };
  char          *bpairs[8] = {
    "NN", "CG", "GC", "GU", "UG", "AU", "UA", "XX"
  };

  unsigned int  i, si, sj, enc_unmod, enc_pp, siu, sju, pair_enc;
  int           e, (*dG)[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET],
  (*dH)[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET];
  double        tempf;
  vrna_md_t     *md;

  md        = &(P->model_details);
  tempf     = (md->temperature + K0) / (37. + K0);
  enc_unmod = params->fallback_encoding;
  dG        = &(params->mismatch_dG);
  dH        = &(params->mismatch_dH);
  nt[5]     = params->one_letter_code;

  if (params->available & MOD_PARAMS_MISMATCH_dG) {
    /*  go through all enclosing base pair types, including those that
     *  arise from the modification. Here we assume that any pair that
     *  involves a modified base is present in both directions, hence
     *  'params->num_ptypes' == 2 * num_pairing partners of the modified
     *  base.
     */
    for (i = 1; i <= NBPAIRS + params->num_ptypes; i++) {
      if (i <= NBPAIRS) {
        /* 'regular' enclosing pairs */
        pair_enc = i;

        bp[0] = bpairs[i][0];
        bp[1] = bpairs[i][1];
      } else {
        /* an enclosing pair with a modification */
        enc_pp = params->pairing_partners_encoding[(i - NBPAIRS - 1) / 2];
        /* pair type of unmodified version as encoded in RNAlib */
        if ((i - NBPAIRS - 1) % 2) {
          pair_enc = md->pair[enc_unmod][enc_pp];

          bp[1] = nt[5];
          bp[0] = nt[enc_pp];
        } else {
          pair_enc = md->pair[enc_pp][enc_unmod];

          bp[0] = nt[5];
          bp[1] = nt[enc_pp];
        }
        if (pair_enc == 0)
          pair_enc = 7;
      }

      for (si = 1; si < MAX_ALPHABET; si++) {
        for (sj = 1; sj < MAX_ALPHABET; sj++) {
          if (si == 5)
            siu = enc_unmod;
          else
            siu = si;

          if (sj == 5)
            sju = enc_unmod;
          else
            sju = sj;

          if ((*dG)[i][si][sj] != INF) {
            if (params->available & MOD_PARAMS_MISMATCH_dH)
              e = (*dH)[i][si][sj] - ((*dH)[i][si][sj] - (*dG)[i][si][sj]) * tempf;
            else
              e = (*dG)[i][si][sj];

            /*
             * take mismatch energies from multiloops as they do not contain anything
             * but mismatch contributions. Also note, that multiloop mismatches are
             * stored such that unpaired bases are outside of pair, in contrast to
             * considering the unpaired bases enclosed. Thus, the pair and 5' and 3'
             * unpaired bases must be rotatated
             */
            diffs->mismatch_diff[i][si][sj] = e - P->mismatchM[pair_enc][sju][siu];

            vrna_log_debug("d_mm(%c%c, %c, %c) = %d = %d - %d",
                   bp[0],
                   bp[1],
                   nt[si],
                   nt[sj],
                   diffs->mismatch_diff[i][si][sj],
                   e,
                   P->mismatchM[pair_enc][sju][siu]);
          }
        }
      }
    }
  }
}


PRIVATE INLINE void
init_dangles(struct vrna_sc_mod_param_s *params,
             energy_corrections         *diffs,
             vrna_param_t               *P)
{
  char          nt[MAX_ALPHABET] = {
    '\0', 'A', 'C', 'G', 'U', 'M'
  };

  char          bp[3] = {
    0
  };
  char          *bpairs[7] = {
    "NN", "CG", "GC", "GU", "UG", "AU", "UA"
  };
  unsigned int  i, si, enc_unmod, enc_pp, siu, pair_enc;
  int           e, (*dG5)[MAX_PAIRS][MAX_ALPHABET], (*dH5)[MAX_PAIRS][MAX_ALPHABET],
  (*dG3)[MAX_PAIRS][MAX_ALPHABET], (*dH3)[MAX_PAIRS][MAX_ALPHABET];
  double        tempf;
  vrna_md_t     *md;

  md        = &(P->model_details);
  tempf     = (md->temperature + K0) / (37. + K0);
  enc_unmod = params->fallback_encoding;
  dG5       = &(params->dangle5_dG);
  dH5       = &(params->dangle5_dH);
  dG3       = &(params->dangle3_dG);
  dH3       = &(params->dangle3_dH);
  nt[5]     = params->one_letter_code;

  if (params->available & MOD_PARAMS_DANGLES_dG) {
    /* process all closing pairs without modified bases */
    for (i = 1; i <= NBPAIRS + params->num_ptypes; i++) {
      if (i <= NBPAIRS) {
        /* 'regular' enclosing pairs */
        pair_enc = i;

        bp[0] = bpairs[i][0];
        bp[1] = bpairs[i][1];
      } else {
        /* an enclosing pair with a modification */
        enc_pp = params->pairing_partners_encoding[(i - NBPAIRS - 1) / 2];
        /* pair type of unmodified version as encoded in RNAlib */
        if ((i - NBPAIRS - 1) % 2) {
          pair_enc = md->pair[enc_unmod][enc_pp];

          bp[1] = nt[5];
          bp[0] = nt[enc_pp];
        } else {
          pair_enc = md->pair[enc_pp][enc_unmod];

          bp[0] = nt[5];
          bp[1] = nt[enc_pp];
        }
        if (pair_enc == 0)
          pair_enc = 7;
      }

      for (si = 1; si < MAX_ALPHABET; si++) {
        if (si == 5)
          siu = enc_unmod;
        else
          siu = si;

        if ((*dG5)[i][si] != INF) {
          if (params->available & MOD_PARAMS_DANGLES_dH)
            e = (*dH5)[i][si] - ((*dH5)[i][si] - (*dG5)[i][si]) * tempf;
          else
            e = (*dG5)[i][si];

          /*
           * Note, that dangling ends arestored such that unpaired bases are
           * outside of pair, in contrast to considering the unpaired bases
           * enclosed. Thus, the pair must be rotatated
           */
          diffs->dangle5_diff[i][si] = e - P->dangle5[pair_enc][siu];

          vrna_log_debug("d_d5(%c%c, %c) = %d = %d - %d",
                 bp[0],
                 bp[1],
                 nt[si],
                 diffs->dangle5_diff[i][si],
                 e,
                 P->dangle5[pair_enc][siu]);
        }

        if ((*dG3)[i][si] != INF) {
          if (params->available & MOD_PARAMS_DANGLES_dH)
            e = (*dH3)[i][si] - ((*dH3)[i][si] - (*dG3)[i][si]) * tempf;
          else
            e = (*dG3)[i][si];

          /*
           * Note, that dangling ends arestored such that unpaired bases are
           * outside of pair, in contrast to considering the unpaired bases
           * enclosed. Thus, the pair must be rotatated
           */
          diffs->dangle3_diff[i][si] = e - P->dangle3[pair_enc][siu];

          vrna_log_debug("d_d3(%c%c, %c) = %d = %d - %d",
                 bp[0],
                 bp[1],
                 nt[si],
                 diffs->dangle3_diff[i][si],
                 e,
                 P->dangle3[pair_enc][siu]);
        }
      }
    }
  }
}


PRIVATE INLINE void
init_terminal(struct vrna_sc_mod_param_s  *params,
              energy_corrections          *diffs,
              vrna_param_t                *P)
{
  char          nt[MAX_ALPHABET] = {
    '\0', 'A', 'C', 'G', 'U', 'M'
  };
  unsigned int  i, enc_unmod, enc_pp, tt;
  int           e, (*dG)[MAX_PAIRS], (*dH)[MAX_PAIRS], Terminal_unmod;
  double        tempf;
  vrna_md_t     *md;

  md        = &(P->model_details);
  tempf     = (md->temperature + K0) / (37. + K0);
  enc_unmod = params->fallback_encoding;
  dG        = &(params->terminal_dG);
  dH        = &(params->terminal_dH);
  nt[5]     = params->one_letter_code;

  if (params->available & MOD_PARAMS_TERMINAL_dG) {
    for (i = 1; i <= params->num_ptypes; i += 2) {
      enc_pp          = params->pairing_partners_encoding[(i - 1) / 2];
      tt              = md->pair[enc_unmod][enc_pp];
      Terminal_unmod  = (tt > 2) ? P->TerminalAU : 0; /* terminal AU only for AU and GU pairs */

      if ((*dG)[i] != INF) {
        if (params->available & MOD_PARAMS_TERMINAL_dH)
          e = (*dH)[i] - ((*dH)[i] - (*dG)[i]) * tempf;
        else
          e = (*dG)[i];

        diffs->terminal_diff[i] = e - Terminal_unmod;

        vrna_log_debug("d_term(%c%c) = %d = %d - %d",
               nt[5],
               nt[enc_pp],
               diffs->terminal_diff[i],
               e,
               Terminal_unmod);
      }

      if ((*dG)[i + 1] != INF) {
        if (params->available & MOD_PARAMS_TERMINAL_dH)
          e = (*dH)[i + 1] - ((*dH)[i + 1] - (*dG)[i + 1]) * tempf;
        else
          e = (*dG)[i + 1];

        diffs->terminal_diff[i + 1] = e - Terminal_unmod;

        vrna_log_debug("d_term(%c%c) = %d = %d - %d",
               nt[enc_pp],
               nt[5],
               diffs->terminal_diff[i + 1],
               e,
               Terminal_unmod);
      }
    }
  }
}


PRIVATE INLINE int
mismatch(vrna_fold_compound_t *fc,
         unsigned int         i,
         unsigned int         j,
         energy_corrections   *data)
{
  short         *enc  = data->enc;
  unsigned int  tt    = data->ptypes[enc[i]][enc[j]];
  vrna_md_t     *md   = &(fc->params->model_details);

  if (tt == 0)
    /* if we don't know the base pair, it must be canonical */
    tt = md->pair[enc[i]][enc[j]];
  else
    tt += NBPAIRS;

  if (j > 1) {
    if (i < fc->length)
      return data->mismatch_diff[tt][enc[i + 1]][enc[j - 1]];
    else
      return data->dangle5_diff[tt][enc[j - 1]];
  } else if (i < fc->length) {
    return data->dangle3_diff[tt][enc[i + 1]];
  }

  return 0;
}


PRIVATE INLINE int
terminal(unsigned int       i,
         unsigned int       j,
         energy_corrections *data)
{
  short         *enc  = data->enc;
  unsigned int  tt    = data->ptypes[enc[i]][enc[j]];

  return data->terminal_diff[tt];
}


/* hairpin loop correction including terminalAU terms */
PRIVATE INLINE int
sc_PAIR_HP_terminal(vrna_fold_compound_t  *fc VRNA_UNUSED,
                    int                   i,
                    int                   j,
                    int                   k VRNA_UNUSED,
                    int                   l VRNA_UNUSED,
                    void                  *d)
{
  return terminal(i, j,
                  (energy_corrections *)d);
}


/* hairpin loop correction including mismatch terms */
PRIVATE INLINE int
sc_PAIR_HP_mismatch(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   k VRNA_UNUSED,
                    int                   l VRNA_UNUSED,
                    void                  *d)
{
  return mismatch(fc,
                  i, j,
                  (energy_corrections *)d);
}


/* hairpin loop correction including mismatch and terminalAU terms */
PRIVATE int
sc_PAIR_HP(vrna_fold_compound_t *fc,
           int                  i,
           int                  j,
           int                  k,
           int                  l,
           void                 *d)
{
  return sc_PAIR_HP_mismatch(fc, i, j, k, l, d) +
         sc_PAIR_HP_terminal(fc, i, j, k, l, d);
}


PRIVATE INLINE int
sc_PAIR_IL_stack(vrna_fold_compound_t *fc VRNA_UNUSED,
                 int                  i,
                 int                  j,
                 int                  k,
                 int                  l,
                 void                 *d)
{
  if ((i + 1 == k) &&
      (l == j - 1)) {
    energy_corrections  *data = (energy_corrections *)d;
    short               *enc  = data->enc;

    unsigned int        enc_i = enc[i];
    unsigned int        enc_j = enc[j];
    unsigned int        enc_k = enc[k];
    unsigned int        enc_l = enc[l];
    unsigned int        t1    = data->ptypes[enc_i][enc_j];
    unsigned int        t2    = data->ptypes[enc_l][enc_k];
    if (t1 != 0)
      return data->stack_diff[t1][enc_l][enc_k];
    else if (t2 != 0)
      return data->stack_diff[t2][enc_i][enc_j];
  }

  return 0;
}


PRIVATE INLINE int
sc_PAIR_IL_terminal(vrna_fold_compound_t  *fc VRNA_UNUSED,
                    int                   i,
                    int                   j,
                    int                   k,
                    int                   l,
                    void                  *d)
{
  if ((i + 1 < k) ||
      (l + 1 < j))
    return terminal(i, j, (energy_corrections *)d) +
           terminal(l, k, (energy_corrections *)d);

  return 0;
}


PRIVATE INLINE int
sc_PAIR_IL_mismatch(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   k,
                    int                   l,
                    void                  *d)
{
  if (((k - i) > 1) &&
      ((j - l) > 1))
    return mismatch(fc, i, j, (energy_corrections *)d) +
           mismatch(fc, l, k, (energy_corrections *)d);

  return 0;
}


PRIVATE INLINE int
sc_PAIR_IL_mismatch_terminal(vrna_fold_compound_t *fc,
                             int                  i,
                             int                  j,
                             int                  k,
                             int                  l,
                             void                 *d)
{
  return sc_PAIR_IL_mismatch(fc, i, j, k, l, d) +
         sc_PAIR_IL_terminal(fc, i, j, k, l, d);
}


PRIVATE INLINE int
sc_PAIR_IL_stack_terminal(vrna_fold_compound_t  *fc,
                          int                   i,
                          int                   j,
                          int                   k,
                          int                   l,
                          void                  *d)
{
  return sc_PAIR_IL_stack(fc, i, j, k, l, d) +
         sc_PAIR_IL_terminal(fc, i, j, k, l, d);
}


PRIVATE INLINE int
sc_PAIR_IL_stack_mismatch(vrna_fold_compound_t  *fc,
                          int                   i,
                          int                   j,
                          int                   k,
                          int                   l,
                          void                  *d)
{
  return sc_PAIR_IL_stack(fc, i, j, k, l, d) +
         sc_PAIR_IL_mismatch(fc, i, j, k, l, d);
}


PRIVATE INLINE int
sc_PAIR_IL(vrna_fold_compound_t *fc,
           int                  i,
           int                  j,
           int                  k,
           int                  l,
           void                 *d)
{
  return sc_PAIR_IL_stack(fc, i, j, k, l, d) +
         sc_PAIR_IL_mismatch(fc, i, j, k, l, d) +
         sc_PAIR_IL_terminal(fc, i, j, k, l, d);
}


PRIVATE INLINE int
sc_PAIR_ML_terminal(vrna_fold_compound_t  *fc VRNA_UNUSED,
                    int                   i,
                    int                   j,
                    int                   k VRNA_UNUSED,
                    int                   l VRNA_UNUSED,
                    void                  *d)
{
  return terminal(i, j,
                  (energy_corrections *)d);
}


PRIVATE INLINE int
sc_PAIR_ML_mismatch(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   k VRNA_UNUSED,
                    int                   l VRNA_UNUSED,
                    void                  *d)
{
  return mismatch(fc,
                  i, j,
                  (energy_corrections *)d);
}


PRIVATE INLINE int
sc_PAIR_ML(vrna_fold_compound_t *fc,
           int                  i,
           int                  j,
           int                  k,
           int                  l,
           void                 *d)
{
  return sc_PAIR_ML_mismatch(fc, i, j, k, l, d) +
         sc_PAIR_ML_terminal(fc, i, j, k, l, d);
}


PRIVATE INLINE int
sc_STEM_terminal(vrna_fold_compound_t *fc VRNA_UNUSED,
                 int                  i VRNA_UNUSED,
                 int                  j VRNA_UNUSED,
                 int                  k,
                 int                  l,
                 void                 *d)
{
  return terminal(l, k,
                  (energy_corrections *)d);
}


PRIVATE INLINE int
sc_STEM_mismatch(vrna_fold_compound_t *fc,
                 int                  i VRNA_UNUSED,
                 int                  j VRNA_UNUSED,
                 int                  k,
                 int                  l,
                 void                 *d)
{
  return mismatch(fc,
                  l, k,
                  (energy_corrections *)d);
}


PRIVATE INLINE int
sc_STEM(vrna_fold_compound_t  *fc,
        int                   i,
        int                   j,
        int                   k,
        int                   l,
        void                  *d)
{
  return sc_STEM_mismatch(fc, i, j, k, l, d) +
         sc_STEM_terminal(fc, i, j, k, l, d);
}


PRIVATE INLINE int
sc_EXT_STEM_EXT_terminal(vrna_fold_compound_t *fc VRNA_UNUSED,
                         int                  i,
                         int                  j VRNA_UNUSED,
                         int                  k,
                         int                  l VRNA_UNUSED,
                         void                 *d)
{
  return terminal(k, i,
                  (energy_corrections *)d);
}


PRIVATE INLINE int
sc_EXT_STEM_EXT_mismatch(vrna_fold_compound_t *fc,
                         int                  i,
                         int                  j VRNA_UNUSED,
                         int                  k,
                         int                  l VRNA_UNUSED,
                         void                 *d)
{
  return mismatch(fc,
                  k, i,
                  (energy_corrections *)d);
}


PRIVATE INLINE int
sc_EXT_STEM_EXT(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j,
                int                   k,
                int                   l,
                void                  *d)
{
  return sc_EXT_STEM_EXT_mismatch(fc, i, j, k, l, d) +
         sc_EXT_STEM_EXT_terminal(fc, i, j, k, l, d);
}


PRIVATE INLINE int
sc_EXT_EXT_STEM_terminal(vrna_fold_compound_t *fc VRNA_UNUSED,
                         int                  i VRNA_UNUSED,
                         int                  j,
                         int                  k VRNA_UNUSED,
                         int                  l,
                         void                 *d)
{
  return terminal(j, l,
                  (energy_corrections *)d);
}


PRIVATE INLINE int
sc_EXT_EXT_STEM_mismatch(vrna_fold_compound_t *fc,
                         int                  i VRNA_UNUSED,
                         int                  j,
                         int                  k VRNA_UNUSED,
                         int                  l,
                         void                 *d)
{
  return mismatch(fc,
                  j, l,
                  (energy_corrections *)d);
}


PRIVATE INLINE int
sc_EXT_EXT_STEM(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j,
                int                   k,
                int                   l,
                void                  *d)
{
  return sc_EXT_EXT_STEM_mismatch(fc, i, j, k, l, d) +
         sc_EXT_EXT_STEM_terminal(fc, i, j, k, l, d);
}


PRIVATE INLINE int
sc_EXT_STEM_OUTSIDE_terminal(vrna_fold_compound_t *fc VRNA_UNUSED,
                             int                  i VRNA_UNUSED,
                             int                  j VRNA_UNUSED,
                             int                  k,
                             int                  l,
                             void                 *d)
{
  return terminal(l, k,
                  (energy_corrections *)d);
}


PRIVATE INLINE int
sc_EXT_STEM_OUTSIDE_mismatch(vrna_fold_compound_t *fc,
                             int                  i VRNA_UNUSED,
                             int                  j VRNA_UNUSED,
                             int                  k,
                             int                  l,
                             void                 *d)
{
  return mismatch(fc,
                  l, k,
                  (energy_corrections *)d);
}


PRIVATE INLINE int
sc_EXT_STEM_OUTSIDE(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   k,
                    int                   l,
                    void                  *d)
{
  return sc_EXT_STEM_OUTSIDE_mismatch(fc, i, j, k, l, d) +
         sc_EXT_STEM_OUTSIDE_terminal(fc, i, j, k, l, d);
}


PRIVATE INLINE int
sc_ML_ML_STEM_terminal(vrna_fold_compound_t *fc VRNA_UNUSED,
                       int                  i VRNA_UNUSED,
                       int                  j,
                       int                  k VRNA_UNUSED,
                       int                  l,
                       void                 *d)
{
  return terminal(j, l,
                  (energy_corrections *)d);
}


PRIVATE INLINE int
sc_ML_ML_STEM_mismatch(vrna_fold_compound_t *fc,
                       int                  i VRNA_UNUSED,
                       int                  j,
                       int                  k VRNA_UNUSED,
                       int                  l,
                       void                 *d)
{
  return mismatch(fc,
                  j, l,
                  (energy_corrections *)d);
}


PRIVATE INLINE int
sc_ML_ML_STEM(vrna_fold_compound_t  *fc,
              int                   i,
              int                   j,
              int                   k,
              int                   l,
              void                  *d)
{
  return sc_ML_ML_STEM_mismatch(fc, i, j, k, l, d) +
         sc_ML_ML_STEM_terminal(fc, i, j, k, l, d);
}


PRIVATE int
prepare_mod_data(vrna_fold_compound_t *fc,
                 void                 *data,
                 unsigned int         event,
                 void                 *event_data)
{
  int ret = 0;

  energy_corrections *diff = (energy_corrections *)data;

  /*  Update the encoding array.
   *
   *  In case we are about to perform sliding-window predictions, only
   *  update the array once at the beginning of the computations
   */
  if ((!(event & VRNA_OPTION_WINDOW)) ||
      ((event & VRNA_OPTION_F3) && ((*(unsigned int *)event_data) == fc->length)) ||
      ((event & VRNA_OPTION_F5) && ((*(unsigned int *)event_data) == 1)) ||
      (diff->enc == NULL)) {
    unsigned int  *ss, *so, strand;
    size_t        i, j, k;

    ss = fc->strand_start;
    so = fc->strand_order;

    free(diff->enc);
    diff->enc = (short *)vrna_alloc(sizeof(short) * (fc->length + 2));

    if (!diff->enc)
      return 1;

    memcpy(diff->enc, fc->sequence_encoding, sizeof(short) * (fc->length + 1));

    /* correct for all known modification sites */
    for (i = 0; i < fc->strands; i++) {
      strand = so[i];

      if (strand > vrna_array_size(diff->modification_sites))
        return 1; /* return with non-zero value to indicate error */

      for (j = 0; j < vrna_array_size(diff->modification_sites[strand]); j++) {
        k = diff->modification_sites[strand][j] + ss[strand] - 1;
        diff->enc[k] = 5;
      }
    }
  }

  return ret;
}


PRIVATE void
free_energy_corrections(void *d)
{
  energy_corrections *diff = (energy_corrections *)d;

  for (size_t i = 0; i < vrna_array_size(diff->modification_sites); i++)
    vrna_array_free(diff->modification_sites[i]);
  vrna_array_free(diff->modification_sites);

  free(diff->enc);
  free(diff);
}
