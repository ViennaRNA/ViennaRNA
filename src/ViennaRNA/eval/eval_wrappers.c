/** \file eval.c */


/*
 *                Free energy evaluation
 *
 *                c Ivo Hofacker, Chrisoph Flamm
 *                original implementation by
 *                Walter Fontana
 *
 *                ViennaRNA Package >= v2.0 by Ronny Lorenz
 *
 *                Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#ifdef _WIN32
#ifdef __MINGW32__
#include <unistd.h>
#else
#include "ViennaRNA/intern/unistd_win.h"
#endif
#else
#include <unistd.h>
#endif

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/eval/structures.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */
#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

PUBLIC int cut_point = -1;  /* set to first pos of second seq for cofolding */
PUBLIC int eos_debug = 0;   /* verbose info from energy_of_struct */

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */
PRIVATE vrna_fold_compound_t *backward_compat_compound = NULL;

#ifdef _OPENMP

#pragma omp threadprivate(backward_compat_compound)

#endif

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE vrna_param_t *
get_updated_params(vrna_param_t *parameters,
                   int          compat);


PRIVATE vrna_fold_compound_t *
recycle_last_call(const char    *string,
                  vrna_param_t  *P);


#endif

PRIVATE INLINE float
eval_structure_simple_v(const char  *string,
                        const char  *structure,
                        int         verbosity_level,
                        int         gquad,
                        int         circular,
                        FILE        *file);


PRIVATE INLINE float
eval_consensus_structure_simple_v(const char  **alignment,
                                  const char  *structure,
                                  int         verbosity_level,
                                  int         gquad,
                                  int         circular,
                                  FILE        *file);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC float
vrna_eval_structure(vrna_fold_compound_t  *fc,
                    const char            *structure)
{
  return vrna_eval_structure_v(fc,
                               structure,
                               VRNA_VERBOSITY_QUIET,
                               NULL);
}


PUBLIC float
vrna_eval_structure_verbose(vrna_fold_compound_t  *fc,
                            const char            *structure,
                            FILE                  *file)
{
  return vrna_eval_structure_v(fc,
                               structure,
                               VRNA_VERBOSITY_DEFAULT,
                               file);
}


PUBLIC float
vrna_eval_structure_simple(const char *string,
                           const char *structure)
{
  return eval_structure_simple_v(string,
                                 structure,
                                 VRNA_VERBOSITY_QUIET,
                                 0,
                                 0,
                                 NULL);
}


PUBLIC float
vrna_eval_consensus_structure_simple(const char **alignment,
                                     const char *structure)
{
  return eval_consensus_structure_simple_v(alignment,
                                           structure,
                                           VRNA_VERBOSITY_QUIET,
                                           0,
                                           0,
                                           NULL);
}


PUBLIC float
vrna_eval_gquad_structure(const char  *string,
                          const char  *structure)
{
  return eval_structure_simple_v(string,
                                 structure,
                                 VRNA_VERBOSITY_QUIET,
                                 1,
                                 0,
                                 NULL);
}


PUBLIC float
vrna_eval_gquad_consensus_structure(const char  **alignment,
                                    const char  *structure)
{
  return eval_consensus_structure_simple_v(alignment,
                                           structure,
                                           VRNA_VERBOSITY_QUIET,
                                           1,
                                           0,
                                           NULL);
}


PUBLIC float
vrna_eval_circ_structure(const char *string,
                         const char *structure)
{
  return eval_structure_simple_v(string,
                                 structure,
                                 VRNA_VERBOSITY_QUIET,
                                 0,
                                 1,
                                 NULL);
}


PUBLIC float
vrna_eval_circ_consensus_structure(const char **alignment,
                                   const char *structure)
{
  return eval_consensus_structure_simple_v(alignment,
                                           structure,
                                           VRNA_VERBOSITY_QUIET,
                                           0,
                                           1,
                                           NULL);
}


PUBLIC float
vrna_eval_circ_gquad_structure(const char *string,
                               const char *structure)
{
  return eval_structure_simple_v(string,
                                 structure,
                                 VRNA_VERBOSITY_QUIET,
                                 1,
                                 1,
                                 NULL);
}


PUBLIC float
vrna_eval_circ_gquad_consensus_structure(const char **alignment,
                                         const char *structure)
{
  return eval_consensus_structure_simple_v(alignment,
                                           structure,
                                           VRNA_VERBOSITY_QUIET,
                                           1,
                                           1,
                                           NULL);
}


PUBLIC float
vrna_eval_structure_simple_verbose(const char *string,
                                   const char *structure,
                                   FILE       *file)
{
  return eval_structure_simple_v(string,
                                 structure,
                                 VRNA_VERBOSITY_DEFAULT,
                                 0,
                                 0,
                                 file);
}


PUBLIC float
vrna_eval_consensus_structure_simple_verbose(const char **alignment,
                                             const char *structure,
                                             FILE       *file)
{
  return eval_consensus_structure_simple_v(alignment,
                                           structure,
                                           VRNA_VERBOSITY_DEFAULT,
                                           0,
                                           0,
                                           file);
}


PUBLIC float
vrna_eval_structure_simple_v(const char *string,
                             const char *structure,
                             int        verbosity_level,
                             FILE       *file)
{
  return eval_structure_simple_v(string,
                                 structure,
                                 verbosity_level,
                                 0,
                                 0,
                                 file);
}


PUBLIC float
vrna_eval_consensus_structure_simple_v(const char **alignment,
                                       const char *structure,
                                       int        verbosity_level,
                                       FILE       *file)
{
  return eval_consensus_structure_simple_v(alignment,
                                           structure,
                                           verbosity_level,
                                           0,
                                           0,
                                           file);
}


PUBLIC float
vrna_eval_circ_structure_v(const char *string,
                           const char *structure,
                           int        verbosity_level,
                           FILE       *file)
{
  return eval_structure_simple_v(string,
                                 structure,
                                 verbosity_level,
                                 0,
                                 1,
                                 file);
}


PUBLIC float
vrna_eval_circ_consensus_structure_v(const char **alignment,
                                     const char *structure,
                                     int        verbosity_level,
                                     FILE       *file)
{
  return eval_consensus_structure_simple_v(alignment,
                                           structure,
                                           verbosity_level,
                                           0,
                                           1,
                                           file);
}


PUBLIC float
vrna_eval_gquad_structure_v(const char  *string,
                            const char  *structure,
                            int         verbosity_level,
                            FILE        *file)
{
  return eval_structure_simple_v(string,
                                 structure,
                                 verbosity_level,
                                 1,
                                 0,
                                 file);
}


PUBLIC float
vrna_eval_gquad_consensus_structure_v(const char  **alignment,
                                      const char  *structure,
                                      int         verbosity_level,
                                      FILE        *file)
{
  return eval_consensus_structure_simple_v(alignment,
                                           structure,
                                           verbosity_level,
                                           1,
                                           0,
                                           file);
}


PUBLIC float
vrna_eval_circ_gquad_structure_v(const char *string,
                                 const char *structure,
                                 int        verbosity_level,
                                 FILE       *file)
{
  return eval_structure_simple_v(string,
                                 structure,
                                 verbosity_level,
                                 1,
                                 1,
                                 file);
}


PUBLIC float
vrna_eval_circ_gquad_consensus_structure_v(const char **alignment,
                                           const char *structure,
                                           int        verbosity_level,
                                           FILE       *file)
{
  return eval_consensus_structure_simple_v(alignment,
                                           structure,
                                           verbosity_level,
                                           1,
                                           1,
                                           file);
}


PUBLIC int
vrna_eval_structure_pt(vrna_fold_compound_t *fc,
                       const short          *pt)
{
  return vrna_eval_structure_pt_v(fc,
                                  pt,
                                  VRNA_VERBOSITY_QUIET,
                                  NULL);
}


PUBLIC int
vrna_eval_structure_pt_verbose(vrna_fold_compound_t *fc,
                               const short          *pt,
                               FILE                 *file)
{
  return vrna_eval_structure_pt_v(fc,
                                  pt,
                                  VRNA_VERBOSITY_DEFAULT,
                                  file);
}


PUBLIC int
vrna_eval_structure_pt_simple(const char  *string,
                              const short *pt)
{
  return vrna_eval_structure_pt_simple_v(string,
                                         pt,
                                         VRNA_VERBOSITY_QUIET,
                                         NULL);
}


PUBLIC int
vrna_eval_consensus_structure_pt_simple(const char  **alignment,
                                        const short *pt)
{
  return vrna_eval_consensus_structure_pt_simple_v(alignment,
                                                   pt,
                                                   VRNA_VERBOSITY_QUIET,
                                                   NULL);
}


PUBLIC int
vrna_eval_structure_pt_simple_verbose(const char  *string,
                                      const short *pt,
                                      FILE        *file)
{
  return vrna_eval_structure_pt_simple_v(string,
                                         pt,
                                         VRNA_VERBOSITY_DEFAULT,
                                         file);
}


PUBLIC int
vrna_eval_consensus_structure_pt_simple_verbose(const char  **alignment,
                                                const short *pt,
                                                FILE        *file)
{
  return vrna_eval_consensus_structure_pt_simple_v(alignment,
                                                   pt,
                                                   VRNA_VERBOSITY_DEFAULT,
                                                   file);
}


PUBLIC int
vrna_eval_structure_pt_simple_v(const char  *string,
                                const short *pt,
                                int         verbosity_level,
                                FILE        *file)
{
  int                   e;
  vrna_fold_compound_t  *fc;

  e = INF;

  if ((string) &&
      (pt)) {
    /* create fold_compound with default parameters and without DP matrices */
    fc = vrna_fold_compound(string, NULL, VRNA_OPTION_EVAL_ONLY);

    /* evaluate structure */
    e = vrna_eval_structure_pt_v(fc, pt, verbosity_level, file);

    /* free fold_compound */
    vrna_fold_compound_free(fc);
  }

  return e;
}


PUBLIC int
vrna_eval_consensus_structure_pt_simple_v(const char  **alignment,
                                          const short *pt,
                                          int         verbosity_level,
                                          FILE        *file)
{
  int                   e;
  vrna_fold_compound_t  *fc;

  e = INF;

  if ((alignment) &&
      (pt)) {
    /* create fold_compound with default parameters and without DP matrices */
    fc = vrna_fold_compound_comparative(alignment, NULL, VRNA_OPTION_DEFAULT);

    /* evaluate structure */
    e = vrna_eval_structure_pt_v(fc, pt, verbosity_level, file);

    /* free fold_compound */
    vrna_fold_compound_free(fc);
  }

  return e;
}


PUBLIC int
vrna_eval_move_pt_simple(const char *string,
                         short      *pt,
                         int        m1,
                         int        m2)
{
  int                   e;
  vrna_fold_compound_t  *fc;

  e = INF;

  if ((string) &&
      (pt)) {
    /* create fold_compound with default parameters and without DP matrices */
    fc = vrna_fold_compound(string, NULL, VRNA_OPTION_EVAL_ONLY);

    /* evaluate structure */
    e = vrna_eval_move_pt(fc, pt, m1, m2);

    /* free fold_compound */
    vrna_fold_compound_free(fc);
  }

  return e;
}


PUBLIC int
vrna_eval_move_shift_pt(vrna_fold_compound_t  *fc,
                        vrna_move_t           *m,
                        short                 *structure)
{
  short       *tmpS;
  int         e, unchangedPosition, insertedPosition, d1, d2, i1, i2;
  vrna_move_t deletion, insertion;

  e = INF;

  if ((fc) &&
      (m) &&
      (structure)) {
    if ((m->pos_5 < 0 && m->pos_3 > 0) ||
        (m->pos_5 > 0 && m->pos_3 < 0)) {
      /* split shift move */
      unchangedPosition = m->pos_5 > 0 ? m->pos_5 : m->pos_3;
      insertedPosition  = m->pos_5 < 0 ? -m->pos_5 : -m->pos_3;
      d1                = -structure[unchangedPosition];
      d2                = -unchangedPosition;

      if (d1 < d2)
        deletion = vrna_move_init(d2, d1);
      else
        deletion = vrna_move_init(d1, d2);

      i1  = unchangedPosition;
      i2  = insertedPosition;

      if (i1 > i2)
        insertion = vrna_move_init(i2, i1);
      else
        insertion = vrna_move_init(i1, i2);

      e     = vrna_eval_move_pt(fc, structure, deletion.pos_5, deletion.pos_3);
      tmpS  = vrna_ptable_copy(structure);
      vrna_move_apply(tmpS, &deletion);
      e += vrna_eval_move_pt(fc, tmpS, insertion.pos_5, insertion.pos_3);
      free(tmpS);
    } else {
      e = vrna_eval_move_pt(fc, structure, m->pos_5, m->pos_3);
    }
  }

  return e;
}


PUBLIC float
vrna_eval_move(vrna_fold_compound_t *fc,
               const char           *structure,
               int                  m1,
               int                  m2)
{
  short *pt;
  int   en;

  en = INF;

  if ((fc) &&
      (structure)) {
    if (strlen(structure) == fc->length) {
      pt  = vrna_ptable(structure);
      en  = vrna_eval_move_pt(fc, pt, m1, m2);
      free(pt);
    } else {
      vrna_log_warning("vrna_eval_move: "
                       "sequence and structure have unequal length (%d vs. %d)",
                       fc->length,
                       strlen(structure));
    }
  }

  return (float)en / 100.;
}


PUBLIC int
vrna_eval_loop_pt(vrna_fold_compound_t  *fc,
                  int                   i,
                  const short           *pt)
{
  return vrna_eval_loop_pt_v(fc, i, pt, VRNA_VERBOSITY_QUIET);
}


/*
 #################################
 # STATIC helper functions below #
 #################################
 */
PRIVATE INLINE float
eval_structure_simple_v(const char  *string,
                        const char  *structure,
                        int         verbosity_level,
                        int         gquad,
                        int         circular,
                        FILE        *file)
{
  char                  *str;
  int                   cp;
  float                 e;
  vrna_fold_compound_t  *fc;
  vrna_md_t             md;

  vrna_md_set_default(&md);

  md.circ   = circular;
  md.gquad  = gquad;

  /* create fold_compound with default parameters and without DP matrices */
  fc = vrna_fold_compound(string, &md, VRNA_OPTION_DEFAULT);

  /* splice-out '&' strand break identifier, if present in structure */
  str = vrna_cut_point_remove(structure, &cp);

  /* evaluate structure */
  e = vrna_eval_structure_v(fc, str, verbosity_level, file);

  /* free fold_compound */
  vrna_fold_compound_free(fc);
  free(str);

  return e;
}


PRIVATE INLINE float
eval_consensus_structure_simple_v(const char  **alignment,
                                  const char  *structure,
                                  int         verbosity_level,
                                  int         gquad,
                                  int         circular,
                                  FILE        *file)
{
  char                  *str;
  int                   cp;
  float                 e;
  vrna_fold_compound_t  *fc;
  vrna_md_t             md;

  vrna_md_set_default(&md);

  md.circ   = circular;
  md.gquad  = gquad;

  /* create fold_compound with default parameters and without DP matrices */
  fc = vrna_fold_compound_comparative(alignment, &md, VRNA_OPTION_DEFAULT);

  /* splice-out '&' strand break identifier, if present in structure */
  str = vrna_cut_point_remove(structure, &cp);

  /* evaluate structure */
  e = vrna_eval_structure_v(fc, str, verbosity_level, file);

  /* free fold_compound */
  vrna_fold_compound_free(fc);
  free(str);

  return e;
}


/*
 #################################
 # DEPRECATED functions below    #
 #################################
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

PRIVATE vrna_fold_compound_t *
recycle_last_call(const char    *string,
                  vrna_param_t  *P)
{
  vrna_fold_compound_t  *fc;
  vrna_md_t             *md;
  int                   cleanup;
  char                  *seq;

  fc      = NULL;
  cleanup = 0;

  if (P) {
    md = &(P->model_details);
  } else {
    md = (vrna_md_t *)vrna_alloc(sizeof(vrna_md_t));
    set_model_details(md);
    cleanup = 1;
  }

  if (string) {
    if (backward_compat_compound) {
      if (!strcmp(string, backward_compat_compound->sequence)) {
        /* check if sequence is the same as before */
        md->window_size = (int)backward_compat_compound->length;
        md->max_bp_span = (int)backward_compat_compound->length;
        /* check if model_details are the same as before */
        if (!memcmp(md, &(backward_compat_compound->params->model_details), sizeof(vrna_md_t)))
          /* re-use previous vrna_fold_compound_t */
          fc = backward_compat_compound;
      }
    }
  }

  /* prepare a new global vrna_fold_compound_t with current settings */
  if (!fc) {
    vrna_fold_compound_free(backward_compat_compound);
    seq                       = vrna_cut_point_insert(string, cut_point);
    backward_compat_compound  = fc = vrna_fold_compound(seq, md, VRNA_OPTION_EVAL_ONLY);
    if (P) {
      free(fc->params);
      fc->params = get_updated_params(P, 1);
    }

    free(seq);
  }

  if (cleanup)
    free(md);

  return fc;
}


PUBLIC float
energy_of_struct(const char *string,
                 const char *structure)
{
  float                 en;
  vrna_fold_compound_t  *fc;

  en = (float)INF / 100.;

  if ((string) &&
      (structure)) {
    fc = recycle_last_call(string, NULL);

    if (eos_debug > 0)
      en = vrna_eval_structure_verbose(fc, structure, NULL);
    else
      en = vrna_eval_structure(fc, structure);
  }

  return en;
}


PUBLIC int
energy_of_struct_pt(const char  *string,
                    short       *pt,
                    short       *s VRNA_UNUSED,
                    short       *s1 VRNA_UNUSED)
{
  int                   en;
  vrna_fold_compound_t  *fc;

  en = INF;

  if ((string) &&
      (pt)) {
    if (pt[0] == (short)strlen(string)) {
      fc  = recycle_last_call(string, NULL);
      en  = vrna_eval_structure_pt_v(fc, pt, eos_debug, NULL);
    } else {
      vrna_log_warning("energy_of_struct_pt: "
                       "string and structure have unequal length (%d vs. %d)",
                       strlen(string),
                       pt[0]);
    }
  }

  return en;
}


PUBLIC float
energy_of_circ_struct(const char  *string,
                      const char  *structure)
{
  float                 en;
  vrna_fold_compound_t  *fc;
  vrna_md_t             *md;

  en = (float)INF / 100.;

  if ((string) &&
      (structure)) {
    fc        = recycle_last_call(string, NULL);
    md        = &(fc->params->model_details);
    md->circ  = 1;

    if (eos_debug > 0)
      en = vrna_eval_structure_verbose(fc, structure, NULL);
    else
      en = vrna_eval_structure(fc, structure);
  }

  return en;
}


PUBLIC float
energy_of_structure(const char  *string,
                    const char  *structure,
                    int         verbosity_level)
{
  float                 en;
  vrna_fold_compound_t  *fc;

  en = (float)INF / 100.;

  if ((string) &&
      (structure)) {
    fc  = recycle_last_call(string, NULL);
    en  = vrna_eval_structure_v(fc, structure, verbosity_level, NULL);
  }

  return en;
}


PUBLIC float
energy_of_struct_par(const char   *string,
                     const char   *structure,
                     vrna_param_t *parameters,
                     int          verbosity_level)
{
  float                 en;
  vrna_fold_compound_t  *fc;

  en = (float)INF / 100.;

  if ((string) &&
      (structure)) {
    fc  = recycle_last_call(string, parameters);
    en  = vrna_eval_structure_v(fc, structure, verbosity_level, NULL);
  }

  return en;
}


PUBLIC float
energy_of_gquad_structure(const char  *string,
                          const char  *structure,
                          int         verbosity_level)
{
  float                 en;
  vrna_fold_compound_t  *fc;
  vrna_md_t             *md;

  en = (float)INF / 100.;

  if ((string) &&
      (structure)) {
    fc        = recycle_last_call(string, NULL);
    md        = &(fc->params->model_details);
    md->gquad = 1;
    en        = vrna_eval_structure_v(fc,
                                      structure,
                                      verbosity_level,
                                      NULL);
  }

  return en;
}


PUBLIC float
energy_of_gquad_struct_par(const char   *string,
                           const char   *structure,
                           vrna_param_t *parameters,
                           int          verbosity_level)
{
  float                 en;
  vrna_fold_compound_t  *fc;
  vrna_md_t             *md;

  en = (float)INF / 100.;

  if ((string) &&
      (structure)) {
    fc        = recycle_last_call(string, parameters);
    md        = &(fc->params->model_details);
    md->gquad = 1;
    en        = vrna_eval_structure_v(fc,
                                      structure,
                                      verbosity_level,
                                      NULL);
  }

  return en;
}


PUBLIC int
energy_of_structure_pt(const char *string,
                       short      *pt,
                       short      *s VRNA_UNUSED,
                       short      *s1 VRNA_UNUSED,
                       int        verbosity_level)
{
  int                   en;
  vrna_fold_compound_t  *fc;

  en = INF;

  if ((string) &&
      (pt)) {
    if (pt[0] == (short)strlen(string)) {
      fc  = recycle_last_call(string, NULL);
      en  = vrna_eval_structure_pt_v(fc, pt, verbosity_level, NULL);
    } else {
      vrna_log_warning("energy_of_structure_pt: "
                       "string and structure have unequal length (%d vs. %d)",
                       strlen(string),
                       pt[0]);
    }
  }

  return en;
}


PUBLIC int
energy_of_struct_pt_par(const char    *string,
                        short         *pt,
                        short         *s VRNA_UNUSED,
                        short         *s1 VRNA_UNUSED,
                        vrna_param_t  *parameters,
                        int           verbosity_level)
{
  int                   en;
  vrna_fold_compound_t  *fc;

  en = INF;

  if ((string) &&
      (pt)) {
    if (pt[0] == (short)strlen(string)) {
      fc  = recycle_last_call(string, parameters);
      en  = vrna_eval_structure_pt_v(fc, pt, verbosity_level, NULL);
    } else {
      vrna_log_warning("energy_of_struct_pt_par: "
                       "string and structure have unequal length (%d vs. %d)",
                       strlen(string),
                       pt[0]);
    }
  }

  return en;
}


PUBLIC float
energy_of_circ_structure(const char *string,
                         const char *structure,
                         int        verbosity_level)
{
  float                 en;
  vrna_fold_compound_t  *fc;
  vrna_md_t             *md;

  en = (float)INF / 100.;

  if ((string) &&
      (structure)) {
    fc        = recycle_last_call(string, NULL);
    md        = &(fc->params->model_details);
    md->circ  = 1;
    en        = vrna_eval_structure_v(fc,
                                      structure,
                                      verbosity_level,
                                      NULL);
  }

  return en;
}


PUBLIC float
energy_of_circ_struct_par(const char    *string,
                          const char    *structure,
                          vrna_param_t  *parameters,
                          int           verbosity_level)
{
  float                 en;
  vrna_fold_compound_t  *fc;
  vrna_md_t             *md;

  en = (float)INF / 100.;

  if ((string) &&
      (structure)) {
    fc        = recycle_last_call(string, parameters);
    md        = &(fc->params->model_details);
    md->circ  = 1;
    en        = vrna_eval_structure_v(fc,
                                      structure,
                                      verbosity_level,
                                      NULL);
  }

  return en;
}


PUBLIC int
loop_energy(short *pt,
            short *s,
            short *s1 VRNA_UNUSED,
            int   i)
{
  char                  *seq;
  int                   en, u;
  vrna_md_t             md;
  vrna_fold_compound_t  *fc;

  en = INF;

  if ((pt) &&
      (s)) {
    set_model_details(&md);

    /* convert encoded sequence back to actual string */
    seq = (char *)vrna_alloc(sizeof(char) * (s[0] + 1));
    for (u = 1; u <= s[0]; u++)
      seq[u - 1] = vrna_nucleotide_decode(s[u], &md);
    seq[u - 1] = '\0';

    fc  = recycle_last_call(seq, NULL);
    en  = vrna_eval_loop_pt_v(fc, i, pt, eos_debug);

    free(seq);
  }

  return en;
}


PUBLIC float
energy_of_move(const char *string,
               const char *structure,
               int        m1,
               int        m2)
{
  float                 en;
  vrna_fold_compound_t  *fc;

  en = (float)INF / 100.;

  if ((string) &&
      (structure)) {
    fc  = recycle_last_call(string, NULL);
    en  = vrna_eval_move(fc, structure, m1, m2);
  }

  return en;
}


PUBLIC int
energy_of_move_pt(short *pt,
                  short *s,
                  short *s1 VRNA_UNUSED,
                  int   m1,
                  int   m2)
{
  int                   en, u;
  char                  *seq;
  vrna_md_t             md;
  vrna_fold_compound_t  *fc;

  en = INF;

  if ((pt) &&
      (s)) {
    set_model_details(&md);

    /* convert encoded sequence back to actual string */
    seq = (char *)vrna_alloc(sizeof(char) * (s[0] + 1));
    for (u = 1; u <= s[0]; u++)
      seq[u - 1] = vrna_nucleotide_decode(s[u], &md);
    seq[u - 1] = '\0';

    fc  = recycle_last_call(seq, NULL);
    en  = vrna_eval_move_pt(fc, pt, m1, m2);

    free(seq);
  }

  return en;
}


PRIVATE vrna_param_t *
get_updated_params(vrna_param_t *parameters,
                   int          compat)
{
  vrna_param_t *P = NULL;

  if (parameters) {
    P = vrna_params_copy(parameters);
  } else {
    vrna_md_t md;
    if (compat)
      set_model_details(&md);
    else
      vrna_md_set_default(&md);

    md.temperature  = temperature;
    P               = vrna_params(&md);
  }

  vrna_md_update(&(P->model_details));
  return P;
}


#endif
