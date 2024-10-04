/* constraints handling */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/params/default.h"
#include "ViennaRNA/params/constants.h" /* defines MINPSCORE */
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/sequences/alignments.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/io/file_formats.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/probing/SHAPE.h"
#include "ViennaRNA/constraints/soft.h"


#ifndef INLINE
# ifdef __GNUC__
#   define INLINE inline
# else
#   define INLINE
# endif
#endif

#define STATE_CLEAN         (unsigned char)0
#define STATE_DIRTY_UP_MFE  (unsigned char)1
#define STATE_DIRTY_UP_PF   (unsigned char)2
#define STATE_DIRTY_BP_MFE  (unsigned char)4
#define STATE_DIRTY_BP_PF   (unsigned char)8

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
sc_reset_up(vrna_fold_compound_t  *fc,
            const FLT_OR_DBL      *constraints,
            unsigned int          options);


PRIVATE int
sc_reset_up_comparative(vrna_fold_compound_t  *fc,
                        const FLT_OR_DBL      **constraints,
                        unsigned int          options);


PRIVATE void
sc_reset_bp(vrna_fold_compound_t  *fc,
            const FLT_OR_DBL      **constraints,
            unsigned int          options);


PRIVATE int
sc_reset_bp_comparative(vrna_fold_compound_t  *fc,
                        const FLT_OR_DBL      ***constraints,
                        unsigned int          options);


PRIVATE void
sc_add_up(vrna_fold_compound_t  *fc,
          unsigned int          i,
          FLT_OR_DBL            energy,
          unsigned int          options);


PRIVATE int
sc_add_up_comparative(vrna_fold_compound_t  *fc,
                      unsigned int          *is,
                      const FLT_OR_DBL      *energies,
                      unsigned int          options);


PRIVATE void
sc_add_bp(vrna_fold_compound_t  *fc,
          unsigned int          i,
          unsigned int          j,
          FLT_OR_DBL            energy,
          unsigned int          options);


PRIVATE int
sc_add_bp_comparative(vrna_fold_compound_t  *fc,
                      unsigned int          *is,
                      unsigned int          *js,
                      const FLT_OR_DBL      *energies,
                      unsigned int          options);


PRIVATE void
prepare_sc_up_mfe(vrna_sc_t     *sc,
                  unsigned int  n,
                  unsigned int  options);


PRIVATE void
prepare_sc_bp_mfe(vrna_sc_t     *sc,
                  unsigned int  n,
                  int           *idx,
                  unsigned int  options);


PRIVATE void
prepare_sc_up_pf(vrna_sc_t    *sc,
                 unsigned int n,
                 double       kT,
                 unsigned int options);


PRIVATE void
prepare_sc_bp_pf(vrna_sc_t    *sc,
                 unsigned int n,
                 int          *idx,
                 double       kT,
                 unsigned int options);


PRIVATE void
prepare_sc_stack_pf(vrna_sc_t     *fc,
                    unsigned int  n,
                    double        kT);


PRIVATE int
prepare_sc_user_cb(vrna_fold_compound_t *fc,
                   vrna_sc_t            *sc,
                   unsigned int         options);


PRIVATE INLINE void
sc_init_up_storage(vrna_sc_t *sc);


PRIVATE INLINE void
populate_sc_up_mfe(vrna_sc_t    *sc,
                   unsigned int i,
                   unsigned int n);


PRIVATE INLINE void
populate_sc_up_pf(vrna_sc_t     *sc,
                  unsigned int  i,
                  unsigned int  n,
                  double        kT);


PRIVATE INLINE void
sc_init_bp_storage(vrna_sc_t *sc);


PRIVATE INLINE void
sc_store_bp(vrna_sc_bp_storage_t  **container,
            unsigned int          i,
            unsigned int          start,
            unsigned int          end,
            int                   e);


PRIVATE INLINE int
get_stored_bp_contributions(vrna_sc_bp_storage_t  *container,
                            unsigned int          j);


PRIVATE INLINE void
populate_sc_bp_mfe(vrna_sc_t    *sc,
                   unsigned int i,
                   unsigned int maxdist,
                   unsigned int n,
                   int          *idx);


PRIVATE INLINE void
populate_sc_bp_pf(vrna_sc_t     *sc,
                  unsigned int  i,
                  unsigned int  maxdist,
                  unsigned int  n,
                  int           *idx,
                  double        kT);


PRIVATE INLINE void
free_sc_up(vrna_sc_t *sc);


PRIVATE INLINE void
free_sc_bp(vrna_sc_t *sc);


PRIVATE vrna_sc_t *
init_sc_default(unsigned int n);


PRIVATE vrna_sc_t *
init_sc_window(unsigned int n);


PRIVATE void
nullify(vrna_sc_t *sc);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC void
vrna_sc_init(vrna_fold_compound_t *fc)
{
  unsigned int s, n, N;

  if (fc) {
    vrna_sc_remove(fc);

    n = fc->length;

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        fc->sc = init_sc_default(n);
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        N       = fc->n_seq;
        fc->scs = (vrna_sc_t **)vrna_alloc(sizeof(vrna_sc_t *) * (N + 1));

        for (s = 0; s < N; s++)
          /*  except for base pairs, constraints for individual sequences
           *  in the MSA are stored per actual nucleotides, so we should
           *  actually initialize with the 'gap-free' size of each sequence
           *  However, for now we initialize with alignment columns until
           *  usage of base pair constraints has been adjusted accordingly
           *
           *  This introduces some waste of memory for gap positions in
           *  unpaired and stack constraints, though
           */
          fc->scs[s] = init_sc_default(n);

        break;
    }
  }
}


PUBLIC void
vrna_sc_init_window(vrna_fold_compound_t *fc)
{
  unsigned int s, n, N;

  if (fc) {
    vrna_sc_remove(fc);

    n = fc->length;

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        fc->sc = init_sc_window(n);
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        N       = fc->n_seq;
        fc->scs = (vrna_sc_t **)vrna_alloc(sizeof(vrna_sc_t *) * (N + 1));

        for (s = 0; s < N; s++)
          /*  except for base pairs, constraints for individual sequences
           *  in the MSA are stored per actual nucleotides, so we should
           *  actually initialize with the 'gap-free' size of each sequence
           *  However, for now we initialize with alignment columns until
           *  usage of base pair constraints has been adjusted accordingly
           *
           *  This introduces some waste of memory for gap positions in
           *  unpaired and stack constraints, though
           */
          fc->scs[s] = init_sc_window(n);

        break;
    }
  }
}


PUBLIC int
vrna_sc_prepare(vrna_fold_compound_t  *fc,
                unsigned int          options)
{
  int ret = 0;

  if (fc) {
    unsigned int  n, s;

    n = fc->length;
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (options & VRNA_OPTION_MFE) {
          prepare_sc_up_mfe(fc->sc, n, options);
          prepare_sc_bp_mfe(fc->sc, n, fc->jindx, options);
        }

        if (options & VRNA_OPTION_PF) {
          prepare_sc_up_pf(fc->sc, n, fc->exp_params->kT, options);
          prepare_sc_bp_pf(fc->sc, n, fc->jindx, fc->exp_params->kT, options);
          prepare_sc_stack_pf(fc->sc, n, fc->exp_params->kT);
        }

        ret &= prepare_sc_user_cb(fc, fc->sc, options);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        if (fc->scs) {
          if (options & VRNA_OPTION_MFE) {
            for (s = 0; s < fc->n_seq; s++) {
              prepare_sc_up_mfe(fc->scs[s], fc->a2s[s][n], options);
              prepare_sc_bp_mfe(fc->scs[s], n, fc->jindx, options);
            }
          }

          if (options & VRNA_OPTION_PF) {
            for (s = 0; s < fc->n_seq; s++) {
              prepare_sc_up_pf(fc->scs[s], fc->a2s[s][n], fc->exp_params->kT, options);
              prepare_sc_bp_pf(fc->scs[s], n, fc->jindx, fc->exp_params->kT, options);
              prepare_sc_stack_pf(fc->scs[s], fc->a2s[s][n], fc->exp_params->kT);
              ret &= prepare_sc_user_cb(fc, fc->scs[s], options);
            }
          }
        }

        break; /* do nothing for now */
    }

  }

  return ret;
}


PUBLIC int
vrna_sc_update(vrna_fold_compound_t *fc,
               unsigned int         i,
               unsigned int         options)
{
  unsigned int  n, maxdist, s, **a2s;
  int           ret = 0;
  vrna_sc_t     *sc, **scs;

  if (fc) {
    n       = fc->length;
    maxdist = (unsigned int)fc->window_size;

    if (i > n) {
      vrna_log_warning("vrna_sc_update(): Position %u out of range!"
                           " (Sequence length: %u)",
                           i, n);
    } else if (i > 0) {
      maxdist = MIN2(maxdist, n - i + 1);
      ret = 1;

      switch (fc->type) {
        case VRNA_FC_TYPE_SINGLE:
          sc = fc->sc;

          if ((sc) &&
              (options & VRNA_OPTION_WINDOW)) {
            /* sliding-window mode, i.e. local structure prediction */
            if (sc->up_storage) {
              if (options & VRNA_OPTION_MFE)
                populate_sc_up_mfe(sc, i, maxdist);

              if (options & VRNA_OPTION_PF)
                populate_sc_up_pf(sc, i, maxdist, fc->exp_params->kT);
            }

            if (sc->bp_storage) {
              if (options & VRNA_OPTION_MFE)
                populate_sc_bp_mfe(fc->sc, i, maxdist, n, fc->jindx);

              if (options & VRNA_OPTION_PF)
                populate_sc_bp_pf(fc->sc, i, maxdist, n, fc->jindx, fc->exp_params->kT);
            }

            if ((sc->data) && (sc->prepare_data))
              ret &= sc->prepare_data(fc, sc->data, options, (void *)&i);
          }

          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          scs = fc->scs;
          a2s = fc->a2s;

          for (s = 0; s < fc->n_seq; s++) {
            if ((scs[s]) &&
                (scs[s]->type == VRNA_SC_WINDOW)) {
              if (scs[s]->up_storage) {
                if (options & VRNA_OPTION_MFE)
                  populate_sc_up_mfe(scs[s], a2s[s][i], a2s[s][maxdist]);

                if (options & VRNA_OPTION_PF)
                  populate_sc_up_pf(scs[s], a2s[s][i], a2s[s][maxdist], fc->exp_params->kT);
              }

              if (scs[s]->bp_storage) {
                if (options & VRNA_OPTION_MFE)
                  populate_sc_bp_mfe(fc->scs[s], i, maxdist, n, fc->jindx);

                if (options & VRNA_OPTION_PF)
                  populate_sc_bp_pf(fc->scs[s], i, maxdist, n, fc->jindx, fc->exp_params->kT);
              }

              if ((scs[s]->data) && (scs[s]->prepare_data))
                ret &= scs[s]->prepare_data(fc, scs[s]->data, options, (void *)&i);
              }
          }

          break;
      }
    }
  }

  return ret;
}


PUBLIC void
vrna_sc_remove(vrna_fold_compound_t *fc)
{
  unsigned int s;

  if (fc) {
    switch (fc->type) {
      case  VRNA_FC_TYPE_SINGLE:
        vrna_sc_free(fc->sc);
        fc->sc = NULL;
        break;

      case  VRNA_FC_TYPE_COMPARATIVE:
        if (fc->scs) {
          for (s = 0; s < fc->n_seq; s++)
            vrna_sc_free(fc->scs[s]);
          free(fc->scs);
        }

        fc->scs = NULL;
        break;
    }
  }
}


PUBLIC void
vrna_sc_free(vrna_sc_t *sc)
{
  if (sc) {
    free_sc_up(sc);
    free_sc_bp(sc);

    free(sc->energy_stack);
    free(sc->exp_energy_stack);

    if (sc->free_data)
      sc->free_data(sc->data);

    free(sc);
  }
}


PUBLIC int
vrna_sc_set_bp(vrna_fold_compound_t *fc,
               const FLT_OR_DBL     **constraints,
               unsigned int         options)
{
  if (fc && (fc->type == VRNA_FC_TYPE_SINGLE)) {
    sc_reset_bp(fc, constraints, options);

    if (options & VRNA_OPTION_MFE)
      prepare_sc_bp_mfe(fc->sc, fc->length, fc->jindx, options);

    if (options & VRNA_OPTION_PF) /* prepare Boltzmann factors for the BP soft constraints */
      prepare_sc_bp_pf(fc->sc, fc->length, fc->jindx, fc->exp_params->kT, options);

    return 1;
  }

  return 0;
}


PUBLIC int
vrna_sc_set_bp_comparative(vrna_fold_compound_t *fc,
                           const FLT_OR_DBL     ***constraints,
                           unsigned int         options)
{
  unsigned int  s;
  int           ret;

  ret = 0;

  if ((fc) &&
      (fc->type == VRNA_FC_TYPE_COMPARATIVE) &&
      (constraints)) {
    ret = sc_reset_bp_comparative(fc, constraints, options);

    if (options & VRNA_OPTION_MFE)
      for (s = 0; s < fc->n_seq; s++)
        prepare_sc_bp_mfe(fc->scs[s], fc->length, fc->jindx, options);

    if (options & VRNA_OPTION_PF) /* prepare Boltzmann factors for the BP soft constraints */
      for (s = 0; s < fc->n_seq; s++)
        prepare_sc_bp_pf(fc->scs[s], fc->length, fc->jindx, fc->exp_params->kT, options);
  }

  return ret;
}


PUBLIC int
vrna_sc_set_bp_comparative_seq(vrna_fold_compound_t *fc,
                               unsigned int         s,
                               const FLT_OR_DBL     **constraints,
                               unsigned int         options)
{
  int               ret;
  const FLT_OR_DBL  ***cs;

  ret = 0;

  if ((fc) &&
      (fc->type == VRNA_FC_TYPE_COMPARATIVE) &&
      (s < fc->n_seq) &&
      (constraints)) {
    cs    = (const FLT_OR_DBL ***)vrna_alloc(sizeof(FLT_OR_DBL **) * fc->n_seq);
    cs[s] = constraints;
    ret   = vrna_sc_set_bp_comparative(fc, cs, options);

    free(cs);
  }
  
  return ret;
}


PUBLIC int
vrna_sc_add_bp(vrna_fold_compound_t *fc,
               unsigned int         i,
               unsigned int         j,
               FLT_OR_DBL           energy,
               unsigned int         options)
{
  if (fc && (fc->type == VRNA_FC_TYPE_SINGLE)) {
    if ((i < 1) ||
        (i > fc->length) ||
        (j < i) ||
        (j > fc->length)) {
      vrna_log_warning("vrna_sc_add_bp(): Base pair (%d, %d) out of range!"
                           " (Sequence length: %d)",
                           i, j, fc->length);
    } else {
      sc_add_bp(fc, i, j, energy, options);

      if (options & VRNA_OPTION_MFE)
        prepare_sc_bp_mfe(fc->sc, fc->length, fc->jindx, options);

      if (options & VRNA_OPTION_PF) /* prepare Boltzmann factors for the BP soft constraints */
        prepare_sc_bp_pf(fc->sc, fc->length, fc->jindx, fc->exp_params->kT, options);

      return 1;
    }
  }

  return 0;
}


PUBLIC int
vrna_sc_add_bp_comparative(vrna_fold_compound_t *fc,
                           unsigned int         *is,
                           unsigned int         *js,
                           const FLT_OR_DBL     *energies,
                           unsigned int         options)
{
  unsigned int  s;
  int           ret;

  ret = 0;

  if ((fc) &&
      (fc->type == VRNA_FC_TYPE_COMPARATIVE) &&
      (is) &&
      (js) &&
      (energies)) {
    for (s = 0; s < fc->n_seq; s++) {
      if (is[s] != 0) {
        if ((is[s] > fc->length) ||
            (js[s] < is[s]) ||
            (js[s] > fc->length)) {
          vrna_log_warning("vrna_sc_add_bp_comparative*(): Base pair (%d, %d) out of range for sequence %d!"
                           " (Alignment length: %d)"
                           "Omitting data!",
                           is[s], js[s], s, fc->length);
          is[s] = 0;
        }
      }
    }

    ret = sc_add_bp_comparative(fc, is, js, energies, options);

    if (options & VRNA_OPTION_MFE)
      for (s = 0; s < fc->n_seq; s++)
        prepare_sc_bp_mfe(fc->scs[s], fc->length, fc->jindx, options);

    if (options & VRNA_OPTION_PF) /* prepare Boltzmann factors for the BP soft constraints */
      for (s = 0; s < fc->n_seq; s++)
        prepare_sc_bp_pf(fc->scs[s], fc->length, fc->jindx, fc->exp_params->kT, options);
  }

  return ret;
}


PUBLIC int
vrna_sc_add_bp_comparative_seq(vrna_fold_compound_t *fc,
                               unsigned int         s,
                               unsigned int         i,
                               unsigned int         j,
                               FLT_OR_DBL           energy,
                               unsigned int         options)
{
  unsigned int  *is, *js;
  int           ret;
  FLT_OR_DBL    *es;

  ret = 0;

  if ((fc) &&
      (fc->type == VRNA_FC_TYPE_COMPARATIVE) &&
      (s < fc->n_seq) &&
      (i != 0) &&
      (i < j)) {
    is    = (unsigned int *)vrna_alloc(sizeof(unsigned int) * fc->n_seq);
    js    = (unsigned int *)vrna_alloc(sizeof(unsigned int) * fc->n_seq);
    es    = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * fc->n_seq);
    is[s] = i;
    js[s] = j;
    es[s] = energy;

    ret = vrna_sc_add_bp_comparative(fc, is, js, es, options);

    free(is);
    free(js);
    free(es);
  }

  return ret;
}


PUBLIC int
vrna_sc_set_up(vrna_fold_compound_t *fc,
               const FLT_OR_DBL     *constraints,
               unsigned int         options)
{
  if (fc && (fc->type == VRNA_FC_TYPE_SINGLE)) {
    sc_reset_up(fc, constraints, options);

    if (options & VRNA_OPTION_MFE)
      prepare_sc_up_mfe(fc->sc, fc->length, options);

    if (options & VRNA_OPTION_PF)
      prepare_sc_up_pf(fc->sc, fc->length, fc->exp_params->kT, options);

    return 1;
  }

  return 0;
}


PUBLIC int
vrna_sc_set_up_comparative(vrna_fold_compound_t *fc,
                           const FLT_OR_DBL     **constraints,
                           unsigned int         options)
{
  unsigned int  s;
  int           ret;

  ret = 0;

  if ((fc) &&
      (fc->type == VRNA_FC_TYPE_COMPARATIVE) &&
      (constraints)) {
    ret = sc_reset_up_comparative(fc, constraints, options);

    if (options & VRNA_OPTION_MFE)
      for (s = 0; s < fc->n_seq; s++)
        prepare_sc_up_mfe(fc->scs[s], fc->a2s[s][fc->length], options);

    if (options & VRNA_OPTION_PF)
      for (s = 0; s < fc->n_seq; s++)
        prepare_sc_up_pf(fc->scs[s], fc->a2s[s][fc->length], fc->exp_params->kT, options);

  }

  return ret;
}


PUBLIC int
vrna_sc_set_up_comparative_seq(vrna_fold_compound_t *fc,
                               unsigned int         s,
                               const FLT_OR_DBL     *constraints,
                               unsigned int         options)
{
  const FLT_OR_DBL  **c;
  int               ret;

  ret = 0;

  if ((fc) &&
      (fc->type == VRNA_FC_TYPE_COMPARATIVE) &&
      (s < fc->n_seq)) {
    c     = (const FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * fc->n_seq);
    c[s]  = constraints;
    ret   = vrna_sc_set_up_comparative(fc, c, options);

    free(c);
  }

  return ret;
}



PUBLIC int
vrna_sc_add_up(vrna_fold_compound_t *fc,
               unsigned int         i,
               FLT_OR_DBL           energy,
               unsigned int         options)
{
  if (fc && (fc->type == VRNA_FC_TYPE_SINGLE)) {
    if ((i < 1) || (i > fc->length)) {
      vrna_log_warning("vrna_sc_add_up(): Nucleotide position %d out of range!"
                           " (Sequence length: %d)",
                           i, fc->length);
    } else {
      sc_add_up(fc, (unsigned int)i, energy, options);

      if (options & VRNA_OPTION_MFE)
        prepare_sc_up_mfe(fc->sc, fc->length, options);

      if (options & VRNA_OPTION_PF)
        prepare_sc_up_pf(fc->sc, fc->length, fc->exp_params->kT, options);

      return 1;
    }
  }

  return 0;
}


PUBLIC int
vrna_sc_add_up_comparative(vrna_fold_compound_t *fc,
                           unsigned int         *is,
                           const FLT_OR_DBL     *energies,
                           unsigned int         options)
{
  unsigned int  s;
  int           ret;


  ret = 0;

  if ((fc) &&
      (fc->type == VRNA_FC_TYPE_COMPARATIVE) &&
      (is) &&
      (energies)) {
    for (s = 0; s < fc->n_seq; s++) {
      if ((is[s] != 0) &&
          (is[s] > fc->a2s[s][fc->length])) {
        vrna_log_warning("vrna_sc_add_up_comparative*(): Nucleotide position %d out of range for sequence %d!"
                         " (Sequence length: %d)\n"
                         "Omitting data!",
                         is[s], s, fc->a2s[s][fc->length]);
        is[s] = 0;
      }
    }

    ret = sc_add_up_comparative(fc, is, energies, options);

    if (options & VRNA_OPTION_MFE)
      for (s = 0; s < fc->n_seq; s++)
        prepare_sc_up_mfe(fc->scs[s], fc->a2s[s][fc->length], options);

    if (options & VRNA_OPTION_PF)
      for (s = 0; s < fc->n_seq; s++)
        prepare_sc_up_pf(fc->scs[s], fc->a2s[s][fc->length], fc->exp_params->kT, options);

  }

  return ret;
}


PUBLIC int
vrna_sc_add_up_comparative_seq(vrna_fold_compound_t *fc,
                               unsigned int         s,
                               unsigned int         i,
                               const FLT_OR_DBL     energy,
                               unsigned int         options)
{
  unsigned int  *is;
  int           ret;
  FLT_OR_DBL    *es;

  ret = 0;

  if ((fc) &&
      (fc->type == VRNA_FC_TYPE_COMPARATIVE) &&
      (s < fc->n_seq)) {
    if ((i == 0) ||
        (i > fc->a2s[s][fc->length])) {
      vrna_log_warning("vrna_sc_add_up_comparative*(): Nucleotide position %d out of range for sequence %d!"
                       " (Sequence length: %d)\n"
                       "Omitting data!",
                       i, s, fc->a2s[s][fc->length]);
    } else {
      is    = (unsigned int *)vrna_alloc(sizeof(unsigned int) * fc->n_seq);
      es    = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * fc->n_seq);
      is[s] = i;
      es[s] = energy;

      ret = vrna_sc_add_up_comparative(fc, is, (const FLT_OR_DBL *)es, options);

      free(is);
      free(es);
    }
  }

  return ret;
}


PUBLIC int
vrna_sc_set_stack(vrna_fold_compound_t  *fc,
                  const FLT_OR_DBL      *constraints,
                  unsigned int          options)
{
  unsigned int i;

  if ((fc) &&
      (constraints) &&
      (fc->type == VRNA_FC_TYPE_SINGLE)) {
    if (!fc->sc) {
      if (options & VRNA_OPTION_WINDOW)
        vrna_sc_init_window(fc);
      else
        vrna_sc_init(fc);
    }

    free(fc->sc->energy_stack);
    fc->sc->energy_stack = (int *)vrna_alloc(sizeof(int) * (fc->length + 1));

    for (i = 1; i <= fc->length; ++i)
      fc->sc->energy_stack[i] = (int)roundf(constraints[i] * 100.);

    return 1;
  }

  return 0;
}


PUBLIC int
vrna_sc_set_stack_comparative(vrna_fold_compound_t  *fc,
                              const FLT_OR_DBL      **constraints,
                              unsigned int          options)
{
  unsigned int  i, s;
  int           ret;

  ret = 0;

  if ((fc) &&
      (constraints) &&
      (fc->type == VRNA_FC_TYPE_COMPARATIVE)) {
    if (!fc->scs) {
      if (options & VRNA_OPTION_WINDOW)
        vrna_sc_init_window(fc);
      else
        vrna_sc_init(fc);
    }

    for (s = 0; s < fc->n_seq; s++) {
      free(fc->scs[s]->energy_stack);
      fc->scs[s]->energy_stack = NULL;

      if (constraints[s]) {
        fc->scs[s]->energy_stack = (int *)vrna_alloc(sizeof(int) * (fc->a2s[s][fc->length] + 1));

        for (i = 1; i <= fc->a2s[s][fc->length]; ++i)
          fc->scs[s]->energy_stack[i] = (int)roundf(constraints[s][i] * 100.);

        ret++;
      }
    }
  }

  return ret;
}


PUBLIC int
vrna_sc_set_stack_comparative_seq(vrna_fold_compound_t  *fc,
                                  unsigned int          s,
                                  const FLT_OR_DBL      *constraints,
                                  unsigned int          options)
{
  int               ret;
  const FLT_OR_DBL  **c;

  ret = 0;

  if ((fc) &&
      (fc->type == VRNA_FC_TYPE_COMPARATIVE) &&
      (s < fc->n_seq) &&
      (constraints)) {
    c     = (const FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * fc->n_seq);
    c[s]  = constraints;
    ret   = vrna_sc_set_stack_comparative(fc, c, options);

    free(c);
  }

  return ret;
}


PUBLIC int
vrna_sc_add_stack(vrna_fold_compound_t  *fc,
                  unsigned int          i,
                  FLT_OR_DBL            energy,
                  unsigned int          options)
{
  if ((fc) &&
      (fc->type == VRNA_FC_TYPE_SINGLE)) {
    if ((i < 1) || (i > fc->length)) {
      vrna_log_warning("vrna_sc_add_stack*(): Nucleotide position %d out of range!"
                           " (Sequence length: %d)",
                           i, fc->length);
    } else {
      if (!fc->sc) {
        if (options & VRNA_OPTION_WINDOW)
          vrna_sc_init_window(fc);
        else
          vrna_sc_init(fc);
      }

      if (!fc->sc->energy_stack)
        fc->sc->energy_stack = (int *)vrna_alloc(sizeof(int) * (fc->length + 1));

      fc->sc->energy_stack[i] += (int)roundf(energy * 100.);

      return 1;
    }
  }

  return 0;
}


PUBLIC int
vrna_sc_add_stack_comparative(vrna_fold_compound_t  *fc,
                              unsigned int          *is,
                              const FLT_OR_DBL      *energies,
                              unsigned int          options)
{
  unsigned int  s;
  int           ret;

  ret = 0;

  if ((fc) &&
      (fc->type == VRNA_FC_TYPE_COMPARATIVE) &&
      (energies)) {
    for (s = 0; s < fc->n_seq; s++) {
      if ((is[s] != 0) &&
          (is[s] > fc->a2s[s][fc->length])) {
        vrna_log_warning("vrna_sc_add_stack_comparative*(): Nucleotide position %d out of range for sequence %d!"
                           " (Sequence length: %d)\n"
                           "Omitting data!",
                           is[s], s, fc->a2s[s][fc->length]);
        is[s] = 0;
      }
    }

    if (!fc->scs) {
      if (options & VRNA_OPTION_WINDOW)
        vrna_sc_init_window(fc);
      else
        vrna_sc_init(fc);
    }

    for (s = 0; s < fc->n_seq; s++) {
      if (is[s] != 0) {
        if (fc->scs[s]->energy_stack == NULL)
          fc->scs[s]->energy_stack = (int *)vrna_alloc(sizeof(int) * (fc->a2s[s][fc->length] + 1));

        fc->scs[s]->energy_stack[is[s]] += (int)roundf(energies[s] * 100.);
        ret++;
      }
    }
  }

  return ret;
}


PUBLIC int
vrna_sc_add_stack_comparative_seq(vrna_fold_compound_t  *fc,
                                  unsigned int          s,
                                  unsigned int          i,
                                  FLT_OR_DBL            energy,
                                  unsigned int          options)
{
  unsigned int  *is;
  int           ret;
  FLT_OR_DBL    *es;

  ret = 0;

  if ((fc) &&
      (fc->type == VRNA_FC_TYPE_COMPARATIVE) &&
      (s < fc->n_seq)) {
    is    = (unsigned int *)vrna_alloc(sizeof(unsigned int) * fc->n_seq);
    es    = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * fc->n_seq);
    is[s] = i;
    es[s] = energy;

    ret   = vrna_sc_add_stack_comparative(fc, is, es, options);
    
    free(is);
    free(es);
  }
  
  return ret;
}


PUBLIC int
vrna_sc_add_data(vrna_fold_compound_t *fc,
                 void                 *data,
                 vrna_auxdata_free_f  free_data)
{
  return vrna_sc_add_auxdata(fc, data, NULL, free_data);
}


PUBLIC int
vrna_sc_add_auxdata(vrna_fold_compound_t    *fc,
                    void                    *data,
                    vrna_auxdata_prepare_f  prepare_cb,
                    vrna_auxdata_free_f     free_cb)
{
  if ((fc) &&
      (fc->type == VRNA_FC_TYPE_SINGLE)) {
    if (!fc->sc)
      vrna_sc_init(fc);

    vrna_sc_t *sc = fc->sc;

    if (sc->free_data)
      sc->free_data(sc->data);

    sc->data          = data;
    sc->free_data     = free_cb;
    sc->prepare_data  = prepare_cb;

    return 1;
  }

  return 0;
}


PUBLIC int
vrna_sc_set_data_comparative(vrna_fold_compound_t *fc,
                             void                 **data,
                             vrna_auxdata_free_f  *free_data,
                             unsigned int         options)
{
  return vrna_sc_set_auxdata_comparative(fc, data, NULL, free_data, options);
}


PUBLIC int
vrna_sc_set_data_comparative_seq(vrna_fold_compound_t    *fc,
                                 unsigned int            s,
                                 void                    *data,
                                 vrna_auxdata_free_f     free_data,
                                 unsigned int            options)
{
  return vrna_sc_set_auxdata_comparative_seq(fc, s, data, NULL, free_data, options);
}



PUBLIC int
vrna_sc_set_auxdata_comparative(vrna_fold_compound_t    *fc,
                                void                    **data,
                                vrna_auxdata_prepare_f  *prepare_cbs,
                                vrna_auxdata_free_f     *free_data,
                                unsigned int            options)
{
  unsigned int  s;
  int           ret;

  ret = 0;

  if ((fc) &&
      (fc->type == VRNA_FC_TYPE_COMPARATIVE) &&
      (data)) {
    if (!fc->scs) {
      if (options & VRNA_OPTION_WINDOW)
        vrna_sc_init_window(fc);
      else
        vrna_sc_init(fc);
    }

    for (s = 0; s < fc->n_seq; s++) {
      /* release previous data entries */
      if (fc->scs[s]->free_data)
        fc->scs[s]->free_data(fc->scs[s]->data);

      fc->scs[s]->free_data     = NULL;
      fc->scs[s]->prepare_data  = NULL;
      fc->scs[s]->data          = NULL;

      /* set current data */
      if (data[s]) {
        fc->scs[s]->data = data[s];
        ret++;
      }
    }

    if (prepare_cbs)
      for (s = 0; s < fc->n_seq; s++)
        if (prepare_cbs[s])
          fc->scs[s]->prepare_data = prepare_cbs[s];

    if (free_data)
      for (s = 0; s < fc->n_seq; s++)
        if (free_data[s])
          fc->scs[s]->free_data = free_data[s];
  }

  return ret;
}


PUBLIC int
vrna_sc_set_auxdata_comparative_seq(vrna_fold_compound_t    *fc,
                                    unsigned int            s,
                                    void                    *data,
                                    vrna_auxdata_prepare_f  prepare_cb,
                                    vrna_auxdata_free_f     free_data,
                                    unsigned int            options)
{
  int                     ret;
  void                    **ds;
  vrna_auxdata_prepare_f  *ps;
  vrna_auxdata_free_f     *fs;

  ret = 0;

  if ((fc) &&
      (fc->type == VRNA_FC_TYPE_COMPARATIVE) &&
      (s < fc->n_seq) &&
      (data)) {
    ds    = (void **)vrna_alloc(sizeof(void *) * fc->n_seq);
    ps    = (vrna_auxdata_prepare_f *)vrna_alloc(sizeof(vrna_auxdata_prepare_f) * fc->n_seq);
    fs    = (vrna_auxdata_free_f *)vrna_alloc(sizeof(vrna_auxdata_free_f) * fc->n_seq);
    ds[s] = data;
    ps[s] = prepare_cb;
    fs[s] = free_data;

    ret   = vrna_sc_set_auxdata_comparative(fc, ds, ps, fs, options);
    
    free(ds);
    free(ps);
    free(fs);
  }

  return ret;
}


PUBLIC int
vrna_sc_add_f(vrna_fold_compound_t  *fc,
              vrna_sc_f             f)
{
  if ((fc) &&
      (f) &&
      (fc->type == VRNA_FC_TYPE_SINGLE)) {
    if (!fc->sc)
      vrna_sc_init(fc);

    fc->sc->f = f;
    return 1;
  }

  return 0;
}


PUBLIC int
vrna_sc_set_f_comparative(vrna_fold_compound_t  *fc,
                          vrna_sc_f             *f,
                          unsigned int          options)
{
  unsigned int  s;
  int           ret;


  ret = 0;

  if ((fc) &&
      (f) &&
      (fc->type == VRNA_FC_TYPE_COMPARATIVE)) {
    if (!fc->scs) {
      if (options & VRNA_OPTION_WINDOW)
        vrna_sc_init_window(fc);
      else
        vrna_sc_init(fc);
    }

    for (s = 0; s < fc->n_seq; s++) {
      fc->scs[s]->f = f[s];

      if (f[s] != NULL)
        ret++;
    }
  }

  return ret;
}


PUBLIC int
vrna_sc_add_bt(vrna_fold_compound_t *fc,
               vrna_sc_bt_f         f)
{
  if ((fc) &&
      (f) &&
      (fc->type == VRNA_FC_TYPE_SINGLE)) {
    if (!fc->sc)
      vrna_sc_init(fc);

    fc->sc->bt = f;
    return 1;
  }

  return 0;
}


PUBLIC int
vrna_sc_add_exp_f(vrna_fold_compound_t  *fc,
                  vrna_sc_exp_f         exp_f)
{
  if ((fc) &&
      (exp_f) &&
      (fc->type == VRNA_FC_TYPE_SINGLE)) {
    if (!fc->sc)
      vrna_sc_init(fc);

    fc->sc->exp_f = exp_f;
    return 1;
  }

  return 0;
}


PUBLIC int
vrna_sc_set_exp_f_comparative(vrna_fold_compound_t  *fc,
                              vrna_sc_exp_f         *exp_f,
                              unsigned int          options)
{
  unsigned int  s;
  int           ret;

  ret = 0;

  if ((fc) &&
      (exp_f) &&
      (fc->type == VRNA_FC_TYPE_COMPARATIVE)) {
    if (!fc->scs) {
      if (options & VRNA_OPTION_WINDOW)
        vrna_sc_init_window(fc);
      else
        vrna_sc_init(fc);
    }

    for (s = 0; s < fc->n_seq; s++) {
      fc->scs[s]->exp_f = exp_f[s];

      if (exp_f[s] != NULL)
        ret++;
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
sc_init_up_storage(vrna_sc_t *sc)
{
  if (!sc->up_storage)
    sc->up_storage = (int *)vrna_alloc(sizeof(int) * (sc->n + 2));
}


/* pupulate sc->energy_up array at position i from sc->up_storage data */
PRIVATE INLINE void
populate_sc_up_mfe(vrna_sc_t    *sc,
                   unsigned int i,
                   unsigned int n)
{
  unsigned int  j;

  sc->energy_up[i][0] = 0;
  for (j = 1; j <= n; j++)
    sc->energy_up[i][j] = sc->energy_up[i][j - 1] +
                          sc->up_storage[i + j - 1];
}


PRIVATE INLINE void
populate_sc_up_pf(vrna_sc_t     *sc,
                  unsigned int  i,
                  unsigned int  n,
                  double        kT)
{
  unsigned int  j;
  double        GT;

  sc->exp_energy_up[i][0] = 1.;

  for (j = 1; j <= n; j++) {
    GT                      = (double)sc->up_storage[i + j - 1] * 10.; /* convert deka-cal/mol to cal/mol */
    sc->exp_energy_up[i][j] = sc->exp_energy_up[i][j - 1]
                              * (FLT_OR_DBL)exp(-GT / kT);
  }
}


PRIVATE INLINE void
sc_init_bp_storage(vrna_sc_t *sc)
{
  unsigned int i;

  if (sc->bp_storage == NULL) {
    sc->bp_storage = (vrna_sc_bp_storage_t **)vrna_alloc(
      sizeof(vrna_sc_bp_storage_t *) * (sc->n + 2));

    for (i = 1; i <= sc->n; i++)
      /* by default we do not change energy contributions for any base pairs */
      sc->bp_storage[i] = NULL;
  }
}


PRIVATE INLINE void
sc_store_bp(vrna_sc_bp_storage_t  **container,
            unsigned int          i,
            unsigned int          start,
            unsigned int          end,
            int                   e)
{
  unsigned int size, cnt = 0;

  if (!container[i]) {
    container[i] = (vrna_sc_bp_storage_t *)vrna_alloc(sizeof(vrna_sc_bp_storage_t) * 2);
  } else {
    /* find out total size of container */
    for (size = 0; container[i][size].interval_start != 0; size++);

    /* find position where we want to insert the new constraint */
    for (cnt = 0; cnt < size; cnt++) {
      if (container[i][cnt].interval_start > start)
        break; /* want to insert before current constraint */

      if (container[i][cnt].interval_end < end)
        continue; /* want to insert after current constraint */
    }
    /* increase memory for bp constraints */
    container[i] = (vrna_sc_bp_storage_t *)vrna_realloc(container[i],
                                                        sizeof(vrna_sc_bp_storage_t) * (size + 2));
    /* shift trailing constraints by 1 entry */
    memmove(container[i] + cnt + 1, container[i] + cnt,
            sizeof(vrna_sc_bp_storage_t) * (size - cnt + 1));
  }

  container[i][cnt].interval_start  = start;
  container[i][cnt].interval_end    = end;
  container[i][cnt].e               = e;
}


PRIVATE INLINE int
get_stored_bp_contributions(vrna_sc_bp_storage_t  *container,
                            unsigned int          j)
{
  unsigned int  cnt;
  int           e;

  e = 0;

  /* go through list of constraints for current position i */
  for (cnt = 0; container[cnt].interval_start != 0; cnt++) {
    if (container[cnt].interval_start > j)
      break; /* only constraints for pairs (i,q) with q > j left */

    if (container[cnt].interval_end < j)
      continue; /* constraint for pairs (i,q) with q < j */

    /* constraint has interval [p,q] with p <= j <= q */
    e += container[cnt].e;
  }

  return e;
}


PRIVATE INLINE void
populate_sc_bp_mfe(vrna_sc_t    *sc,
                   unsigned int i,
                   unsigned int maxdist,
                   unsigned int n,
                   int          *idx)
{
  unsigned int  j, k;
  int           e;

  if (sc->bp_storage[i]) {
    for (k = 1; k < maxdist; k++) {
      j = i + k;

      if (j > n)
        break;

      e = get_stored_bp_contributions(sc->bp_storage[i], j);

      switch (sc->type) {
        case VRNA_SC_DEFAULT:
          sc->energy_bp[idx[j] + i] = e;
          break;

        case VRNA_SC_WINDOW:
          sc->energy_bp_local[i][j - i] = e;
          break;
      }
    }
  } else {
    for (k = 1; k < maxdist; k++) {
      j = i + k;
      if (j > n)
        break;

      switch (sc->type) {
        case VRNA_SC_DEFAULT:
          sc->energy_bp[idx[j] + i] = 0;
          break;

        case VRNA_SC_WINDOW:
          sc->energy_bp_local[i][j - i] = 0;
          break;
      }
    }
  }
}


PRIVATE INLINE void
populate_sc_bp_pf(vrna_sc_t     *sc,
                  unsigned int  i,
                  unsigned int  maxdist,
                  unsigned int  n,
                  int           *idx,
                  double        kT)
{
  unsigned int      j, k;
  int               e;
  FLT_OR_DBL        q;
  double            GT;

  if (sc->bp_storage[i]) {
    for (k = 1; k < maxdist; k++) {
      j = i + k;

      if (j > n)
        break;

      e = get_stored_bp_contributions(sc->bp_storage[i], j);

      GT  = e * 10.;
      q   = (FLT_OR_DBL)exp(-GT / kT);

      switch (sc->type) {
        case VRNA_SC_DEFAULT:
          sc->exp_energy_bp[idx[j] + i] = q;
          break;

        case VRNA_SC_WINDOW:
          sc->exp_energy_bp_local[i][j - i] = q;
          break;
      }
    }
  } else {
    for (k = 1; k < maxdist; k++) {
      j = i + k;
      if (j > n)
        break;

      switch (sc->type) {
        case VRNA_SC_DEFAULT:
          sc->exp_energy_bp[idx[j] + i] = 1.;
          break;

        case VRNA_SC_WINDOW:
          sc->exp_energy_bp_local[i][j - i] = 1.;
          break;
      }
    }
  }
}


PRIVATE void
sc_add_bp(vrna_fold_compound_t  *fc,
          unsigned int          i,
          unsigned int          j,
          FLT_OR_DBL            energy,
          unsigned int          options)
{
  vrna_sc_t *sc;

  if ((options & VRNA_OPTION_WINDOW) && (!fc->sc))
    vrna_sc_init_window(fc);
  else if (!fc->sc)
    vrna_sc_init(fc);

  sc = fc->sc;
  sc_init_bp_storage(sc);
  sc_store_bp(sc->bp_storage, i, j, j, (int)roundf(energy * 100.));
  sc->state |= STATE_DIRTY_BP_MFE | STATE_DIRTY_BP_PF;
}


PRIVATE int
sc_add_bp_comparative(vrna_fold_compound_t  *fc,
                      unsigned int          *is,
                      unsigned int          *js,
                      const FLT_OR_DBL      *energies,
                      unsigned int          options)
{
  unsigned int  s;
  int           ret;
  vrna_sc_t     *sc;

  ret = 0;
  if ((options & VRNA_OPTION_WINDOW) && (!fc->sc))
    vrna_sc_init_window(fc);
  else if (!fc->sc)
    vrna_sc_init(fc);

  for (s = 0; s < fc->n_seq; s++) {
    if (is[s] != 0) {
      sc = fc->scs[s];
      sc_init_bp_storage(sc);
      sc_store_bp(sc->bp_storage, is[s], js[s], js[s], (int)roundf(energies[s] * 100.));
      sc->state |= STATE_DIRTY_BP_MFE | STATE_DIRTY_BP_PF;

      ret++;
    }
  }

  return ret;
}


PRIVATE INLINE void
free_sc_up(vrna_sc_t *sc)
{
  unsigned int i;

  if (sc) {
    free(sc->up_storage);

    sc->up_storage = NULL;

    if (sc->type == VRNA_SC_DEFAULT) {
      if (sc->energy_up)
        for (i = 0; i <= sc->n + 1; i++)
          free(sc->energy_up[i]);

      if (sc->exp_energy_up)
        for (i = 0; i <= sc->n + 1; i++)
          free(sc->exp_energy_up[i]);
    }

    free(sc->energy_up);
    sc->energy_up = NULL;

    free(sc->exp_energy_up);
    sc->exp_energy_up = NULL;

    sc->state &= ~(STATE_DIRTY_UP_MFE | STATE_DIRTY_UP_PF);
  }
}


PRIVATE INLINE void
free_sc_bp(vrna_sc_t *sc)
{
  unsigned int i;

  if (sc) {
    if (sc->bp_storage) {
      for (i = 1; i <= sc->n; i++)
        free(sc->bp_storage[i]);
      free(sc->bp_storage);
      sc->bp_storage = NULL;
    }

    switch (sc->type) {
      case VRNA_SC_DEFAULT:
        free(sc->energy_bp);
        sc->energy_bp = NULL;

        free(sc->exp_energy_bp);
        sc->energy_bp = NULL;

        break;

      case VRNA_SC_WINDOW:
        free(sc->energy_bp_local);
        sc->energy_bp_local = NULL;

        free(sc->exp_energy_bp_local);
        sc->exp_energy_bp_local = NULL;

        break;
    }
    sc->state &= ~(STATE_DIRTY_BP_MFE | STATE_DIRTY_BP_PF);
  }
}


PRIVATE void
sc_reset_up(vrna_fold_compound_t  *fc,
            const FLT_OR_DBL      *constraints,
            unsigned int          options)
{
  unsigned int  i, n;
  vrna_sc_t     *sc;

  n = fc->length;

  if (!fc->sc) {
    if (options & VRNA_OPTION_WINDOW)
      vrna_sc_init_window(fc);
    else
      vrna_sc_init(fc);
  }

  sc = fc->sc;

  if (constraints) {
    free_sc_up(sc);

    /* initialize container for unpaired probabilities */
    sc_init_up_storage(sc);

    /* add contributions to storage container */
    for (i = 1; i <= n; i++)
      sc->up_storage[i] = (int)roundf(constraints[i] * 100.); /* convert to 10kal/mol */

    sc->state |= STATE_DIRTY_UP_MFE | STATE_DIRTY_UP_PF;
  } else {
    free_sc_up(sc);
  }
}


PRIVATE int
sc_reset_up_comparative(vrna_fold_compound_t  *fc,
                        const FLT_OR_DBL      **constraints,
                        unsigned int          options)
{
  unsigned int  i, s, n;
  int           ret;
  vrna_sc_t     *sc;

  ret = 0;

  if (!fc->scs) {
    if (options & VRNA_OPTION_WINDOW)
      vrna_sc_init_window(fc);
    else
      vrna_sc_init(fc);
  }

  for (s = 0; s < fc->n_seq; s++) {
    sc  = fc->scs[s];

    if (constraints[s]) {
      n   = fc->a2s[s][fc->length];

      free_sc_up(sc);

      /* initialize container for unpaired probabilities */
      sc_init_up_storage(sc);

      /* add contributions to storage container */
      for (i = 1; i <= n; i++)
        sc->up_storage[i] = (int)roundf(constraints[s][i] * 100.); /* convert to 10kal/mol */

      sc->state |= STATE_DIRTY_UP_MFE | STATE_DIRTY_UP_PF;

      ret++;
    } else {
      free_sc_up(sc);
    }
  }

  return ret;
}


PRIVATE void
sc_reset_bp(vrna_fold_compound_t  *fc,
            const FLT_OR_DBL      **constraints,
            unsigned int          options)
{
  unsigned int  i, j, n;
  vrna_sc_t     *sc;

  n = fc->length;

  if (!fc->sc) {
    if (options & VRNA_OPTION_WINDOW)
      vrna_sc_init_window(fc);
    else
      vrna_sc_init(fc);
  }

  sc = fc->sc;

  if (constraints) {
    free_sc_bp(sc);

    /* initialize container for base pair constraints */
    sc_init_bp_storage(sc);

    /* add contributions to storage container */
    for (i = 1; i < n; i++)
      for (j = i + 1; j <= n; j++)
        sc_store_bp(sc->bp_storage, i, j, j, (int)roundf(constraints[i][j] * 100.));

    sc->state |= STATE_DIRTY_BP_MFE | STATE_DIRTY_BP_PF;
  } else {
    free_sc_bp(sc);
  }
}


PRIVATE int
sc_reset_bp_comparative(vrna_fold_compound_t  *fc,
                        const FLT_OR_DBL      ***constraints,
                        unsigned int          options)
{
  unsigned int  i, j, s, n;
  int           ret;
  vrna_sc_t     *sc;

  ret = 0;

  if (!fc->scs) {
    if (options & VRNA_OPTION_WINDOW)
      vrna_sc_init_window(fc);
    else
      vrna_sc_init(fc);
  }

  n = fc->length;

  for (s = 0; s < fc->n_seq; s++) {
    sc  = fc->scs[s];

    if (constraints[s]) {
      free_sc_bp(sc);

      /* initialize container for base pair constraints */
      sc_init_bp_storage(sc);

      /* add contributions to storage container */
      for (i = 1; i < n; i++)
        for (j = i + 1; j <= n; j++)
          sc_store_bp(sc->bp_storage, i, j, j, (int)roundf(constraints[s][i][j] * 100.));

      sc->state |= STATE_DIRTY_BP_MFE | STATE_DIRTY_BP_PF;

      ret++;
    } else {
      free_sc_bp(sc);
    }
  }

  return ret;
}


PRIVATE void
sc_add_up(vrna_fold_compound_t  *fc,
          unsigned int          i,
          FLT_OR_DBL            energy,
          unsigned int          options)
{
  vrna_sc_t *sc;

  if ((options & VRNA_OPTION_WINDOW) && (!fc->sc))
    vrna_sc_init_window(fc);
  else if (!fc->sc)
    vrna_sc_init(fc);

  sc = fc->sc;
  sc_init_up_storage(sc);
  sc->up_storage[i] += (int)roundf(energy * 100.);
  sc->state         |= STATE_DIRTY_UP_MFE | STATE_DIRTY_UP_PF;
}


PRIVATE int
sc_add_up_comparative(vrna_fold_compound_t  *fc,
                      unsigned int          *is,
                      const FLT_OR_DBL      *energies,
                      unsigned int          options)
{
  unsigned int  s;
  int           ret;
  vrna_sc_t     *sc;

  ret = 0;

  if (!fc->scs) {
    if (options & VRNA_OPTION_WINDOW)
      vrna_sc_init_window(fc);
    else
      vrna_sc_init(fc);
  }

  for (s = 0; s < fc->n_seq; s++) {
    if (is[s] != 0) {
      sc = fc->scs[s];
      sc_init_up_storage(sc);
      sc->up_storage[is[s]] += (int)roundf(energies[s] * 100.);
      sc->state             |= STATE_DIRTY_UP_MFE | STATE_DIRTY_UP_PF;

      ret++;
    }
  }

  return ret;
}


/* populate sc->energy_up arrays for usage in MFE computations */
PRIVATE void
prepare_sc_up_mfe(vrna_sc_t     *sc,
                  unsigned int  n,
                  unsigned int  options)
{
  unsigned int i;

  /* prepare sc for unpaired nucleotides only if we actually have some to apply */
  if (sc) {
    if (sc->up_storage) {
      if (sc->state & STATE_DIRTY_UP_MFE) {
        /*  allocate memory such that we can access the soft constraint
         *  energies of a subsequence of length j starting at position i
         *  via sc->energy_up[i][j]
         */
        sc->energy_up = (int **)vrna_realloc(sc->energy_up, sizeof(int *) * (n + 2));

        if (options & VRNA_OPTION_WINDOW) {
          /*
           *  simply init with NULL pointers, since the sliding-window implementation must take
           *  care of allocating the required memory and filling in appropriate energy contributions
           */
          for (i = 0; i <= n + 1; i++)
            sc->energy_up[i] = NULL;
        } else {
          for (i = 1; i <= n; i++)
            sc->energy_up[i] = (int *)vrna_realloc(sc->energy_up[i], sizeof(int) * (n - i + 2));


          sc->energy_up[0]      = (int *)vrna_realloc(sc->energy_up[0], sizeof(int));
          sc->energy_up[n + 1]  = (int *)vrna_realloc(sc->energy_up[n + 1], sizeof(int));

          /* now add soft constraints as stored in container for unpaired sc */
          for (i = 1; i <= n; i++)
            populate_sc_up_mfe(sc, i, n - i + 1);

          sc->energy_up[0][0]     = 0;
          sc->energy_up[n + 1][0] = 0;
        }

        sc->state &= ~STATE_DIRTY_UP_MFE;
      }
    } else if (sc->energy_up) {
      /* remove any unpaired sc if storage container is empty */
      free_sc_up(sc);
    }
  }
}


/* populate sc->exp_energy_up arrays for usage in partition function computations */
PRIVATE void
prepare_sc_up_pf(vrna_sc_t    *sc,
                 unsigned int n,
                 double       kT,
                 unsigned int options)
{
  unsigned int  i;

  if (sc) {
    /* prepare sc for unpaired nucleotides only if we actually have some to apply */
    if (sc->up_storage) {
      if (sc->state & STATE_DIRTY_UP_PF) {
        /*  allocate memory such that we can access the soft constraint
         *  energies of a subsequence of length j starting at position i
         *  via sc->exp_energy_up[i][j]
         */
        sc->exp_energy_up = (FLT_OR_DBL **)vrna_realloc(sc->exp_energy_up,
                                                        sizeof(FLT_OR_DBL *) * (n + 2));

        if (options & VRNA_OPTION_WINDOW) {
          /*
           *  simply init with NULL pointers, since the sliding-window implementation must take
           *  care of allocating the required memory and filling in appropriate Boltzmann factors
           */
          for (i = 0; i <= n + 1; i++)
            sc->exp_energy_up[i] = NULL;
        } else {
          for (i = 1; i <= n; i++)
            sc->exp_energy_up[i] =
              (FLT_OR_DBL *)vrna_realloc(sc->exp_energy_up[i],
                                         sizeof(FLT_OR_DBL) * (n - i + 2));

          sc->exp_energy_up[0] =
            (FLT_OR_DBL *)vrna_realloc(sc->exp_energy_up[0], sizeof(FLT_OR_DBL));
          sc->exp_energy_up[n + 1] = (FLT_OR_DBL *)vrna_realloc(sc->exp_energy_up[n + 1],
                                                                sizeof(FLT_OR_DBL));

          for (i = 1; i <= n; i++)
            populate_sc_up_pf(sc, i, n - i + 1, kT);

          sc->exp_energy_up[0][0]     = 1.;
          sc->exp_energy_up[n + 1][0] = 1.;
        }

        sc->state &= ~STATE_DIRTY_UP_PF;
      }
    }
  }
}


/* populate sc->energy_bp arrays for usage in MFE computations */
PRIVATE void
prepare_sc_bp_mfe(vrna_sc_t     *sc,
                  unsigned int  n,
                  int           *idx,
                  unsigned int  options)
{
  unsigned int  i;

  if (sc) {
    /* prepare sc for base paired positions only if we actually have some to apply */
    if (sc->bp_storage) {
      if (sc->state & STATE_DIRTY_BP_MFE) {
        if (options & VRNA_OPTION_WINDOW) {
          sc->energy_bp_local =
            (int **)vrna_realloc(sc->energy_bp_local, sizeof(int *) * (n + 2));
        } else {
          sc->energy_bp =
            (int *)vrna_realloc(sc->energy_bp, sizeof(int) * (((n + 1) * (n + 2)) / 2));

          for (i = 1; i < n; i++)
            populate_sc_bp_mfe(sc, i, n, n, idx);
        }

        sc->state &= ~STATE_DIRTY_BP_MFE;
      }
    } else {
      free_sc_bp(sc);
    }
  }
}


/* populate sc->exp_energy_bp arrays for usage in partition function computations */
PRIVATE void
prepare_sc_bp_pf(vrna_sc_t    *sc,
                 unsigned int n,
                 int          *idx,
                 double       kT,
                 unsigned int options)
{
  unsigned int  i;

  if (sc) {
    /* prepare sc for base paired positions only if we actually have some to apply */
    if (sc->bp_storage) {
      if (sc->state & STATE_DIRTY_BP_PF) {
        if (options & VRNA_OPTION_WINDOW) {
          sc->exp_energy_bp_local =
            (FLT_OR_DBL **)vrna_realloc(sc->exp_energy_bp_local,
                                        sizeof(FLT_OR_DBL *) * (n + 2));
        } else {
          sc->exp_energy_bp =
            (FLT_OR_DBL *)vrna_realloc(sc->exp_energy_bp,
                                       sizeof(FLT_OR_DBL) * (((n + 1) * (n + 2)) / 2));

          for (i = 1; i < n; i++)
            populate_sc_bp_pf(sc, i, n, n, idx, kT);
        }

        sc->state &= ~STATE_DIRTY_BP_PF;
      }
    }
  }
}


PRIVATE void
prepare_sc_stack_pf(vrna_sc_t     *sc,
                    unsigned int  n,
                    double        kT)
{
  unsigned int           i;

  if (sc) {
    if(sc->energy_stack) {
      if (!sc->exp_energy_stack) {
        sc->exp_energy_stack = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
        for (i = 0; i <= n; ++i)
          sc->exp_energy_stack[i] = 1.;
      }

      for (i = 1; i <= n; ++i)
        sc->exp_energy_stack[i] = (FLT_OR_DBL)exp(-(sc->energy_stack[i] * 10.) / kT);
    }  
  }
}


PRIVATE int
prepare_sc_user_cb(vrna_fold_compound_t *fc,
                   vrna_sc_t            *sc,
                   unsigned int         options)
{
  int ret = 0;

  if ((sc) &&
      (sc->data) &&
      (sc->prepare_data) &&
      (sc->type == VRNA_SC_DEFAULT))
    ret = sc->prepare_data(fc, sc->data, options, NULL);

  return ret;
}


PRIVATE vrna_sc_t *
init_sc_default(unsigned int n)
{
  vrna_sc_t *sc, init = {
    .type = VRNA_SC_DEFAULT
  };

  sc = vrna_alloc(sizeof(vrna_sc_t));

  if (sc) {
    memcpy(sc, &init, sizeof(vrna_sc_t));
    nullify(sc);
    sc->n = n;
  }

  return sc;
}


PRIVATE vrna_sc_t *
init_sc_window(unsigned int n)
{
  vrna_sc_t *sc, init = {
    .type = VRNA_SC_WINDOW
  };

  sc = vrna_alloc(sizeof(vrna_sc_t));

  if (sc) {
    memcpy(sc, &init, sizeof(vrna_sc_t));
    nullify(sc);
    sc->n = n;
  }

  return sc;
}


PRIVATE void
nullify(vrna_sc_t *sc)
{
  if (sc) {
    sc->state             = STATE_CLEAN;
    sc->up_storage        = NULL;
    sc->bp_storage        = NULL;
    sc->energy_up         = NULL;
    sc->energy_stack      = NULL;
    sc->exp_energy_stack  = NULL;
    sc->exp_energy_up     = NULL;

    sc->f             = NULL;
    sc->exp_f         = NULL;
    sc->data          = NULL;
    sc->prepare_data  = NULL;
    sc->free_data     = NULL;

    switch (sc->type) {
      case VRNA_SC_DEFAULT:
        sc->energy_bp     = NULL;
        sc->exp_energy_bp = NULL;
        break;

      case VRNA_SC_WINDOW:
        sc->energy_bp_local     = NULL;
        sc->exp_energy_bp_local = NULL;
        break;
    }
  }
}
