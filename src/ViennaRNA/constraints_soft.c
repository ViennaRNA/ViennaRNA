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

#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/energy_const.h" /* defines MINPSCORE */
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/aln_util.h"
#include "ViennaRNA/file_formats.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/constraints_soft.h"


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
sc_reset_up(vrna_fold_compound_t  *vc,
            const FLT_OR_DBL      *constraints);


PRIVATE void
sc_reset_bp(vrna_fold_compound_t  *vc,
            const FLT_OR_DBL      **constraints);


PRIVATE void
sc_add_up(vrna_fold_compound_t  *vc,
          int                   i,
          FLT_OR_DBL            energy,
          unsigned int          options);


PRIVATE void
sc_add_bp(vrna_fold_compound_t  *vc,
          int                   i,
          int                   j,
          FLT_OR_DBL            energy,
          unsigned int          options);


PRIVATE void
prepare_sc_up_mfe(vrna_fold_compound_t  *vc,
                  unsigned int          options);


PRIVATE void
prepare_sc_bp_mfe(vrna_fold_compound_t  *vc,
                  unsigned int          options);


PRIVATE void
prepare_sc_up_pf(vrna_fold_compound_t *vc,
                 unsigned int         options);


PRIVATE void
prepare_sc_bp_pf(vrna_fold_compound_t *vc,
                 unsigned int         options);


PRIVATE INLINE void
sc_init_up_storage(vrna_sc_t *sc);


PRIVATE INLINE void
sc_init_bp_storage(vrna_sc_t *sc);


PRIVATE INLINE void
sc_store_bp(vrna_sc_bp_storage_t  **container,
            int                   i,
            int                   start,
            int                   end,
            int                   e);


PRIVATE INLINE void
free_sc_up(vrna_sc_t *sc);


PRIVATE INLINE void
free_sc_bp(vrna_sc_t *sc);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC void
vrna_sc_init(vrna_fold_compound_t *vc)
{
  unsigned int  s;
  vrna_sc_t     *sc;

  if (vc) {
    vrna_sc_remove(vc);

    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        sc                    = (vrna_sc_t *)vrna_alloc(sizeof(vrna_sc_t));
        sc->type              = VRNA_SC_DEFAULT;
        sc->n                 = vc->length;
        sc->state             = STATE_CLEAN;
        sc->up_storage        = NULL;
        sc->bp_storage        = NULL;
        sc->energy_up         = NULL;
        sc->energy_bp         = NULL;
        sc->energy_stack      = NULL;
        sc->exp_energy_stack  = NULL;
        sc->exp_energy_up     = NULL;
        sc->exp_energy_bp     = NULL;
        sc->f                 = NULL;
        sc->exp_f             = NULL;
        sc->data              = NULL;
        sc->free_data         = NULL;

        vc->sc = sc;
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        vc->scs = (vrna_sc_t **)vrna_alloc(sizeof(vrna_sc_t *) * (vc->n_seq + 1));
        for (s = 0; s < vc->n_seq; s++) {
          sc                    = (vrna_sc_t *)vrna_alloc(sizeof(vrna_sc_t));
          sc->type              = VRNA_SC_DEFAULT;
          sc->n                 = vc->length;
          sc->state             = STATE_CLEAN;
          sc->up_storage        = NULL;
          sc->bp_storage        = NULL;
          sc->energy_up         = NULL;
          sc->energy_bp         = NULL;
          sc->energy_stack      = NULL;
          sc->exp_energy_stack  = NULL;
          sc->exp_energy_up     = NULL;
          sc->exp_energy_bp     = NULL;
          sc->f                 = NULL;
          sc->exp_f             = NULL;
          sc->data              = NULL;
          sc->free_data         = NULL;

          vc->scs[s] = sc;
        }
        break;
      default:                      /* do nothing */
        break;
    }
  }
}


PUBLIC void
vrna_sc_init_window(vrna_fold_compound_t *vc)
{
  unsigned int  s;
  vrna_sc_t     *sc;

  if (vc) {
    vrna_sc_remove(vc);

    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        sc                      = (vrna_sc_t *)vrna_alloc(sizeof(vrna_sc_t));
        sc->type                = VRNA_SC_WINDOW;
        sc->n                   = vc->length;
        sc->state               = STATE_CLEAN;
        sc->up_storage          = NULL;
        sc->bp_storage          = NULL;
        sc->energy_up           = NULL;
        sc->energy_bp_local     = NULL;
        sc->energy_stack        = NULL;
        sc->exp_energy_stack    = NULL;
        sc->exp_energy_up       = NULL;
        sc->exp_energy_bp_local = NULL;
        sc->f                   = NULL;
        sc->exp_f               = NULL;
        sc->data                = NULL;
        sc->free_data           = NULL;

        vc->sc = sc;
        break;

      default:                      /* do nothing */
        break;
    }
  }
}


PUBLIC void
vrna_sc_prepare(vrna_fold_compound_t  *vc,
                unsigned int          options)
{
  if (vc) {
    if (options & VRNA_OPTION_MFE) {
      prepare_sc_up_mfe(vc, options);
      prepare_sc_bp_mfe(vc, options);
    }

    if (options & VRNA_OPTION_PF) {
      prepare_sc_up_pf(vc, options);
      prepare_sc_bp_pf(vc, options);
      if ((vc->sc) && (!vc->sc->exp_energy_stack))
        vrna_sc_add_SHAPE_deigan(vc, NULL, 0, 0, VRNA_OPTION_PF);
    }
  }
}


PUBLIC void
vrna_sc_update(vrna_fold_compound_t *vc,
               int                  i,
               unsigned int         options)
{
  int       j, k, n, maxdist, turn, cnt, e;
  vrna_sc_t *sc;

  if (vc) {
    sc = vc->sc;

    if (options & VRNA_OPTION_WINDOW) {
      /* sliding-window mode, i.e. local structure prediction */
      if (sc && (i > 0)) {
        if (sc->up_storage) {
          maxdist = MIN2(vc->window_size, vc->length - i + 1);

          /* additional energy contributions per unpaired nucleotide */
          sc->energy_up[i][0] = 0;
          for (k = 1; k <= maxdist; k++)
            sc->energy_up[i][k] = sc->energy_up[i][k - 1] +
                                  sc->up_storage[i + k - 1];
        }

        if (sc->bp_storage) {
          /* additional energy contributions per base pair */
          n       = (int)vc->length;
          maxdist = vc->window_size;
          turn    = vc->params->model_details.min_loop_size;

          switch (vc->type) {
            case VRNA_FC_TYPE_SINGLE:
              for (k = turn + 1; k < maxdist; k++) {
                j = i + k;
                e = 0;
                if (j > n)
                  break;

                /* check whether we have constraints on any pairing partner i or j */
                if (sc->bp_storage[i]) {
                  /* go through list of constraints on position i */
                  for (cnt = 0; sc->bp_storage[i][cnt].interval_start != 0; cnt++) {
                    if (sc->bp_storage[i][cnt].interval_start > j)
                      break; /* only constraints for pairs (i,q) with q > j left */

                    if (sc->bp_storage[i][cnt].interval_end < j)
                      continue; /* constraint for pairs (i,q) with q < j */

                    /* constraint has interval [p,q] with p <= j <= q */
                    e += sc->bp_storage[i][cnt].e;
                  }
                }

                sc->energy_bp_local[i][j - i] = e;
              }

              break;
          }
        }
      }
    } else {
      /* do nothing here, until we know what a reasonable action for global folding would be */
    }
  }
}


PUBLIC void
vrna_sc_remove(vrna_fold_compound_t *vc)
{
  int s;

  if (vc) {
    switch (vc->type) {
      case  VRNA_FC_TYPE_SINGLE:
        vrna_sc_free(vc->sc);
        vc->sc = NULL;
        break;

      case  VRNA_FC_TYPE_COMPARATIVE:
        if (vc->scs) {
          for (s = 0; s < vc->n_seq; s++)
            vrna_sc_free(vc->scs[s]);
          free(vc->scs);
        }

        vc->scs = NULL;
        break;
      default:  /* do nothing */
        break;
    }
  }
}


PUBLIC void
vrna_sc_free(vrna_sc_t *sc)
{
  unsigned int i;

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


PUBLIC void
vrna_sc_set_bp(vrna_fold_compound_t *vc,
               const FLT_OR_DBL     **constraints,
               unsigned int         options)
{
  if (vc && (vc->type == VRNA_FC_TYPE_SINGLE)) {
    sc_reset_bp(vc, constraints);

    if (options & VRNA_OPTION_MFE)
      prepare_sc_bp_mfe(vc, options);

    if (options & VRNA_OPTION_PF) /* prepare Boltzmann factors for the BP soft constraints */
      prepare_sc_bp_pf(vc, options);
  }
}


PUBLIC void
vrna_sc_add_bp(vrna_fold_compound_t *vc,
               int                  i,
               int                  j,
               FLT_OR_DBL           energy,
               unsigned int         options)
{
  if (vc && (vc->type == VRNA_FC_TYPE_SINGLE)) {
    if ((i < 1) || (i > vc->length) || (j < i) || (j > vc->length)) {
      vrna_message_warning("vrna_sc_add_bp(): Base pair (%d, %d) out of range!"
                           " (Sequence length: %d)",
                           i, j, vc->length);
    } else {
      sc_add_bp(vc, i, j, energy, options);

      if (options & VRNA_OPTION_MFE)
        prepare_sc_bp_mfe(vc, options);

      if (options & VRNA_OPTION_PF) /* prepare Boltzmann factors for the BP soft constraints */
        prepare_sc_bp_pf(vc, options);
    }
  }
}


PUBLIC void
vrna_sc_set_up(vrna_fold_compound_t *vc,
               const FLT_OR_DBL     *constraints,
               unsigned int         options)
{
  if (vc && (vc->type == VRNA_FC_TYPE_SINGLE)) {
    sc_reset_up(vc, constraints);

    if (options & VRNA_OPTION_MFE)
      prepare_sc_up_mfe(vc, options);

    if (options & VRNA_OPTION_PF)
      prepare_sc_up_pf(vc, options);
  }
}


PUBLIC void
vrna_sc_add_up(vrna_fold_compound_t *vc,
               int                  i,
               FLT_OR_DBL           energy,
               unsigned int         options)
{
  if (vc && (vc->type == VRNA_FC_TYPE_SINGLE)) {
    if ((i < 1) || (i > vc->length)) {
      vrna_message_warning("vrna_sc_add_up(): Nucleotide position %d out of range!"
                           " (Sequence length: %d)",
                           i, vc->length);
    } else {
      sc_add_up(vc, i, energy, options);

      if (options & VRNA_OPTION_MFE)
        prepare_sc_up_mfe(vc, options);

      if (options & VRNA_OPTION_PF)
        prepare_sc_up_pf(vc, options);
    }
  }
}


PUBLIC void
vrna_sc_set_stack(vrna_fold_compound_t  *vc,
                  const FLT_OR_DBL      *constraints,
                  unsigned int          options)
{
  int i;

  if ((vc) && (constraints)) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (options & VRNA_OPTION_WINDOW) {
          if (!vc->sc)
            vrna_sc_init_window(vc);
        } else {
          if (!vc->sc)
            vrna_sc_init(vc);
        }

        free(vc->sc->energy_stack);
        vc->sc->energy_stack = (int *)vrna_alloc(sizeof(int) * (vc->length + 1));

        for (i = 1; i <= vc->length; ++i)
          vc->sc->energy_stack[i] = (int)roundf(constraints[i] * 100.);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        break;
    }
  }
}


PUBLIC void
vrna_sc_add_stack(vrna_fold_compound_t  *vc,
                  int                   i,
                  FLT_OR_DBL            energy,
                  unsigned int          options)
{
  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if ((i < 1) || (i > vc->length)) {
          vrna_message_warning("vrna_sc_add_stack(): Nucleotide position %d out of range!"
                               " (Sequence length: %d)",
                               i, vc->length);
        } else {
          if (options & VRNA_OPTION_WINDOW) {
            if (!vc->sc)
              vrna_sc_init_window(vc);
          } else {
            if (!vc->sc)
              vrna_sc_init(vc);
          }

          if (!vc->sc->energy_stack)
            vc->sc->energy_stack = (int *)vrna_alloc(sizeof(int) * (vc->length + 1));

          vc->sc->energy_stack[i] += (int)roundf(energy * 100.);
        }

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        break;
    }
  }
}


PUBLIC void
vrna_sc_add_data(vrna_fold_compound_t       *vc,
                 void                       *data,
                 vrna_callback_free_auxdata *free_data)
{
  if (vc) {
    if (vc->type == VRNA_FC_TYPE_SINGLE) {
      if (!vc->sc)
        vrna_sc_init(vc);

      vc->sc->data      = data;
      vc->sc->free_data = free_data;
    }
  }
}


PUBLIC void
vrna_sc_add_f(vrna_fold_compound_t    *vc,
              vrna_callback_sc_energy *f)
{
  if (vc && f) {
    if (vc->type == VRNA_FC_TYPE_SINGLE) {
      if (!vc->sc)
        vrna_sc_init(vc);

      vc->sc->f = f;
    }
  }
}


PUBLIC void
vrna_sc_add_bt(vrna_fold_compound_t       *vc,
               vrna_callback_sc_backtrack *f)
{
  if (vc && f) {
    if (vc->type == VRNA_FC_TYPE_SINGLE) {
      if (!vc->sc)
        vrna_sc_init(vc);

      vc->sc->bt = f;
    }
  }
}


PUBLIC void
vrna_sc_add_exp_f(vrna_fold_compound_t        *vc,
                  vrna_callback_sc_exp_energy *exp_f)
{
  if (vc && exp_f) {
    if (vc->type == VRNA_FC_TYPE_SINGLE) {
      if (!vc->sc)
        vrna_sc_init(vc);

      vc->sc->exp_f = exp_f;
    }
  }
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
            int                   i,
            int                   start,
            int                   end,
            int                   e)
{
  int size, cnt = 0;

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


PRIVATE void
sc_add_bp(vrna_fold_compound_t  *vc,
          int                   i,
          int                   j,
          FLT_OR_DBL            energy,
          unsigned int          options)
{
  unsigned int  n;
  vrna_sc_t     *sc;
  int           *idx;

  n = vc->length;

  if ((options & VRNA_OPTION_WINDOW) && (!vc->sc))
    vrna_sc_init_window(vc);
  else if (!vc->sc)
    vrna_sc_init(vc);

  sc = vc->sc;
  sc_init_bp_storage(sc);
  sc_store_bp(sc->bp_storage, i, j, j, (int)roundf(energy * 100.));
  sc->state |= STATE_DIRTY_BP_MFE | STATE_DIRTY_BP_PF;
}


PRIVATE INLINE void
free_sc_up(vrna_sc_t *sc)
{
  unsigned int i;

  free(sc->up_storage);

  sc->up_storage = NULL;

  if (sc->energy_up) {
    for (i = 0; i <= sc->n; i++)
      free(sc->energy_up[i]);

    free(sc->energy_up);
    sc->energy_up = NULL;
  }

  if (sc->exp_energy_up) {
    for (i = 0; i <= sc->n; i++)
      free(sc->exp_energy_up[i]);
    free(sc->exp_energy_up);
    sc->exp_energy_up = NULL;
  }

  sc->state &= ~(STATE_DIRTY_UP_MFE | STATE_DIRTY_UP_PF);
}


PRIVATE INLINE void
free_sc_bp(vrna_sc_t *sc)
{
  unsigned int i;

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


PRIVATE void
sc_reset_up(vrna_fold_compound_t  *vc,
            const FLT_OR_DBL      *constraints)
{
  unsigned int  i, j, n;
  vrna_sc_t     *sc;

  n = vc->length;

  if (!vc->sc)
    vrna_sc_init(vc);

  sc = vc->sc;

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


PRIVATE void
sc_reset_bp(vrna_fold_compound_t  *vc,
            const FLT_OR_DBL      **constraints)
{
  unsigned int  i, j, n;
  vrna_sc_t     *sc;
  int           *idx;

  n = vc->length;

  if (!vc->sc)
    vrna_sc_init(vc);

  sc = vc->sc;

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


PRIVATE void
sc_add_up(vrna_fold_compound_t  *vc,
          int                   i,
          FLT_OR_DBL            energy,
          unsigned int          options)
{
  unsigned int  j, u, n;
  vrna_sc_t     *sc;

  n = vc->length;

  if ((options & VRNA_OPTION_WINDOW) && (!vc->sc))
    vrna_sc_init_window(vc);
  else if (!vc->sc)
    vrna_sc_init(vc);

  sc = vc->sc;
  sc_init_up_storage(sc);
  sc->up_storage[i] += (int)roundf(energy * 100.);
  sc->state         |= STATE_DIRTY_UP_MFE | STATE_DIRTY_UP_PF;
}


/* populate sc->energy_up arrays for usage in MFE computations */
PRIVATE void
prepare_sc_up_mfe(vrna_fold_compound_t  *vc,
                  unsigned int          options)
{
  unsigned int  i, j, n;
  vrna_sc_t     *sc;

  n = vc->length;
  switch (vc->type) {
    case VRNA_FC_TYPE_SINGLE:
      sc = vc->sc;
      if (sc) {
        /* prepare sc for unpaired nucleotides only if we actually have some to apply */
        if (sc->up_storage) {
          if (sc->state & STATE_DIRTY_UP_MFE) {
            /*  allocate memory such that we can access the soft constraint
             *  energies of a subsequence of length j starting at position i
             *  via sc->energy_up[i][j]
             */
            sc->energy_up = (int **)vrna_realloc(sc->energy_up, sizeof(int *) * (n + 2));

            if (options & VRNA_OPTION_WINDOW) {
              /*
               *  simply init with NULL pointers, since the sliding-window implementation should take
               *  care of allocating the required memory and filling in appropriate energy contributions
               */
              for (i = 0; i <= n + 1; i++)
                sc->energy_up[i] = NULL;
            } else {
              for (i = 1; i <= n; i++)
                sc->energy_up[i] = (int *)vrna_realloc(sc->energy_up[i], sizeof(int) * (n - i + 2));

              sc->energy_up[0]      = NULL;
              sc->energy_up[n + 1]  = NULL;

              /* now add soft constraints as stored in container for unpaired sc */
              for (i = 1; i <= n; i++) {
                sc->energy_up[i][0] = 0;
                for (j = 1; j <= (n - i + 1); j++)
                  sc->energy_up[i][j] = sc->energy_up[i][j - 1]
                                        + sc->up_storage[i + j - 1];
              }
            }

            sc->state &= ~STATE_DIRTY_UP_MFE;
          }
        } else if (sc->energy_up) {
          /* remove any unpaired sc if storage container is empty */
          free_sc_up(sc);
        }
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      break; /* do nothing for now */
  }
}


/* populate sc->exp_energy_up arrays for usage in partition function computations */
PRIVATE void
prepare_sc_up_pf(vrna_fold_compound_t *vc,
                 unsigned int         options)
{
  unsigned int      i, j, n;
  vrna_sc_t         *sc;

  n = vc->length;

  vrna_exp_param_t  *exp_params = vc->exp_params;
  double            GT          = 0.;
  double            temperature = exp_params->temperature;
  double            kT          = exp_params->kT;
  double            TT          = (temperature + K0) / (Tmeasure);

  switch (vc->type) {
    case VRNA_FC_TYPE_SINGLE:
      sc = vc->sc;
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
               *  simply init with NULL pointers, since the sliding-window implementation should take
               *  care of allocating the required memory and filling in appropriate Boltzmann factors
               */
              for (i = 0; i <= n + 1; i++)
                sc->exp_energy_up[i] = NULL;
            } else {
              for (i = 1; i <= n; i++)
                sc->exp_energy_up[i] =
                  (FLT_OR_DBL *)vrna_realloc(sc->exp_energy_up[i],
                                             sizeof(FLT_OR_DBL) * (n - i + 2));

              sc->exp_energy_up[0]      = NULL;
              sc->exp_energy_up[n + 1]  = NULL;

              for (i = 1; i <= n; i++) {
                sc->exp_energy_up[i][0] = 1.;

                for (j = 1; j <= (n - i + 1); j++) {
                  GT                      = (double)sc->up_storage[i + j - 1] * 10.; /* convert deka-cal/mol to cal/mol */
                  sc->exp_energy_up[i][j] = sc->exp_energy_up[i][j - 1]
                                            * (FLT_OR_DBL)exp(-GT / kT);
                }
              }
            }

            sc->state &= ~STATE_DIRTY_UP_PF;
          }
        }
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      break; /* do nothing for now */
  }
}


/* populate sc->energy_bp arrays for usage in MFE computations */
PRIVATE void
prepare_sc_bp_mfe(vrna_fold_compound_t  *vc,
                  unsigned int          options)
{
  unsigned int  i, j, k, n, cnt;
  int           *idx, turn, e;
  vrna_sc_t     *sc;

  n     = vc->length;
  turn  = vc->params->model_details.min_loop_size;
  idx   = vc->jindx;

  switch (vc->type) {
    case VRNA_FC_TYPE_SINGLE:
      sc = vc->sc;
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

              for (i = 1; i < n; i++) {
                /* check whether we have constraints on any pairing partner i or j */
                if (!(sc->bp_storage[i]))
                  continue;

                for (k = turn + 1; k < n; k++) {
                  j = i + k;
                  e = 0;
                  if (j > n)
                    break;

                  /* go through list of constraints on position i */
                  for (cnt = 0; sc->bp_storage[i][cnt].interval_start != 0; cnt++) {
                    if (sc->bp_storage[i][cnt].interval_start > j)
                      break; /* only constraints for pairs (i,q) with q > j left */

                    if (sc->bp_storage[i][cnt].interval_end < j)
                      continue; /* constraint for pairs (i,q) with q < j */

                    /* constraint has interval [p,q] with p <= j <= q */
                    e += sc->bp_storage[i][cnt].e;
                  }

                  sc->energy_bp[idx[j] + i] = e;
                }
              }
            }

            sc->state &= ~STATE_DIRTY_BP_MFE;
          }
        } else {
          free_sc_bp(sc);
        }
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      break; /* do nothing for now */
  }
}


/* populate sc->exp_energy_bp arrays for usage in partition function computations */
PRIVATE void
prepare_sc_bp_pf(vrna_fold_compound_t *vc,
                 unsigned int         options)
{
  unsigned int      i, j, k, n, cnt;
  vrna_sc_t         *sc;
  int               *idx, *jidx, turn, e;

  n = vc->length;

  vrna_exp_param_t  *exp_params = vc->exp_params;
  double            GT          = 0.;
  double            temperature = exp_params->temperature;
  double            kT          = exp_params->kT;
  double            TT          = (temperature + K0) / (Tmeasure);
  turn  = vc->params->model_details.min_loop_size;
  idx   = vc->iindx;
  jidx  = vc->jindx;

  switch (vc->type) {
    case VRNA_FC_TYPE_SINGLE:
      sc = vc->sc;
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

              for (i = 1; i < n; i++) {
                /* check whether we have constraints on any pairing partner i or j */
                if (!(sc->bp_storage[i]))
                  continue;

                for (k = turn + 1; k < n; k++) {
                  j = i + k;
                  e = 0;
                  if (j > n)
                    break;

                  /* go through list of constraints on position i */
                  for (cnt = 0; sc->bp_storage[i][cnt].interval_start != 0; cnt++) {
                    if (sc->bp_storage[i][cnt].interval_start > j)
                      break; /* only constraints for pairs (i,q) with q > j left */

                    if (sc->bp_storage[i][cnt].interval_end < j)
                      continue; /* constraint for pairs (i,q) with q < j */

                    /* constraint has interval [p,q] with p <= j <= q */
                    e += sc->bp_storage[i][cnt].e;
                  }

                  GT                            = e * 10.;
                  sc->exp_energy_bp[idx[i] - j] = (FLT_OR_DBL)exp(-GT / kT);
                }
              }
            }

            sc->state &= ~STATE_DIRTY_BP_PF;
          }
        }
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      break; /* do nothing for now */
  }
}
