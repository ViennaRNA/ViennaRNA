/* constraints handling */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
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
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/constraints_hard.h"


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
hc_init_up_storage(vrna_hc_t *hc);


PRIVATE void
hc_init_bp_storage(vrna_hc_t *hc);


PRIVATE void
hc_store_bp(vrna_hc_bp_storage_t  **container,
            int                   i,
            int                   start,
            int                   end,
            unsigned char         loop_type);


PRIVATE void
hc_add_up(vrna_fold_compound_t  *vc,
          int                   i,
          unsigned char         option);


PRIVATE INLINE void
hc_cant_pair(unsigned int   i,
             unsigned char  c_option,
             unsigned char  *hc,
             unsigned int   length,
             unsigned int   min_loop_size,
             int            *index);


PRIVATE INLINE void
hc_must_pair(unsigned int   i,
             unsigned char  c_option,
             unsigned char  *hc,
             int            *index);


PRIVATE INLINE void
hc_pairs_upstream(unsigned int  i,
                  unsigned char c_option,
                  unsigned char *hc,
                  unsigned int  length,
                  int           *index);


PRIVATE INLINE void
hc_pairs_downstream(unsigned int  i,
                    unsigned char c_option,
                    unsigned char *hc,
                    unsigned int  length,
                    int           *index);


PRIVATE INLINE void
hc_allow_pair(unsigned int  i,
              unsigned int  j,
              unsigned char c_option,
              unsigned char *hc,
              int           *index);


PRIVATE INLINE void
hc_weak_enforce_pair(unsigned int   i,
                     unsigned int   j,
                     unsigned char  c_option,
                     unsigned char  *hc,
                     unsigned int   length,
                     unsigned int   min_loop_size,
                     int            *index);


PRIVATE INLINE void
hc_enforce_pair(unsigned int  i,
                unsigned int  j,
                unsigned char c_option,
                unsigned char *hc,
                unsigned int  length,
                unsigned int  min_loop_size,
                int           *index);


PRIVATE INLINE void
hc_intramolecular_only(unsigned int   i,
                       unsigned char  c_option,
                       unsigned char  *hc,
                       unsigned int   length,
                       unsigned int   min_loop_size,
                       int            cut,
                       int            *index);


PRIVATE INLINE void
hc_intermolecular_only(unsigned int   i,
                       unsigned char  c_option,
                       unsigned char  *hc,
                       unsigned int   length,
                       unsigned int   min_loop_size,
                       int            cut,
                       int            *index);


PRIVATE void
apply_DB_constraint(vrna_fold_compound_t  *vc,
                    const char            *constraint,
                    unsigned int          options);


PRIVATE void
hc_reset_to_default(vrna_fold_compound_t *vc);


PRIVATE void
hc_update_up(vrna_fold_compound_t *vc);


PRIVATE void
hc_update_up_window(vrna_fold_compound_t  *vc,
                    int                   i);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC void
vrna_message_constraint_options_all(void)
{
  vrna_message_constraint_options(VRNA_CONSTRAINT_DB_PIPE
                                  | VRNA_CONSTRAINT_DB_DOT
                                  | VRNA_CONSTRAINT_DB_X
                                  | VRNA_CONSTRAINT_DB_ANG_BRACK
                                  | VRNA_CONSTRAINT_DB_RND_BRACK);
}


PUBLIC void
vrna_message_constraint_options(unsigned int option)
{
  printf("Input structure constraints using the following notation:\n");
  if (option & VRNA_CONSTRAINT_DB_PIPE)
    printf("| : paired with another base\n");

  if (option & VRNA_CONSTRAINT_DB_DOT)
    printf(". : no constraint at all\n");

  if (option & VRNA_CONSTRAINT_DB_X)
    printf("x : base must not pair\n");

  if (option & VRNA_CONSTRAINT_DB_ANG_BRACK)
    printf("< : base i is paired with a base j<i\n> : base i is paired with a base j>i\n");

  if (option & VRNA_CONSTRAINT_DB_RND_BRACK)
    printf("matching brackets ( ): base i pairs base j\n");
}


PUBLIC void
vrna_hc_init(vrna_fold_compound_t *vc)
{
  unsigned int  n;
  vrna_hc_t     *hc;

  n = vc->length;

  /* free previous hard constraints */
  vrna_hc_free(vc->hc);

  /* allocate memory new hard constraints data structure */
  hc          = (vrna_hc_t *)vrna_alloc(sizeof(vrna_hc_t));
  hc->type    = VRNA_HC_DEFAULT;
  hc->n       = n;
  hc->matrix  = (unsigned char *)vrna_alloc(sizeof(unsigned char) * ((n * (n + 1)) / 2 + 2));
  hc->up_ext  = (int *)vrna_alloc(sizeof(int) * (n + 2));
  hc->up_hp   = (int *)vrna_alloc(sizeof(int) * (n + 2));
  hc->up_int  = (int *)vrna_alloc(sizeof(int) * (n + 2));
  hc->up_ml   = (int *)vrna_alloc(sizeof(int) * (n + 2));

  /* set new hard constraints */
  vc->hc = hc;

  /* prefill default values  */
  hc_reset_to_default(vc);

  /* add null pointers for the generalized hard constraint feature */
  hc->f         = NULL;
  hc->data      = NULL;
  hc->free_data = NULL;

  /* update */
  hc_update_up(vc);
}


PUBLIC void
vrna_hc_init_window(vrna_fold_compound_t *vc)
{
  unsigned int  i, n, window_size;
  vrna_hc_t     *hc;

  n           = vc->length;
  window_size = vc->window_size;

  /* free previous hard constraints */
  vrna_hc_free(vc->hc);

  /* allocate memory new hard constraints data structure */
  hc                = (vrna_hc_t *)vrna_alloc(sizeof(vrna_hc_t));
  hc->type          = VRNA_HC_WINDOW;
  hc->n             = n;
  hc->matrix_local  = (unsigned char **)vrna_alloc(sizeof(unsigned char *) * (n + 2));
  hc->up_storage    = NULL;
  hc->bp_storage    = NULL;
  hc->up_ext        = NULL;
  hc->up_hp         = NULL;
  hc->up_int        = NULL;
  hc->up_ml         = NULL;

  /* set new hard constraints */
  vc->hc = hc;

  /* add null pointers for the generalized hard constraint feature */
  hc->f         = NULL;
  hc->data      = NULL;
  hc->free_data = NULL;
}


PUBLIC void
vrna_hc_update(vrna_fold_compound_t *vc,
               int                  i)
{
  int       j, k, type, n, maxdist, pairSize, turn, noLP;
  short     *S;
  char      **ptype;
  vrna_md_t *md;
  vrna_hc_t *hc;

  n         = (int)vc->length;
  S         = vc->sequence_encoding2;
  hc        = vc->hc;
  ptype     = vc->ptype_local;
  maxdist   = vc->window_size;
  md        = &(vc->params->model_details);
  turn      = md->min_loop_size;
  noLP      = md->noLP;
  pairSize  = md->max_bp_span;

  /* init up_xx arrays if necessary */
  if (!hc->up_ext) {
    hc->up_ext  = (int *)vrna_alloc(sizeof(int) * (n + 2));
    hc->up_hp   = (int *)vrna_alloc(sizeof(int) * (n + 2));
    hc->up_int  = (int *)vrna_alloc(sizeof(int) * (n + 2));
    hc->up_ml   = (int *)vrna_alloc(sizeof(int) * (n + 2));

    hc_update_up(vc);
  }

  /* ######################### */
  /* fill with default values  */
  /* ######################### */

  if (hc->up_storage) {
    /* We use user-defined constraints for unpaired nucleotides */
    hc->matrix_local[i][0] = hc->up_storage[i];
  } else {
    /* ... or simply allow unpaired nucleotides in all contexts */
    hc->matrix_local[i][0] = VRNA_CONSTRAINT_CONTEXT_EXT_LOOP
                             | VRNA_CONSTRAINT_CONTEXT_HP_LOOP
                             | VRNA_CONSTRAINT_CONTEXT_INT_LOOP
                             | VRNA_CONSTRAINT_CONTEXT_MB_LOOP;
  }

  /* 2. add default base pairing rules */
  switch (vc->type) {
    case VRNA_FC_TYPE_SINGLE:
      for (k = turn + 1; k < maxdist; k++) {
        j = i + k;
        if (j > n)
          break;

        unsigned char opt = (unsigned char)0;
        if ((j - i + 1) <= pairSize) {
          type = md->pair[S[i]][S[j]];
          switch (type) {
            case 0:
              break;
            case 3:
            /* fallthrough */
            case 4:
              if (md->noGU) {
                break;
              } else if (md->noGUclosure) {
                opt = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
                opt &= ~(VRNA_CONSTRAINT_CONTEXT_HP_LOOP | VRNA_CONSTRAINT_CONTEXT_MB_LOOP);
                break;
              }

            /* else fallthrough */
            default:
              opt = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
              break;
          }
        }

        /* check whether we have constraints on any pairing partner i or j */
        if ((hc->bp_storage) && (hc->bp_storage[i])) {
          unsigned char constraint = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
          int           cnt, set;
          /* go through list of constraints on position i */
          for (set = cnt = 0; hc->bp_storage[i][cnt].interval_start != 0; cnt++) {
            if (hc->bp_storage[i][cnt].interval_start > j)
              break; /* only constraints for pairs (i,q) with q > j left */

            if (hc->bp_storage[i][cnt].interval_end < j)
              continue; /* constraint for pairs (i,q) with q < j */

            /* constraint has interval [p,q] with p <= j <= q */
            constraint  &= hc->bp_storage[i][cnt].loop_type;
            set         = 1;
          }
          if (set && (opt == (unsigned char)0))
            /* overwrite if bp is non-canonical */
            opt = constraint;
          else
            /* apply constraint to canonical bp */
            opt &= constraint;
        }

        hc->matrix_local[i][j - i] = opt;
      }

      /* here space for noLP option */
      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      for (k = turn + 1; k <= maxdist; k++) {
        j = i + k;
        if (j > n)
          break;

        unsigned char opt = (unsigned char)0;
        if ((j - i + 1) <= pairSize)
          if (vc->pscore_local[i][j - i] >= md->cv_fact * MINPSCORE)
            opt = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

        /* check whether we have constraints on any pairing partner i or j */
        if ((hc->bp_storage) && (hc->bp_storage[i])) {
          unsigned char constraint = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
          int           cnt, set;
          /* go through list of constraints on position i */
          for (set = cnt = 0; hc->bp_storage[i][cnt].interval_start != 0; cnt++) {
            if (hc->bp_storage[i][cnt].interval_start > j)
              break; /* only constraints for pairs (i,q) with q > j left */

            if (hc->bp_storage[i][cnt].interval_end < j)
              continue; /* constraint for pairs (i,q) with q < j */

            /* constraint has interval [p,q] with p <= j <= q */
            constraint  &= hc->bp_storage[i][cnt].loop_type;
            set         = 1;
          }
          if (set && (opt == (unsigned char)0))
            /* overwrite if bp is non-canonical */
            opt = constraint;
          else
            /* apply constraint to canonical bp */
            opt &= constraint;
        }

        hc->matrix_local[i][j - i] = opt;
      }
      break;

    default:
      break;
  }

  hc_update_up_window(vc, i);
}


PUBLIC void
vrna_hc_add_up(vrna_fold_compound_t *vc,
               int                  i,
               unsigned char        option)
{
  int j;

  if (vc) {
    if (vc->hc) {
      if ((i <= 0) || (i > vc->length)) {
        vrna_message_warning("vrna_hc_add_up: position out of range, not doing anything");
        return;
      }

      hc_add_up(vc, i, option);

      if (vc->hc->type != VRNA_HC_WINDOW)
        hc_update_up(vc);
    }
  }
}


PUBLIC int
vrna_hc_add_up_batch(vrna_fold_compound_t *vc,
                     vrna_hc_up_t         *constraints)
{
  int i, ret;

  ret = 0; /* failure */

  if (vc) {
    if (vc->hc && constraints) {
      for (i = 0; constraints[i].position != 0; i++) {
        int           pos     = constraints[i].position;
        unsigned char options = constraints[i].options;
        if ((pos <= 0) || (pos > vc->length)) {
          vrna_message_warning(
            "vrna_hc_add_up_batch: position out of range, application of hard constraints stops here!");
          return ret;
        }

        hc_add_up(vc, pos, options);
      }

      if (vc->hc->type != VRNA_HC_WINDOW)
        hc_update_up(vc);

      ret = 1; /* success */
    }
  }

  return ret;
}


PRIVATE void
hc_init_up_storage(vrna_hc_t *hc)
{
  unsigned int i;

  if (hc->up_storage == NULL) {
    free(hc->up_storage);
    hc->up_storage = (unsigned char *)vrna_alloc(sizeof(unsigned char) * (hc->n + 2));

    for (i = 1; i <= hc->n; i++) {
      /* by default unpaired nucleotides are allowed in all contexts */
      hc->up_storage[i] = VRNA_CONSTRAINT_CONTEXT_EXT_LOOP
                          | VRNA_CONSTRAINT_CONTEXT_HP_LOOP
                          | VRNA_CONSTRAINT_CONTEXT_INT_LOOP
                          | VRNA_CONSTRAINT_CONTEXT_MB_LOOP;
    }
  }
}


PRIVATE void
hc_init_bp_storage(vrna_hc_t *hc)
{
  unsigned int i;

  if (hc->bp_storage == NULL) {
    hc->bp_storage = (vrna_hc_bp_storage_t **)vrna_alloc(
      sizeof(vrna_hc_bp_storage_t *) * (hc->n + 2));

    for (i = 1; i <= hc->n; i++)
      /* by default we do not limit base pairs to any context */
      hc->bp_storage[i] = NULL;
  }
}


PRIVATE void
hc_store_bp(vrna_hc_bp_storage_t  **container,
            int                   i,
            int                   start,
            int                   end,
            unsigned char         loop_type)
{
  int size, cnt = 0;

  if (!container[i]) {
    container[i] = (vrna_hc_bp_storage_t *)vrna_alloc(sizeof(vrna_hc_bp_storage_t) * 2);
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
    container[i] = (vrna_hc_bp_storage_t *)vrna_realloc(container[i],
                                                        sizeof(vrna_hc_bp_storage_t) * (size + 2));
    /* shift trailing constraints by 1 entry */
    memmove(container[i] + cnt + 1, container[i] + cnt,
            sizeof(vrna_hc_bp_storage_t) * (size - cnt + 1));
  }

  container[i][cnt].interval_start  = start;
  container[i][cnt].interval_end    = end;
  container[i][cnt].loop_type       = loop_type;
}


PRIVATE void
hc_add_up(vrna_fold_compound_t  *vc,
          int                   i,
          unsigned char         option)
{
  int           j;
  unsigned char type = (unsigned char)0;

  if (vc->hc->type == VRNA_HC_WINDOW) {
    if (option & VRNA_CONSTRAINT_CONTEXT_ENFORCE) {
      hc_init_up_storage(vc->hc);
      type = option & (unsigned char)(VRNA_CONSTRAINT_CONTEXT_EXT_LOOP
                                      | VRNA_CONSTRAINT_CONTEXT_HP_LOOP
                                      | VRNA_CONSTRAINT_CONTEXT_INT_LOOP
                                      | VRNA_CONSTRAINT_CONTEXT_MB_LOOP);

      vc->hc->up_storage[i] = type;

      if (!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
        hc_init_bp_storage(vc->hc);
        /* add constraints for all pairs (j, i) with j < i */
        for (j = 1; j < i; j++)
          hc_store_bp(vc->hc->bp_storage, j, i, i, (unsigned char)0);

        /* add constraints for all pairs (i, j) with i < j */
        hc_store_bp(vc->hc->bp_storage, i, i + 1, vc->length, (unsigned char)0);
      }
    } else {
      hc_init_up_storage(vc->hc);
      vc->hc->up_storage[i] = (unsigned char)(VRNA_CONSTRAINT_CONTEXT_EXT_LOOP
                                              | VRNA_CONSTRAINT_CONTEXT_HP_LOOP
                                              | VRNA_CONSTRAINT_CONTEXT_INT_LOOP
                                              | VRNA_CONSTRAINT_CONTEXT_MB_LOOP);

      if (!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
        hc_init_bp_storage(vc->hc);
        /* add constraints for all pairs (j, i) with j < i */
        for (j = 1; j < i; j++)
          hc_store_bp(vc->hc->bp_storage, j, i, i, ~type);

        /* add constraints for all pairs (i, j) with i < j */
        hc_store_bp(vc->hc->bp_storage, i, i + 1, vc->length, ~type);
      }
    }
  } else {
    if (option & VRNA_CONSTRAINT_CONTEXT_ENFORCE) {
      /* force nucleotide to appear unpaired within a certain type of loop */
      /* do not allow i to be paired with any other nucleotide */
      if (!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
        for (j = 1; j < i; j++)
          vc->hc->matrix[vc->jindx[i] + j] = (unsigned char)0;
        for (j = i + 1; j <= vc->length; j++)
          vc->hc->matrix[vc->jindx[j] + i] = (unsigned char)0;
      }

      type = option & (unsigned char)(VRNA_CONSTRAINT_CONTEXT_EXT_LOOP
                                      | VRNA_CONSTRAINT_CONTEXT_HP_LOOP
                                      | VRNA_CONSTRAINT_CONTEXT_INT_LOOP
                                      | VRNA_CONSTRAINT_CONTEXT_MB_LOOP);

      vc->hc->matrix[vc->jindx[i] + i] = type;
    } else {
      type = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

      /* do not allow i to be paired with any other nucleotide (in context type) */
      if (!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
        for (j = 1; j < i; j++)
          vc->hc->matrix[vc->jindx[i] + j] &= ~type;
        for (j = i + 1; j <= vc->length; j++)
          vc->hc->matrix[vc->jindx[j] + i] &= ~type;
      }

      vc->hc->matrix[vc->jindx[i] + i] = (unsigned char)(VRNA_CONSTRAINT_CONTEXT_EXT_LOOP
                                                         | VRNA_CONSTRAINT_CONTEXT_HP_LOOP
                                                         | VRNA_CONSTRAINT_CONTEXT_INT_LOOP
                                                         | VRNA_CONSTRAINT_CONTEXT_MB_LOOP);
    }
  }
}


PUBLIC void
vrna_hc_add_bp_nonspecific(vrna_fold_compound_t *vc,
                           int                  i,
                           int                  d,
                           unsigned char        option)
{
  int           p;
  unsigned char type, t1, t2;

  if (vc) {
    if (vc->hc) {
      if ((i <= 0) || (i > vc->length)) {
        vrna_message_warning("vrna_hc_add_bp_nonspecific: position out of range, not doing anything");
        return;
      }

      /* position i may pair in provided contexts */
      type = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
      /* acknowledge pairing direction */
      t1  = (d <= 0) ? type : (unsigned char)0;
      t2  = (d >= 0) ? type : (unsigned char)0;

      if (vc->hc->type == VRNA_HC_WINDOW) {
        /* nucleotide mustn't be unpaired */
        hc_init_up_storage(vc->hc);
        vc->hc->up_storage[i] = (unsigned char)0;

        /* force pairing direction */
        hc_init_bp_storage(vc->hc);
        for (p = 1; p < i; p++)
          hc_store_bp(vc->hc->bp_storage, p, i, i, t1);

        hc_store_bp(vc->hc->bp_storage, i, i + 1, vc->length, t2);
      } else {
        if (option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE) {
          /* only allow for possibly non-canonical pairs, do not enforce them */
          for (p = 1; p < i; p++)
            vc->hc->matrix[vc->jindx[i] + p] |= t1;
          for (p = i + 1; p <= vc->length; p++)
            vc->hc->matrix[vc->jindx[p] + i] |= t2;
        } else {
          /* force pairing direction */
          for (p = 1; p < i; p++)
            vc->hc->matrix[vc->jindx[i] + p] &= t1;
          for (p = i + 1; p <= vc->length; p++)
            vc->hc->matrix[vc->jindx[p] + i] &= t2;
          /* nucleotide mustn't be unpaired */
          vc->hc->matrix[vc->jindx[i] + i] = (unsigned char)0;
        }

        hc_update_up(vc);
      }
    }
  }
}


PUBLIC void
vrna_hc_add_bp(vrna_fold_compound_t *vc,
               int                  i,
               int                  j,
               unsigned char        option)
{
  int           k, l;
  unsigned char type;

  if (vc) {
    if (vc->hc) {
      if ((i <= 0) || (j <= i) || (j > vc->length)) {
        vrna_message_warning("vrna_hc_add_bp: position out of range, not doing anything");
        return;
      }

      if (vc->hc->type == VRNA_HC_WINDOW) {
        hc_init_bp_storage(vc->hc);
        hc_store_bp(vc->hc->bp_storage, i, j, j, option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS);

        if (!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
          /*
           * remove all conflicting base pairs, i.e. do not allow i or j to pair
           * with any other nucleotide k
           */
          for (k = 1; k < i; k++)
            hc_store_bp(vc->hc->bp_storage, k, i, j, (unsigned char)0);             /* (k, i), (k, i + 1), ..., (k, j) with 1 <= k < i */

          hc_store_bp(vc->hc->bp_storage, i, i + 1, j - 1, (unsigned char)0);       /* (i, k), i < k < j */

          for (k = i + 1; k < j; k++)
            hc_store_bp(vc->hc->bp_storage, k, j, vc->length, (unsigned char)0);    /* (i + 1, k), (i + 1, k), ..., (j - 1, k) with (j < k <= n */

          hc_store_bp(vc->hc->bp_storage, i, j + 1, vc->length, (unsigned char)0);  /* (i, k), j < k <= n */
          hc_store_bp(vc->hc->bp_storage, j, j + 1, vc->length, (unsigned char)0);  /* (j, k), j < k <= n */
        }

        if (option & VRNA_CONSTRAINT_CONTEXT_ENFORCE) {
          /* do not allow i,j to be unpaired */
          hc_init_up_storage(vc->hc);
          vc->hc->up_storage[i] = (unsigned char)0;
          vc->hc->up_storage[j] = (unsigned char)0;
        }
      } else {
        /* reset ptype in case (i,j) is a non-canonical pair */
        if (option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS) {
          if (vc->hc->matrix[vc->jindx[j] + i])
            if (vc->ptype[vc->jindx[j] + i] == 0)
              vc->ptype[vc->jindx[j] + i] = 7;
        }

        vc->hc->matrix[vc->jindx[j] + i] = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

        if (!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
          /*
           * remove all conflicting base pairs, i.e. do not allow i,j to pair
           * with any other nucleotide k
           */
          for (k = 1; k < i; k++) {
            vc->hc->matrix[vc->jindx[i] + k]  = (unsigned char)0;
            vc->hc->matrix[vc->jindx[j] + k]  = (unsigned char)0;
            for (l = i + 1; l < j; l++)
              vc->hc->matrix[vc->jindx[l] + k] = (unsigned char)0;
          }
          for (k = i + 1; k < j; k++) {
            vc->hc->matrix[vc->jindx[k] + i]  = (unsigned char)0;
            vc->hc->matrix[vc->jindx[j] + k]  = (unsigned char)0;
            for (l = j + 1; l <= vc->length; l++)
              vc->hc->matrix[vc->jindx[l] + k] = (unsigned char)0;
          }
          for (k = j + 1; k <= vc->length; k++) {
            vc->hc->matrix[vc->jindx[k] + i]  = (unsigned char)0;
            vc->hc->matrix[vc->jindx[k] + j]  = (unsigned char)0;
          }
        }

        if (option & VRNA_CONSTRAINT_CONTEXT_ENFORCE) {
          /* do not allow i,j to be unpaired */
          vc->hc->matrix[vc->jindx[i] + i]  = (unsigned char)0;
          vc->hc->matrix[vc->jindx[j] + j]  = (unsigned char)0;

          hc_update_up(vc);
        }
      }
    }
  }
}


PUBLIC void
vrna_hc_free(vrna_hc_t *hc)
{
  if (hc) {
    if (hc->type == VRNA_HC_DEFAULT) {
      free(hc->matrix);
    } else if (hc->type == VRNA_HC_WINDOW) {
      unsigned int i;
      free(hc->matrix_local);
      free(hc->up_storage);
      if (hc->bp_storage) {
        for (i = 1; i <= hc->n; i++)
          free(hc->bp_storage[i]);
        free(hc->bp_storage);
      }
    }

    free(hc->up_ext);
    free(hc->up_hp);
    free(hc->up_int);
    free(hc->up_ml);

    if (hc->free_data)
      hc->free_data(hc->data);

    free(hc);
  }
}


PUBLIC void
vrna_hc_add_f(vrna_fold_compound_t      *vc,
              vrna_callback_hc_evaluate *f)
{
  if (vc && f) {
    if (vc->type == VRNA_FC_TYPE_SINGLE) {
      if (!vc->hc)
        vrna_hc_init(vc);

      vc->hc->f = f;
    }
  }
}


PUBLIC void
vrna_hc_add_data(vrna_fold_compound_t       *vc,
                 void                       *data,
                 vrna_callback_free_auxdata *f)
{
  if (vc && data) {
    if (vc->type == VRNA_FC_TYPE_SINGLE) {
      if (!vc->hc)
        vrna_hc_init(vc);

      vc->hc->data      = data;
      vc->hc->free_data = f;
    }
  }
}


PUBLIC int
vrna_hc_add_from_db(vrna_fold_compound_t  *vc,
                    const char            *constraint,
                    unsigned int          options)
{
  const char  *structure_constraint;
  char        *tmp;
  int         i, d, ret;
  vrna_md_t   *md;

  ret = 0; /* Failure */

  if (vc) {
    tmp = NULL;
    if (vc->params)
      md = &(vc->params->model_details);
    else if (vc->exp_params)
      md = &(vc->exp_params->model_details);
    else
      return ret;

    if (!vc->hc)
      vrna_hc_init(vc);

    if (options & VRNA_CONSTRAINT_DB_WUSS) {
      tmp                   = vrna_db_from_WUSS(constraint);
      structure_constraint  = (const char *)tmp;
    } else {
      structure_constraint = constraint;
    }

    /* apply hard constraints from dot-bracket notation */
    apply_DB_constraint(vc, structure_constraint, options);
    hc_update_up(vc);

    ret = 1; /* Success */

    free(tmp);
  }

  return ret;
}


PRIVATE void
apply_DB_constraint(vrna_fold_compound_t  *vc,
                    const char            *constraint,
                    unsigned int          options)
{
  char          c_option, *hc, *sequence;
  short         *S;
  unsigned int  length, min_loop_size;
  int           n, i, j, hx, *stack, *index, cut;
  vrna_md_t     *md;

  if (constraint == NULL)
    return;

  sequence      = vc->sequence;
  length        = (int)vc->length;
  S             = vc->sequence_encoding2;
  hc            = vc->hc->matrix;
  md            = &(vc->params->model_details);
  min_loop_size = md->min_loop_size;
  cut           = vc->cutpoint;
  n             = (int)strlen(constraint);
  stack         = (int *)vrna_alloc(sizeof(int) * (n + 1));
  index         = vrna_idx_col_wise(length);
  c_option      = VRNA_CONSTRAINT_CONTEXT_EXT_LOOP
                  | VRNA_CONSTRAINT_CONTEXT_HP_LOOP
                  | VRNA_CONSTRAINT_CONTEXT_INT_LOOP
                  | VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC
                  | VRNA_CONSTRAINT_CONTEXT_MB_LOOP
                  | VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC;

  for (hx = 0, j = 1; j <= n; j++) {
    switch (constraint[j - 1]) {
      /* can't pair */
      case 'x':
        if (options & VRNA_CONSTRAINT_DB_X)
          hc_cant_pair(j, c_option, hc, length, min_loop_size, index);

        break;

      /* must pair, i.e. may not be unpaired */
      case '|':
        if (options & VRNA_CONSTRAINT_DB_PIPE)
          if (options & VRNA_CONSTRAINT_DB_ENFORCE_BP)
            hc_must_pair(j, c_option, hc, index);

        break;

      /* weak enforced pair 'open' */
      case '(':
        if (options & VRNA_CONSTRAINT_DB_RND_BRACK)
          stack[hx++] = j;

        break;

      /* weak enforced pair 'close' */
      case ')':
        if (options & VRNA_CONSTRAINT_DB_RND_BRACK) {
          if (hx <= 0)
            vrna_message_error("%s\nunbalanced brackets in constraints", constraint);

          i = stack[--hx];

          if (options & VRNA_CONSTRAINT_DB_CANONICAL_BP) {
            /* check whether this pair forms a non-canoncial base pair */
            if (md->pair[S[i]][S[j]] == 0) {
              vrna_message_warning("Removing non-canonical base pair %c%c (%d,%d) from constraint",
                                   sequence[i - 1], sequence[j - 1],
                                   i, j);
              break;
            }
          }

          if (options & VRNA_CONSTRAINT_DB_ENFORCE_BP)
            hc_enforce_pair(i, j, c_option, hc, length, min_loop_size, index);
          else
            hc_weak_enforce_pair(i, j, c_option, hc, length, min_loop_size, index);
        }

        break;

      /* pairs upstream */
      case '<':
        if (options & VRNA_CONSTRAINT_DB_ANG_BRACK) {
          hc_pairs_downstream(j, c_option, hc, length, index);
          if (options & VRNA_CONSTRAINT_DB_ENFORCE_BP)
            hc_must_pair(j, c_option, hc, index);
        }

        break;

      /* pairs downstream */
      case '>':
        if (options & VRNA_CONSTRAINT_DB_ANG_BRACK) {
          hc_pairs_upstream(j, c_option, hc, length, index);
          if (options & VRNA_CONSTRAINT_DB_ENFORCE_BP)
            hc_must_pair(j, c_option, hc, index);
        }

        break;

      /* only intramolecular basepairing */
      case 'l':
        if (options & VRNA_CONSTRAINT_DB_INTRAMOL)
          hc_intramolecular_only(j, c_option, hc, length, min_loop_size, cut, index);

        break;

      /* only intermolecular bp */
      case 'e':
        if (options & VRNA_CONSTRAINT_DB_INTERMOL)
          hc_intermolecular_only(j, c_option, hc, length, min_loop_size, cut, index);

        break;

      case '.':
        break;

      default:
        vrna_message_warning(
          "Unrecognized character '%c' in pseudo dot-bracket notation constraint string",
          constraint[j - 1]);
        break;
    }
  }

  if (hx != 0)
    vrna_message_error("%s\nunbalanced brackets in constraint string", constraint);

  /* clean up */
  free(index);
  free(stack);
}


PRIVATE INLINE void
hc_intramolecular_only(unsigned int   i,
                       unsigned char  c_option,
                       unsigned char  *hc,
                       unsigned int   length,
                       unsigned int   min_loop_size,
                       int            cut,
                       int            *index)
{
  unsigned int l;

  if (cut > 1) {
    if (i < cut)
      for (l = MAX2(i + min_loop_size, cut); l <= length; l++)
        hc[index[l] + i] &= ~c_option;
    else
      for (l = 1; l < MIN2(cut, i - min_loop_size); l++)
        hc[index[i] + l] &= ~c_option;
  }
}


PRIVATE INLINE void
hc_intermolecular_only(unsigned int   i,
                       unsigned char  c_option,
                       unsigned char  *hc,
                       unsigned int   length,
                       unsigned int   min_loop_size,
                       int            cut,
                       int            *index)
{
  unsigned int l;

  if (cut > 1) {
    if (i < cut) {
      for (l = 1; l < i; l++)
        hc[index[i] + l] &= ~c_option;
      for (l = i + 1; l < cut; l++)
        hc[index[l] + i] &= ~c_option;
    } else {
      for (l = cut; l < i; l++)
        hc[index[i] + l] &= ~c_option;
      for (l = i + 1; l <= length; l++)
        hc[index[l] + i] &= ~c_option;
    }
  }
}


PRIVATE INLINE void
hc_cant_pair(unsigned int   i,
             unsigned char  c_option,
             unsigned char  *hc,
             unsigned int   length,
             unsigned int   min_loop_size,
             int            *index)
{
  hc_pairs_upstream(i, c_option, hc, length, index);
  hc_pairs_downstream(i, c_option, hc, length, index);
}


PRIVATE INLINE void
hc_must_pair(unsigned int   i,
             unsigned char  c_option,
             unsigned char  *hc,
             int            *index)
{
  hc[index[i] + i] &= ~c_option;
}


PRIVATE INLINE void
hc_pairs_upstream(unsigned int  i,
                  unsigned char c_option,
                  unsigned char *hc,
                  unsigned int  length,
                  int           *index)
{
  unsigned int l;

  /* prohibit downstream pairs */
  for (l = length; l > i; l--)
    hc[index[l] + i] = (unsigned char)0;
  /* allow upstream pairs of given type */
  for (l = i - 1; l >= 1; l--)
    hc[index[i] + l] &= c_option;
}


PRIVATE INLINE void
hc_pairs_downstream(unsigned int  i,
                    unsigned char c_option,
                    unsigned char *hc,
                    unsigned int  length,
                    int           *index)
{
  unsigned int l;

  /* allow downstream pairs of given type */
  for (l = length; l > i; l--)
    hc[index[l] + i] &= c_option;
  /* forbid upstream pairs */
  for (l = i - 1; l >= 1; l--)
    hc[index[i] + l] = (unsigned char)0;
}


PRIVATE INLINE void
hc_allow_pair(unsigned int  i,
              unsigned int  j,
              unsigned char c_option,
              unsigned char *hc,
              int           *index)
{
  hc[index[j] + i] |= c_option;
}


PRIVATE INLINE void
hc_weak_enforce_pair(unsigned int   i,
                     unsigned int   j,
                     unsigned char  c_option,
                     unsigned char  *hc,
                     unsigned int   length,
                     unsigned int   min_loop_size,
                     int            *index)
{
  unsigned int k, l;

  /* don't allow pairs (k,i) 1 <= k < i */
  /* don't allow pairs (i,k) i < k <= n */
  hc_pairs_upstream(i, (unsigned char)0, hc, length, index);
  /* don't allow pairs (k,j) 1 <= k < j */
  /* don't allow pairs (j,k) j < k <= n */
  hc_pairs_upstream(j, (unsigned char)0, hc, length, index);

  /* don't allow pairs i < k < j < l */
  for (k = i + 1; k < j; k++)
    for (l = j + 1; l <= length; l++)
      hc[index[l] + k] = 0;

  /* don't allow pairs k<i<l<j */
  for (k = 1; k < i; k++)
    for (l = i + 1; l < j; l++)
      hc[index[l] + k] = 0;

  /* allow base pair (i,j) */
  hc[index[j] + i] |= c_option;
}


PRIVATE INLINE void
hc_enforce_pair(unsigned int  i,
                unsigned int  j,
                unsigned char c_option,
                unsigned char *hc,
                unsigned int  length,
                unsigned int  min_loop_size,
                int           *index)
{
  hc_weak_enforce_pair(i,
                       j,
                       c_option,
                       hc,
                       length,
                       min_loop_size,
                       index);

  /* forbid i and j to be unpaired */
  hc[index[i] + i]  = 0;
  hc[index[j] + j]  = 0;
}


PRIVATE void
hc_reset_to_default(vrna_fold_compound_t *vc)
{
  unsigned int  i, j, ij, min_loop_size, n;
  int           max_span, *idx;
  vrna_md_t     *md;
  vrna_hc_t     *hc;
  short         *S;

  md  = NULL;
  n   = vc->length;
  hc  = vc->hc;
  idx = vc->jindx;
  S   = vc->sequence_encoding2;

  if (vc->params)
    md = &(vc->params->model_details);
  else if (vc->exp_params)
    md = &(vc->exp_params->model_details);
  else
    vrna_message_error("missing model_details in fold_compound");

  min_loop_size = md->min_loop_size;
  max_span      = md->max_bp_span;

  if ((max_span < 5) || (max_span > n))
    max_span = n;

  /* ######################### */
  /* fill with default values  */
  /* ######################### */

  /* 1. unpaired nucleotides are allowed in all contexts */
  for (i = 1; i <= n; i++)
    hc->matrix[idx[i] + i] = VRNA_CONSTRAINT_CONTEXT_EXT_LOOP
                             | VRNA_CONSTRAINT_CONTEXT_HP_LOOP
                             | VRNA_CONSTRAINT_CONTEXT_INT_LOOP
                             | VRNA_CONSTRAINT_CONTEXT_MB_LOOP;

  /* 2. all base pairs with pscore above threshold are allowed in all contexts */
  switch (vc->type) {
    case VRNA_FC_TYPE_COMPARATIVE:
      for (j = n; j > min_loop_size + 1; j--) {
        ij = idx[j] + 1;
        for (i = 1; i < j - min_loop_size; i++, ij++) {
          unsigned char opt = (unsigned char)0;
          if ((j - i + 1) <= max_span)
            if (vc->pscore[idx[j] + i] >= md->cv_fact * MINPSCORE)
              opt = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

          hc->matrix[ij] = opt;
        }
      }
      break;

    case VRNA_FC_TYPE_SINGLE:
      for (j = n; j > min_loop_size + 1; j--) {
        ij = idx[j] + 1;
        for (i = 1; i < j - min_loop_size; i++, ij++) {
          unsigned char opt = (unsigned char)0;
          if ((j - i + 1) <= max_span) {
            int t = md->pair[S[i]][S[j]];
            switch (t) {
              case 0:
                break;
              case 3:                               /* fallthrough */
              case 4:
                if (md->noGU) {
                  break;
                } else if (md->noGUclosure) {
                  opt = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
                  opt &= ~(VRNA_CONSTRAINT_CONTEXT_HP_LOOP | VRNA_CONSTRAINT_CONTEXT_MB_LOOP);
                  break;
                }                                     /* else fallthrough */

              default:
                opt = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
                break;
            }
          }

          hc->matrix[ij] = opt;
        }
      }

      /* correct for no lonely pairs (assuming that ptypes already incorporate noLP status) */
      /* this should be fixed such that ij loses its hard constraint type if it does not
       * allow for enclosing an interior loop, etc.
       */
      /*  ???????
       *  Is this necessary? We could leave the noLP option somewhere else, i.e. do not enforce it
       *  on the level of ptype/constraints, but an the level of recursions...
       *  ???????
       */
      if (md->noLP) {
        if (!vc->ptype)
          vc->ptype = vrna_ptypes(vc->sequence_encoding2, md);

        for (i = 1; i < n; i++)
          for (j = i + min_loop_size + 1; j <= n; j++) {
            if (hc->matrix[idx[j] + i])
              if (!vc->ptype[idx[j] + i])
                hc->matrix[idx[j] + i] = (unsigned char)0;
          }
      }

      break;

    default:
      break;
  }

  /* should we reset the generalized hard constraint feature here? */
  if (hc->f || hc->data) {
    if (hc->free_data)
      hc->free_data(hc->data);

    hc->f         = NULL;
    hc->data      = NULL;
    hc->free_data = NULL;
  }
}


PRIVATE void
hc_update_up(vrna_fold_compound_t *vc)
{
  unsigned int  i, n, u;
  int           *idx;
  vrna_hc_t     *hc;

  n   = vc->length;
  idx = vc->jindx;
  hc  = vc->hc;

  if (hc->type == VRNA_HC_WINDOW) {
    /* do we actually have any constraints on unpaired positions? */
    if (hc->up_storage) {
      for (hc->up_ext[n + 1] = 0, i = n; i > 0; i--) /* unpaired stretch in exterior loop */
        hc->up_ext[i] = (hc->up_storage[i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) ? 1 +
                        hc->up_ext[i + 1] : 0;

      for (hc->up_hp[n + 1] = 0, i = n; i > 0; i--)  /* unpaired stretch in hairpin loop */
        hc->up_hp[i] = (hc->up_storage[i] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) ? 1 +
                       hc->up_hp[i + 1] : 0;

      for (hc->up_int[n + 1] = 0, i = n; i > 0; i--) /* unpaired stretch in interior loop */
        hc->up_int[i] = (hc->up_storage[i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) ? 1 +
                        hc->up_int[i + 1] : 0;

      for (hc->up_ml[n + 1] = 0, i = n; i > 0; i--)  /* unpaired stretch in multibranch loop */
        hc->up_ml[i] = (hc->up_storage[i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) ? 1 +
                       hc->up_ml[i + 1] : 0;
    } else {
      /* no constraints on unpaired positions */
      for (u = n, i = 1; i <= n; i++, u--)
        hc->up_ext[i] = hc->up_hp[i] = hc->up_int[i] = hc->up_ml[i] = u;
    }
  } else {
    for (hc->up_ext[n + 1] = 0, i = n; i > 0; i--) /* unpaired stretch in exterior loop */
      hc->up_ext[i] = (hc->matrix[idx[i] + i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) ? 1 +
                      hc->up_ext[i + 1] : 0;

    for (hc->up_hp[n + 1] = 0, i = n; i > 0; i--)  /* unpaired stretch in hairpin loop */
      hc->up_hp[i] = (hc->matrix[idx[i] + i] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) ? 1 +
                     hc->up_hp[i + 1] : 0;

    for (hc->up_int[n + 1] = 0, i = n; i > 0; i--) /* unpaired stretch in interior loop */
      hc->up_int[i] = (hc->matrix[idx[i] + i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) ? 1 +
                      hc->up_int[i + 1] : 0;

    for (hc->up_ml[n + 1] = 0, i = n; i > 0; i--)  /* unpaired stretch in multibranch loop */
      hc->up_ml[i] = (hc->matrix[idx[i] + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) ? 1 +
                     hc->up_ml[i + 1] : 0;

    /*
     *  loop arround once more until we find a nucleotide that mustn't
     *  be unpaired (needed for circular folding)
     */

    if (hc->matrix[idx[1] + 1] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
      hc->up_ext[n + 1] = hc->up_ext[1];
      for (i = n; i > 0; i--) {
        if (hc->matrix[idx[i] + i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP)
          hc->up_ext[i] = MIN2(n, 1 + hc->up_ext[i + 1]);
        else
          break;
      }
    }

    if (hc->matrix[idx[1] + 1] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) {
      hc->up_hp[n + 1] = hc->up_hp[1];
      for (i = n; i > 0; i--) {
        if (hc->matrix[idx[i] + i] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP)
          hc->up_hp[i] = MIN2(n, 1 + hc->up_hp[i + 1]);
        else
          break;
      }
    }

    if (hc->matrix[idx[1] + 1] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
      hc->up_int[n + 1] = hc->up_int[1];
      for (i = n; i > 0; i--) {
        if (hc->matrix[idx[i] + i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP)
          hc->up_int[i] = MIN2(n, 1 + hc->up_int[i + 1]);
        else
          break;
      }
    }

    if (hc->matrix[idx[1] + 1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
      hc->up_ml[n + 1] = hc->up_ml[1];
      for (i = n; i > 0; i--) {
        if (hc->matrix[idx[i] + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP)
          hc->up_ml[i] = MIN2(n, 1 + hc->up_ml[i + 1]);
        else
          break;
      }
    }
  }
}


PRIVATE void
hc_update_up_window(vrna_fold_compound_t  *vc,
                    int                   i)
{
  vrna_hc_t *hc;

  hc = vc->hc;

  hc->up_ext[i] = (hc->matrix_local[i][0] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) ?
                  1 + hc->up_ext[i + 1] :
                  0;
  hc->up_hp[i] = (hc->matrix_local[i][0] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) ?
                 1 + hc->up_hp[i + 1] :
                 0;
  hc->up_int[i] = (hc->matrix_local[i][0] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) ?
                  1 + hc->up_int[i + 1] :
                  0;
  hc->up_ml[i] = (hc->matrix_local[i][0] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) ?
                 1 + hc->up_ml[i + 1] :
                 0;
}


#ifdef  VRNA_BACKWARD_COMPAT

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

PUBLIC void
print_tty_constraint_full(void)
{
  vrna_message_constraint_options_all();
}


PUBLIC void
print_tty_constraint(unsigned int option)
{
  vrna_message_constraint_options(option);
}


PUBLIC void
constrain_ptypes(const char   *constraint,
                 unsigned int length,
                 char         *ptype,
                 int          *BP,
                 int          min_loop_size,
                 unsigned int idx_type)
{
  int   n, i, j, k, l;
  int   hx, *stack;
  char  type;
  int   *index;

  if (constraint == NULL)
    return;

  n = (int)strlen(constraint);

  stack = vrna_alloc(sizeof(int) * (n + 1));

  if (!idx_type) {
    /* index allows access in energy matrices at pos (i,j) via index[j]+i */
    index = vrna_idx_col_wise(length);

    for (hx = 0, j = 1; j <= n; j++) {
      switch (constraint[j - 1]) {
        case '|':
          if (BP)
            BP[j] = -1;

          break;
        case 'x':   /* can't pair */
          for (l = 1; l < j - min_loop_size; l++)
            ptype[index[j] + l] = 0;
          for (l = j + min_loop_size + 1; l <= (int)length; l++)
            ptype[index[l] + j] = 0;
          break;
        case '(':
          stack[hx++] = j;
        /* fallthrough */
        case '<':   /* pairs upstream */
          for (l = 1; l < j - min_loop_size; l++)
            ptype[index[j] + l] = 0;
          break;
        case ')':
          if (hx <= 0)
            vrna_message_error("%s\nunbalanced brackets in constraint", constraint);

          i     = stack[--hx];
          type  = ptype[index[j] + i];
          for (k = i + 1; k <= (int)length; k++)
            ptype[index[k] + i] = 0;
          /* don't allow pairs i<k<j<l */
          for (l = j; l <= (int)length; l++)
            for (k = i + 1; k <= j; k++)
              ptype[index[l] + k] = 0;
          /* don't allow pairs k<i<l<j */
          for (l = i; l <= j; l++)
            for (k = 1; k <= i; k++)
              ptype[index[l] + k] = 0;
          for (k = 1; k < j; k++)
            ptype[index[j] + k] = 0;
          ptype[index[j] + i] = (type == 0) ? 7 : type;
        /* fallthrough */
        case '>':   /* pairs downstream */
          for (l = j + min_loop_size + 1; l <= (int)length; l++)
            ptype[index[l] + j] = 0;
          break;
      }
    }
  } else {
    /* index allows access in energy matrices at pos (i,j) via index[i]-j */
    index = vrna_idx_row_wise(length);

    for (hx = 0, j = 1; j <= n; j++) {
      switch (constraint[j - 1]) {
        case 'x':   /* can't pair */
          for (l = 1; l < j - min_loop_size; l++)
            ptype[index[l] - j] = 0;
          for (l = j + min_loop_size + 1; l <= (int)length; l++)
            ptype[index[j] - l] = 0;
          break;
        case '(':
          stack[hx++] = j;
        /* fallthrough */
        case '<':   /* pairs upstream */
          for (l = 1; l < j - min_loop_size; l++)
            ptype[index[l] - j] = 0;
          break;
        case ')':
          if (hx <= 0)
            vrna_message_error("%s\nunbalanced brackets in constraints", constraint);

          i     = stack[--hx];
          type  = ptype[index[i] - j];
          /* don't allow pairs i<k<j<l */
          for (k = i; k <= j; k++)
            for (l = j; l <= (int)length; l++)
              ptype[index[k] - l] = 0;
          /* don't allow pairs k<i<l<j */
          for (k = 1; k <= i; k++)
            for (l = i; l <= j; l++)
              ptype[index[k] - l] = 0;
          ptype[index[i] - j] = (type == 0) ? 7 : type;
        /* fallthrough */
        case '>':   /* pairs downstream */
          for (l = j + min_loop_size + 1; l <= (int)length; l++)
            ptype[index[j] - l] = 0;
          break;
      }
    }
  }

  if (hx != 0)
    vrna_message_error("%s\nunbalanced brackets in constraint string", constraint);

  free(index);
  free(stack);
}


#endif
