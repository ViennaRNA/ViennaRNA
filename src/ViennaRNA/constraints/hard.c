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

#include "ViennaRNA/params/default.h"
#include "ViennaRNA/params/constants.h" /* defines MINPSCORE */
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/alignments.h"
#include "ViennaRNA/io/file_formats.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/constraints/basic.h"
#include "ViennaRNA/constraints/hard.h"


#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#define FULL_HC_MX  1


struct hc_basepair {
  unsigned int  strand;
  unsigned int  start;
  unsigned int  end;
  unsigned char context;
  unsigned char replace;
};


struct vrna_hc_depot_s {
  unsigned char       **up; /* position-wise unpaired-context constraints (0th dimension: strand, 1st: i */
  struct hc_basepair  ***bp;  /* list of base pair constraints (0th dimension: strand, 1st: i, 2nd: entries for each interval [j:k] */
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
PRIVATE unsigned char
default_pair_constraint(vrna_fold_compound_t  *fc,
                        int                   i,
                        int                   j);


PRIVATE void
hc_init_up_storage(vrna_hc_t *hc);


PRIVATE void
hc_init_bp_storage(vrna_hc_t *hc);


PRIVATE void
hc_store_bp(vrna_hc_bp_storage_t  **container,
            int                   i,
            int                   start,
            int                   end,
            unsigned char         loop_type,
            unsigned char         replace);


PRIVATE void
hc_store_bp_override(vrna_hc_bp_storage_t **container,
                     int                  i,
                     int                  start,
                     int                  end,
                     unsigned char        loop_type);


PRIVATE void
hc_store_bp_add(vrna_hc_bp_storage_t  **container,
                int                   i,
                int                   start,
                int                   end,
                unsigned char         loop_type);


PRIVATE INLINE void
apply_stored_bp_hc(unsigned char        *current,
                   vrna_hc_bp_storage_t *container,
                   unsigned int         j);


PRIVATE INLINE void
populate_hc_up(vrna_fold_compound_t *fc,
               unsigned int         i);


PRIVATE INLINE void
populate_hc_bp(vrna_fold_compound_t *fc,
               unsigned int         i,
               unsigned int         maxdist);


PRIVATE void
hc_add_up(vrna_fold_compound_t  *vc,
          int                   i,
          unsigned char         option);


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


PRIVATE void
hc_depot_init(vrna_fold_compound_t *fc);


PRIVATE void
hc_depot_store_up(vrna_fold_compound_t  *fc,
                  unsigned int          i,
                  unsigned int          strand,
                  unsigned char         context,
                  unsigned char         replace);


PRIVATE void
hc_depot_store_bp(vrna_fold_compound_t  *fc,
                  unsigned int          i,
                  unsigned int          strand,
                  unsigned int          start,
                  unsigned int          end,
                  unsigned int          target_strand,
                  unsigned int          context,
                  unsigned char         replace);


PRIVATE void
hc_depot_free(vrna_fold_compound_t *fc);


PRIVATE void
prepare_hc_up(vrna_fold_compound_t *fc,
              unsigned int         options);


PRIVATE void
prepare_hc_bp(vrna_fold_compound_t *fc,
              unsigned int         options);

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
    printf("< : base i is paired downstream with a base i < j\n"
           "> : base i is paired upstream with a base j < i\n");

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
  hc->mx      = (unsigned char *)vrna_alloc(sizeof(unsigned char) * ((n + 1) * (n + 1)));
  hc->up_ext  = (int *)vrna_alloc(sizeof(int) * (n + 2));
  hc->up_hp   = (int *)vrna_alloc(sizeof(int) * (n + 2));
  hc->up_int  = (int *)vrna_alloc(sizeof(int) * (n + 2));
  hc->up_ml   = (int *)vrna_alloc(sizeof(int) * (n + 2));
  hc->depot   = NULL;

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
  unsigned int  n;
  vrna_hc_t     *hc;

  n = vc->length;

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
  hc->depot         = NULL;

  /* set new hard constraints */
  vc->hc = hc;

  /* add null pointers for the generalized hard constraint feature */
  hc->f         = NULL;
  hc->data      = NULL;
  hc->free_data = NULL;
}


PUBLIC void
vrna_hc_update(vrna_fold_compound_t *fc,
               unsigned int         i)
{
  unsigned int  n, maxdist;
  vrna_hc_t     *hc;

  if (fc) {
    n       = fc->length;
    maxdist = fc->window_size;
    hc      = fc->hc;

    if (i > n) {
      vrna_message_warning("vrna_hc_update(): Position %u out of range!",
                           " (Sequence length: %u)",
                           i, n);
    } else {
      /* init up_xx arrays if necessary */
      if (!hc->up_ext) {
        hc->up_ext  = (int *)vrna_alloc(sizeof(int) * (n + 2));
        hc->up_hp   = (int *)vrna_alloc(sizeof(int) * (n + 2));
        hc->up_int  = (int *)vrna_alloc(sizeof(int) * (n + 2));
        hc->up_ml   = (int *)vrna_alloc(sizeof(int) * (n + 2));

        hc_update_up(fc);
      }

      populate_hc_up(fc, i);
      populate_hc_bp(fc, i, maxdist);
    }
  }
}


PUBLIC int
vrna_hc_prepare(vrna_fold_compound_t *fc,
                unsigned int         options)
{
  int ret = 0;

  if (fc) {
    if (options & VRNA_OPTION_WINDOW) {
      /* check for minimal hard constraints structure */
      if ((!fc->hc) || (fc->hc->type != VRNA_HC_WINDOW) || (!fc->hc->matrix_local))
        vrna_hc_init_window(fc);
    } else {
      prepare_hc_up(fc, options);
      prepare_hc_bp(fc, options);
    }

    ret = 1;
  }

  return ret;
}


PUBLIC void
vrna_hc_add_up(vrna_fold_compound_t *fc,
               int                  i,
               unsigned char        option)
{
  unsigned char type;
  unsigned int  actual_i, strand_i;

  if (fc) {
    if (fc->hc) {
      if ((i <= 0) || (i > fc->length)) {
        vrna_message_warning("vrna_hc_add_up: position out of range, not doing anything");
        return;
      }

      strand_i = fc->strand_number[i];
      actual_i = (unsigned int)i - fc->strand_start[strand_i] + 1;

      type = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

      hc_depot_store_up(fc,
                        actual_i,
                        strand_i,
                        type,
                        0);

      hc_add_up(fc, i, option);

      if (fc->hc->type != VRNA_HC_WINDOW)
        hc_update_up(fc);
    }
  }
}


PUBLIC int
vrna_hc_add_up_strand(vrna_fold_compound_t *fc,
                      unsigned int         i,
                      unsigned int         strand,
                      unsigned char        option)
{
  int ret = 0;

  if ((fc) &&
      (fc->hc) &&
      (strand < fc->strands) &&
      (i > 0)) {
    unsigned int n_i = (fc->type == VRNA_FC_TYPE_SINGLE) ?
                        fc->nucleotides[strand].length :
                        fc->alignment[strand].sequences[0].length;

    if (i > n_i)
      return ret;

    hc_depot_store_up(fc,
                      i,
                      strand,
                      option,
                      0);

    ret = 1;
  }

  return ret;
}


PUBLIC int
vrna_hc_add_up_batch(vrna_fold_compound_t *fc,
                     vrna_hc_up_t         *constraints)
{
  unsigned char options;
  unsigned int  strand_i, actual_i, *ss, *sn;
  int           i, ret, pos;

  ret = 0; /* failure */

  if ((fc) &&
      (constraints)) {
    if (fc->hc) {
      sn = fc->strand_number;
      ss = fc->strand_start;

      for (i = 0; constraints[i].position != 0; i++) {
        pos     = constraints[i].position;
        options = constraints[i].options;

        if ((pos <= 0) || (pos > fc->length))
          return ret;

        strand_i = sn[pos];
        actual_i = pos - ss[strand_i] + 1;

        hc_depot_store_up(fc,
                          actual_i,
                          strand_i,
                          options,
                          0);

        hc_add_up(fc, pos, options);
        ret++;
      }

      if (fc->hc->type != VRNA_HC_WINDOW)
        hc_update_up(fc);
    }
  }

  return ret;
}


PUBLIC int
vrna_hc_add_up_strand_batch(vrna_fold_compound_t *fc,
                            vrna_hc_up_t         *constraints)
{
  unsigned char options;
  unsigned int  i, strand, pos, n_pos, *ss, *sn;
  int           ret;

  ret = 0; /* failure */

  if ((fc) &&
      (constraints)) {
    if (fc->hc) {
      sn = fc->strand_number;
      ss = fc->strand_start;

      for (i = 0; constraints[i].position != 0; i++) {
        pos     = (unsigned int)constraints[i].position;
        strand  = (unsigned int)constraints[i].strand;
        options = constraints[i].options;

        if (strand < fc->strands) {
          n_pos = (fc->type == VRNA_FC_TYPE_SINGLE) ?
                  fc->nucleotides[strand].length :
                  fc->alignment[strand].sequences[0].length;

          if (pos > n_pos)
            break;

          hc_depot_store_up(fc,
                            pos,
                            strand,
                            options,
                            0);

          ret++;
        } else {
          break;
        }
      }

      if (fc->hc->type != VRNA_HC_WINDOW)
        hc_update_up(fc);
    }
  }

  return ret;
}


PUBLIC void
vrna_hc_add_bp_nonspecific(vrna_fold_compound_t *vc,
                           int                  i,
                           int                  d,
                           unsigned char        option)
{
  unsigned char type, t1, t2;
  unsigned int  n, strand, actual_i, *sn, *ss, *se;
  int           p;
  vrna_hc_t     *hc;

  if (vc) {

    if (vc->hc) {
      if ((i <= 0) || (i > vc->length)) {
        vrna_message_warning("vrna_hc_add_bp_nonspecific: position out of range, not doing anything");
        return;
      }

      hc        = vc->hc;
      n         = hc->n;
      sn        = vc->strand_number;
      ss        = vc->strand_start;
      se        = vc->strand_end;
      strand    = sn[i];
      actual_i  = i - ss[strand] + 1;

      /* position i may pair in provided contexts */
      type = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
      /* acknowledge pairing direction */
      t1  = (d <= 0) ? type : VRNA_CONSTRAINT_CONTEXT_NONE;
      t2  = (d >= 0) ? type : VRNA_CONSTRAINT_CONTEXT_NONE;

      if (!hc->depot)
        hc_depot_init(vc);

      /* mark this nucleotide to not be unpaired */
      hc_depot_store_up(vc,
                        actual_i,
                        strand,
                        VRNA_CONSTRAINT_CONTEXT_NONE,
                        1);

      /* enforce pairing direction (if any) */
      if (d != 0) {
        /*
            upstream/downstream pairing makes no sense in multi-fold
            context, so we limit the contraint to intra-molecular
            pairing partners
        */
        for (p = 1; p < actual_i; p++)
          hc_depot_store_bp(vc,
                            p,
                            strand,
                            actual_i,
                            actual_i,
                            strand,
                            t1,
                            1);

        hc_depot_store_bp(vc,
                          actual_i,
                          strand,
                          actual_i + 1,
                          vc->nucleotides[strand].length,
                          strand,
                          t2,
                          1);
      }

      if (hc->type == VRNA_HC_WINDOW) {
        /* nucleotide mustn't be unpaired */
        hc_init_up_storage(hc);
        hc->up_storage[i] = VRNA_CONSTRAINT_CONTEXT_NONE;

        /* force pairing direction */
        hc_init_bp_storage(hc);
        for (p = 1; p < i; p++)
          hc_store_bp_add(hc->bp_storage, p, i, i, t1);

        hc_store_bp_add(hc->bp_storage, i, i + 1, n, t2);
      } else {
        if (option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE) {
          /* only allow for possibly non-canonical pairs, do not enforce them */
          for (p = 1; p < i; p++) {
            hc->matrix[vc->jindx[i] + p]  |= t1;
            hc->mx[n * i + p]             |= t1;
            hc->mx[n * p + i]             |= t1;
          }
          for (p = i + 1; p <= vc->length; p++) {
            hc->matrix[vc->jindx[p] + i]  |= t2;
            hc->mx[n * i + p]             |= t2;
            hc->mx[n * p + i]             |= t2;
          }
        } else {
          /* force pairing direction */
          for (p = 1; p < i; p++) {
            hc->matrix[vc->jindx[i] + p]  &= t1;
            hc->mx[n * i + p]             &= t1;
            hc->mx[n * p + i]             &= t1;
          }
          for (p = i + 1; p <= vc->length; p++) {
            hc->matrix[vc->jindx[p] + i]  &= t2;
            hc->mx[n * i + p]             &= t2;
            hc->mx[n * p + i]             &= t2;
          }
          /* nucleotide mustn't be unpaired */
          hc->matrix[vc->jindx[i] + i]  = VRNA_CONSTRAINT_CONTEXT_NONE;
          hc->mx[n * i + i]             = VRNA_CONSTRAINT_CONTEXT_NONE;
        }

        hc_update_up(vc);
      }
    }
  }
}


PUBLIC int
vrna_hc_add_bp_strand(vrna_fold_compound_t *fc,
                      unsigned int         i,
                      unsigned int         strand_i,
                      unsigned int         j,
                      unsigned int         strand_j,
                      unsigned char        option)
{
  unsigned int  n, *sn, *se, *ss, n_i, n_j, k, l;
  int           ret, turn;
  vrna_hc_t     *hc;

  ret = 0;

  if ((fc) &&
      (fc->hc) &&
      (strand_i < fc->strands) &&
      (strand_j < fc->strands) &&
      (i > 0) &&
      (j > 0)) {
    sn  = fc->strand_number;
    ss  = fc->strand_start;
    se  = fc->strand_end;
    hc        = fc->hc;
    n_i       = (fc->type == VRNA_FC_TYPE_SINGLE) ?
                fc->nucleotides[strand_i].length :
                fc->alignment[strand_i].sequences[0].length;
    n_j       = (fc->type == VRNA_FC_TYPE_SINGLE) ?
                fc->nucleotides[strand_j].length :
                fc->alignment[strand_j].sequences[0].length;
    turn      = fc->params->model_details.min_loop_size;

    /* check for strand length ranges */
    if ((i > n_i) || (j > n_j))
      return ret;

    /* check for minimum hairpin loop length for intra-molecular pairs */
    if ((strand_i == strand_j) &&
        ((j - i - 1) < turn))
      return ret;

    if (!hc->depot)
      hc_depot_init(fc);

    /* store the constraint */
    hc_depot_store_bp(fc,
                      i,
                      strand_i,
                      j,
                      j,
                      strand_j,
                      option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS,
                      0); /* add instead of replace */

    /* remove conflicting pairs if necessary */
    if (!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
      /*
          this doesn't make too much sense in multi-fold context,
          so we only remove possible pairs if i and j are on the
          same strand
      */
      if (strand_i == strand_j) {
        for (k = 1; k < i; k++)
          hc_depot_store_bp(fc,
                            k,
                            strand_i,
                            i,
                            j,
                            strand_i,
                            VRNA_CONSTRAINT_CONTEXT_NONE,
                            1);

        hc_depot_store_bp(fc,
                          i,
                          strand_i,
                          i + 1,
                          j - 1,
                          strand_i,
                          VRNA_CONSTRAINT_CONTEXT_NONE,
                          1);

        for (k = i + 1; k < j; k++)
          hc_depot_store_bp(fc,
                            k,
                            strand_i,
                            j,
                            n_i,
                            strand_i,
                            VRNA_CONSTRAINT_CONTEXT_NONE,
                            1);
      }
    }

    /* enforce base pairing if necessary */
    if (option & VRNA_CONSTRAINT_CONTEXT_ENFORCE) {
      hc_depot_store_up(fc,
                        i,
                        strand_i,
                        VRNA_CONSTRAINT_CONTEXT_NONE,
                        1);
      hc_depot_store_up(fc,
                        j,
                        strand_j,
                        VRNA_CONSTRAINT_CONTEXT_NONE,
                        1);
    }

    ret = 1;
  }

  return ret;
}


PUBLIC void
vrna_hc_add_bp(vrna_fold_compound_t *vc,
               int                  i,
               int                  j,
               unsigned char        option)
{
  unsigned int  n, *sn, *se, *ss, strand_i, strand_j, actual_i, actual_j;
  int           k, l;
  vrna_hc_t     *hc;

  if (vc) {
    sn  = vc->strand_number;
    ss  = vc->strand_start;
    se  = vc->strand_end;

    if (vc->hc) {
      if ((i <= 0) || (j <= i) || (j > vc->length)) {
        vrna_message_warning("vrna_hc_add_bp: position out of range, omitting constraint");
        return;
      } else if ((sn[i] == sn[j]) && ((j - i - 1) < vc->params->model_details.min_loop_size)) {
        vrna_message_warning(
          "vrna_hc_add_bp: Pairing partners (%d, %d) violate minimum loop size settings of %dnt, omitting constraint",
          i,
          j,
          vc->params->model_details.min_loop_size);
        return;
      }

      hc        = vc->hc;
      n         = hc->n;

      /*
          determine the corresponding strand numbers and the actual
          position (relative to the strand)
      */
      strand_i  = sn[i];
      strand_j  = sn[j];
      actual_i  = i - ss[strand_i] + 1;
      actual_j  = j - ss[strand_j] + 1;

      if(!vrna_hc_add_bp_strand(vc,
                                actual_i,
                                strand_i,
                                actual_j,
                                strand_j,
                                option))
        return;

      if (hc->type == VRNA_HC_WINDOW) {
        hc_init_bp_storage(hc);
        hc_store_bp_override(hc->bp_storage, i, j, j,
                             option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS);

        if (!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
          /*
           * remove all conflicting base pairs, i.e. do not allow i or j to pair
           * with any other nucleotide k
           */
          for (k = 1; k < i; k++)
            hc_store_bp_add(hc->bp_storage, k, i, j, VRNA_CONSTRAINT_CONTEXT_NONE);             /* (k, i), (k, i + 1), ..., (k, j) with 1 <= k < i */

          hc_store_bp_add(hc->bp_storage, i, i + 1, j - 1, VRNA_CONSTRAINT_CONTEXT_NONE);       /* (i, k), i < k < j */

          for (k = i + 1; k < j; k++)
            hc_store_bp_add(hc->bp_storage, k, j, vc->length, VRNA_CONSTRAINT_CONTEXT_NONE);    /* (i + 1, k), (i + 1, k), ..., (j - 1, k) with (j < k <= n */

          hc_store_bp_add(hc->bp_storage, i, j + 1, vc->length, VRNA_CONSTRAINT_CONTEXT_NONE);  /* (i, k), j < k <= n */
          hc_store_bp_add(hc->bp_storage, j, j + 1, vc->length, VRNA_CONSTRAINT_CONTEXT_NONE);  /* (j, k), j < k <= n */
        }

        if (option & VRNA_CONSTRAINT_CONTEXT_ENFORCE) {
          /* do not allow i,j to be unpaired */
          hc_init_up_storage(hc);
          hc->up_storage[i] = VRNA_CONSTRAINT_CONTEXT_NONE;
          hc->up_storage[j] = VRNA_CONSTRAINT_CONTEXT_NONE;
        }
      } else {
        /* reset ptype in case (i,j) is a non-canonical pair */
        if ((vc->type == VRNA_FC_TYPE_SINGLE) && (option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS)) {
          if (hc->matrix[vc->jindx[j] + i])
            if (vc->ptype[vc->jindx[j] + i] == 0)
              vc->ptype[vc->jindx[j] + i] = 7;

          if (hc->mx[n * i + j])
            if (vc->ptype[vc->jindx[j] + i] == 0)
              vc->ptype[vc->jindx[j] + i] = 7;
        }

        hc->matrix[vc->jindx[j] + i]  = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
        hc->mx[n * i + j]             = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
        hc->mx[n * j + i]             = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

        if (!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
          /*
           * remove all conflicting base pairs, i.e. do not allow i,j to pair
           * with any other nucleotide k
           */
          for (k = 1; k < i; k++) {
            hc->matrix[vc->jindx[i] + k]  = VRNA_CONSTRAINT_CONTEXT_NONE;
            hc->matrix[vc->jindx[j] + k]  = VRNA_CONSTRAINT_CONTEXT_NONE;

            hc->mx[n * i + k] = VRNA_CONSTRAINT_CONTEXT_NONE;
            hc->mx[n * k + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
            hc->mx[n * j + k] = VRNA_CONSTRAINT_CONTEXT_NONE;
            hc->mx[n * k + j] = VRNA_CONSTRAINT_CONTEXT_NONE;

            for (l = i + 1; l < j; l++) {
              hc->matrix[vc->jindx[l] + k] = VRNA_CONSTRAINT_CONTEXT_NONE;

              hc->mx[n * k + l] = VRNA_CONSTRAINT_CONTEXT_NONE;
              hc->mx[n * l + k] = VRNA_CONSTRAINT_CONTEXT_NONE;
            }
          }
          for (k = i + 1; k < j; k++) {
            hc->matrix[vc->jindx[k] + i]  = VRNA_CONSTRAINT_CONTEXT_NONE;
            hc->matrix[vc->jindx[j] + k]  = VRNA_CONSTRAINT_CONTEXT_NONE;

            hc->mx[n * i + k] = VRNA_CONSTRAINT_CONTEXT_NONE;
            hc->mx[n * k + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
            hc->mx[n * j + k] = VRNA_CONSTRAINT_CONTEXT_NONE;
            hc->mx[n * k + j] = VRNA_CONSTRAINT_CONTEXT_NONE;

            for (l = j + 1; l <= vc->length; l++) {
              hc->matrix[vc->jindx[l] + k] = VRNA_CONSTRAINT_CONTEXT_NONE;

              hc->mx[n * k + l] = VRNA_CONSTRAINT_CONTEXT_NONE;
              hc->mx[n * l + k] = VRNA_CONSTRAINT_CONTEXT_NONE;
            }
          }
          for (k = j + 1; k <= vc->length; k++) {
            hc->matrix[vc->jindx[k] + i]  = VRNA_CONSTRAINT_CONTEXT_NONE;
            hc->matrix[vc->jindx[k] + j]  = VRNA_CONSTRAINT_CONTEXT_NONE;

            hc->mx[n * i + k] = VRNA_CONSTRAINT_CONTEXT_NONE;
            hc->mx[n * k + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
            hc->mx[n * j + k] = VRNA_CONSTRAINT_CONTEXT_NONE;
            hc->mx[n * k + j] = VRNA_CONSTRAINT_CONTEXT_NONE;
          }
        }

        if (option & VRNA_CONSTRAINT_CONTEXT_ENFORCE) {
          /* do not allow i,j to be unpaired */
          hc->matrix[vc->jindx[i] + i]  = VRNA_CONSTRAINT_CONTEXT_NONE;
          hc->matrix[vc->jindx[j] + j]  = VRNA_CONSTRAINT_CONTEXT_NONE;

          hc->mx[n * i + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
          hc->mx[n * j + j] = VRNA_CONSTRAINT_CONTEXT_NONE;

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
      free(hc->mx);
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
  int         ret;

  ret = 0; /* Failure */

  if (vc) {
    tmp = NULL;
    if ((!vc->params) && (!vc->exp_params))
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
    if (vc->hc->type != VRNA_HC_WINDOW)
      hc_update_up(vc);

    ret = 1; /* Success */

    free(tmp);
  }

  return ret;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE unsigned char
default_pair_constraint(vrna_fold_compound_t  *fc,
                        int                   i,
                        int                   j)
{
  unsigned char constraint, can_stack;
  short         *S;
  unsigned int  *sn;
  int           type;
  vrna_md_t     *md;

  sn          = fc->strand_number;
  md          = &(fc->params->model_details);
  constraint  = VRNA_CONSTRAINT_CONTEXT_NONE;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      S = fc->sequence_encoding2;
      if ((sn[i] != sn[j]) ||
          (((j - i) < md->max_bp_span) && ((j - i) > md->min_loop_size))) {
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
              constraint  = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
              constraint  &= ~(VRNA_CONSTRAINT_CONTEXT_HP_LOOP | VRNA_CONSTRAINT_CONTEXT_MB_LOOP);
              break;
            }

          /* else fallthrough */
          default:
            constraint = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
            break;
        }

        if (md->noLP) {
          /* check, whether this nucleotide can actually stack with anything or only forms isolated pairs */
          can_stack = VRNA_CONSTRAINT_CONTEXT_NONE;

          /* can it be enclosed by another base pair? */
          if ((i > 1) &&
              (j < fc->length) &&
              (((j - i + 2) < md->max_bp_span) || (sn[i - 1] != sn[j + 1])) &&
              (md->pair[S[i - 1]][S[j + 1]]))
            can_stack = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

          /* can it enclose another base pair? */
          if ((i + 2 < j) &&
              (((j - i - 2) > md->min_loop_size) || (sn[i + 1] != sn[j - 1])) &&
              (md->pair[S[i + 1]][S[j - 1]]))
            can_stack = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

          constraint &= can_stack;
        }
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      if ((sn[i] != sn[j]) ||
          (((j - i + 1) <= md->max_bp_span) && ((j - i - 1) >= md->min_loop_size))) {
        int min_score = md->cv_fact * MINPSCORE;
        int act_score = (fc->hc->type == VRNA_HC_WINDOW) ?
                        fc->pscore_local[i][j - i] :
                        fc->pscore[fc->jindx[j] + i];
        if (act_score >= min_score)
          constraint = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

        if (md->noLP) {
          /* check, whether this nucleotide can actually stack with anything or only forms isolated pairs */
          can_stack = VRNA_CONSTRAINT_CONTEXT_NONE;

          /* can it be enclosed by another base pair? */
          if ((i > 1) &&
              (j < fc->length) &&
              (((j - i + 2) < md->max_bp_span) || (sn[i - 1] != sn[j + 1]))) {
            int outer_pscore = (fc->hc->type == VRNA_HC_WINDOW) ?
                               fc->pscore_local[i - 1][j - i + 2] :
                               fc->pscore[fc->jindx[j + 1] + i - 1];
            if (outer_pscore >= min_score)
              can_stack = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
          }

          /* can it enclose another base pair? */
          if ((i + 2 < j) &&
              (((j - i - 2) > md->min_loop_size) || (sn[i + 1] != sn[j - 1]))) {
            int inner_pscore = (fc->hc->type == VRNA_HC_WINDOW) ?
                               fc->pscore_local[i + 1][j - i - 2] :
                               fc->pscore[fc->jindx[j - 1] + i + 1];
            if (inner_pscore >= min_score)
              can_stack = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
          }

          constraint &= can_stack;
        }
      }

      break;
  }

  return constraint;
}


PRIVATE void
hc_init_up_storage(vrna_hc_t *hc)
{
  unsigned int i;

  if (hc->up_storage == NULL) {
    free(hc->up_storage);
    hc->up_storage = (unsigned char *)vrna_alloc(sizeof(unsigned char) * (hc->n + 2));

    for (i = 1; i <= hc->n; i++)
      /* by default unpaired nucleotides are allowed in all contexts */
      hc->up_storage[i] = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
  }
}


PRIVATE INLINE void
populate_hc_up(vrna_fold_compound_t *fc,
               unsigned int         i)
{
  vrna_hc_t *hc = fc->hc;

  if (hc->type == VRNA_HC_WINDOW) {
    if (hc->up_storage)
      /* We use user-defined constraints for unpaired nucleotides */
      hc->matrix_local[i][0] = hc->up_storage[i];
    else
      /* ... or simply allow unpaired nucleotides in all contexts */
      hc->matrix_local[i][0] = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

    hc_update_up_window(fc, i);
  } else {
    /* do something reasonable here... */
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
hc_depot_init(vrna_fold_compound_t *fc)
{
  vrna_hc_t *hc = fc->hc;
  if (!hc->depot) {
    hc->depot = (vrna_hc_depot_t *)vrna_alloc(sizeof(vrna_hc_depot_t));

    /*
        by default, we only allocate memory for potential constraints for
        each strand. Missing entries, i.e. NULL pointers for any of the nt
        a constraint may be specified for are considered unconstrained, i.e
        default rules apply
    */
    if (fc->strands > 0) {
      hc->depot->up = (unsigned char **)vrna_alloc(sizeof(unsigned char *) * fc->strands);
      hc->depot->bp = (struct hc_basepair ***)vrna_alloc(sizeof(struct hc_basepair **) * fc->strands);
    } else {
      hc->depot->up = NULL;
      hc->depot->bp = NULL;
    }
  }
}


PRIVATE void
hc_depot_store_up(vrna_fold_compound_t  *fc,
                  unsigned int          i,
                  unsigned int          strand,
                  unsigned char         context,
                  unsigned char         replace)
{
  vrna_hc_t *hc = fc->hc;

  hc_depot_init(fc);

  if (!hc->depot->up[strand]) {
    unsigned int length = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->nucleotides[strand].length : fc->alignment[strand].sequences[0].length;
    hc->depot->up[strand] = (unsigned char *)vrna_alloc(sizeof(unsigned char) * (length + 1));
    memset(hc->depot->up[strand] + 1, (int)VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS, sizeof(unsigned char) * length);
  }

  printf("storing hc up %d (strand %d) = %u [%d]\n", i, strand, context, replace);

  if (replace)
    hc->depot->up[strand][i] = context;
  else
    hc->depot->up[strand][i] |= context;
}


PRIVATE void
hc_depot_store_bp(vrna_fold_compound_t  *fc,
                  unsigned int          i,
                  unsigned int          strand,
                  unsigned int          start,
                  unsigned int          end,
                  unsigned int          target_strand,
                  unsigned int          context,
                  unsigned char         replace)
{
  size_t              list_size;
  vrna_hc_t           *hc;
  struct hc_basepair  *list, *ptr, *ptr_last;

  hc_depot_init(fc);

  hc = fc->hc;

  if (!hc->depot->bp[strand]) {
    unsigned int length = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->nucleotides[strand].length : fc->alignment[strand].sequences[0].length;
    hc->depot->bp[strand] = (struct hc_basepair **)vrna_alloc(sizeof(struct hc_basepair *) * (length + 1));
  }

  if (!hc->depot->bp[strand][i]) {
    hc->depot->bp[strand][i] = (struct hc_basepair *)vrna_alloc(sizeof(struct hc_basepair) * 2);
    ptr = hc->depot->bp[strand][i];
  } else {
    /* find out total size of depot */
    for (ptr = hc->depot->bp[strand][i]; ptr->end != 0; ptr++);

    list_size = ptr - hc->depot->bp[strand][i];

    /* increase memory for bp constraints */
    hc->depot->bp[strand][i] = (struct hc_basepair *)vrna_realloc(hc->depot->bp[strand][i],
                                                                  sizeof(struct hc_basepair) * (list_size + 2));

    list = hc->depot->bp[strand][i];

    /* find position where we want to insert the new constraint */
    for (ptr = list; ptr->end != 0; ptr++) {
      if (ptr->strand < target_strand)
        continue;
      else if (ptr->strand > target_strand)
        break;

      /* we only perform the following check for ptr->strand == strand */
      if (ptr->start > start)
        break; /* want to insert before current constraint */

      if (ptr->end < end)
        continue; /* want to insert after current constraint */
    }

    size_t  offset = ptr - list;

    /* shift trailing constraints by 1 entry */
    memmove(list + offset + 1, list + offset,
            sizeof(struct hc_basepair) * (list_size - offset + 1));
  }

  ptr->start    = start;
  ptr->end      = end;
  ptr->strand   = target_strand;
  ptr->context  = context;
  ptr->replace  = replace;
}


PRIVATE void
hc_depot_free(vrna_fold_compound_t *fc)
{
  unsigned int    s, i, n;
  vrna_hc_depot_t *depot = fc->hc->depot;

  if (depot) {
    if (depot->up) {
      for (s = 0; s < fc->strands; s++)
        free(depot->up[s]);

      free(depot->up);
    }

    if (depot->bp) {
      for (s = 0; s < fc->strands; s++) {
        if (depot->bp[s]) {
          n = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->nucleotides[s].length : fc->alignment[s].sequences[0].length;
          for (i = 1; i <= n; i++)
            free(depot->bp[s][i]);

          free(depot->bp[s]);
        }
      }

      free(depot->bp);
    }

    free(depot);
  }
  
  fc->hc->depot = NULL;
}


PRIVATE void
prepare_hc_up(vrna_fold_compound_t *fc,
              unsigned int         options)
{
  unsigned int i, j, n, s, start, end, *sn, *ss, *se;
  int          *idx;
  vrna_hc_t    *hc;

  hc = fc->hc;
  sn = fc->strand_number;
  ss = fc->strand_start;
  se = fc->strand_end;

  if (options & VRNA_OPTION_WINDOW) {
  
  } else {
    n   = fc->length;
    idx = fc->jindx;

    for (i = 1; i <= n; i++) {
      hc->matrix[idx[i] + i] = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
      hc->mx[n * i + i]      = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
    }

    /* 2. apply constraints stored in depot */
    if ((hc->depot) && (hc->depot->up)) {
      for (s = 0; s < fc->strands; s++) {
        if (hc->depot->up[s]) {
          start = ss[s];
          end   = se[s];
          for (j = 1, i = start; i <= end; i++, j++) {
            hc->matrix[idx[i] + i] = hc->depot->up[s][j];
            hc->mx[n * i + i]      = hc->depot->up[s][j];
          }
        }
      }
    }

    hc_update_up(fc);
  }
}


PRIVATE void
prepare_hc_bp(vrna_fold_compound_t *fc,
              unsigned int         options)
{
  unsigned char constraint;
  unsigned int  i, j, k, n, start, end, t_start, t_end, s, *sn, *ss, *se;
  int           *idx, ij;
  vrna_hc_t    *hc;
  struct hc_basepair  *list, *ptr, *ptr_last;

  hc = fc->hc;
  sn = fc->strand_number;
  ss = fc->strand_start;
  se = fc->strand_end;

  if (options & VRNA_OPTION_WINDOW) {
  
  } else {
    n   = fc->length;
    idx = fc->jindx;

    for (j = n; j > 1; j--) {
      ij = idx[j] + 1;
      for (i = 1; i < j; i++, ij++) {
        hc->matrix[ij] = default_pair_constraint(fc, i, j);

        hc->mx[n * i + j] = default_pair_constraint(fc, i, j);
        hc->mx[n * j + i] = hc->mx[n * i + j];
      }
    }

    /* 2. apply constraints stored in depot */
    if ((hc->depot) && (hc->depot->bp)) {
      for (s = 0; s < fc->strands; s++) {
        if (hc->depot->bp[s]) {
          start = ss[s];
          end   = se[s];
          for (j = 1, i = start; i <= end; i++, j++) {
            if (hc->depot->bp[s][j]) {
              /* we actually have some constraints on the jth nucleotide of strand s */
              list = hc->depot->bp[s][j];
              for (ptr = list; ptr->end != 0; ptr++) {
                /* convert local to global positions */
                t_start = ptr->start + ss[ptr->strand] - 1;
                t_end   = ptr->end + ss[ptr->strand] - 1;
                for (k = t_start; k <= t_end; k++) {
                  constraint = hc->mx[n * j + k];
                  if (ptr->replace)
                    constraint = ptr->context;
                  else
                    constraint |= ptr->context;

                  ij = (j > k) ? idx[j] + k : idx[k] + j;
                  hc->matrix[ij] = constraint;

                  hc->mx[n * j + k] = hc->mx[n * k + j] = constraint;
                }
              }
            }
          }
        }
      }
    }
  }
}


PRIVATE void
hc_store_bp_override(vrna_hc_bp_storage_t **container,
                     int                  i,
                     int                  start,
                     int                  end,
                     unsigned char        loop_type)
{
  hc_store_bp(container, i, start, end, loop_type, 1);
}


PRIVATE void
hc_store_bp_add(vrna_hc_bp_storage_t  **container,
                int                   i,
                int                   start,
                int                   end,
                unsigned char         loop_type)
{
  hc_store_bp(container, i, start, end, loop_type, 0);
}


PRIVATE void
hc_store_bp(vrna_hc_bp_storage_t  **container,
            int                   i,
            int                   start,
            int                   end,
            unsigned char         loop_type,
            unsigned char         replace)
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
  container[i][cnt].replace         = replace ? 1 : 0;
}


PRIVATE INLINE void
apply_stored_bp_hc(unsigned char        *current,
                   vrna_hc_bp_storage_t *container,
                   unsigned int         j)
{
  unsigned int  cnt, replace;
  unsigned char constraint = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

  replace = 0;
  /* go through list of constraints for current position i */
  for (cnt = 0; container[cnt].interval_start != 0; cnt++) {
    if (container[cnt].interval_start > j)
      break; /* only constraints for pairs (i,q) with q > j left */

    if (container[cnt].interval_end < j)
      continue; /* constraint for pairs (i,q) with q < j */

    /* constraint has interval [p,q] with p <= j <= q */
    constraint &= container[cnt].loop_type;

    /* is this a replacement or addition constraint? */
    replace = (container[cnt].replace) ? 1 : 0;
  }

  if (replace)
    /* overwrite current constraint */
    *current = constraint;
  else
    /* apply constraint to current (canonical) bp */
    *current &= constraint;
}


PRIVATE INLINE void
populate_hc_bp(vrna_fold_compound_t *fc,
               unsigned int         i,
               unsigned int         maxdist)
{
  unsigned char constraint;
  unsigned int  j, k, n, turn;
  vrna_hc_t     *hc;

  n     = fc->length;
  turn  = fc->params->model_details.min_loop_size;
  hc    = fc->hc;

  /* 2. add default base pairing rules */
  for (k = turn + 1; k < maxdist; k++) {
    j = i + k;
    if (j > n)
      break;

    constraint = default_pair_constraint(fc, i, j);

    /* check whether we have constraints on any pairing partner i or j */
    if ((hc->bp_storage) && (hc->bp_storage[i]))
      apply_stored_bp_hc(&constraint, hc->bp_storage[i], j);

    if (hc->type == VRNA_HC_WINDOW) {
      hc->matrix_local[i][j - i] = constraint;
    } else {
      hc->matrix[fc->jindx[j] + i] = constraint;

      hc->mx[n * i + j] = constraint;
      hc->mx[n * j + i] = constraint;
    }
  }
}


PRIVATE void
hc_add_up(vrna_fold_compound_t  *vc,
          int                   i,
          unsigned char         option)
{
  unsigned char type = VRNA_CONSTRAINT_CONTEXT_NONE;
  unsigned int  n;
  int           j;
  vrna_hc_t     *hc;

  hc  = vc->hc;
  n   = hc->n;

  if (hc->type == VRNA_HC_WINDOW) {
    if (option & VRNA_CONSTRAINT_CONTEXT_ENFORCE) {
      hc_init_up_storage(hc);
      type = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

      hc->up_storage[i] = type;

      if (!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
        hc_init_bp_storage(hc);
        /* add constraints for all pairs (j, i) with j < i */
        for (j = 1; j < i; j++)
          hc_store_bp_add(hc->bp_storage, j, i, i, VRNA_CONSTRAINT_CONTEXT_NONE);

        /* add constraints for all pairs (i, j) with i < j */
        hc_store_bp_add(hc->bp_storage, i, i + 1, n, VRNA_CONSTRAINT_CONTEXT_NONE);
      }
    } else {
      type = ~option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
      hc_init_up_storage(hc);
      hc->up_storage[i] = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

      if (!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
        hc_init_bp_storage(hc);
        /* add constraints for all pairs (j, i) with j < i */
        for (j = 1; j < i; j++)
          hc_store_bp_add(hc->bp_storage, j, i, i, type);

        /* add constraints for all pairs (i, j) with i < j */
        hc_store_bp_add(hc->bp_storage, i, i + 1, n, type);
      }
    }
  } else {
    if (option & VRNA_CONSTRAINT_CONTEXT_ENFORCE) {
      /* force nucleotide to appear unpaired within a certain type of loop */
      /* do not allow i to be paired with any other nucleotide */
      if (!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
        for (j = 1; j < i; j++) {
          hc->matrix[vc->jindx[i] + j] = VRNA_CONSTRAINT_CONTEXT_NONE;

          hc->mx[n * i + j] = VRNA_CONSTRAINT_CONTEXT_NONE;
          hc->mx[n * j + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
        }
        for (j = i + 1; j <= n; j++) {
          hc->matrix[vc->jindx[j] + i] = VRNA_CONSTRAINT_CONTEXT_NONE;

          hc->mx[n * i + j] = VRNA_CONSTRAINT_CONTEXT_NONE;
          hc->mx[n * j + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
        }
      }

      type = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

      hc->matrix[vc->jindx[i] + i] = type;

      hc->mx[n * i + i] = type;
    } else {
      type = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

      /* do not allow i to be paired with any other nucleotide (in context type) */
      if (!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
        for (j = 1; j < i; j++) {
          hc->matrix[vc->jindx[i] + j] &= ~type;

          hc->mx[n * i + j] &= ~type;
          hc->mx[n * j + i] &= ~type;
        }
        for (j = i + 1; j <= n; j++) {
          hc->matrix[vc->jindx[j] + i] &= ~type;

          hc->mx[n * i + j] &= ~type;
          hc->mx[n * j + i] &= ~type;
        }
      }

      hc->matrix[vc->jindx[i] + i] = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

      hc->mx[n * i + i] = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
    }
  }
}


struct hc_bp {
  int           i;
  int           j;
  unsigned char options;
};


PRIVATE void
apply_DB_constraint(vrna_fold_compound_t  *vc,
                    const char            *constraint,
                    unsigned int          options)
{
  char          *sequence;
  short         *S;
  unsigned int  length, min_loop_size, num_up, num_bp, num_bp_unspecific,
                size_up, size_bp, size_bp_unspecific;
  int           n, i, j, hx, *stack, cut;
  vrna_md_t     *md;
  vrna_hc_up_t  *up;
  struct hc_bp  *bp, *bp_unspecific;

  if (constraint == NULL)
    return;

  sequence      = vc->sequence;
  length        = (int)vc->length;
  S             = vc->sequence_encoding2;
  md            = &(vc->params->model_details);
  min_loop_size = md->min_loop_size;
  cut           = vc->cutpoint;
  n             = (int)strlen(constraint);
  stack         = (int *)vrna_alloc(sizeof(int) * (n + 1));
  size_up       = size_bp = size_bp_unspecific = 10;
  num_up        = num_bp = num_bp_unspecific = 0;
  up            = (vrna_hc_up_t *)vrna_alloc(sizeof(vrna_hc_up_t) * size_up);
  bp            = (struct hc_bp *)vrna_alloc(sizeof(struct hc_bp) * size_bp);
  bp_unspecific = (struct hc_bp *)vrna_alloc(sizeof(struct hc_bp) * size_bp_unspecific);


  for (hx = 0, j = 1; j <= n; j++) {
    switch (constraint[j - 1]) {
      /* can't pair */
      case 'x':
        if (options & VRNA_CONSTRAINT_DB_X) {
          up[num_up].position = j;
          up[num_up].options  = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

          num_up++;

          if (num_up == size_up) {
            size_up *= 1.4;
            up      = (vrna_hc_up_t *)vrna_realloc(up, sizeof(vrna_hc_up_t) * size_up);
          }
        }

        break;

      /* must pair, i.e. may not be unpaired */
      case '|':
        if (options & VRNA_CONSTRAINT_DB_PIPE) {
          bp_unspecific[num_bp_unspecific].i        = j;  /* position */
          bp_unspecific[num_bp_unspecific].j        = 0;  /* direction */
          bp_unspecific[num_bp_unspecific].options  = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

          num_bp_unspecific++;

          if (num_bp_unspecific == size_bp_unspecific) {
            size_bp_unspecific  *= 1.4;
            bp_unspecific       = (struct hc_bp *)vrna_realloc(bp_unspecific,
                                                               sizeof(struct hc_bp) *
                                                               size_bp_unspecific);
          }
        }

        break;

      /* weak enforced pair 'open' */
      case '(':
        if (options & VRNA_CONSTRAINT_DB_RND_BRACK)
          stack[hx++] = j;

        break;

      /* weak enforced pair 'close' */
      case ')':
        if (options & VRNA_CONSTRAINT_DB_RND_BRACK) {
          if (hx <= 0) {
            vrna_message_warning("vrna_hc_add_from_db: "
                                 "Unbalanced brackets in constraint string\n%s\n"
                                 "No constraints will be applied!",
                                 constraint);
            goto db_constraints_exit;
          }

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

          if ((j - i - 1) < md->min_loop_size) {
            vrna_message_warning("vrna_hc_add_from_db: "
                                 "Pairing partners (%d, %d) violate minimum loop size settings of %dnt, omitting constraint",
                                 i,
                                 j,
                                 md->min_loop_size);
            break;
          }

          bp[num_bp].i        = i;
          bp[num_bp].j        = j;
          bp[num_bp].options  = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

          if (options & VRNA_CONSTRAINT_DB_ENFORCE_BP)
            bp[num_bp].options |= VRNA_CONSTRAINT_CONTEXT_ENFORCE;

          num_bp++;

          if (num_bp == size_bp) {
            size_bp *= 1.4;
            bp      = (struct hc_bp *)vrna_realloc(bp, sizeof(struct hc_bp) * size_bp);
          }
        }

        break;

      /* pairs downstream */
      case '<':
        if (options & VRNA_CONSTRAINT_DB_ANG_BRACK) {
          bp_unspecific[num_bp_unspecific].i        = j;
          bp_unspecific[num_bp_unspecific].j        = 1;
          bp_unspecific[num_bp_unspecific].options  = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

          num_bp_unspecific++;

          if (num_bp_unspecific == size_bp_unspecific) {
            size_bp_unspecific  *= 1.4;
            bp_unspecific       = (struct hc_bp *)vrna_realloc(bp_unspecific,
                                                               sizeof(struct hc_bp) *
                                                               size_bp_unspecific);
          }

          if (!(options & VRNA_CONSTRAINT_DB_ENFORCE_BP)) {
            /* (re-)allow this nucleotide to stay unpaired for nostalgic reasons */
            up[num_up].position = j;
            up[num_up].options  = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS |
                                  VRNA_CONSTRAINT_CONTEXT_NO_REMOVE;

            num_up++;

            if (num_up == size_up) {
              size_up *= 1.4;
              up      = (vrna_hc_up_t *)vrna_realloc(up, sizeof(vrna_hc_up_t) * size_up);
            }
          }
        }

        break;

      /* pairs upstream */
      case '>':
        if (options & VRNA_CONSTRAINT_DB_ANG_BRACK) {
          bp_unspecific[num_bp_unspecific].i        = j;
          bp_unspecific[num_bp_unspecific].j        = -1;
          bp_unspecific[num_bp_unspecific].options  = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

          num_bp_unspecific++;

          if (num_bp_unspecific == size_bp_unspecific) {
            size_bp_unspecific  *= 1.4;
            bp_unspecific       = (struct hc_bp *)vrna_realloc(bp_unspecific,
                                                               sizeof(struct hc_bp) *
                                                               size_bp_unspecific);
          }

          if (!(options & VRNA_CONSTRAINT_DB_ENFORCE_BP)) {
            /* (re-)allow this nucleotide to stay unpaired for nostalgic reasons */
            up[num_up].position = j;
            up[num_up].options  = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS |
                                  VRNA_CONSTRAINT_CONTEXT_NO_REMOVE;

            num_up++;

            if (num_up == size_up) {
              size_up *= 1.4;
              up      = (vrna_hc_up_t *)vrna_realloc(up, sizeof(vrna_hc_up_t) * size_up);
            }
          }
        }

        break;

      /* only intramolecular basepairing */
      case 'l':
        if (options & VRNA_CONSTRAINT_DB_INTRAMOL) {
          unsigned int l;
          if (cut > 1) {
            if (j < cut) {
              for (l = MAX2(j + min_loop_size, cut); l <= length; l++) {
                bp[num_bp].i        = j;
                bp[num_bp].j        = l;
                bp[num_bp].options  = VRNA_CONSTRAINT_CONTEXT_NONE |
                                      VRNA_CONSTRAINT_CONTEXT_NO_REMOVE;

                num_bp++;

                if (num_bp == size_bp) {
                  size_bp *= 1.4;
                  bp      = (struct hc_bp *)vrna_realloc(bp, sizeof(struct hc_bp) * size_bp);
                }
              }
            } else {
              for (l = 1; l < MIN2(cut, j - min_loop_size); l++) {
                bp[num_bp].i        = l;
                bp[num_bp].j        = j;
                bp[num_bp].options  = VRNA_CONSTRAINT_CONTEXT_NONE |
                                      VRNA_CONSTRAINT_CONTEXT_NO_REMOVE;

                num_bp++;

                if (num_bp == size_bp) {
                  size_bp *= 1.4;
                  bp      = (struct hc_bp *)vrna_realloc(bp, sizeof(struct hc_bp) * size_bp);
                }
              }
            }
          }
        }

        break;

      /* only intermolecular bp */
      case 'e':
        if (options & VRNA_CONSTRAINT_DB_INTERMOL) {
          unsigned int l;
          if (cut > 1) {
            if (j < cut) {
              for (l = 1; l < j; l++) {
                bp[num_bp].i        = l;
                bp[num_bp].j        = j;
                bp[num_bp].options  = VRNA_CONSTRAINT_CONTEXT_NONE |
                                      VRNA_CONSTRAINT_CONTEXT_NO_REMOVE;

                num_bp++;

                if (num_bp == size_bp) {
                  size_bp *= 1.4;
                  bp      = (struct hc_bp *)vrna_realloc(bp, sizeof(struct hc_bp) * size_bp);
                }
              }

              for (l = j + 1; l < cut; l++) {
                bp[num_bp].i        = j;
                bp[num_bp].j        = l;
                bp[num_bp].options  = VRNA_CONSTRAINT_CONTEXT_NONE |
                                      VRNA_CONSTRAINT_CONTEXT_NO_REMOVE;

                num_bp++;

                if (num_bp == size_bp) {
                  size_bp *= 1.4;
                  bp      = (struct hc_bp *)vrna_realloc(bp, sizeof(struct hc_bp) * size_bp);
                }
              }
            } else {
              for (l = cut; l < j; l++) {
                bp[num_bp].i        = l;
                bp[num_bp].j        = j;
                bp[num_bp].options  = VRNA_CONSTRAINT_CONTEXT_NONE |
                                      VRNA_CONSTRAINT_CONTEXT_NO_REMOVE;

                num_bp++;

                if (num_bp == size_bp) {
                  size_bp *= 1.4;
                  bp      = (struct hc_bp *)vrna_realloc(bp, sizeof(struct hc_bp) * size_bp);
                }
              }

              for (l = j + 1; l <= length; l++) {
                bp[num_bp].i        = j;
                bp[num_bp].j        = l;
                bp[num_bp].options  = VRNA_CONSTRAINT_CONTEXT_NONE |
                                      VRNA_CONSTRAINT_CONTEXT_NO_REMOVE;

                num_bp++;

                if (num_bp == size_bp) {
                  size_bp *= 1.4;
                  bp      = (struct hc_bp *)vrna_realloc(bp, sizeof(struct hc_bp) * size_bp);
                }
              }
            }
          }
        }

        break;

      case '.':
        break;

      default:
        vrna_message_warning(
          "vrna_hc_add_from_db: "
          "Unrecognized character '%c' in constraint string",
          constraint[j - 1]);
        break;
    }
  }

  if (hx != 0) {
    vrna_message_warning("vrna_hc_add_from_db: "
                         "Unbalanced brackets in constraint string\n%s\n"
                         "No constraints will be applied!",
                         constraint);
    goto db_constraints_exit;
  }

  /* finally, apply constraints */

  /* 1st, unspecific pairing states */
  for (i = 0; i < num_bp_unspecific; i++)
    vrna_hc_add_bp_nonspecific(vc,
                               bp_unspecific[i].i,  /* nucleotide position */
                               bp_unspecific[i].j,  /* pairing direction */
                               bp_unspecific[i].options);

  /* 2nd, specific base pairs */
  for (i = 0; i < num_bp; i++)
    vrna_hc_add_bp(vc,
                   bp[i].i,
                   bp[i].j,
                   bp[i].options);

  /* 3rd, unpaired constraints */
  if (num_up > 0) {
    up[num_up].position = 0;  /* end of list marker */
    vrna_hc_add_up_batch(vc, up);
  }

db_constraints_exit:

  /* clean up */
  free(up);
  free(bp);
  free(bp_unspecific);
  free(stack);
}


PRIVATE void
hc_reset_to_default(vrna_fold_compound_t *vc)
{
  unsigned int  i, j, ij, n, s, strands;
  int           *idx;
  vrna_hc_t     *hc;

  n   = vc->length;
  hc  = vc->hc;
  idx = vc->jindx;

  /* ######################### */
  /* fill with default values  */
  /* ######################### */

  /* 1. unpaired nucleotides are allowed in all contexts */
  for (i = 1; i <= n; i++) {
    hc->matrix[idx[i] + i] = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

    hc->mx[n * i + i] = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
  }

  /* 2. base pairs follow default rules, i.e. canonical pairs, maybe without isolated pairs (if noLP) */
  for (j = n; j > 1; j--) {
    ij = idx[j] + 1;
    for (i = 1; i < j; i++, ij++) {
      hc->matrix[ij] = default_pair_constraint(vc, i, j);

      hc->mx[n * i + j] = default_pair_constraint(vc, i, j);
      hc->mx[n * j + i] = hc->mx[n * i + j];
    }
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
  vrna_hc_t     *hc;

  n   = vc->length;
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
      hc->up_ext[i] = (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) ? 1 +
                      hc->up_ext[i + 1] : 0;

    for (hc->up_hp[n + 1] = 0, i = n; i > 0; i--)  /* unpaired stretch in hairpin loop */
      hc->up_hp[i] = (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) ? 1 +
                     hc->up_hp[i + 1] : 0;

    for (hc->up_int[n + 1] = 0, i = n; i > 0; i--) /* unpaired stretch in interior loop */
      hc->up_int[i] = (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) ? 1 +
                      hc->up_int[i + 1] : 0;

    for (hc->up_ml[n + 1] = 0, i = n; i > 0; i--)  /* unpaired stretch in multibranch loop */
      hc->up_ml[i] = (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) ? 1 +
                     hc->up_ml[i + 1] : 0;

    /*
     *  loop arround once more until we find a nucleotide that mustn't
     *  be unpaired (needed for circular folding)
     */

    if (hc->mx[n + 1] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
      hc->up_ext[n + 1] = hc->up_ext[1];
      for (i = n; i > 0; i--) {
        if (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP)
          hc->up_ext[i] = MIN2(n, 1 + hc->up_ext[i + 1]);
        else
          break;
      }
    }

    if (hc->mx[n + 1] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) {
      hc->up_hp[n + 1] = hc->up_hp[1];
      for (i = n; i > 0; i--) {
        if (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP)
          hc->up_hp[i] = MIN2(n, 1 + hc->up_hp[i + 1]);
        else
          break;
      }
    }

    if (hc->mx[n + 1] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
      hc->up_int[n + 1] = hc->up_int[1];
      for (i = n; i > 0; i--) {
        if (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP)
          hc->up_int[i] = MIN2(n, 1 + hc->up_int[i + 1]);
        else
          break;
      }
    }

    if (hc->mx[n + 1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
      hc->up_ml[n + 1] = hc->up_ml[1];
      for (i = n; i > 0; i--) {
        if (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP)
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


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

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
