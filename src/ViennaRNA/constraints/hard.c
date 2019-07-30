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

#define STATE_CLEAN         (unsigned char)0
#define STATE_DIRTY_UP      (unsigned char)1
#define STATE_DIRTY_BP      (unsigned char)2
#define STATE_UNINITIALIZED (unsigned char)4

#include "hc_depot.inc"

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
               unsigned int         maxdist,
               unsigned int         options);


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
default_hc_up(vrna_fold_compound_t *fc,
              unsigned int         options);

PRIVATE void
default_hc_bp(vrna_fold_compound_t *fc,
              unsigned int         options);

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
  hc->state   = STATE_UNINITIALIZED;

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
  hc->state         = STATE_UNINITIALIZED;

  /* set new hard constraints */
  vc->hc = hc;

  /* add null pointers for the generalized hard constraint feature */
  hc->f         = NULL;
  hc->data      = NULL;
  hc->free_data = NULL;
}


PUBLIC void
vrna_hc_update(vrna_fold_compound_t *fc,
               unsigned int         i,
               unsigned int         options)
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
      populate_hc_bp(fc, i, maxdist, options);
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
      if (fc->hc->state & STATE_UNINITIALIZED) {
        default_hc_up(fc, options);
        default_hc_bp(fc, options);
      }

      if (fc->hc->state & STATE_DIRTY_UP)
        prepare_hc_up(fc, options);

      if (fc->hc->state & STATE_DIRTY_BP)
        prepare_hc_bp(fc, options);

      if (fc->hc->state & ~STATE_CLEAN)
        hc_update_up(fc);
    }

    fc->hc->state = STATE_CLEAN;
    ret = 1;
  }

  return ret;
}


PUBLIC void
vrna_hc_add_up(vrna_fold_compound_t *fc,
               int                  i,
               unsigned char        option)
{
  unsigned int  actual_i, strand_i;

  if (fc) {
    if (fc->hc) {
      if ((i <= 0) || (i > fc->length)) {
        vrna_message_warning("vrna_hc_add_up: position out of range, not doing anything");
        return;
      }

      strand_i = fc->strand_number[i];
      actual_i = (unsigned int)i - fc->strand_start[strand_i] + 1;

      hc_depot_store_up(fc,
                        actual_i,
                        strand_i,
                        option);

      fc->hc->state |= STATE_DIRTY_UP;
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
                      option);

    fc->hc->state |= STATE_DIRTY_UP;

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
          break;

        strand_i = sn[pos];
        actual_i = pos - ss[strand_i] + 1;

        hc_depot_store_up(fc,
                          actual_i,
                          strand_i,
                          options);

        ret++;
      }
    }
  }

  if (ret)
    fc->hc->state |= STATE_DIRTY_UP;

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
                            options);

          ret++;
        } else {
          break;
        }
      }
    }
  }

  if (ret)
    fc->hc->state |= STATE_DIRTY_UP;

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

      hc_depot_store_nonspec(vc,
                             actual_i,
                             strand,
                             d,
                             option);

      hc->state |= STATE_DIRTY_UP;

      if (hc->type == VRNA_HC_WINDOW) {
        /* nucleotide mustn't be unpaired */
        hc_init_up_storage(hc);
        hc->up_storage[i] = VRNA_CONSTRAINT_CONTEXT_NONE;

        /* force pairing direction */
        hc_init_bp_storage(hc);
        for (p = 1; p < i; p++)
          hc_store_bp_add(hc->bp_storage, p, i, i, t1);

        hc_store_bp_add(hc->bp_storage, i, i + 1, n, t2);
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

    /* store the constraint */
    hc_depot_store_bp(fc,
                      i,
                      strand_i,
                      j,
                      strand_j,
                      option);

    hc->state |= STATE_DIRTY_BP;

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

      hc->state |= STATE_DIRTY_BP;

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

    hc_depot_free(hc);

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
      if (((j - i) < md->max_bp_span) &&
          ((sn[i] != sn[j]) || ((j - i) > md->min_loop_size))) {
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
  unsigned char context;
  unsigned int actual_i, strand;
  vrna_hc_t *hc = fc->hc;

  if (hc->type == VRNA_HC_WINDOW) {
    strand   = fc->strand_number[i];
    actual_i = fc->strand_start[strand] + i - 1;

    if ((hc->depot) &&
        (hc->depot->up) &&
        (hc->depot->up_size[strand] >= i)) {
      /* We use user-defined constraints for unpaired nucleotides */
      context = hc->depot->up[strand][i].context;

      if (hc->depot->up[strand][i].nonspec) {
        /* this nucleotide must pair */
        hc->matrix_local[actual_i][0] = VRNA_CONSTRAINT_CONTEXT_NONE;

        /* pairing direction will be enforced by populate_hc_bp() */
      } else if (context & VRNA_CONSTRAINT_CONTEXT_ENFORCE) {
        hc->matrix_local[i][0] = context & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
        /* removal of base pairs involving i will be handled by populate_hc_bp() */
      } else {
        hc->matrix_local[i][0] = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
        /* restriction of base pairing contexts involving i will be handled by populate_hc_bp() */
      }
    } else if (hc->up_storage) {
      /* We use user-defined constraints for unpaired nucleotides */
      hc->matrix_local[i][0] = hc->up_storage[i];
    } else {
      /* ... or simply allow unpaired nucleotides in all contexts */
      hc->matrix_local[i][0] = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
    }

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
default_hc_up(vrna_fold_compound_t *fc,
              unsigned int         options)
{
  unsigned int    i, n;
  int             *idx;
  vrna_hc_t       *hc;
  vrna_hc_depot_t *depot;

  hc = fc->hc;

  if (options & VRNA_OPTION_WINDOW) {
  
  } else {
    n   = fc->length;
    idx = fc->jindx;

    for (i = 1; i <= n; i++) {
      hc->matrix[idx[i] + i] = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
      hc->mx[n * i + i]      = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
    }
  }
}

PRIVATE void
prepare_hc_up(vrna_fold_compound_t *fc,
              unsigned int         options)
{
  unsigned char   option, type, t1, t2;
  unsigned int    i, j, k, n, s, *ss;
  int             *idx;
  vrna_hc_t       *hc;
  vrna_hc_depot_t *depot;

  hc = fc->hc;

  if (options & VRNA_OPTION_WINDOW) {
  
  } else {
    n     = fc->length;
    idx   = fc->jindx;
    ss    = fc->strand_start;
    depot = hc->depot;

    /* 2. apply constraints as stored in depot */
    if ((depot) && (depot->up)) {
      for (s = 0; s < depot->strands; s++) {
        for (k = 1; k <= depot->up_size[s]; k++) {
          /* process nucleotide-specific constraint */
          option = depot->up[s][k].context;
          i      = ss[s] + k - 1; /* constraint position in current strand order */

          if (depot->up[s][k].nonspec) {
            /* this is actually a must-pair constraint */
            type = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
            /* acknowledge pairing direction */
            t1  = (depot->up[s][k].direction <= 0) ? type : VRNA_CONSTRAINT_CONTEXT_NONE;
            t2  = (depot->up[s][k].direction >= 0) ? type : VRNA_CONSTRAINT_CONTEXT_NONE;

            if (option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE) {
              /* only allow for possibly non-canonical pairs, do not enforce them */
              for (j = 1; j < i; j++) {
                hc->matrix[idx[i] + j]  |= t1;
                hc->mx[n * i + j]       |= t1;
                hc->mx[n * j + i]       |= t1;
              }
              for (j = i + 1; j <= n; j++) {
                hc->matrix[idx[j] + i]  |= t2;
                hc->mx[n * i + j]       |= t2;
                hc->mx[n * j + i]       |= t2;
              }
            } else {
              /* force pairing direction */
              for (j = 1; j < i; j++) {
                hc->matrix[idx[i] + j]  &= t1;
                hc->mx[n * i + j]       &= t1;
                hc->mx[n * j + i]       &= t1;
              }
              for (j = i + 1; j <= n; j++) {
                hc->matrix[idx[j] + i]  &= t2;
                hc->mx[n * i + j]       &= t2;
                hc->mx[n * j + i]       &= t2;
              }
              /* nucleotide mustn't be unpaired */
              hc->matrix[idx[i] + i]  = VRNA_CONSTRAINT_CONTEXT_NONE;
              hc->mx[n * i + i]       = VRNA_CONSTRAINT_CONTEXT_NONE;
            }
          } else {
            /* 'regular' nucleotide-specific constraint */
            if (option & VRNA_CONSTRAINT_CONTEXT_ENFORCE) {
              /* force nucleotide to appear unpaired within a certain type of loop */
              /* do not allow i to be paired with any other nucleotide */
              if (!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
                for (j = 1; j < i; j++) {
                  hc->matrix[idx[i] + j] = VRNA_CONSTRAINT_CONTEXT_NONE;

                  hc->mx[n * i + j] = VRNA_CONSTRAINT_CONTEXT_NONE;
                  hc->mx[n * j + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
                }
                for (j = i + 1; j <= n; j++) {
                  hc->matrix[idx[j] + i] = VRNA_CONSTRAINT_CONTEXT_NONE;

                  hc->mx[n * i + j] = VRNA_CONSTRAINT_CONTEXT_NONE;
                  hc->mx[n * j + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
                }
              }

              type = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

              hc->matrix[idx[i] + i] = type;

              hc->mx[n * i + i] = type;
            } else {
              type = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

              /* do not allow i to be paired with any other nucleotide (in context type) */
              if (!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
                for (j = 1; j < i; j++) {
                  hc->matrix[idx[i] + j] &= ~type;

                  hc->mx[n * i + j] &= ~type;
                  hc->mx[n * j + i] &= ~type;
                }
                for (j = i + 1; j <= n; j++) {
                  hc->matrix[idx[j] + i] &= ~type;

                  hc->mx[n * i + j] &= ~type;
                  hc->mx[n * j + i] &= ~type;
                }
              }

              hc->matrix[idx[i] + i]  = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
              hc->mx[n * i + i]       = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
            }
          }
        }
      }
    }
  }
}


PRIVATE void
default_hc_bp(vrna_fold_compound_t *fc,
              unsigned int         options)
{
  unsigned int  i, j, n;
  int           *idx, ij;
  vrna_hc_t    *hc;

  hc = fc->hc;

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
  }
}


PRIVATE void
prepare_hc_bp(vrna_fold_compound_t *fc,
              unsigned int         options)
{
  unsigned char   option, type;
  unsigned int    i, j, k, p, q, n, strand_j, start, end, t_start, t_end, s, *sn, *ss, *se;
  int             *idx, ij;
  vrna_hc_t       *hc;
  vrna_hc_depot_t *depot;

  hc    = fc->hc;
  depot = hc->depot;
  sn = fc->strand_number;
  ss = fc->strand_start;
  se = fc->strand_end;

  if ((!depot) || (!depot->bp))
    return;

  if (options & VRNA_OPTION_WINDOW) {
  
  } else {
    n   = fc->length;
    idx = fc->jindx;

    /* 2. apply constraints stored in depot */
    for (s = 0; s < depot->strands; s++) {
      for (k = 0; k < depot->bp_num[s]; k++) {
        /* process base pair specific constraint */
        option    = depot->bp[s][k].context;
        i         = ss[s] + depot->bp[s][k].i - 1;  /* constraint position in current strand order */
        strand_j  = depot->bp[s][k].strand_j;
        j         = ss[strand_j] + depot->bp[s][k].j - 1;  /* constraint position in current strand order */

        /* actually apply the constraint */
        hc->matrix[idx[j] + i]  = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
        hc->mx[n * i + j]       = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
        hc->mx[n * j + i]       = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

        /* is the ptype reset actually required??? */
        if (option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS) {
          /* reset ptype in case (i,j) is a non-canonical pair */
          if (fc->ptype[idx[j] + i] == 0)
            fc->ptype[idx[j] + i] = 7;
        }

        if (!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
          /*
           * remove all conflicting base pairs, i.e. do not allow i,j to pair
           * with any other nucleotide k
           */
          for (p = 1; p < i; p++) {
            hc->matrix[idx[i] + p]  = VRNA_CONSTRAINT_CONTEXT_NONE;
            hc->matrix[idx[j] + p]  = VRNA_CONSTRAINT_CONTEXT_NONE;

            hc->mx[n * i + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
            hc->mx[n * p + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
            hc->mx[n * j + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
            hc->mx[n * p + j] = VRNA_CONSTRAINT_CONTEXT_NONE;

            for (q = i + 1; q < j; q++) {
              hc->matrix[idx[q] + p] = VRNA_CONSTRAINT_CONTEXT_NONE;

              hc->mx[n * p + q] = VRNA_CONSTRAINT_CONTEXT_NONE;
              hc->mx[n * q + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
            }
          }
          for (p = i + 1; p < j; p++) {
            hc->matrix[idx[p] + i]  = VRNA_CONSTRAINT_CONTEXT_NONE;
            hc->matrix[idx[j] + p]  = VRNA_CONSTRAINT_CONTEXT_NONE;

            hc->mx[n * i + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
            hc->mx[n * p + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
            hc->mx[n * j + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
            hc->mx[n * p + j] = VRNA_CONSTRAINT_CONTEXT_NONE;

            for (q = j + 1; q <= n; q++) {
              hc->matrix[idx[q] + p] = VRNA_CONSTRAINT_CONTEXT_NONE;

              hc->mx[n * p + q] = VRNA_CONSTRAINT_CONTEXT_NONE;
              hc->mx[n * q + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
            }
          }
          for (p = j + 1; p <= n; p++) {
            hc->matrix[idx[p] + i]  = VRNA_CONSTRAINT_CONTEXT_NONE;
            hc->matrix[idx[p] + j]  = VRNA_CONSTRAINT_CONTEXT_NONE;

            hc->mx[n * i + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
            hc->mx[n * p + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
            hc->mx[n * j + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
            hc->mx[n * p + j] = VRNA_CONSTRAINT_CONTEXT_NONE;
          }
        }

        if (option & VRNA_CONSTRAINT_CONTEXT_ENFORCE) {
          /* do not allow i,j to be unpaired */
          hc->matrix[idx[i] + i]  = VRNA_CONSTRAINT_CONTEXT_NONE;
          hc->matrix[idx[j] + j]  = VRNA_CONSTRAINT_CONTEXT_NONE;

          hc->mx[n * i + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
          hc->mx[n * j + j] = VRNA_CONSTRAINT_CONTEXT_NONE;
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
               unsigned int         maxdist,
               unsigned int         options)
{
  unsigned char constraint;
  unsigned int  j, k, p, n, turn, strand, sj, actual_i, actual_j, *sn, *ss;
  vrna_hc_t     *hc;

  n         = fc->length;
  sn        = fc->strand_number;
  ss        = fc->strand_start;
  strand    = sn[i];
  actual_i  = i - ss[strand] + 1;
  turn      = fc->params->model_details.min_loop_size;
  hc        = fc->hc;

  if (options & VRNA_CONSTRAINT_WINDOW_UPDATE_3) {
    /* the sliding window moves from 3' to 5' side */

    /* apply default constraints first */
    for (k = turn + 1; k < maxdist; k++) {
      j = i + k;
      if (j > n)
        break;

      constraint = default_pair_constraint(fc, i, j);

      hc->matrix_local[i][j - i] = constraint;
    }

    if (hc->depot) {
      /* apply remainder of (partly) applied nucleotide-specific constraints */
      if (hc->depot->up) {

        /*
            apply remainder of (partly) applied nucleotide-specific
            constraints for i
        */
        if ((hc->depot->up[strand]) &&
            (hc->depot->up_size[strand] >= actual_i)) {
          constraint = hc->depot->up[strand][actual_i].context;

          if (hc->depot->up[strand][actual_i].nonspec != 0) {
            /* i must be paired, check preferred direction */
            if (hc->depot->up[sj][actual_j].direction < 0) {
              /* i is only allowed to pair upstream */
              /* remove all base pairs (i, j) with i < j < i + maxist */
              for (p = i + 1; p < MIN2(i + maxdist, n + 1); p++)
                hc->matrix_local[i][p - i] = VRNA_CONSTRAINT_CONTEXT_NONE;
            } else if (hc->depot->up[sj][actual_j].direction > 0) {
              /* nothing to do */
            }
          } else if (!(constraint & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
            if (constraint & VRNA_CONSTRAINT_CONTEXT_ENFORCE) {
              /* remove all base pairs (i, j) with i < j < i + maxdist */
              for (p = i + 1; p < MIN2(i + maxdist, n + 1); p++)
                hc->matrix_local[i][p - i] = VRNA_CONSTRAINT_CONTEXT_NONE;
            } else {
              constraint = ~constraint & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
              for (p = i + 1; p < MIN2(i + maxdist, n + 1); p++)
                hc->matrix_local[i][p - i] = constraint;
            }
          }
        }

        /*
            apply remainder of (partly) applied nucleotide-specific
            constraints for j with i < j < i + maxdist
        */
        for (k = 1; k < maxdist; k++) {
          j         = i + k;
          sj        = sn[j];
          actual_j  = j - ss[sj] + 1;

          if (j > n)
            break;

          if ((hc->depot->up[sj]) &&
              (hc->depot->up_size[sj] >= actual_j)) {
            constraint = hc->depot->up[sj][actual_j].context;

            if (hc->depot->up[sj][actual_j].nonspec != 0) {
              /* j must be paired, check preferred direction */
              if (hc->depot->up[sj][actual_j].direction < 0) {
                /* nothing to do */
              } else if (hc->depot->up[sj][actual_j].direction > 0) {
                /* remove (i, j) base pair */
                hc->matrix_local[i][j - i] = VRNA_CONSTRAINT_CONTEXT_NONE;
              }
            } else if (!(constraint & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
              if (constraint & VRNA_CONSTRAINT_CONTEXT_ENFORCE) {
                /* remove (i, j) base pair */
                hc->matrix_local[i][j - i] = VRNA_CONSTRAINT_CONTEXT_NONE;
              } else {
                constraint = ~constraint & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
                hc->matrix_local[i][j - i] = constraint;
              }
            }
          }
        }
      }

      /* now for the actual base pair constraints */
      if (hc->depot->bp) {

      }
    }

    for (k = turn + 1; k < maxdist; k++) {
      j = i + k;
      if (j > n)
        break;

      /* check whether we have constraints on any pairing partner i or j */
      if ((hc->bp_storage) && (hc->bp_storage[i]))
        apply_stored_bp_hc(&constraint, hc->bp_storage[i], j);
    }
  } else if (options & VRNA_CONSTRAINT_WINDOW_UPDATE_5) {
    /* the sliding window moves from 5' to 3' side (i is 3' nucleotide) */

    /* apply default constraints first */
    unsigned int j_start = 1;
    unsigned int j_stop  = i;

    if (i > maxdist)
      j_start = i - maxdist + 1;

    if (i > turn + 1)
      j_stop = i - turn - 1;

    for (j = j_start; j <= j_stop; j++)
      hc->matrix_local[j][i - j] = default_pair_constraint(fc, j, i);

    /* next are user-defined constraints */
    if (hc->depot) {
      /* apply remainder of (partly) applied nucleotide-specific constraints */
      if (hc->depot->up) {

        /*
            apply remainder of (partly) applied nucleotide-specific
            constraints for i
        */
        if ((hc->depot->up[strand]) &&
            (hc->depot->up_size[strand] >= actual_i)) {
          constraint = hc->depot->up[strand][actual_i].context;

          if (hc->depot->up[strand][actual_i].nonspec != 0) {
            /* i must be paired, check preferred direction */
            if (hc->depot->up[sj][actual_j].direction < 0) {
              /* nothing to do */
            } else if (hc->depot->up[sj][actual_j].direction > 0) {
              /* i is only allowed to pair downstream */
              /* remove all base pairs (j, i) with i - maxdist < j < i */
              for (j = j_start; j <= j_stop; j++)
                hc->matrix_local[j][i - j] = VRNA_CONSTRAINT_CONTEXT_NONE;
            }
          } else if (!(constraint & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
            if (constraint & VRNA_CONSTRAINT_CONTEXT_ENFORCE) {
              /* remove all base pairs (j, i) with i - maxdist < j < i */
              for (j = j_start; j <= j_stop; j++)
                hc->matrix_local[j][i - j] = VRNA_CONSTRAINT_CONTEXT_NONE;
            } else {
              constraint = ~constraint & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
              for (j = j_start; j <= j_stop; j++)
                hc->matrix_local[j][i - j] = constraint;
            }
          }
        }

        /*
            apply remainder of (partly) applied nucleotide-specific
            constraints for j with i - maxdist < j < i
        */
        for (j = j_start; j <= j_stop; j++) {
          sj        = sn[j];
          actual_j  = j - ss[sj] + 1;

          if ((hc->depot->up[sj]) &&
              (hc->depot->up_size[sj] >= actual_j)) {
            constraint = hc->depot->up[sj][actual_j].context;

            if (hc->depot->up[sj][actual_j].nonspec != 0) {
              /* j must be paired, check preferred direction */
              if (hc->depot->up[sj][actual_j].direction < 0) {
                /* remove (j, i) base pair */
                hc->matrix_local[j][i - j] = VRNA_CONSTRAINT_CONTEXT_NONE;
              } else if (hc->depot->up[sj][actual_j].direction > 0) {
                /* nothing to do */
              }
            } else if (!(constraint & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
              if (constraint & VRNA_CONSTRAINT_CONTEXT_ENFORCE) {
                /* remove (j, i) base pair */
                hc->matrix_local[j][i - j] = VRNA_CONSTRAINT_CONTEXT_NONE;
              } else {
                constraint = ~constraint & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
                hc->matrix_local[j][i - j] = constraint;
              }
            }
          }
        }
      }

      /* now for the actual base pair constraints */
      if (hc->depot->bp) {
      
      }
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
  vrna_hc_t *hc = vc->hc;

  /* ######################### */
  /* fill with default values  */
  /* ######################### */

  default_hc_up(vc, 0);

  default_hc_bp(vc, 0);

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
     *  Note, circular fold is only possible for single strand predictions
     */
    if (vc->strands < 2) {
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
