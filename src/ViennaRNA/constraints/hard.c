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
#include "ViennaRNA/sequences/alignments.h"
#include "ViennaRNA/utils/log.h"
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
                        unsigned int          i,
                        unsigned int          j);


PRIVATE INLINE void
populate_hc_up(vrna_fold_compound_t *fc,
               unsigned int         i,
               unsigned int         options);


PRIVATE INLINE void
populate_hc_bp(vrna_fold_compound_t *fc,
               unsigned int         i,
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
                    unsigned int          i,
                    unsigned int          options);


PRIVATE void
default_hc_up(vrna_fold_compound_t  *fc,
              unsigned int          options);


PRIVATE void
default_hc_bp(vrna_fold_compound_t  *fc,
              unsigned int          options);


PRIVATE void
prepare_hc_up(vrna_fold_compound_t  *fc,
              unsigned int          options);


PRIVATE void
prepare_hc_bp(vrna_fold_compound_t  *fc,
              unsigned int          options);


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
  hc->mx      = (unsigned char *)vrna_alloc(sizeof(unsigned char) * ((n + 1) * (n + 1) + 1));
  hc->up_ext  = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n + 2));
  hc->up_hp   = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n + 2));
  hc->up_int  = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n + 2));
  hc->up_ml   = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n + 2));
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
  unsigned int  n;
  vrna_hc_t     *hc;

  if (fc) {
    n   = fc->length;
    hc  = fc->hc;

    if (i > n) {
      vrna_log_warning("vrna_hc_update(): Position %u out of range!",
                           " (Sequence length: %u)",
                           i, n);
    } else {
      /* init up_xx arrays if necessary */
      if (!hc->up_ext) {
        hc->up_ext  = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n + 2));
        hc->up_hp   = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n + 2));
        hc->up_int  = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n + 2));
        hc->up_ml   = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n + 2));

        hc_update_up(fc);
      }

      populate_hc_up(fc, i, options);
      populate_hc_bp(fc, i, options);
    }
  }
}


PUBLIC int
vrna_hc_prepare(vrna_fold_compound_t  *fc,
                unsigned int          options)
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
    ret           = 1;
  }

  return ret;
}


PUBLIC void
vrna_hc_add_up(vrna_fold_compound_t *fc,
               unsigned int         i,
               unsigned char        option)
{
  (void)vrna_hc_add_up_strand(fc, i, -1, option);
}


PUBLIC int
vrna_hc_add_up_strand(vrna_fold_compound_t  *fc,
                      unsigned int          i,
                      int                   strand_indicator,
                      unsigned char         option)
{
  unsigned int  strand, n_i;
  int           ret = 0;

  if ((fc) &&
      (fc->hc) &&
      (strand_indicator < (int)fc->strands) &&
      (i > 0)) {

    /* auto-detect strand? */
    if (strand_indicator < 0) {
      strand  = fc->strand_number[i];
      i       = i - fc->strand_start[strand] + 1;
    } else {
      strand = (unsigned int)strand_indicator;
    }

    n_i = (fc->type == VRNA_FC_TYPE_SINGLE) ?
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
      sn  = fc->strand_number;
      ss  = fc->strand_start;

      for (i = 0; constraints[i].position != 0; i++) {
        pos     = constraints[i].position;
        options = constraints[i].options;

        if ((pos <= 0) || ((unsigned int)pos > fc->length))
          break;

        strand_i  = sn[pos];
        actual_i  = pos - ss[strand_i] + 1;

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
vrna_hc_add_up_strand_batch(vrna_fold_compound_t  *fc,
                            vrna_hc_up_t          *constraints)
{
  unsigned char options;
  unsigned int  i, strand, pos, n_pos;
  int           ret;

  ret = 0; /* failure */

  if ((fc) &&
      (constraints)) {
    if (fc->hc) {
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
                           unsigned int         i,
                           int                  d,
                           unsigned char        option)
{
  unsigned int  strand, actual_i, *sn, *ss;
  vrna_hc_t     *hc;

  if (vc) {
    if (vc->hc) {
      if ((i <= 0) || (i > vc->length)) {
        vrna_log_warning("vrna_hc_add_bp_nonspecific: position out of range, not doing anything");
        return;
      }

      hc        = vc->hc;
      sn        = vc->strand_number;
      ss        = vc->strand_start;
      strand    = sn[i];
      actual_i  = i - ss[strand] + 1;

      hc_depot_store_nonspec(vc,
                             actual_i,
                             strand,
                             d,
                             option);

      hc->state |= STATE_DIRTY_UP;
    }
  }
}


PUBLIC int
vrna_hc_add_bp_strand(vrna_fold_compound_t  *fc,
                      unsigned int          i,
                      unsigned int          j,
                      int                   strand_i,
                      int                   strand_j,
                      unsigned char         option)
{
  unsigned int  n_i, n_j, *sn, *ss, turn;
  int           ret;
  vrna_hc_t     *hc;

  ret = 0;

  if ((fc) &&
      (fc->hc) &&
      (strand_i < (int)fc->strands) &&
      (strand_j < (int)fc->strands) &&
      (i > 0) &&
      (j > 0)) {
    sn  = fc->strand_number;
    ss  = fc->strand_start;
    hc  = fc->hc;

    /* check whether we need to autodetect strand(s) */
    if (strand_i < 0) {
      strand_i  = sn[i];
      i  = i - ss[strand_i] + 1;
    }

    if (strand_j < 0) {
      strand_j  = sn[j];
      j         = j - ss[strand_j] + 1;
    }

    n_i = (fc->type == VRNA_FC_TYPE_SINGLE) ?
          fc->nucleotides[strand_i].length :
          fc->alignment[strand_i].sequences[0].length;
    n_j = (fc->type == VRNA_FC_TYPE_SINGLE) ?
          fc->nucleotides[strand_j].length :
          fc->alignment[strand_j].sequences[0].length;
    turn = fc->params->model_details.min_loop_size;

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


PUBLIC int
vrna_hc_add_bp(vrna_fold_compound_t *vc,
               unsigned int         i,
               unsigned int         j,
               unsigned char        option)
{
  unsigned int  *sn, *ss, strand_i, strand_j, actual_i, actual_j;
  int           ret;

  ret = 0;

  if (vc) {
    sn  = vc->strand_number;
    ss  = vc->strand_start;

    if (vc->hc) {
      if ((i <= 0) ||
          (j <= i) ||
          (j > vc->length)) {
        vrna_log_warning("vrna_hc_add_bp: position out of range, omitting constraint");
      } else if ((sn[i] == sn[j]) &&
                 ((j - i - 1) < (unsigned int)vc->params->model_details.min_loop_size)) {
        vrna_log_warning(
          "vrna_hc_add_bp: Pairing partners (%d, %d) violate minimum loop size settings of %dnt, omitting constraint",
          i,
          j,
          vc->params->model_details.min_loop_size);
      } else {
        /*
         *  determine the corresponding strand numbers and the actual
         *  position (relative to the strand)
         */
        strand_i  = sn[i];
        strand_j  = sn[j];
        actual_i  = i - ss[strand_i] + 1;
        actual_j  = j - ss[strand_j] + 1;

        ret = vrna_hc_add_bp_strand(vc,
                                    actual_i,
                                    actual_j,
                                    strand_i,
                                    strand_j,
                                    option);
      }
    }
  }

  return ret;
}


PUBLIC void
vrna_hc_free(vrna_hc_t *hc)
{
  if (hc) {
    if (hc->type == VRNA_HC_DEFAULT)
      free(hc->mx);
    else if (hc->type == VRNA_HC_WINDOW)
      free(hc->matrix_local);

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
vrna_hc_add_f(vrna_fold_compound_t  *vc,
              vrna_hc_eval_f        f)
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
vrna_hc_add_data(vrna_fold_compound_t *vc,
                 void                 *data,
                 vrna_auxdata_free_f  f)
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
                        unsigned int          i,
                        unsigned int          j)
{
  unsigned char constraint, can_stack;
  short         *S;
  unsigned int  *sn;
  unsigned int  type;
  vrna_md_t     *md;

  sn          = fc->strand_number;
  md          = &(fc->params->model_details);
  constraint  = VRNA_CONSTRAINT_CONTEXT_NONE;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      S = fc->sequence_encoding2;
      if (((j - i) < (unsigned int)md->max_bp_span) &&
          ((sn[i] != sn[j]) || ((j - i) > (unsigned int)md->min_loop_size))) {
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
              (((j - i + 2) < (unsigned int)md->max_bp_span) || (sn[i - 1] != sn[j + 1])) &&
              (md->pair[S[i - 1]][S[j + 1]]))
            can_stack = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

          /* can it enclose another base pair? */
          if ((i + 2 < j) &&
              (((j - i - 2) > (unsigned int)md->min_loop_size) || (sn[i + 1] != sn[j - 1])) &&
              (md->pair[S[i + 1]][S[j - 1]]))
            can_stack = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

          constraint &= can_stack;
        }
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      if ((sn[i] != sn[j]) ||
          (((j - i + 1) <= (unsigned int)md->max_bp_span) && ((j - i - 1) >= (unsigned int)md->min_loop_size))) {
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
              (((j - i + 2) < (unsigned int)md->max_bp_span) || (sn[i - 1] != sn[j + 1]))) {
            int outer_pscore = (fc->hc->type == VRNA_HC_WINDOW) ?
                               fc->pscore_local[i - 1][j - i + 2] :
                               fc->pscore[fc->jindx[j + 1] + i - 1];
            if (outer_pscore >= min_score)
              can_stack = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
          }

          /* can it enclose another base pair? */
          if ((i + 2 < j) &&
              (((j - i - 2) > (unsigned int)md->min_loop_size) || (sn[i + 1] != sn[j - 1]))) {
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


PRIVATE INLINE void
populate_hc_up(vrna_fold_compound_t *fc,
               unsigned int         i,
               unsigned int         options VRNA_UNUSED)
{
  unsigned char context;
  unsigned int  actual_i, strand;

  vrna_hc_t     *hc = fc->hc;

  if (hc->type == VRNA_HC_WINDOW) {
    strand    = fc->strand_number[i];
    actual_i  = i - fc->strand_start[strand] + 1;

    if ((hc->depot) &&
        (hc->depot->up) &&
        (hc->depot->up_size[strand] >= actual_i)) {
      /* We use user-defined constraints for unpaired nucleotides */
      context = hc->depot->up[strand][actual_i].context;

      if (hc->depot->up[strand][actual_i].nonspec) {
        if (context & VRNA_CONSTRAINT_CONTEXT_ENFORCE)
          /* this nucleotide must not stay unpaired */
          hc->matrix_local[i][0] = VRNA_CONSTRAINT_CONTEXT_NONE;
        /* the actual pairing direction will be enforced by populate_hc_bp() */
      } else if (context & VRNA_CONSTRAINT_CONTEXT_ENFORCE) {
        hc->matrix_local[i][0] = context & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
        /* restriction of base pairing contexts and removal of possible pairing will be handled by populate_hc_bp() */
      } else {
        hc->matrix_local[i][0] = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
        /* restriction of base pairing contexts and removal of possible pairing will be handled by populate_hc_bp() */
      }
    } else {
      /* if no constraints are available, we simply allow it to be unpaired in all contexts */
      hc->matrix_local[i][0] = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
    }
  } else {
    /* do something reasonable here... */
  }
}


PRIVATE void
default_hc_up(vrna_fold_compound_t  *fc,
              unsigned int          options)
{
  unsigned int  i, n;
  vrna_hc_t     *hc;

  hc = fc->hc;

  if (options & VRNA_OPTION_WINDOW) {
  } else {
    n = fc->length;

    for (i = 1; i <= n; i++)
      hc->mx[n * i + i] = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
  }
}


PRIVATE void
prepare_hc_up(vrna_fold_compound_t  *fc,
              unsigned int          options)
{
  unsigned char   option, type, t1, t2;
  unsigned int    i, j, k, n, s, *ss;
  vrna_hc_t       *hc;
  vrna_hc_depot_t *depot;

  hc = fc->hc;

  if (options & VRNA_OPTION_WINDOW) {
  } else {
    n     = fc->length;
    ss    = fc->strand_start;
    depot = hc->depot;

    /* 2. apply constraints as stored in depot */
    if ((depot) && (depot->up)) {
      for (s = 0; s < depot->strands; s++) {
        for (k = 1; k <= depot->up_size[s]; k++) {
          /* process nucleotide-specific constraint */
          option  = depot->up[s][k].context;
          i       = ss[s] + k - 1; /* constraint position in current strand order */

          if (depot->up[s][k].nonspec) {
            /* this is actually a must-pair constraint */
            type = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
            /* acknowledge pairing direction */
            t1  = (depot->up[s][k].direction <= 0) ? type : VRNA_CONSTRAINT_CONTEXT_NONE;
            t2  = (depot->up[s][k].direction >= 0) ? type : VRNA_CONSTRAINT_CONTEXT_NONE;

            if (option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE) {
              /* only allow for possibly non-canonical pairs, do not enforce them */
              for (j = 1; j < i; j++) {
                hc->mx[n * i + j] |= t1;
                hc->mx[n * j + i] |= t1;
              }
              for (j = i + 1; j <= n; j++) {
                hc->mx[n * i + j] |= t2;
                hc->mx[n * j + i] |= t2;
              }
            } else {
              /* force pairing direction */
              for (j = 1; j < i; j++) {
                hc->mx[n * i + j] &= t1;
                hc->mx[n * j + i] &= t1;
              }
              for (j = i + 1; j <= n; j++) {
                hc->mx[n * i + j] &= t2;
                hc->mx[n * j + i] &= t2;
              }
            }

            /* nucleotide mustn't be unpaired */
            if (option & VRNA_CONSTRAINT_CONTEXT_ENFORCE)
              hc->mx[n * i + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
          } else {
            /* 'regular' nucleotide-specific constraint */
            if (option & VRNA_CONSTRAINT_CONTEXT_ENFORCE) {
              /*
               * force nucleotide to appear unpaired within a certain type of loop
               * do not allow i to be paired with any other nucleotide
               */
              if (!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
                for (j = 1; j < i; j++) {
                  hc->mx[n * i + j] = VRNA_CONSTRAINT_CONTEXT_NONE;
                  hc->mx[n * j + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
                }
                for (j = i + 1; j <= n; j++) {
                  hc->mx[n * i + j] = VRNA_CONSTRAINT_CONTEXT_NONE;
                  hc->mx[n * j + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
                }
              }

              type = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

              hc->mx[n * i + i] = type;
            } else {
              type = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

              /* do not allow i to be paired with any other nucleotide (in context type) */
              if (!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
                for (j = 1; j < i; j++) {
                  hc->mx[n * i + j] &= ~type;
                  hc->mx[n * j + i] &= ~type;
                }
                for (j = i + 1; j <= n; j++) {
                  hc->mx[n * i + j] &= ~type;
                  hc->mx[n * j + i] &= ~type;
                }
              }

              hc->mx[n * i + i] = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
            }
          }
        }
      }
    }
  }
}


PRIVATE void
default_hc_bp(vrna_fold_compound_t  *fc,
              unsigned int          options)
{
  unsigned int  i, j, n;
  vrna_hc_t     *hc;

  hc = fc->hc;

  if (options & VRNA_OPTION_WINDOW) {
  } else {
    n = fc->length;

    for (j = n; j > 1; j--) {
      for (i = 1; i < j; i++) {
        hc->mx[n * i + j] = default_pair_constraint(fc, i, j);
        hc->mx[n * j + i] = hc->mx[n * i + j];
      }
    }
  }
}


PRIVATE void
prepare_hc_bp(vrna_fold_compound_t  *fc,
              unsigned int          options)
{
  unsigned char   option;
  unsigned int    i, j, k, p, q, n, actual_i, actual_j, strand_j, s, *ss;
  int             *idx;
  vrna_hc_t       *hc;
  vrna_hc_depot_t *depot;

  hc    = fc->hc;
  depot = hc->depot;
  ss    = fc->strand_start;

  if ((!depot) || (!depot->bp))
    return;

  if (options & VRNA_OPTION_WINDOW) {
  } else {
    n   = fc->length;
    idx = fc->jindx;

    /* 2. apply constraints stored in depot */
    for (s = 0; s < depot->strands; s++) {
      for (actual_i = 1; actual_i <= depot->bp_size[s]; actual_i++) {
        for (k = 0; k < depot->bp[s][actual_i].list_size; k++) {
          option    = depot->bp[s][actual_i].context[k];
          strand_j  = depot->bp[s][actual_i].strand_j[k];
          actual_j  = depot->bp[s][actual_i].j[k];
          i         = ss[s] + actual_i - 1;         /* constraint position in current strand order */
          j         = ss[strand_j] + actual_j - 1;  /* constraint position in current strand order */

          if (i < j) {
            /* apply the constraint */
            hc->mx[n * i + j] = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
            hc->mx[n * j + i] = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

            /* is the ptype reset actually required??? */
            if ((fc->type == VRNA_FC_TYPE_SINGLE) &&
                (option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS)) {
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
                hc->mx[n * i + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
                hc->mx[n * p + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
                hc->mx[n * j + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
                hc->mx[n * p + j] = VRNA_CONSTRAINT_CONTEXT_NONE;

                for (q = i + 1; q < j; q++) {
                  hc->mx[n * p + q] = VRNA_CONSTRAINT_CONTEXT_NONE;
                  hc->mx[n * q + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
                }
              }
              for (p = i + 1; p < j; p++) {
                hc->mx[n * i + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
                hc->mx[n * p + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
                hc->mx[n * j + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
                hc->mx[n * p + j] = VRNA_CONSTRAINT_CONTEXT_NONE;

                for (q = j + 1; q <= n; q++) {
                  hc->mx[n * p + q] = VRNA_CONSTRAINT_CONTEXT_NONE;
                  hc->mx[n * q + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
                }
              }
              for (p = j + 1; p <= n; p++) {
                hc->mx[n * i + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
                hc->mx[n * p + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
                hc->mx[n * j + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
                hc->mx[n * p + j] = VRNA_CONSTRAINT_CONTEXT_NONE;
              }
            }

            if (option & VRNA_CONSTRAINT_CONTEXT_ENFORCE) {
              /* do not allow i,j to be unpaired */
              hc->mx[n * i + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
              hc->mx[n * j + j] = VRNA_CONSTRAINT_CONTEXT_NONE;
            }
          }
        }
      }
    }
  }
}


PRIVATE INLINE void
populate_hc_bp(vrna_fold_compound_t *fc,
               unsigned int         i,
               unsigned int         options)
{
  unsigned char constraint, type, t_down, t_up;
  unsigned int  max_span, maxdist, j, k, p, n, l, strand, sj, sl, actual_i, actual_j,
                actual_l, *sn, *ss;
  vrna_hc_t     *hc;

  n         = fc->length;
  maxdist   = fc->window_size;
  max_span  = fc->params->model_details.max_bp_span;
  sn        = fc->strand_number;
  ss        = fc->strand_start;
  strand    = sn[i];
  actual_i  = i - ss[strand] + 1;
  hc        = fc->hc;

  if (options & VRNA_OPTION_F3) {
    /* the sliding window moves from 3' to 5' side */

    /* apply default constraints first */
    for (k = 1; k < max_span; k++) {
      j = i + k;
      if (j > n)
        break;

      constraint = default_pair_constraint(fc, i, j);

      hc->matrix_local[i][j - i] = constraint;
    }

    if (hc->depot) {
      /* 1. apply remainder of (partly) applied nucleotide-specific constraints */
      if (hc->depot->up) {
        /* 1.a apply nucleotide specific constraints for i */
        if ((hc->depot->up[strand]) &&
            (hc->depot->up_size[strand] >= actual_i)) {
          constraint = hc->depot->up[strand][actual_i].context;
          type = constraint & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

          if (hc->depot->up[strand][actual_i].nonspec) {
            /* handle unspecific pairing contraint */
            t_down  = (hc->depot->up[strand][actual_i].direction >= 0) ? type : VRNA_CONSTRAINT_CONTEXT_NONE;

            if (constraint & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)
              /* only allow for possibly non-canonical pairs, do not enforce them */
              for (k = i + 1; k < MIN2(i + max_span, n + 1); k++)
                hc->matrix_local[i][k - i] |= t_down;
            else
              /* remove downstream pairs if necessary */
              for (k = i + 1; k < MIN2(i + max_span, n + 1); k++)
                hc->matrix_local[i][k - i] &= t_down;
          } else {
            /* handle 'regular' unpairedness constraint */

            if (constraint & VRNA_CONSTRAINT_CONTEXT_ENFORCE) {
              if (!(constraint & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
                /* do not allow i to be paired with any other nucleotide */
                for (k = i + 1; k < MIN2(i + max_span, n + 1); k++)
                  hc->matrix_local[i][k - i] = VRNA_CONSTRAINT_CONTEXT_NONE;
              }
            } else {
              if (!(constraint & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
                /* only allow i to be paired for particular types */
                for (k = i + 1; k < MIN2(i + max_span, n + 1); k++)
                  hc->matrix_local[i][k - i] &= ~type;
              }
            }
          }
        }

        /*
         *  1.b apply remainder of (partly) applied nucleotide-specific
         *  constraints for j with i < j < i + maxdist
         */
        for (k = 1; k < max_span; k++) {
          j         = i + k;
          sj        = sn[j];
          actual_j  = j - ss[sj] + 1;

          if (j > n)
            break;

          if ((hc->depot->up[sj]) &&
              (hc->depot->up_size[sj] >= actual_j)) {
            constraint = hc->depot->up[sj][actual_j].context;
            type = constraint & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

            if (hc->depot->up[sj][actual_j].nonspec) {
              t_up  = (hc->depot->up[sj][actual_j].direction <= 0) ? type : VRNA_CONSTRAINT_CONTEXT_NONE;
              if (constraint & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)
                hc->matrix_local[i][j - i] |= t_up;
              else
                hc->matrix_local[i][j - i] &= t_up;
            } else if (!(constraint & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
              if (constraint & VRNA_CONSTRAINT_CONTEXT_ENFORCE) {
                /* remove (i, j) base pair */
                hc->matrix_local[i][j - i] = VRNA_CONSTRAINT_CONTEXT_NONE;
              } else {
                hc->matrix_local[i][j - i]  &= ~type;
              }
            }
          }
        }
      }

      /* 2. now for the actual base pair constraints */
      if (hc->depot->bp) {
        /* 2.a apply base pair specific constraint for nucleotide i */
        if ((hc->depot->bp[strand]) &&
            (hc->depot->bp_size[strand] >= actual_i) &&
            (hc->depot->bp[strand][actual_i].list_size > 0)) {
          /* go through list of all constraints for this nucleotide */
          for (size_t cnt = 0; cnt < hc->depot->bp[strand][actual_i].list_size; cnt++) {
            constraint  = hc->depot->bp[strand][actual_i].context[cnt];
            actual_j    = hc->depot->bp[strand][actual_i].j[cnt];
            sj          = hc->depot->bp[strand][actual_i].strand_j[cnt];
            j           = ss[sj] + actual_j - 1;

            /* apply the constraint */

            /* do not allow i to stay unpaired if necessary */
            if (constraint & VRNA_CONSTRAINT_CONTEXT_ENFORCE)
              hc->matrix_local[i][0] = VRNA_CONSTRAINT_CONTEXT_NONE;

            if (j > i) {
              /* j is 3' base of base pair (i, j) */
              if (i + max_span > j)
                hc->matrix_local[i][j - i] = constraint & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

              /* remove other base pairs violating the constraint */
              if (!(constraint & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
                /* remove base pairs (i, k) with k != j */
                for (k = i + 1; k < j; k++)
                  hc->matrix_local[i][k - i] = VRNA_CONSTRAINT_CONTEXT_NONE;

                for (k = j + 1; k < MIN2(i + max_span, n + 1); k++)
                  hc->matrix_local[i][k - i] = VRNA_CONSTRAINT_CONTEXT_NONE;
              }
            } else {
              /* j is 5' base of base pair (j, i) */

              if (!(constraint & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
                /* this base pairs upstream, so remove all downstream pairing partners */
                for (k = 1; k < max_span; k++) {
                  j = i + k;
                  if (j > n)
                    break;

                  hc->matrix_local[i][j - i] = VRNA_CONSTRAINT_CONTEXT_NONE;
                }
              }
            }
          }
        }

        /* 2.b acknowledge base pair constraints that involve j with i < j < i + max_span */
        for (k = 1; k < max_span; k++) {
          j         = i + k;
          sj        = sn[j];
          actual_j  = j - ss[sj] + 1;

          if (j > n)
            break;

          if ((hc->depot->bp[sj]) &&
              (hc->depot->bp_size[sj] >= actual_j) &&
              (hc->depot->bp[sj][actual_j].list_size > 0)) {
            /* go through list of all constraints for this nucleotide */
            for (size_t cnt = 0; cnt < hc->depot->bp[sj][actual_j].list_size; cnt++) {
              constraint  = hc->depot->bp[sj][actual_j].context[cnt];
              actual_l    = hc->depot->bp[sj][actual_j].j[cnt];
              sl          = hc->depot->bp[sj][actual_j].strand_j[cnt];
              l           = ss[sl] + actual_l - 1;

              if (!(constraint & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
                if (l > j) {
                  /* remove base pairs (i, p) with i < j <= p <= l */
                  for (p = j; p < MIN2(i + max_span, l + 1); p++)
                    hc->matrix_local[i][p - i] = VRNA_CONSTRAINT_CONTEXT_NONE;
                } else if (i < l) {
                  /* remove base pairs (i, p) with i < l <= p < j */
                  for (p = l; p <= j; p++)
                    hc->matrix_local[i][p - i] = VRNA_CONSTRAINT_CONTEXT_NONE;
                } else if (l < i) {
                  /* remove base pairs (i, p) with i < j <= p */
                  for (p = j; p < MIN2(i + max_span, n + 1); p++)
                    hc->matrix_local[i][p - i] = VRNA_CONSTRAINT_CONTEXT_NONE;
                }
              }
            }
          }
        }
      }
    }

    hc_update_up_window(fc, i, options);
  } else if (options & VRNA_OPTION_F5) {
    /* the sliding window moves from 5' to 3' side (i is 3' nucleotide) */

    j         = i;
    sj        = strand;
    actual_j  = actual_i;

    if ((j <= n) && (j > 0)) {
      unsigned int start_i, min_i;

      start_i = (j > max_span) ? j - max_span + 1 : 1;
      min_i   = (j > maxdist) ? j - maxdist + 1 : 1;

      /* apply default constraints first */
      for (i = start_i; i < j; i++)
        hc->matrix_local[i][j - i] = default_pair_constraint(fc, i, j);

      /* next are user-defined constraints */
      if (hc->depot) {
        /* 1. apply remainder of (partly) applied nucleotide-specific constraints */
        if (hc->depot->up) {
          /* 1.a apply nucleotide-specific constraints for j */
          if ((hc->depot->up[sj]) &&
              (hc->depot->up_size[sj] >= actual_j)) {
            constraint = hc->depot->up[sj][actual_j].context;
            type  = constraint & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

            if (hc->depot->up[sj][actual_j].nonspec) {
              /* handle unspecific pairing constraint */
              t_up  = (hc->depot->up[sj][actual_j].direction <= 0) ? type : VRNA_CONSTRAINT_CONTEXT_NONE;

              if (constraint & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)
                for (i = start_i; i < j; i++)
                  hc->matrix_local[i][j - i] |= t_up;
              else
                for (i = start_i; i < j; i++)
                  hc->matrix_local[i][j - i] &= t_up;
            } else {
              /* handle 'regular' unpairedness constraint */

              if (constraint & VRNA_CONSTRAINT_CONTEXT_ENFORCE) {
                if (!(constraint & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
                  /* do not allow j to be paired with any other nucleotide */
                  for (i = start_i; i < j; i++)
                    hc->matrix_local[i][j - i] = VRNA_CONSTRAINT_CONTEXT_NONE;
                }
              } else {
                if (!(constraint & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
                  /* only allow j to be paired for particular types */
                  for (i = start_i; i < j; i++)
                    hc->matrix_local[i][j - i] &= ~type;
                }
              }
            }
          }

          /* 1.b apply remainder of (partly) applied nucleotide-specific
           *  constraints for i with j - maxdist < i < j
           */

          for (i = start_i; i < j; i++) {
            strand    = sn[i];
            actual_i  = i - ss[strand] + 1;

            if ((hc->depot->up[strand]) &&
                (hc->depot->up_size[strand] >= actual_i)) {
              constraint = hc->depot->up[strand][actual_i].context;
              type = constraint & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

              if (hc->depot->up[strand][actual_i].nonspec) {
                t_down  = (hc->depot->up[sj][actual_j].direction >= 0) ? type : VRNA_CONSTRAINT_CONTEXT_NONE;
                if (constraint & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)
                  hc->matrix_local[i][j - i] |= t_down;
                else
                  hc->matrix_local[i][j - i] &= t_down;
              } else if (!(constraint & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
                if (constraint & VRNA_CONSTRAINT_CONTEXT_ENFORCE) {
                  /* remove (i, j) base pair */
                  hc->matrix_local[i][j - i] = VRNA_CONSTRAINT_CONTEXT_NONE;
                } else {
                  hc->matrix_local[i][j - i]  &= ~type;
                }
              }
            }
          }
        }

        /* 2. now for the actual base pair constraints */
        if (hc->depot->bp) {
          /* 2.a apply base pair specific constraint for nucleotide j */
          if ((hc->depot->bp[sj]) &&
              (hc->depot->bp_size[sj] >= actual_j) &&
              (hc->depot->bp[sj][actual_j].list_size > 0)) {
            for (size_t cnt = 0; cnt < hc->depot->bp[sj][actual_j].list_size; cnt++) {
              constraint  = hc->depot->bp[sj][actual_j].context[cnt];
              actual_i    = hc->depot->bp[sj][actual_j].j[cnt];
              strand      = hc->depot->bp[sj][actual_j].strand_j[cnt];
              i           = ss[strand] + actual_i - 1;

              /*
               * apply the constraint
               * do not allow j to stay unpaired if necessary
               */
              if (constraint & VRNA_CONSTRAINT_CONTEXT_ENFORCE)
                hc->matrix_local[j][0] = VRNA_CONSTRAINT_CONTEXT_NONE;

              if (i < j) {
                /* j is 3' base of base pair (i, j) */
                if (i + max_span > j)
                  hc->matrix_local[i][j - i] = constraint & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

                /* remove other base pairs violating the constraint */
                if (!(constraint & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
                  /* remove base pairs (k, j) with k != i */
                  for (k = start_i; k < i; k++)
                    hc->matrix_local[k][j - k] = VRNA_CONSTRAINT_CONTEXT_NONE;

                  for (k = i + 1; k < j; k++)
                    hc->matrix_local[k][j - k] = VRNA_CONSTRAINT_CONTEXT_NONE;
                }
              } else {
                /* j is 5' base of base pair (j, i) */

                if (!(constraint & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
                  /* this base pair pairs downstream, so remove all upstream pairing partners */
                  for (k = start_i; k < j; k++)
                    hc->matrix_local[k][j - k] = VRNA_CONSTRAINT_CONTEXT_NONE;
                }
              }
            }
          }

          /* 2.b acknowledge base pair constraints that involve i with j - maxdist < i < j */
          for (i = min_i; i < j; i++) {
            strand    = sn[i];
            actual_i  = i - ss[strand] + 1;

            if ((hc->depot->bp[strand]) &&
                (hc->depot->bp_size[strand] >= actual_i) &&
                (hc->depot->bp[strand][actual_i].list_size > 0)) {
              for (size_t cnt = 0; cnt < hc->depot->bp[strand][actual_i].list_size; cnt++) {
                constraint  = hc->depot->bp[strand][actual_i].context[cnt];
                actual_l    = hc->depot->bp[strand][actual_i].j[cnt];
                sl          = hc->depot->bp[strand][actual_i].strand_j[cnt];
                l           = ss[sl] + actual_l - 1;

                if (!(constraint & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
                  if (l < i) {
                    /* remove base pairs (p, j) with l <= p <= i < j */
                    for (p = MAX2(l, start_i); p <= i; p++)
                      hc->matrix_local[p][j - p] = VRNA_CONSTRAINT_CONTEXT_NONE;
                  } else if (l < j) {
                    /* remove base pairs (p, j) with i <= p <= l < j */
                    for (p = MAX2(i, start_i); p <= l; p++)
                      hc->matrix_local[p][j - p] = VRNA_CONSTRAINT_CONTEXT_NONE;
                  } else if (l > j) {
                    /* remove base pairs (p, j) with p <= i */
                    for (p = start_i; p <= i; p++)
                      hc->matrix_local[p][j - p] = VRNA_CONSTRAINT_CONTEXT_NONE;
                  }
                }
              }
            }
          }
        }
      }
    }

    hc_update_up_window(fc, j, options);
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
  unsigned int  n, i, j, spos;
  int           hx, *stack;
  vrna_md_t     *md;

  vrna_array(vrna_hc_up_t)  up;
  vrna_array(struct hc_bp)  bp;
  vrna_array(struct hc_bp)  bp_unspecific;

  if (constraint == NULL)
    return;

  sequence      = vc->sequence;
  S             = vc->sequence_encoding2;
  md            = &(vc->params->model_details);
  n             = strlen(constraint);
  stack         = (int *)vrna_alloc(sizeof(int) * (n + 1));

  vrna_array_init(up);
  vrna_array_init(bp);
  vrna_array_init(bp_unspecific);

  ;
  for (hx = 0, spos = 0, j = 1; spos < n; j++, spos++) {
    switch (constraint[spos]) {
      case '&':
        /* multistrand delimiter */
        j--;
        break;

      /* can't pair */
      case 'x':
        if (options & VRNA_CONSTRAINT_DB_X) {
          vrna_array_append(up,
                            ((vrna_hc_up_t){
                              .position = j,
                              .options = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS
                            }));
        }

        break;

      /* must pair, i.e. may not be unpaired */
      case '|':
        if (options & VRNA_CONSTRAINT_DB_PIPE) {
          vrna_array_append(bp_unspecific,
                            ((struct hc_bp){
                              .i = j,
                              .j = 0,
                              .options = (options & VRNA_CONSTRAINT_DB_ENFORCE_BP) ?
                                         VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS | VRNA_CONSTRAINT_CONTEXT_ENFORCE:
                                         VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS
                            }));
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
            vrna_log_warning("vrna_hc_add_from_db: "
                                 "Unbalanced brackets in constraint string\n%s\n"
                                 "No constraints will be applied!",
                                 constraint);
            goto db_constraints_exit;
          }

          i = stack[--hx];

          if (options & VRNA_CONSTRAINT_DB_CANONICAL_BP) {
            /* check whether this pair forms a non-canoncial base pair */
            if (md->pair[S[i]][S[j]] == 0) {
              vrna_log_warning("Removing non-canonical base pair %c%c (%d,%d) from constraint",
                                   sequence[i - 1], sequence[j - 1],
                                   i, j);
              break;
            }
          }

          if ((j - i - 1) < (unsigned int)md->min_loop_size) {
            vrna_log_warning("vrna_hc_add_from_db: "
                                 "Pairing partners (%d, %d) violate minimum loop size settings of %dnt, omitting constraint",
                                 i,
                                 j,
                                 md->min_loop_size);
            break;
          }

          vrna_array_append(bp,
                            ((struct hc_bp){
                              .i = i,
                              .j = j,
                              .options = (options & VRNA_CONSTRAINT_DB_ENFORCE_BP) ?
                                          VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS | VRNA_CONSTRAINT_CONTEXT_ENFORCE:
                                          VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS
                            }));
        }

        break;

      /* pairs downstream */
      case '<':
        if (options & VRNA_CONSTRAINT_DB_ANG_BRACK) {
          vrna_array_append(bp_unspecific,
                            ((struct hc_bp){
                              .i = j,
                              .j = 1,
                              .options = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS
                            }));
        }

        break;

      /* pairs upstream */
      case '>':
        if (options & VRNA_CONSTRAINT_DB_ANG_BRACK) {
          vrna_array_append(bp_unspecific,
                            ((struct hc_bp){
                              .i        = j,
                              .j        = -1,
                              .options  = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS
                            }));
        }

        break;

      /* only intramolecular basepairing */
      case 'l':
        if (options & VRNA_CONSTRAINT_DB_INTRAMOL) {
          unsigned int l;
          if (vc->strands > 1) {
            for (l = 1; l < vc->strand_start[vc->strand_number[j]]; l++) {
              vrna_array_append(bp,
                                ((struct hc_bp){
                                  .i        = j,
                                  .j        = l,
                                  .options  = VRNA_CONSTRAINT_CONTEXT_NONE |
                                              VRNA_CONSTRAINT_CONTEXT_NO_REMOVE
                                }));
            }

            for (l = vc->strand_end[vc->strand_number[j]] + 1; l <= vc->length; l++) {
              vrna_array_append(bp,
                                ((struct hc_bp){
                                  .i        = j,
                                  .j        = l,
                                  .options  = VRNA_CONSTRAINT_CONTEXT_NONE |
                                              VRNA_CONSTRAINT_CONTEXT_NO_REMOVE
                                }));
            }
          }
        }

        break;

      /* only intermolecular bp */
      case 'e':
        if (options & VRNA_CONSTRAINT_DB_INTERMOL) {
          unsigned int l;
          if (vc->strands > 1) {
            for (l = vc->strand_start[vc->strand_number[j]]; l < j; l++) {
                vrna_array_append(bp,
                                  ((struct hc_bp){
                                    .i        = l,
                                    .j        = j,
                                    .options  = VRNA_CONSTRAINT_CONTEXT_NONE |
                                                VRNA_CONSTRAINT_CONTEXT_NO_REMOVE
                                  }));
            }

            for (l = j + 1; l <= vc->strand_end[vc->strand_number[j]]; l++) {
                vrna_array_append(bp,
                                  ((struct hc_bp){
                                    .i        = j,
                                    .j        = l,
                                    .options  = VRNA_CONSTRAINT_CONTEXT_NONE |
                                                VRNA_CONSTRAINT_CONTEXT_NO_REMOVE
                                  }));
            }
          }
        }

        break;

      case '.':
        break;

      default:
        vrna_log_warning(
          "vrna_hc_add_from_db: "
          "Unrecognized character '%c' in constraint string",
          constraint[spos]);
        break;
    }
  }

  if (hx != 0) {
    vrna_log_warning("vrna_hc_add_from_db: "
                         "Unbalanced brackets in constraint string\n%s\n"
                         "No constraints will be applied!",
                         constraint);
    goto db_constraints_exit;
  }

  /* finally, apply constraints */

  /* 1st, unspecific pairing states */
  for (i = 0; i < vrna_array_size(bp_unspecific); i++)
    vrna_hc_add_bp_nonspecific(vc,
                               bp_unspecific[i].i,  /* nucleotide position */
                               bp_unspecific[i].j,  /* pairing direction */
                               bp_unspecific[i].options);

  /* 2nd, specific base pairs */
  for (i = 0; i < vrna_array_size(bp); i++)
    vrna_hc_add_bp(vc,
                   bp[i].i,
                   bp[i].j,
                   bp[i].options);

  /* 3rd, unpaired constraints */
  if (vrna_array_size(up) > 0) {
    vrna_array_append(up,
                      ((vrna_hc_up_t){
                        .position = 0 /* end of list marker */
                      }));
    vrna_hc_add_up_batch(vc, up);
  }

db_constraints_exit:

  /* clean up */
  vrna_array_free(up);
  vrna_array_free(bp);
  vrna_array_free(bp_unspecific);
  free(stack);
}


PRIVATE void
hc_reset_to_default(vrna_fold_compound_t *vc)
{
  vrna_hc_t *hc = vc->hc;

  /*
   * #########################
   * fill with default values
   * #########################
   */

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
  unsigned int  i, n;
  vrna_hc_t     *hc;

  n   = vc->length;
  hc  = vc->hc;

  if (hc->type == VRNA_HC_WINDOW) {
    /* do nothing for now! */
  } else {
    for (hc->up_ext[n + 1] = 0, i = n; i > 0; i--) /* unpaired stretch in exterior loop */
      hc->up_ext[i] = (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) ? 1 +
                      hc->up_ext[i + 1] : 0;

    for (hc->up_hp[n + 1] = 0, i = n; i > 0; i--)  /* unpaired stretch in hairpin loop */
      hc->up_hp[i] = (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) ? 1 +
                     hc->up_hp[i + 1] : 0;

    for (hc->up_int[n + 1] = 0, i = n; i > 0; i--) /* unpaired stretch in internal loop */
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
                    unsigned int          i,
                    unsigned int          options)
{
  vrna_hc_t     *hc;
  unsigned int  k, kmin, winsize, up_ext, up_hp, up_int, up_ml;

  hc      = vc->hc;
  winsize = vc->window_size;

  if (options & VRNA_OPTION_F5) {
    up_ext = (hc->matrix_local[i][0] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) ?
             1 : 0;
    up_hp = (hc->matrix_local[i][0] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) ?
            1 : 0;
    up_int = (hc->matrix_local[i][0] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) ?
             1 : 0;
    up_ml = (hc->matrix_local[i][0] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) ?
            1 : 0;
  } else {
    up_ext = (hc->matrix_local[i][0] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) ?
             1 + hc->up_ext[i + 1] :
             0;
    up_hp = (hc->matrix_local[i][0] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) ?
            1 + hc->up_hp[i + 1] :
            0;
    up_int = (hc->matrix_local[i][0] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) ?
             1 + hc->up_int[i + 1] :
             0;
    up_ml = (hc->matrix_local[i][0] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) ?
            1 + hc->up_ml[i + 1] :
            0;
  }

  hc->up_ext[i] = up_ext;
  hc->up_hp[i]  = up_hp;
  hc->up_int[i] = up_int;
  hc->up_ml[i]  = up_ml;

  if (options & VRNA_OPTION_F5) {
    /* the sliding window proceeds from 5' to 3' so we update constraints 3' to 5' */
    if (up_ext > 0) {
      kmin = 1;
      if (i > winsize)
        kmin = i - winsize;
      for (k = i - 1; k >= kmin; k--) {
        if (hc->up_ext[k] < 1)
          break;

        hc->up_ext[k] += up_ext;
      }
    }

    if (up_hp > 0) {
      kmin = 1;
      if (i > winsize)
        kmin = i - winsize;
      for (k = i - 1; k >= kmin; k--) {
        if (hc->up_hp[k] < 1)
          break;

        hc->up_hp[k] += up_hp;
      }
    }

    if (up_int > 0) {
      kmin = 1;
      if (i > winsize)
        kmin = i - winsize;
      for (k = i - 1; k >= kmin; k--) {
        if (hc->up_int[k] < 1)
          break;

        hc->up_int[k] += up_int;
      }
    }

    if (up_ml > 0) {
      kmin = 1;
      if (i > winsize)
        kmin = i - winsize;
      for (k = i - 1; k >= kmin; k--) {
        if (hc->up_ml[k] < 1)
          break;

        hc->up_ml[k] += up_ml;
      }
    }
  }
}


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
 * ###########################################
 * # deprecated functions below              #
 *###########################################
 */
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
            vrna_log_error("%s\nunbalanced brackets in constraint", constraint);

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
            vrna_log_error("%s\nunbalanced brackets in constraints", constraint);

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
    vrna_log_error("%s\nunbalanced brackets in constraint string", constraint);

  free(index);
  free(stack);
}


#endif
