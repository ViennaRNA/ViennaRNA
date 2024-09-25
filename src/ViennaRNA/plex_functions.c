/*
 *         compute potentially pseudoknotted structures of an RNA sequence
 *                           Ivo Hofacker
 *                        Vienna RNA package
 */

/*
 * library containing the function used in PKplex
 * it generates pseudoknotted structures by letting the sequence form a duplex structure with itself
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <time.h>

#include "ViennaRNA/params/default.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/eval/hairpin.h"
#include "ViennaRNA/eval/internal.h"
#include "ViennaRNA/eval/multibranch.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/LPfold.h"
#include "ViennaRNA/mfe/global.h"
#include "ViennaRNA/datastructures/heap.h"

#include "ViennaRNA/constraints/exterior_hc.inc"

#include "ViennaRNA/pk_plex.h"

#undef  MAXLOOP
#define MAXLOOP 10

#define DEFAULT_PENALTY             810
#define DEFAULT_DELTA               0
#define DEFAULT_INTERACTION_LENGTH  12
#define DEFAULT_CUTOFF              1e-6

struct vrna_pk_plex_option_s {
  unsigned int                delta;
  unsigned int                max_interaction_length;
  int                         min_penalty;
  vrna_pk_plex_score_f scoring_function;
  void                        *scoring_data;
};

typedef struct {
  int penalty;
} default_data;


typedef struct {
  double  kT;
  double  cutoff;
  int     **open_en;
} access_data;

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
PRIVATE int
default_pk_plex_penalty(const short *pt,
                        int         i,
                        int         j,
                        int         k,
                        int         l,
                        void        *data);


PRIVATE vrna_heap_t
duplexfold_XS(vrna_fold_compound_t        *fc,
              const int                   **access_s1,
              const int                   max_interaction_length,
              vrna_pk_plex_score_f scoring_function,
              void                        *scoring_data);


PRIVATE char *
backtrack_XS(vrna_fold_compound_t *fc,
             int                  kk,
             int                  ll,
             const int            ii,
             const int            jj,
             const int            max_interaction_length,
             int                  ***c3);


PRIVATE int
prepare_PKplex(vrna_fold_compound_t *fc);


PRIVATE int
PKplex_heap_cmp(const void  *a,
                const void  *b,
                void        *data);


PRIVATE void
store_pU_callback(double        *pU,
                  int           size,
                  int           k,
                  int           ulength,
                  unsigned int  type,
                  void          *data);

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_pk_plex_opt_t
vrna_pk_plex_opt_defaults(void)
{
  return vrna_pk_plex_opt(DEFAULT_DELTA,
                          DEFAULT_INTERACTION_LENGTH,
                          DEFAULT_PENALTY);
}


PUBLIC vrna_pk_plex_opt_t
vrna_pk_plex_opt(unsigned int delta,
                 unsigned int max_interaction_length,
                 int          pk_penalty)
{
  vrna_pk_plex_opt_t opt;

  opt = (vrna_pk_plex_opt_t)vrna_alloc(sizeof(struct vrna_pk_plex_option_s));

  opt->delta                  = delta;
  opt->max_interaction_length = max_interaction_length;
  opt->min_penalty            = pk_penalty;
  opt->scoring_function       = NULL;
  opt->scoring_data           = NULL;

  return opt;
}


PUBLIC vrna_pk_plex_opt_t
vrna_pk_plex_opt_fun(unsigned int                 delta,
                     unsigned int                 max_interaction_length,
                     vrna_pk_plex_score_f  scoring_function,
                     void                         *scoring_data)
{
  vrna_pk_plex_opt_t opt = NULL;

  if (scoring_function) {
    opt = (vrna_pk_plex_opt_t)vrna_alloc(sizeof(struct vrna_pk_plex_option_s));

    opt->delta                  = delta;
    opt->max_interaction_length = max_interaction_length;
    opt->scoring_function       = scoring_function;
    opt->scoring_data           = scoring_data;
  }

  return opt;
}


PUBLIC vrna_pk_plex_t *
vrna_pk_plex(vrna_fold_compound_t *fc,
             const int            **accessibility,
             vrna_pk_plex_opt_t   user_options)
{
  size_t              cnt;
  char                *mfe_struct, *constraint;
  int                 i, **opening_energies;
  double              mfe, mfe_pk, constrainedEnergy, pk_penalty, subopts, best_e;
  default_data        scoring_dat;
  vrna_pk_plex_t      *hits, *hit_ptr, *ptr, *mfe_entry;
  vrna_pk_plex_opt_t  options;
  vrna_heap_t         interactions, results;

  hits              = NULL;
  results           = NULL;
  opening_energies  = NULL;

  if (fc) {
    mfe_struct = (char *)vrna_alloc(sizeof(char) * (fc->length + 1));

    mfe = mfe_pk = (double)vrna_mfe(fc, mfe_struct);

    options = (user_options) ? user_options : vrna_pk_plex_opt_defaults();

    /* apply simplified scoring function with constant penalty */
    if (!options->scoring_function) {
      scoring_dat.penalty       = options->min_penalty;
      options->scoring_function = &default_pk_plex_penalty;
      options->scoring_data     = (void *)&scoring_dat;
    }

    /* sanity check for maximum interaction length */
    options->max_interaction_length = MIN2(DEFAULT_INTERACTION_LENGTH, fc->length - 3);

    /* compute opening energies if not passed as argument */
    if (!accessibility) {
      vrna_fold_compound_t *fca = vrna_fold_compound(fc->sequence,
                                                     &(fc->params->model_details),
                                                     VRNA_OPTION_DEFAULT | VRNA_OPTION_WINDOW);

      vrna_exp_params_rescale(fca, &mfe);
      opening_energies = vrna_pk_plex_accessibility(fca,
                                                    options->max_interaction_length,
                                                    DEFAULT_CUTOFF);
      vrna_fold_compound_free(fca);
    }

    /* obtain list of potential PK interactions */
    interactions = duplexfold_XS(fc,
                                 (accessibility) ? accessibility : (const int **)opening_energies,
                                 options->max_interaction_length,
                                 options->scoring_function,
                                 options->scoring_data);

    pk_penalty =
      ((double)options->scoring_function(NULL, 0, 0, 0, 0, options->scoring_data) / 100.);
    subopts = ((double)options->delta / 100.);

    if (vrna_heap_size(interactions) > 0) {
      results = vrna_heap_init(vrna_heap_size(interactions) + 2,
                               PKplex_heap_cmp,
                               NULL,
                               NULL,
                               NULL);

      /*  now we re-evaluate the energies and thereby prune the list of pre-results
       *  such that the re-avaluation process is not done too often.
       */

      constraint = (char *)vrna_alloc(sizeof(char) * (fc->length + 1));

      while ((hit_ptr = vrna_heap_pop(interactions))) {
        /*
         * simple check whether we believe that this structure might achieve
         * a net free energy within our boundaries
         * For that, we add up the interaction energy, the pk-free MFE
         * the PK penalty, and the minimal energy to unfold at least
         * one of the PK interaction sites. This for now seems a good
         * choice for a lower bound of the actual energy
         */
        best_e = hit_ptr->dGint +
                 mfe +
                 pk_penalty +
                 MIN2(hit_ptr->dG1, hit_ptr->dG2);

        if (best_e <= mfe_pk + subopts) {
          /* now for the exact evaluation of the structures energy incl. PKs */

          /* prepare the structure constraint for constrained folding */
          vrna_hc_init(fc);

          for (i = hit_ptr->start_5; i <= hit_ptr->end_5; i++)
            vrna_hc_add_up(fc, i, VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS);
          for (i = hit_ptr->start_3; i <= hit_ptr->end_3; i++)
            vrna_hc_add_up(fc, i, VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS);

          /* energy evaluation */
          constrainedEnergy = vrna_mfe(fc, constraint);

          if (options->scoring_function != &default_pk_plex_penalty) {
            short *pt = vrna_ptable(constraint);
            hit_ptr->dGpk = ((double)options->scoring_function((const short *)pt,
                                                               hit_ptr->start_5,
                                                               hit_ptr->end_5,
                                                               hit_ptr->start_3,
                                                               hit_ptr->end_3,
                                                               options->scoring_data) / 100.);
            free(pt);
          } else {
            hit_ptr->dGpk = pk_penalty;
          }

          /*
           *  compute net free energy, i.e.:
           *  interaction energy +
           *  constrained MFE +
           *  pseudo-knot penalty
           */
          hit_ptr->energy = hit_ptr->dGint +
                            constrainedEnergy +
                            hit_ptr->dGpk;

          /* check if this structure is worth keeping */
          if (hit_ptr->energy <= mfe_pk + subopts) {
            /* add pseudo-knot brackets to the structure */
            for (i = hit_ptr->start_5 - 1; i < hit_ptr->end_5; i++)
              if (hit_ptr->structure[i - hit_ptr->start_5 + 1] == '(')
                constraint[i] = '[';

            for (i = hit_ptr->start_3 - 1; i < hit_ptr->end_3; i++)
              if (hit_ptr->structure[i - hit_ptr->start_3 + 1 + 1 + 1 +
                                     hit_ptr->end_5 - hit_ptr->start_5] == ')')
                constraint[i] = ']';

            if (hit_ptr->energy < mfe_pk)
              mfe_pk = hit_ptr->energy;

            free(hit_ptr->structure);
            hit_ptr->structure = constraint;

            vrna_heap_insert(results, hit_ptr);

            constraint = (char *)vrna_alloc(sizeof(char) * (fc->length + 1));

            continue;
          }
        }

        free(hit_ptr->structure);
        free(hit_ptr);
      }

      free(constraint);
    }

    /* Add the PK-free MFE as potential result */
    mfe_entry = (vrna_pk_plex_t *)vrna_alloc(sizeof(vrna_pk_plex_t));

    mfe_entry->structure  = mfe_struct;
    mfe_entry->energy     = mfe;
    mfe_entry->start_5    = 0;

    if (results == NULL) {
      results = vrna_heap_init(1,
                               PKplex_heap_cmp,
                               NULL,
                               NULL,
                               NULL);
    }

    vrna_heap_insert(results, (void *)mfe_entry);

    /*
     * now go through the active hits again and filter those out that are above
     * the subopt threshold. This is necessary due to the fact that while we've
     * been processing the list once we always compared against the best known
     * pk-MFE. But this value might have changed during processing.
     */
    cnt   = 0;
    hits  = (vrna_pk_plex_t *)vrna_alloc(sizeof(vrna_pk_plex_t) * (vrna_heap_size(results) + 1));
    /* collect all final results */
    while ((ptr = vrna_heap_pop(results))) {
      if (ptr->energy > mfe_pk + subopts)
        break;

      hits[cnt++] = *ptr;
    }

    /* end of list marker */
    hits[cnt].structure = NULL;

    /* remove intermediate hits that didn't surpass the threshold */
    while ((ptr = vrna_heap_pop(results))) {
      free(ptr->structure);
      free(ptr);
    }

    /* cleanup the heap storages */
    vrna_heap_free(interactions);
    vrna_heap_free(results);

    if (opening_energies) {
      i = opening_energies[0][0];
      while (--i > -1)
        free(opening_energies[i]);
      free(opening_energies);
    }

    if (options != user_options)
      free(options);
  }

  return hits;
}


PUBLIC int **
vrna_pk_plex_accessibility(vrna_fold_compound_t *fc,
                           unsigned int unpaired,
                           double       cutoff)
{
  unsigned int          n, i, j;
  int                   **a = NULL, r;
  double                **pup, kT;
  plist                 *dpp = NULL;
  vrna_param_t          *P;
  vrna_md_t             *md;
  access_data           data;

  data.open_en = NULL;

  if (fc) {

    n   = fc->length;
    P   = fc->params;
    md  = &(P->model_details);

    data.kT       = (md->temperature + K0) * GASCONST / 1000.0;
    data.cutoff   = (cutoff > 0) ? cutoff : 0;
    data.open_en  = (int **)vrna_alloc(sizeof(int *) * (unpaired + 2));

    for (j = 0; j < unpaired + 2; j++) {
      data.open_en[j] = (int *)vrna_alloc(sizeof(int) * (n + 1));
      for (i = 0; i <= n; i++)
        data.open_en[j][i] = INF;
    }

    data.open_en[0][0] = unpaired + 2;

    r = vrna_probs_window(fc, unpaired, VRNA_PROBS_WINDOW_UP, &store_pU_callback, (void *)(&data));
  }

  return data.open_en;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE int
default_pk_plex_penalty(const short *pt,
                        int         i,
                        int         j,
                        int         k,
                        int         l,
                        void        *data)
{
  return ((default_data *)data)->penalty;
}


PRIVATE int ***
get_array(unsigned int  n,
          unsigned int  interaction_length)
{
  unsigned int  i, j;
  int           ***c3;

  c3 = (int ***)vrna_alloc(sizeof(int **) * interaction_length);

  for (i = 0; i < interaction_length; i++) {
    c3[i] = (int **)vrna_alloc(sizeof(int *) * n);

    for (j = 0; j < n; j++)
      c3[i][j] = (int *)vrna_alloc(sizeof(int) * interaction_length);
  }

  return c3;
}


PRIVATE void
reset_array(int           ***c3,
            unsigned int  n,
            unsigned int  interaction_length)
{
  unsigned int i, j, k;

  for (i = 0; i < interaction_length; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < interaction_length; k++)
        c3[i][j][k] = INF;
}


PRIVATE void
free_array(int          ***c3,
           unsigned int n,
           unsigned int interaction_length)
{
  unsigned int i, j;

  for (i = 0; i < interaction_length; i++) {
    for (j = 0; j < n; j++)
      free(c3[i][j]);
    free(c3[i]);
  }

  free(c3);
}


PRIVATE int
prepare_PKplex(vrna_fold_compound_t *fc)
{
  /* prepare Boltzmann factors if required */
  vrna_params_prepare(fc, VRNA_OPTION_MFE);

  /* prepare ptype array(s) */
  vrna_ptypes_prepare(fc, VRNA_OPTION_MFE);

  /* prepare hard constraints */
  vrna_hc_prepare(fc, VRNA_OPTION_MFE);

  /* prepare soft constraints data structure, if required */
  vrna_sc_prepare(fc, VRNA_OPTION_MFE);

  return 1;
}


PRIVATE void
store_pU_callback(double        *pU,
                  int           size,
                  int           k,
                  int           ulength,
                  unsigned int  type,
                  void          *data)
{
  int         i;
  double      kT, cutoff;
  access_data *d = (access_data *)data;

  if ((type & VRNA_PROBS_WINDOW_UP) && ((type & VRNA_ANY_LOOP) == VRNA_ANY_LOOP)) {
    kT      = d->kT;
    cutoff  = d->cutoff;

    for (i = 1; i <= size; i++) {
      if (pU[i] > cutoff)
        d->open_en[i][k] = rint(100 * (-log(pU[i])) * kT);
    }
  }
}


PRIVATE vrna_heap_t
duplexfold_XS(vrna_fold_compound_t        *fc,
              const int                   **access_s1,
              const int                   max_interaction_length,
              vrna_pk_plex_score_f scoring_function,
              void                        *scoring_data)
{
  char                      *struc;
  short                     *S, *S1, si, sk, sl, sp, sq;
  unsigned int              n, type, type2, type3;
  int                       ***c3, i, j, k, l, p, q, Emin, l_min, k_min, j_min, E,
                            tempK, i_pos_begin, j_pos_end, dGx, dGy, inter,
                            turn, penalty;

  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_heap_t               storage;
  vrna_pk_plex_t            *entry;

  vrna_hc_eval_f evaluate_ext;
  struct hc_ext_def_dat     hc_dat_local;

  struc   = NULL;
  n       = fc->length;
  S       = fc->sequence_encoding2;
  S1      = fc->sequence_encoding;
  P       = fc->params;
  md      = &(P->model_details);
  turn    = md->min_loop_size;
  hc      = fc->hc;
  penalty = scoring_function(NULL, 0, 0, 0, 0, scoring_data);

  storage = vrna_heap_init(128,
                           PKplex_heap_cmp,
                           NULL,
                           NULL,
                           NULL);

  evaluate_ext = prepare_hc_ext_def(fc, &hc_dat_local);

  c3 = get_array(n, max_interaction_length);

  if (n > turn + 1) {
    for (i = n - turn - 1; i > 0; i--) {
      Emin  = INF;
      j_min = 0;
      l_min = 0;
      k_min = 0;

      reset_array(c3, n, max_interaction_length);

      si = S1[i + 1];

      /* matrix starting values for (i,j)-basepairs */
      for (j = i + turn + 1; j <= n; j++) {
        if (evaluate_ext(i, j, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
          type                                      = md->pair[S[j]][S[i]];
          c3[max_interaction_length - 1][j - 1][0]  = vrna_E_exterior_stem(type,
                                                                      S1[j - 1],
                                                                      si,
                                                                      P);
        }
      }

      i_pos_begin = MAX2(0, i - max_interaction_length);

      /* fill matrix */
      for (k = i - 1; k > i_pos_begin; k--) {
        tempK = max_interaction_length - i + k - 1;
        sk    = S1[k + 1];
        for (l = i + turn + 1; l <= n; l++) {
          if (hc->mx[n * k + l] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
            type2 = md->pair[S[k]][S[l]];
            sl    = S1[l - 1];

            for (p = k + 1; (p <= i) && (p <= k + MAXLOOP + 1); p++) {
              sp = S1[p - 1];

              for (q = l - 1; (q >= i + turn + 1) && (q >= l - MAXLOOP - 1); q--) {
                if (p - k + l - q - 2 > MAXLOOP)
                  break;

                if (hc->mx[n * p + q] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) {
                  type3 = md->pair[S[q]][S[p]];
                  sq    = S1[q + 1];

                  E = vrna_E_internal(p - k - 1,
                                l - q - 1,
                                type2,
                                type3,
                                sk,
                                sl,
                                sp,
                                sq,
                                P);
                  for (j = MAX2(i + turn + 1, l - max_interaction_length + 1); j <= q; j++) {
                    if (hc->mx[n * i + j] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
                      type                    = md->pair[S[i]][S[j]];
                      c3[tempK][j - 1][l - j] =
                        MIN2(c3[tempK][j - 1][l - j],
                             c3[max_interaction_length - i + p - 1][j - 1][q - j] + E);
                    }
                  }
                }
              } /* next j */
            }   /* next q */
          }     /* next p */
        }       /* next l */
      }         /* next k */

      /* read out matrix minimum */
      for (k = i - 1; k > i_pos_begin; k--) {
        if (access_s1[i - k + 1][i] == INF)
          continue;

        int **c3k = c3[max_interaction_length - i + k - 1];
        int a1    = access_s1[i - k + 1][i];
        sk = (k > 1) ? S1[k - 1] : -1;

        for (j = i + turn + 1; j <= n; j++) {
          if (evaluate_ext(i, j, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
            j_pos_end = MIN2(n + 1, j + max_interaction_length);
            int *c3kj = c3k[j - 1];

            for (l = j + 1; l < j_pos_end; l++) {
              if (access_s1[l - j + 1][l] == INF)
                continue;

              if (evaluate_ext(k, l, k, l, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
                type2 = md->pair[S[k]][S[l]];
                sl    = (l < n) ? S1[l + 1] : -1;
                E     = c3kj[l - j] +
                        vrna_E_exterior_stem(type2, sk, sl, P) +
                        a1 +
                        access_s1[l - j + 1][l];

                if (E < Emin) {
                  Emin  = E;
                  k_min = k;
                  l_min = l;
                  j_min = j;
                }
              }
            }
          }
        }
      }

      /*
       *  Only consider hits where the interaction is more
       *  stable than the PK penalty, i.e. the net free
       *  energy is negative
       */
      if (Emin < -penalty) {
        struc = backtrack_XS(fc, k_min, l_min, i, j_min, max_interaction_length, c3);

        dGx   = access_s1[i - k_min + 1][i];
        dGy   = access_s1[l_min - j_min + 1][l_min];
        inter = Emin - dGx - dGy;

        entry             = (vrna_pk_plex_t *)vrna_alloc(sizeof(vrna_pk_plex_t));
        entry->start_5    = k_min;
        entry->end_5      = i;
        entry->start_3    = j_min;
        entry->end_3      = l_min;
        entry->energy     = (double)Emin * 0.01;
        entry->dG1        = (double)dGx * 0.01;
        entry->dG2        = (double)dGy * 0.01;
        entry->dGint      = (double)inter * 0.01;
        entry->structure  = struc;

        vrna_heap_insert(storage, entry);
      }
    }
  }

  free_array(c3, n, max_interaction_length);

  return storage;
}


PRIVATE char *
backtrack_XS(vrna_fold_compound_t *fc,
             int                  k,
             int                  l,
             const int            i,
             const int            j,
             const int            max_interaction_length,
             int                  ***c3)
{
  /* backtrack structure going backwards from i, and forwards from j
   * return structure in bracket notation with & as separator */
  short         *S, *S1;
  int           p, q, type, type2, E, traced, i0, j0, *rtype;
  char          *st1, *st2, *struc;
  vrna_param_t  *P;
  vrna_md_t     *md;

  S               = fc->sequence_encoding2;
  S1              = fc->sequence_encoding;
  P               = fc->params;
  md              = &(P->model_details);
  rtype           = &(md->rtype[0]);
  st1             = (char *)vrna_alloc(sizeof(char) * (i - k + 2));
  st1[i - k + 1]  = '\0';
  st2             = (char *)vrna_alloc(sizeof(char) * (l - j + 2));
  st2[l - j + 1]  = '\0';

  i0  = k;
  j0  = l;
  while (k <= i && l >= j) {
    E           = c3[max_interaction_length - i + k - 1][j - 1][l - j];
    traced      = 0;
    st1[k - i0] = '(';
    st2[l - j]  = ')';

    type = md->pair[S[k]][S[l]];

    if (!type) {
      vrna_log_error("backtrack failed in fold duplex bli");
      free(st1);
      free(st2);
      return NULL;
    }

    for (p = k + 1; p <= i; p++) {
      for (q = l - 1; q >= j; q--) {
        int LE;
        if (p - k + l - q - 2 > MAXLOOP)
          break;

        type2 = md->pair[S[q]][S[p]];

        if (!type2)
          continue;

        LE = vrna_E_internal(p - k - 1,
                       l - q - 1,
                       type,
                       type2,
                       S1[k + 1],
                       S1[l - 1],
                       S1[p - 1],
                       S1[q + 1],
                       P);
        if (E == c3[max_interaction_length - i + p - 1][j - 1][q - j] + LE) {
          traced  = 1;
          k       = p;
          l       = q;
          break;
        }
      }
      if (traced)
        break;
    }
    if (!traced) {
      E -= vrna_E_exterior_stem(rtype[type],
                           S1[l - 1],
                           S1[k + 1],
                           P);

      if (E != 0) {
        vrna_log_error("backtrack failed in fold duplex bal");
        free(st1);
        free(st2);
        return NULL;
      } else {
        break;
      }
    }
  }
  struc = (char *)vrna_alloc(sizeof(char) * (k - i0 + 1 + j0 - l + 1 + 2));

  for (p = 0; p <= i - i0; p++)
    if (!st1[p])
      st1[p] = '.';

  for (p = 0; p <= j0 - j; p++)
    if (!st2[p])
      st2[p] = '.';

  strcpy(struc, st1);
  strcat(struc, "&");
  strcat(struc, st2);
  free(st1);
  free(st2);
  return struc;
}


PRIVATE
int
PKplex_heap_cmp(const void  *a,
                const void  *b,
                void        *data)
{
  vrna_pk_plex_t  *p1 = (vrna_pk_plex_t *)a;
  vrna_pk_plex_t  *p2 = (vrna_pk_plex_t *)b;

  if (p1->energy > p2->energy)
    return 1;
  else if (p1->energy < p2->energy)
    return -1;

  return 0;
}


/*
 * ###########################################
 * # deprecated functions below              #
 *###########################################
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

PUBLIC dupVar *
PKLduplexfold_XS(const char *s1,
                 const int  **access_s1,
                 int        penalty,
                 int        max_interaction_length,
                 int        delta)
{
  vrna_fold_compound_t  *fc;
  dupVar                *hits;
  default_data          scoring_dat;
  size_t                cnt;
  vrna_heap_t           interactions;
  vrna_pk_plex_t        *entry;

  hits = NULL;

  if ((s1) &&
      (access_s1)) {
    fc = vrna_fold_compound(s1, NULL, VRNA_OPTION_DEFAULT);

    prepare_PKplex(fc);

    scoring_dat.penalty = -penalty;

    interactions = duplexfold_XS(fc,
                                 access_s1,
                                 max_interaction_length,
                                 &default_pk_plex_penalty,
                                 (void *)&scoring_dat);

    cnt   = 0;
    hits  = (dupVar *)vrna_alloc(sizeof(dupVar) * (vrna_heap_size(interactions) + 2));
    while ((entry = vrna_heap_pop(interactions))) {
      hits[cnt].structure = entry->structure;
      hits[cnt].tb        = entry->start_5;
      hits[cnt].te        = entry->end_5;
      hits[cnt].qb        = entry->start_3;
      hits[cnt].qe        = entry->end_3;
      hits[cnt].ddG       = entry->energy;
      hits[cnt].dG1       = entry->dG1;
      hits[cnt].dG2       = entry->dG2;
      hits[cnt].energy    = entry->dGint;
      hits[cnt].inactive  = 0;
      hits[cnt].processed = 0;
      free(entry);
      cnt++;
    }

    hits[cnt].inactive  = 1;
    hits[cnt].structure = NULL;

    vrna_heap_free(interactions);
    vrna_fold_compound_free(fc);
  }

  return hits;
}


#endif
