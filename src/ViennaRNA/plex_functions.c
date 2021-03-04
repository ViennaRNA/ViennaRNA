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
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/LPfold.h"
#include "ViennaRNA/mfe.h"
#include "ViennaRNA/datastructures/heap.h"

#include "ViennaRNA/loops/external_hc.inc"

#include "ViennaRNA/pk_plex.h"

#undef  MAXLOOP
#define MAXLOOP 10

#define DEFAULT_PENALTY             8.10
#define DEFAULT_DELTA               0
#define DEFAULT_INTERACTION_LENGTH  12

#define WITH_HEAP

struct vrna_pk_plex_option_s {
  unsigned int                delta;
  unsigned int                max_interaction_length;
  int                         min_penalty;
  vrna_callback_pk_plex_score *scoring_function;
  void                        *scoring_data;
};

typedef struct {
  int penalty;
} default_data;

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


#ifdef WITH_HEAP
PRIVATE vrna_heap_t
#else
PRIVATE vrna_pkplex_t *
#endif
duplexfold_XS(vrna_fold_compound_t        *fc,
              const int                   **access_s1,
              const int                   max_interaction_length,
              vrna_callback_pk_plex_score *scoring_function,
              void                        *scoring_data);


PRIVATE char *
backtrack_XS(vrna_fold_compound_t *fc,
             int        kk,
             int        ll,
             const int  ii,
             const int  jj,
             const int  max_interaction_length,
             int        ***c3);


PRIVATE int
prepare_PKplex(vrna_fold_compound_t *fc);


PRIVATE int
PlexHit_cmp(const void  *c1,
            const void  *c2);


PRIVATE int
PlexHit_cmp_energy(const void *c1,
                   const void *c2);


PRIVATE int
PlexHit_cmp_active_energy(const void *c1,
                          const void *c2);

PRIVATE
int PKplex_heap_cmp(const void *a,
                    const void *b,
                    void       *data);

PRIVATE
int PKplex_heap_result_cmp(const void *a,
                           const void *b,
                           void       *data);

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
                 int      pk_penalty)
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
                     vrna_callback_pk_plex_score  *scoring_function,
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



PUBLIC vrna_pkplex_t *
vrna_pk_plex(vrna_fold_compound_t *fc,
             const int            **accessibility,
             vrna_pk_plex_opt_t   user_options)
{
  char *mfe_struct, *constraint;
  int i;
  double mfe, mfe_pk, constrainedEnergy, pk_penalty, subopts;
  size_t  NumberOfHits;

  vrna_pkplex_t       *hits, *hit_ptr;
  default_data        scoring_dat;
  vrna_pk_plex_opt_t  options;
#ifdef WITH_HEAP
  vrna_heap_t         interactions, results;
#else
  vrna_pkplex_t       *interactions;
#endif

  hits  = NULL;
  results = NULL;


  if ((fc) &&
      (accessibility)) {
    mfe_struct = (char *)vrna_alloc(sizeof(char) * (fc->length + 1));

    mfe = mfe_pk = (double)vrna_mfe(fc, mfe_struct);

    options = (user_options) ? user_options : vrna_pk_plex_opt_defaults();

    /* apply simplified scoring function with constant penalty */
    if (!options->scoring_function) {
      scoring_dat.penalty       = options->min_penalty;
      options->scoring_function = &default_pk_plex_penalty;
      options->scoring_data     = (void *)&scoring_dat;
    }

    /* obtain list of potential PK interactions */
    interactions = duplexfold_XS(fc,
                         accessibility,
                         options->max_interaction_length,
                         options->scoring_function,
                         options->scoring_data);

    NumberOfHits = 0;

    pk_penalty = ((double)options->scoring_function(NULL, 0, 0, 0, 0, options->scoring_data) / 100.);
    subopts     = ((double)options->delta / 100.);

#ifdef WITH_HEAP
    if (vrna_heap_size(interactions) > 0) {
      results = vrna_heap_init(vrna_heap_size(interactions) + 2,
                                           PKplex_heap_result_cmp,
                                            NULL,
                                            NULL,
                                            NULL);
#else
    if (interactions) {
      for (hit_ptr = interactions; hit_ptr->structure; hit_ptr++)
        NumberOfHits++;

      /* first sort all the pre-results descending by their estimated energy */
      qsort(interactions, NumberOfHits, sizeof(vrna_pkplex_t), PlexHit_cmp);
#endif

      /*  now we re-evaluate the energies and thereby prune the list of pre-results
      *  such that the re-avaluation process is not done too often.
      */

      constraint  = (char *)vrna_alloc(sizeof(char) * (fc->length + 1));

#ifdef WITH_HEAP
      while (hit_ptr = vrna_heap_pop(interactions)) {
#else
      for (hit_ptr = interactions; hit_ptr->structure; hit_ptr++) {
#endif
        /*
        * simple check whether we believe that this structure might achieve
        * a net free energy within our boundaries
        * For that, we add up the interaction energy, the pk-free MFE
        * the PK penalty, and the minimal energy to unfold at least
        * one of the PK interaction sites. This for now seems a good
        * choice for a lower bound of the actual energy
        */
        double best_e = hit_ptr->energy +
                        mfe +
                        pk_penalty +
                        MIN2(hit_ptr->dG1, hit_ptr->dG2);

        if (best_e <= mfe_pk + subopts) {
          /* now for the exact evaluation of the structures energy incl. PKs */

          /* prepare the structure constraint for constrained folding */
          vrna_hc_init(fc);

          for (i = hit_ptr->tb; i <= hit_ptr->te; i++)
            vrna_hc_add_up(fc, i, VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS);
          for (i = hit_ptr->qb; i <= hit_ptr->qe; i++)
            vrna_hc_add_up(fc, i, VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS);

          /* energy evaluation */
          constrainedEnergy = vrna_mfe(fc, constraint);

          /* compute net free energy */
          hit_ptr->energy += constrainedEnergy +
                             pk_penalty;

          /* check if this structure is worth keeping */
          if (hit_ptr->energy <= mfe_pk + subopts) {
            /* add pseudo-knot brackets to the structure */
            for (i = hit_ptr->tb - 1; i < hit_ptr->te; i++)
              if (hit_ptr->structure[i - hit_ptr->tb + 1] == '(')
                constraint[i] = '[';

            for (i = hit_ptr->qb - 1; i < hit_ptr->qe; i++)
              if (hit_ptr->structure[i - hit_ptr->qb + 1 + 1 + 1 +
                                              hit_ptr->te - hit_ptr->tb] == ')')
                constraint[i] = ']';

            if (hit_ptr->energy < mfe_pk)
              mfe_pk = hit_ptr->energy;

            free(hit_ptr->structure);
            hit_ptr->structure = constraint;
            hit_ptr->processed = 1;
            hit_ptr->inactive  = 0;

#ifdef WITH_HEAP
            vrna_heap_insert(results, hit_ptr);
#endif
            constraint = (char *)vrna_alloc(sizeof(char) * (fc->length + 1));
            continue;
          }
        }

#ifdef WITH_HEAP
        free(hit_ptr->structure);
        free(hit_ptr);
#else
        hit_ptr->inactive = 1;
#endif
      }

      free(constraint);
    }

#ifdef WITH_HEAP
    vrna_pkplex_t *mfe_entry = (vrna_pkplex_t *)vrna_alloc(sizeof(vrna_pkplex_t));

    mfe_entry->structure = mfe_struct;
    mfe_entry->energy    = mfe;
    mfe_entry->tb        = 0;
    mfe_entry->inactive  = 0;

    if (results == NULL) {
      results = vrna_heap_init(1,
                               PKplex_heap_result_cmp,
                               NULL,
                               NULL,
                               NULL);
    }
    vrna_heap_insert(results, (void *)mfe_entry);

#else
    /* add the MFE structure without PKs to list */
    hits = (vrna_pkplex_t *)vrna_realloc(interactions, sizeof(vrna_pkplex_t) * (NumberOfHits + 2));
    hits[NumberOfHits].structure  = mfe_struct;
    hits[NumberOfHits].energy     = mfe;
    hits[NumberOfHits].tb         = 0;
    hits[NumberOfHits++].inactive = 0;
    hits[NumberOfHits].structure  = NULL;
    hits[NumberOfHits].energy     = (double)INF * 0.01;
    hits[NumberOfHits].inactive   = 1;
#endif

    /*
    * now go through the active hits again and filter those out that are above
    * the subopt threshold. This is necessary due to the fact that while we've
    * been processing the list once we always compared against the best known
    * pk-MFE. But this value might have changed during processing.
    */
#ifdef WITH_HEAP
    size_t  cnt = 0;
    vrna_pkplex_t *ptr;

    hits = (vrna_pkplex_t *)vrna_alloc(sizeof(vrna_pkplex_t) * (vrna_heap_size(results) + 1));
    /* collect all final results */
    while((ptr = vrna_heap_pop(results))) {
      if (ptr->energy > mfe_pk + subopts)
        break;

      hits[cnt++] = *ptr;
    }

    hits[cnt].inactive   = 1;
    hits[cnt].structure  = NULL;

    /* remove intermediate hits that didn't surpass the threshold */
    while((ptr = vrna_heap_pop(results))) {
      free(ptr->structure);
      free(ptr);
    }

    /* cleanup the heap storages */
    vrna_heap_free(interactions);
    vrna_heap_free(results);
#else
    for (hit_ptr = hits; hit_ptr->structure; hit_ptr++) {
      if ((!hit_ptr->inactive) &&
          (hit_ptr->energy > mfe_pk + subopts))
        hit_ptr->inactive = 1;
    }

    /* now sort the actual results again according to their energy */
    qsort(hits, NumberOfHits, sizeof(vrna_pkplex_t), PlexHit_cmp_active_energy);
#endif
    if (options != user_options)
      free(options);
  }

  return hits;
}

PUBLIC vrna_pkplex_t *
PKLduplexfold_XS(const char *s1,
                 const int  **access_s1,
                 int  penalty,
                 int  max_interaction_length,
                 int  delta)
{
  vrna_fold_compound_t  *fc;
  vrna_pkplex_t         *hits;
  default_data          scoring_dat;
#ifdef WITH_HEAP
  size_t                cnt;
  vrna_heap_t           interactions;
  vrna_pkplex_t         *entry;
#else
  vrna_pkplex_t         *interactions;
#endif

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

#ifdef WITH_HEAP
    cnt  = 0;
    hits = (vrna_pkplex_t *)vrna_alloc(sizeof(vrna_pkplex_t) * (vrna_heap_size(interactions) + 2));
    while((entry = vrna_heap_pop(interactions)))
      hits[cnt++] = *entry;

    hits[cnt].inactive   = 1;
    hits[cnt].structure  = NULL;
#else
    hits = interactions;
#endif

    vrna_fold_compound_free(fc);
  }

  return hits;
}


PUBLIC int **
vrna_pk_plex_accessibility(const char    *sequence,
                           unsigned int  unpaired,
                           double        cutoff)
{
  unsigned int  n, i, j;
  int           **a = NULL;
  double        **pup, kT;
  plist         *dpp = NULL;
  vrna_fold_compound_t  *fc;
  vrna_param_t  *P;
  vrna_md_t     *md;

  if (sequence) {
    fc = vrna_fold_compound(sequence, NULL, VRNA_OPTION_DEFAULT | VRNA_OPTION_WINDOW);

    n   = fc->length;
    P   = fc->params;
    md  = &(P->model_details);

    pup       = (double **)vrna_alloc((n + 1) * sizeof(double *));
    pup[0]    = (double *)vrna_alloc(sizeof(double));   /*I only need entry 0*/
    pup[0][0] = (double)unpaired;

    (void)pfl_fold(fc->sequence, n, n, cutoff, pup, &dpp, NULL, NULL);

    kT = (md->temperature + K0) * GASCONST / 1000.0;

    /* prepare the accesibility array */
    a = (int **)vrna_alloc(sizeof(int *) * (unpaired + 2));

    for (i = 0; i < unpaired + 2; i++)
      a[i] = (int *)vrna_alloc(sizeof(int) * (n + 1));

    for (i = 0; i <= n; i++)
      for (j = 0; j < unpaired + 2; j++)
        a[j][i] = INF;

    for (i = 1; i <= n; i++) {
      for (j = 1; j < unpaired + 1; j++)
        if (pup[i][j] > 0)
          a[j][i] = rint(100 * (-log(pup[i][j])) * kT);
    }

    a[0][0] = unpaired + 2;

    vrna_fold_compound_free(fc);
  }

  return a;
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
get_array(unsigned int n,
          unsigned int interaction_length)
{
  unsigned int  i, j, k;
  int           ***c3;

  c3 = (int ***)vrna_alloc(sizeof(int **) * (n));

  for (i = 0; i < n; i++) {
    c3[i] = (int **)vrna_alloc(sizeof(int *) * interaction_length);

    for (j = 0; j < interaction_length; j++)
      c3[i][j] = (int *)vrna_alloc(sizeof(int) * interaction_length);
  }

  return c3;
}


PRIVATE void
reset_array(int ***c3,
            unsigned int n,
            unsigned int interaction_length)
{
  unsigned int  i, j, k;

  for (i = 0; i < n; i++) {
    for (j = 0; j < interaction_length; j++) {
      for (k = 0; k < interaction_length; k++)
        c3[i][j][k] = INF;
    }
  }
}


PRIVATE void
free_array(int ***c3,
           unsigned int n,
           unsigned int interaction_length)
{
  unsigned int i, j;

  for (i = 0; i < n; i++) {
    for (j = 0; j < interaction_length; j++)
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
}


#ifdef WITH_HEAP
PRIVATE vrna_heap_t
#else
PRIVATE vrna_pkplex_t *
#endif
duplexfold_XS(vrna_fold_compound_t        *fc,
              const int                   **access_s1,
              const int                   max_interaction_length,
              vrna_callback_pk_plex_score *scoring_function,
              void                        *scoring_data)
{
  char          *struc;
  short         *S, *S1, si, sk, sl, sp, sq;
  size_t        storage_size, storage_fill;
  unsigned int  n, type, type2, type3;
  int           ***c3, i, j, k, l, p, q, Emin, l_min, k_min, j_min, E,
                tempK, *rtype, i_pos_begin, j_pos_end, dGx, dGy, inter,
                turn, penalty;

  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;
#ifdef WITH_HEAP
  vrna_heap_t   storage;
  vrna_pkplex_t *entry;
#else
  vrna_pkplex_t *storage;
#endif

  vrna_callback_hc_evaluate *evaluate_ext;
  struct default_data       hc_dat_local;

  struc   = NULL;
  n       = fc->length;
  S       = fc->sequence_encoding2;
  S1      = fc->sequence_encoding;
  P       = fc->params;
  md      = &(P->model_details);
  turn    = md->min_loop_size;
  rtype   = &(md->rtype[0]);
  hc      = fc->hc;
  penalty = scoring_function(NULL, 0, 0, 0, 0, scoring_data);

#ifdef WITH_HEAP
  storage = vrna_heap_init(128,
                           PKplex_heap_cmp,
                           NULL,
                           NULL,
                           NULL);
#else
  storage_size  = 64;
  storage_fill  = 0;
  storage       = (vrna_pkplex_t *)vrna_alloc(sizeof(vrna_pkplex_t) * storage_size);
#endif

  evaluate_ext  = prepare_hc_default(fc, &hc_dat_local);

  c3  = get_array(n, max_interaction_length);

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
          type = md->pair[S[j]][S[i]];
          c3[j - 1][max_interaction_length - 1][0] = vrna_E_ext_stem(type,
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
            /* again, why 9 less then the sequence length ? */
            type2 = md->pair[S[k]][S[l]];
            sl    = S1[l - 1];

            for (p = k + 1; (p <= i) && (p <= k + MAXLOOP + 1); p++) {
              sp  = S1[p - 1];

              for (q = l - 1; (q >= i + turn + 1) && (q >= l - MAXLOOP - 1); q--) {
                if (p - k + l - q - 2 > MAXLOOP)
                  break;

                if (hc->mx[n * p + q] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) {
                  type3 = md->pair[S[q]][S[p]];
                  sq    = S1[q + 1];

                  E = E_IntLoop(p - k - 1,
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
                      type = md->pair[S[i]][S[j]];
                      c3[j - 1][tempK][l - j] =
                          MIN2(c3[j - 1][tempK][l - j],
                               c3[j - 1][max_interaction_length - i + p - 1][q - j] + E);
                    }
                  }
                }
              } /* next j */
            }   /* next q */
          }     /* next p */
        }       /* next l */
      }         /* next k */

      /* read out matrix minimum */
      for (j = i + turn + 1; j <= n; j++) {
        if (evaluate_ext(i, j, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
          j_pos_end = MIN2(n + 1, j + max_interaction_length);
          for (k = i - 1; k > i_pos_begin; k--) {
            sk = (k > i_pos_begin + 1) ? S1[k - 1] : -1; /* should actually be k > 10 */
            for (l = j + 1; l < j_pos_end; l++) {
              if (evaluate_ext(k, l, k, l, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
                type2 = md->pair[S[k]][S[l]];
                sl  = (l < j_pos_end - 1) ? S1[l + 1] : -1; /* should actually be l < n - 10 */
                E   = c3[j - 1][max_interaction_length - i + k - 1][l - j] +
                      vrna_E_ext_stem(type2, sk, sl, P) +
                      access_s1[i - k + 1][i] +
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
          Only consider hits where the interaction is more
          stable than the PK penalty, i.e. the net free
          energy is negative
      */
      if (Emin < -penalty) { 
        struc = backtrack_XS(fc, k_min, l_min, i, j_min, max_interaction_length, c3);

        dGx   = access_s1[i - k_min + 1][i];
        dGy   = access_s1[l_min - j_min + 1][l_min];
        inter = Emin - dGx - dGy;

#ifdef WITH_HEAP
        entry = (vrna_pkplex_t *)vrna_alloc(sizeof(vrna_pkplex_t));
        entry->tb        = k_min;
        entry->te        = i;
        entry->qb        = j_min;
        entry->qe        = l_min;
        entry->ddG       = (double)Emin * 0.01;
        entry->dG1       = (double)dGx * 0.01;
        entry->dG2       = (double)dGy * 0.01;
        entry->energy    = (double)inter * 0.01;
        entry->structure = struc;
        entry->inactive  = 0;
        entry->processed = 0;

        vrna_heap_insert(storage, entry);
#else
        storage[storage_fill].tb        = k_min;
        storage[storage_fill].te        = i;
        storage[storage_fill].qb        = j_min;
        storage[storage_fill].qe        = l_min;
        storage[storage_fill].ddG       = (double)Emin * 0.01;
        storage[storage_fill].dG1       = (double)dGx * 0.01;
        storage[storage_fill].dG2       = (double)dGy * 0.01;
        storage[storage_fill].energy    = (double)inter * 0.01;
        storage[storage_fill].structure = struc;
        storage[storage_fill].inactive  = 0;
        storage[storage_fill].processed = 0;

        storage_fill++;

        if (storage_fill == storage_size - 1) {
          storage_size *= 1.4;
          storage       = (vrna_pkplex_t *)vrna_realloc(storage,
                                                        sizeof(vrna_pkplex_t) *
                                                        storage_size);
        }
#endif
      }
    }
  }

  free_array(c3, n, max_interaction_length);

#ifndef WITH_HEAP
  /* resize to space actually required */
  if (storage_fill > 0) {
    storage = (vrna_pkplex_t *)vrna_realloc(storage,
                                            sizeof(vrna_pkplex_t) *
                                            (storage_fill + 1));

    /* add end-of-list identifier */
    storage[storage_fill].structure = NULL;
    storage[storage_fill].inactive  = 1;
  } else {
    free(storage);
    storage = NULL;
  }
#endif

  return storage;
}


PRIVATE char *
backtrack_XS(vrna_fold_compound_t *fc,
             int        k,
             int        l,
             const int  i,
             const int  j,
             const int  max_interaction_length,
             int        ***c3)
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
    E           = c3[j - 1][max_interaction_length - i + k - 1][l - j];
    traced      = 0;
    st1[k - i0] = '(';
    st2[l - j]  = ')';

    type = md->pair[S[k]][S[l]];

    if (!type)
      vrna_message_error("backtrack failed in fold duplex bli");

    for (p = k + 1; p <= i; p++) {
      for (q = l - 1; q >= j; q--) {
        int LE;
        if (p - k + l - q - 2 > MAXLOOP)
          break;

        type2 = md->pair[S[q]][S[p]];

        if (!type2)
          continue;

        LE = E_IntLoop(p - k - 1,
                       l - q - 1,
                       type,
                       type2,
                       S1[k + 1],
                       S1[l - 1],
                       S1[p - 1],
                       S1[q + 1],
                       P);
        if (E == c3[j - 1][max_interaction_length - i + p - 1][q - j] + LE) {
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
      E -= vrna_E_ext_stem(rtype[type],
                           S1[l - 1],
                           S1[k + 1],
                           P);

      if (E != 0)
        vrna_message_error("backtrack failed in fold duplex bal");
      else
        break;
    }
  }
  struc = (char *)vrna_alloc(k - i0 + 1 + j0 - l + 1 + 2);

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


PRIVATE int
PlexHit_cmp(const void  *c1,
            const void  *c2)
{
  dupVar  *p1 = (dupVar *)c1;
  dupVar  *p2 = (dupVar *)c2;

  return p1->ddG >= p2->ddG;
}


PRIVATE int
PlexHit_cmp_energy(const void *c1,
                   const void *c2)
{
  dupVar  *p1 = (dupVar *)c1;
  dupVar  *p2 = (dupVar *)c2;

  if (p1->energy > p2->energy)
    return 1;
  else if (p1->energy < p2->energy)
    return -1;

  return 0;
}

PRIVATE int
PlexHit_cmp_active_energy(const void *c1,
                          const void *c2)
{
  dupVar  *p1 = (dupVar *)c1;
  dupVar  *p2 = (dupVar *)c2;

  if (p1->inactive > p2->inactive) {
    return 1;
  } else if (p1->inactive < p2->inactive) {
    return -1;
  } else if (!p2->inactive) {
    if (p1->energy < p2->energy)
      return -1;
    else if (p1->energy > p2->energy)
      return 1;
  }

  return 0;
}

PRIVATE
int PKplex_heap_cmp(const void *a,
                    const void *b,
                    void       *data)
{
  dupVar  *p1 = (dupVar *)a;
  dupVar  *p2 = (dupVar *)b;

  return p1->ddG >= p2->ddG;
}

PRIVATE
int PKplex_heap_result_cmp(const void *a,
                           const void *b,
                           void       *data)
{
  dupVar  *p1 = (dupVar *)a;
  dupVar  *p2 = (dupVar *)b;

  if (p1->energy > p2->energy)
    return 1;
  else if (p1->energy < p2->energy)
    return -1;

  return 0;
}
