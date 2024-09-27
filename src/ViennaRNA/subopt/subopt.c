/*
 * suboptimal folding - Stefan Wuchty, Walter Fontana & Ivo Hofacker
 *
 *                     Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#ifdef _WIN32
#ifdef __MINGW32__
#include <unistd.h>
#else
#include "ViennaRNA/intern/unistd_win.h"
#endif
#else
#include <unistd.h>
#endif

#include <ctype.h>
#include <string.h>
#include <math.h>
#include "ViennaRNA/fold.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/datastructures/lists.h"
#include "ViennaRNA/eval/structures.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/eval/hairpin.h"
#include "ViennaRNA/eval/internal.h"
#include "ViennaRNA/eval/multibranch.h"
#include "ViennaRNA/eval/gquad.h"
#include "ViennaRNA/cofold.h"
#include "ViennaRNA/subopt/gquad.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/subopt/wuchty.h"

#include "ViennaRNA/constraints/exterior_hc.inc"
#include "ViennaRNA/constraints/hairpin_hc.inc"
#include "ViennaRNA/constraints/internal_hc.inc"
#include "ViennaRNA/constraints/multibranch_hc.inc"

#include "ViennaRNA/constraints/exterior_sc.inc"
#include "ViennaRNA/constraints/hairpin_sc.inc"
#include "ViennaRNA/constraints/internal_sc.inc"
#include "ViennaRNA/constraints/multibranch_sc.inc"

/* hack */
#include "ViennaRNA/intern/color_output.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#define true              1
#define false             0

typedef struct {
  struct hc_ext_def_dat     hc_dat_ext;
  vrna_hc_eval_f hc_eval_ext;

  struct hc_hp_def_dat      hc_dat_hp;
  vrna_hc_eval_f hc_eval_hp;

  struct hc_int_def_dat     hc_dat_int;
  eval_hc                   hc_eval_int;

  struct hc_mb_def_dat      hc_dat_mb;
  vrna_hc_eval_f hc_eval_mb;

  struct sc_f5_dat          sc_dat_ext;
  struct sc_hp_dat          sc_dat_hp;
  struct sc_int_dat         sc_dat_int;
  struct sc_mb_dat          sc_dat_mb;
} constraint_helpers;

/**
 *  @brief  Sequence interval stack element used in subopt.c
 */
typedef struct INTERVAL {
  int           i;
  int           j;
  unsigned int  array_flag;
} INTERVAL;

typedef struct {
  char  *structure;
  LIST  *Intervals;
  int   partial_energy;
  int   is_duplex;
  /* int best_energy;   */ /* best attainable energy */
} STATE;

typedef struct {
  LIST  *Intervals;
  LIST  *Stack;
  int   nopush;
} subopt_env;


struct old_subopt_dat {
  unsigned long           max_sol;
  unsigned long           n_sol;
  vrna_subopt_solution_t  *SolutionList;
  FILE                    *fp;
  unsigned int            strands;
  unsigned int            *strand_start;
};

/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */
PUBLIC int    subopt_sorted = 0;      /* output sorted by energy */
PUBLIC int    density_of_states[MAXDOS + 1];
PUBLIC double print_energy = 9999;    /* printing threshold for use with logML */

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */

/* some backward compatibility stuff */
PRIVATE int                   backward_compat           = 0;
PRIVATE vrna_fold_compound_t  *backward_compat_compound = NULL;

#ifdef _OPENMP

#pragma omp threadprivate(backward_compat_compound, backward_compat)

#endif

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE void
init_constraint_helpers(vrna_fold_compound_t  *fc,
                        constraint_helpers    *d);


PRIVATE void
free_constraint_helpers(constraint_helpers *d);


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

PRIVATE vrna_subopt_solution_t *
wrap_subopt(char          *seq,
            char          *structure,
            vrna_param_t  *parameters,
            int           delta,
            int           is_constrained,
            int           is_circular,
            FILE          *fp);


#endif

PRIVATE void
make_pair(int   i,
          int   j,
          STATE *state);


/* mark a gquadruplex in the resulting dot-bracket structure */
PRIVATE void
make_gquad(unsigned int i,
           unsigned int n,
           unsigned int L,
           unsigned int l[3],
           STATE        *state);


PRIVATE INTERVAL *
make_interval(int           i,
              int           j,
              unsigned int  array_flag);


PRIVATE STATE *
make_state(LIST *Intervals,
           char *structure,
           int  partial_energy,
           int  is_duplex,
           int  length);


PRIVATE STATE *
copy_state(STATE *state);


PRIVATE void
print_state(STATE *state);


PRIVATE void
print_stack(LIST *list) VRNA_UNUSED;


PRIVATE LIST *
make_list(void);


PRIVATE void
push(LIST *list,
     void *data);


PRIVATE void *
pop(LIST *list);


PRIVATE int
best_attainable_energy(vrna_fold_compound_t *fc,
                       STATE                *state);


PRIVATE void
scan_interval(vrna_fold_compound_t  *fc,
              int                   i,
              int                   j,
              unsigned int          array_flag,
              int                   threshold,
              STATE                 *state,
              subopt_env            *env,
              constraint_helpers    *constraints_dat);


PRIVATE INLINE void
scan_mb(vrna_fold_compound_t  *fc,
        int                   i,
        int                   j,
        int                   threshold,
        STATE                 *state,
        subopt_env            *env,
        constraint_helpers    *constraints_dat);


PRIVATE INLINE void
scan_m1(vrna_fold_compound_t  *fc,
        int                   i,
        int                   j,
        unsigned int          array_flag,
        int                   threshold,
        STATE                 *state,
        subopt_env            *env,
        constraint_helpers    *constraints_dat);


PRIVATE INLINE void
scan_m2(vrna_fold_compound_t  *fc,
        unsigned int          i,
        unsigned int          j,
        int                   threshold,
        STATE                 *state,
        subopt_env            *env,
        constraint_helpers    *constraints_dat);


PRIVATE INLINE void
scan_pair(vrna_fold_compound_t  *fc,
          int                   i,
          int                   j,
          int                   threshold,
          STATE                 *state,
          subopt_env            *env,
          constraint_helpers    *constraints_dat);


PRIVATE INLINE void
scan_ext(vrna_fold_compound_t *fc,
         int                  i,
         int                  j,
         int                  threshold,
         STATE                *state,
         subopt_env           *env,
         constraint_helpers   *constraints_dat);


PRIVATE INLINE void
scan_circular(vrna_fold_compound_t  *fc,
              unsigned int          i,
              unsigned int          j,
              int                   threshold,
              STATE                 *state,
              subopt_env            *env,
              constraint_helpers    *constraints_dat);


PRIVATE INLINE void
scan_fms5(vrna_fold_compound_t  *fc,
          unsigned int          i,
          unsigned int          strand,
          int                   threshold,
          STATE                 *state,
          subopt_env            *env,
          constraint_helpers    *constraints_dat);


PRIVATE INLINE void
scan_fms3(vrna_fold_compound_t  *fc,
          unsigned int          i,
          unsigned int          strand,
          int                   threshold,
          STATE                 *state,
          subopt_env            *env,
          constraint_helpers    *constraints_dat);


PRIVATE INLINE void
scan_gquad(vrna_fold_compound_t *fc,
           int                  i,
           int                  j,
           int                  threshold,
           STATE                *state,
           subopt_env           *env,
           constraint_helpers   *constraints_dat);


PRIVATE void
free_interval_node(INTERVAL *node);


PRIVATE void
free_state_node(STATE *node);


PRIVATE void
push_back(LIST  *Stack,
          STATE *state);


PRIVATE char *
get_structure(STATE *state);


PRIVATE int
compare(const void  *a,
        const void  *b);


PRIVATE int
compare_en(const void *a,
           const void *b);


PRIVATE void
make_output(vrna_subopt_solution_t  *SL,
            unsigned int            strands,
            const unsigned int      *strand_start,
            int                     compressed,
            FILE                    *fp);


PRIVATE void
repeat(vrna_fold_compound_t *fc,
       unsigned int         i,
       unsigned int         j,
       STATE                *state,
       int                  part_energy,
       int                  temp_energy,
       int                  best_energy,
       int                  threshold,
       subopt_env           *env,
       constraint_helpers   *constraints_dat);


PRIVATE void
repeat_gquad(vrna_fold_compound_t *fc,
             int                  i,
             int                  j,
             STATE                *state,
             int                  part_energy,
             int                  temp_energy,
             int                  best_energy,
             int                  threshold,
             subopt_env           *env,
             constraint_helpers   *constraints_dat);


PRIVATE void
old_subopt_print(const char *structure,
                 float      energy,
                 void       *data);


PRIVATE void
old_subopt_store(const char *structure,
                 float      energy,
                 void       *data);


PRIVATE void
old_subopt_store_compressed(const char  *structure,
                            float       energy,
                            void        *data);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_subopt_solution_t *
vrna_subopt(vrna_fold_compound_t  *fc,
            int                   delta,
            int                   sorted,
            FILE                  *fp)
{
  struct old_subopt_dat data;
  vrna_subopt_result_f  cb;

  data.SolutionList = NULL;
  data.max_sol      = 128;
  data.n_sol        = 0;
  data.fp           = fp;
  data.strands      = fc->strands;
  data.strand_start = fc->strand_start;

  if (fc) {
    /* SolutionList stores the suboptimal structures found */

    data.SolutionList =
      (vrna_subopt_solution_t *)vrna_alloc(data.max_sol * sizeof(vrna_subopt_solution_t));

    /* end initialize ------------------------------------------------------- */

    if (fp) {
      float min_en;
      char  *SeQ = NULL, *energies = NULL;
      min_en = vrna_mfe(fc, NULL);
      char  *tmp_seq = strdup(fc->sequence);

      if (fc->strands == 1) {
        SeQ = tmp_seq;
      } else {
        for (unsigned int i = 1; i < fc->strands; i++) {
          SeQ = vrna_cut_point_insert(tmp_seq, (int)fc->strand_start[i] + (i - 1));
          free(tmp_seq);
          tmp_seq = SeQ;
        }
      }

      energies  = vrna_strdup_printf(" %6.2f %6.2f", min_en, (float)delta / 100.);
      print_structure(fp, SeQ, energies);
      free(SeQ);
      free(energies);

      vrna_mx_mfe_free(fc);
    }

    cb = old_subopt_store;

    if (fp) {
      if (!sorted)
        cb = old_subopt_print;
      else if (!(fc->params->model_details.gquad))
        cb = old_subopt_store_compressed;
    }

    /* call subopt() */
    vrna_subopt_cb(fc, delta, cb, (void *)&data);

    if (sorted) {
      /* sort structures by energy */
      if (data.n_sol > 0) {
        int (*compare_fun)(const void *a,
                           const void *b);

        switch (sorted) {
          case VRNA_SORT_BY_ENERGY_ASC:
            compare_fun = compare_en;
            break;

          default: /* a.k.a. VRNA_SORT_BY_ENERGY_LEXICOGRAPHIC_ASC */
            compare_fun = compare;
            break;
        }

        qsort(data.SolutionList, data.n_sol - 1, sizeof(vrna_subopt_solution_t), compare_fun);
      }

      if (fp)
        make_output(data.SolutionList,
                    fc->strands,
                    fc->strand_start,
                    !(fc->params->model_details.gquad),
                    fp);
    }

    if (fp) {
      /* we've printed everything -- free solutions */
      vrna_subopt_solution_t *sol;
      for (sol = data.SolutionList; sol->structure != NULL; sol++)
        free(sol->structure);
      free(data.SolutionList);
      data.SolutionList = NULL;
    }
  }

  return data.SolutionList;
}


PUBLIC void
vrna_subopt_cb(vrna_fold_compound_t *fc,
               int                  delta,
               vrna_subopt_result_f cb,
               void                 *data)
{
  subopt_env          *env;
  STATE               *state;
  INTERVAL            *interval;
  int                 maxlevel, count, partial_energy, old_dangles, logML, dangle_model, length,
                      circular,
                      threshold;
  double              structure_energy, min_en, eprint;
  char                *struc, *structure;
  float               correction;
  vrna_param_t        *P;
  vrna_md_t           *md;
  int                 minimal_energy;
  int                 Fc;
  int                 *f5;
  constraint_helpers  constraints_dat;

  vrna_fold_compound_prepare(fc, VRNA_OPTION_MFE);

  length  = fc->length;
  P       = fc->params;
  md      = &(P->model_details);

  /*
   * do mfe folding to get fill arrays and get ground state energy
   * in case dangles is neither 0 or 2, set dangles=2 while folding
   */

  circular    = md->circ;
  logML       = md->logML;
  old_dangles = dangle_model = md->dangles;

  if (md->uniq_ML != 1) /* failsafe mechanism to enforce valid fM1 array */
    md->uniq_ML = 1;

  /* temporarily set dangles to 2 if necessary */
  if ((md->dangles != 0) && (md->dangles != 2))
    md->dangles = 2;

  struc = (char *)vrna_alloc(sizeof(char) * (length + 1));

  min_en = vrna_mfe(fc, struc);

  /* restore dangle model */
  md->dangles = old_dangles;

  /* re-evaluate in case we're using logML etc */
  min_en  = vrna_eval_structure(fc, struc);
  f5      = fc->matrices->f5;
  Fc      = fc->matrices->Fc;

  free(struc);
  eprint = print_energy + min_en;

  correction = (min_en < 0) ? -0.1 : 0.1;

  /* Initialize ------------------------------------------------------------ */
  init_constraint_helpers(fc, &constraints_dat);

  maxlevel        = 0;
  count           = 0;
  partial_energy  = 0;

  /* Initialize the stack ------------------------------------------------- */

  minimal_energy  = (circular) ? Fc : f5[length];
  threshold       = minimal_energy + delta;
  if (threshold >= INF) {
    vrna_log_warning("Energy range too high, limiting to reasonable value");
    threshold = INF - EMAX;
  }

  /* init env data structure */
  env             = (subopt_env *)vrna_alloc(sizeof(subopt_env));
  env->Stack      = NULL;
  env->nopush     = true;
  env->Stack      = make_list();                      /* anchor */
  env->Intervals  = make_list();                      /* initial state: */
  interval        = make_interval(1, length, VRNA_MX_FLAG_F5);      /* interval [1,length, F5 array] */
  push(env->Intervals, interval);
  env->nopush = false;
  state       = make_state(env->Intervals, NULL, partial_energy, 0, length);
  /* state->best_energy = minimal_energy; */
  push(env->Stack, state);
  env->nopush = false;

  /* end initialize ------------------------------------------------------- */


  while (1) {
    /* forever, til nothing remains on stack */

    maxlevel = (env->Stack->count > maxlevel ? env->Stack->count : maxlevel);

    if (LST_EMPTY(env->Stack)) {
      /*
       * we are done! clean up and quit
       * fprintf(stderr, "maxlevel: %d\n", maxlevel);
       */

      lst_kill(env->Stack, free_state_node);

      cb(NULL, 0, data);   /* NULL (last time to call callback function */

      break;
    }

    /* pop the last element ---------------------------------------------- */

    state = pop(env->Stack);                       /* current state to work with */

    if (LST_EMPTY(state->Intervals)) {
      int e;
      /* state has no intervals left: we got a solution */

      count++;
      structure         = get_structure(state);
      structure_energy  = state->partial_energy / 100.;

#ifdef CHECK_ENERGY
      structure_energy = vrna_eval_structure(fc, structure);

      if (!logML) {
        if ((double)(state->partial_energy / 100.) != structure_energy) {
          vrna_log_error("%s %6.2f %6.2f",
                             structure,
                             state->partial_energy / 100.,
                             structure_energy);
          exit(1);
        }
      }

#endif
      if (logML || (dangle_model == 1) || (dangle_model == 3)) /* recalc energy */
        structure_energy = vrna_eval_structure(fc, structure);

      e = (int)((structure_energy - min_en) * 10. - correction); /* avoid rounding errors */
      if (e > MAXDOS)
        e = MAXDOS;

      density_of_states[e]++;

      if (structure_energy <= eprint) {
        char  *outstruct = NULL;
        char  *tmp_struct = strdup(structure);

        if (fc->strands == 1) {
          outstruct = tmp_struct;
        } else {
          for (unsigned int i = 1; i < fc->strands; i++) {
            outstruct = vrna_cut_point_insert(tmp_struct, (int)fc->strand_start[i] + (i - 1));
            free(tmp_struct);
            tmp_struct = outstruct;
          }
        }
        cb((const char *)tmp_struct, structure_energy, data);
        free(tmp_struct);
      }

      free(structure);
    } else {
      /* get (and remove) next interval of state to analyze */

      interval = pop(state->Intervals);
      scan_interval(fc,
                    interval->i,
                    interval->j,
                    interval->array_flag,
                    threshold,
                    state, env,
                    &constraints_dat);

      free_interval_node(interval);        /* free the current interval */
    }

    free_state_node(state);                     /* free the current state */
  } /* end of while (1) */

  /* cleanup memory */
  free_constraint_helpers(&constraints_dat);

  free(env);
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE void
init_constraint_helpers(vrna_fold_compound_t  *fc,
                        constraint_helpers    *d)
{
  /* hard constraints first */
  d->hc_eval_ext  = prepare_hc_ext_def(fc, &(d->hc_dat_ext));
  d->hc_eval_hp   = prepare_hc_hp_def(fc, &(d->hc_dat_hp));
  d->hc_eval_int  = prepare_hc_int_def(fc, &(d->hc_dat_int));
  d->hc_eval_mb   = prepare_hc_mb_def(fc, &(d->hc_dat_mb));

  init_sc_f5(fc, &(d->sc_dat_ext));
  init_sc_hp(fc, &(d->sc_dat_hp));
  init_sc_int(fc, &(d->sc_dat_int));
  init_sc_mb(fc, &(d->sc_dat_mb));
}


PRIVATE void
free_constraint_helpers(constraint_helpers *d)
{
  /* currently only required for comparative folding soft constraints, but here for consistency reasons */
  free_sc_f5(&(d->sc_dat_ext));
  free_sc_hp(&(d->sc_dat_hp));
  free_sc_int(&(d->sc_dat_int));
  free_sc_mb(&(d->sc_dat_mb));
}


/*
 * ---------------------------------------------------------------------------
 * List routines--------------------------------------------------------------
 *---------------------------------------------------------------------------
 */
PRIVATE void
make_pair(int   i,
          int   j,
          STATE *state)
{
  state->structure[i - 1] = '(';
  state->structure[j - 1] = ')';
}


PRIVATE void
make_gquad(unsigned int i,
           unsigned int L,
           unsigned int n,
           unsigned int l[3],
           STATE        *state)
{
  vrna_db_insert_gq(state->structure, i, L, l, n);
}


PRIVATE INTERVAL *
make_interval(int           i,
              int           j,
              unsigned int  array_flag)
{
  INTERVAL *interval;

  interval              = lst_newnode(sizeof(INTERVAL));
  interval->i           = i;
  interval->j           = j;
  interval->array_flag  = array_flag;
  return interval;
}


PRIVATE void
free_interval_node(INTERVAL *node)
{
  lst_freenode(node);
}


PRIVATE void
free_state_node(STATE *node)
{
  free(node->structure);
  if (node->Intervals)
    lst_kill(node->Intervals, lst_freenode);

  lst_freenode(node);
}


PRIVATE STATE *
make_state(LIST *Intervals,
           char *structure,
           int  partial_energy,
           int  is_duplex VRNA_UNUSED,
           int  length)
{
  STATE *state;

  state = lst_newnode(sizeof(STATE));

  if (Intervals)
    state->Intervals = Intervals;
  else
    state->Intervals = lst_init();

  if (structure) {
    state->structure = structure;
  } else {
    int i;
    state->structure = (char *)vrna_alloc(length + 1);
    for (i = 0; i < length; i++)
      state->structure[i] = '.';
  }

  state->partial_energy = partial_energy;

  return state;
}


PRIVATE STATE *
copy_state(STATE *state)
{
  STATE     *new_state;
  void      *after;
  INTERVAL  *new_interval, *next;

  new_state                 = lst_newnode(sizeof(STATE));
  new_state->Intervals      = lst_init();
  new_state->partial_energy = state->partial_energy;
  /* new_state->best_energy = state->best_energy; */

  if (state->Intervals->count) {
    after = LST_HEAD(new_state->Intervals);
    for (next = lst_first(state->Intervals); next; next = lst_next(next)) {
      new_interval  = lst_newnode(sizeof(INTERVAL));
      *new_interval = *next;
      lst_insertafter(new_state->Intervals, new_interval, after);
      after = new_interval;
    }
  }

  new_state->structure = strdup(state->structure);
  if (!new_state->structure) {
    vrna_log_error("out of memory");
    return NULL;
  }

  return new_state;
}


/*@unused @*/ PRIVATE void
print_state(STATE *state)
{
  INTERVAL *next;

  if (state->Intervals->count) {
    printf("%d intervals:\n", state->Intervals->count);
    for (next = lst_first(state->Intervals); next; next = lst_next(next))
      printf("[%d,%d],%u ", next->i, next->j, next->array_flag);
    printf("\n");
  }

  printf("partial structure: %s\n", state->structure);
  printf("\n");
  printf(" partial_energy: %d\n", state->partial_energy);
  /* printf(" best_energy: %d\n", state->best_energy); */
  (void)fflush(stdout);
}


/*@unused @*/ PRIVATE void
print_stack(LIST *list)
{
  void *rec;

  printf("================\n");
  printf("%d states\n", list->count);
  for (rec = lst_first(list); rec; rec = lst_next(rec)) {
    printf("state-----------\n");
    print_state(rec);
  }
  printf("================\n");
}


PRIVATE LIST *
make_list(void)
{
  return lst_init();
}


PRIVATE void
push(LIST *list,
     void *data)
{
  lst_insertafter(list, data, LST_HEAD(list));
}


/*
 * PRIVATE void
 * push_stack(STATE *state) { */ /* keep the stack sorted by energy
 *   STATE *after, *next;
 *   nopush = false;
 *   next = after = LST_HEAD(Stack);
 *   while ( next = lst_next(next)) {
 *     if ( next->best_energy >= state->best_energy ) break;
 *     after = next;
 *   }
 *   lst_insertafter(Stack, state, after);
 * }
 */
PRIVATE void *
pop(LIST *list)
{
  void *data;

  data = lst_deletenext(list, LST_HEAD(list));
  return data;
}


/*
 * ---------------------------------------------------------------------------
 * auxiliary routines---------------------------------------------------------
 *---------------------------------------------------------------------------
 */
PRIVATE int
best_attainable_energy(vrna_fold_compound_t *fc,
                       STATE                *state)
{
  /* evaluation of best possible energy attainable within remaining intervals */

  register int  sum;
  INTERVAL      *next;
  vrna_md_t     *md;
  vrna_mx_mfe_t *matrices;
  int           *indx;

  md        = &(fc->params->model_details);
  matrices  = fc->matrices;
  indx      = fc->jindx;

  sum = state->partial_energy;  /* energy of already found elements */

  for (next = lst_first(state->Intervals); next; next = lst_next(next)) {
    if (next->array_flag == VRNA_MX_FLAG_F5)
      sum += (md->circ) ? matrices->Fc : matrices->f5[next->j];
    else if (next->array_flag == VRNA_MX_FLAG_M)
      sum += matrices->fML[indx[next->j] + next->i];
    else if (next->array_flag == VRNA_MX_FLAG_M2)
      sum += matrices->fM2_real[indx[next->j] + next->i];
    else if (next->array_flag == VRNA_MX_FLAG_C)
      sum += matrices->c[indx[next->j] + next->i];
    else if (next->array_flag == VRNA_MX_FLAG_M1)
      sum += matrices->fM1[indx[next->j] + next->i];
    else if (next->array_flag == VRNA_MX_FLAG_MS5)
      sum += matrices->fms5[next->j][next->i];
    else if (next->array_flag == VRNA_MX_FLAG_MS3)
      sum += matrices->fms3[next->j][next->i];
    else if (next->array_flag == VRNA_MX_FLAG_G)
#ifndef VRNA_DISABLE_C11_FEATURES
      sum += vrna_smx_csr_get(matrices->c_gq, next->i, next->j, INF);
#else
      sum += vrna_smx_csr_int_get(matrices->c_gq, next->i, next->j, INF);
#endif
  }

  return sum;
}


PRIVATE void
push_back(LIST  *Stack,
          STATE *state)
{
  push(Stack, copy_state(state));
  return;
}


PRIVATE char *
get_structure(STATE *state)
{
  char *structure;

  structure = strdup(state->structure);
  return structure;
}


PRIVATE int
compare(const void  *a,
        const void  *b)
{
  if (((vrna_subopt_solution_t *)a)->energy > ((vrna_subopt_solution_t *)b)->energy)
    return 1;

  if (((vrna_subopt_solution_t *)a)->energy < ((vrna_subopt_solution_t *)b)->energy)
    return -1;

  return strcmp(((vrna_subopt_solution_t *)a)->structure,
                ((vrna_subopt_solution_t *)b)->structure);
}


PRIVATE int
compare_en(const void *a,
           const void *b)
{
  if (((vrna_subopt_solution_t *)a)->energy > ((vrna_subopt_solution_t *)b)->energy)
    return 1;

  if (((vrna_subopt_solution_t *)a)->energy < ((vrna_subopt_solution_t *)b)->energy)
    return -1;

  return 0;
}


PRIVATE void
make_output(vrna_subopt_solution_t  *SL,
            unsigned int            strands,
            const unsigned int      *strand_start,
            int                     compressed,
            FILE                    *fp)                  /* prints stuff */
{
  vrna_subopt_solution_t *sol;

  for (sol = SL; sol->structure != NULL; sol++) {
    char  *e_string = vrna_strdup_printf(" %6.2f", sol->energy);
    char  *ss       = (compressed) ? vrna_db_unpack(sol->structure) : strdup(sol->structure);
    char  *s        = ss;

    if (strands > 1) {
      for (unsigned int i = 1; i < strands; i++) {
        s = vrna_cut_point_insert(ss, (int)strand_start[i] + (i - 1));
        free(ss);
        ss = s;
      }
    }

    print_structure(fp, s, e_string);
    free(s);
    free(e_string);
  }
}


PRIVATE STATE *
derive_new_state(int          i,
                 int          j,
                 STATE        *s,
                 int          e,
                 unsigned int flag)
{
  STATE     *s_new  = copy_state(s);
  INTERVAL  *ival   = make_interval(i, j, flag);

  push(s_new->Intervals, ival);

  s_new->partial_energy += e;

  return s_new;
}


PRIVATE void
fork_state(int          i,
           int          j,
           STATE        *s,
           int          e,
           unsigned int flag,
           subopt_env   *env)
{
  STATE *s_new = derive_new_state(i, j, s, e, flag);

  push(env->Stack, s_new);
  env->nopush = false;
}


PRIVATE void
fork_int_state(int        i,
               int        j,
               int        p,
               int        q,
               STATE      *s,
               int        e,
               subopt_env *env)
{
  STATE *s_new = derive_new_state(p, q, s, e, VRNA_MX_FLAG_C);

  make_pair(i, j, s_new);
  make_pair(p, q, s_new);
  push(env->Stack, s_new);
  env->nopush = false;
}


PRIVATE void
fork_state_pair(int         i,
                int         j,
                STATE       *s,
                int         e,
                subopt_env  *env)
{
  STATE *new_state;

  new_state = copy_state(s);
  make_pair(i, j, new_state);
  new_state->partial_energy += e;
  push(env->Stack, new_state);
  env->nopush = false;
}


PRIVATE void
fork_two_states_pair(int          i,
                     int          j,
                     int          k,
                     STATE        *s,
                     int          e,
                     unsigned int flag1,
                     unsigned int flag2,
                     subopt_env   *env)
{
  INTERVAL  *interval1, *interval2;
  STATE     *new_state;

  new_state = copy_state(s);
  interval1 = make_interval(i + 1, k - 1, flag1);
  interval2 = make_interval(k, j - 1, flag2);
  if (k - i < j - k) {
    /* push larger interval first */
    push(new_state->Intervals, interval1);
    push(new_state->Intervals, interval2);
  } else {
    push(new_state->Intervals, interval2);
    push(new_state->Intervals, interval1);
  }

  make_pair(i, j, new_state);
  new_state->partial_energy += e;

  push(env->Stack, new_state);
  env->nopush = false;
}


PRIVATE void
fork_state_pair_interval(int          i,
                         int          j,
                         int          k,
                         int          l,
                         STATE        *s,
                         int          e,
                         unsigned int flag,
                         subopt_env   *env)
{
  INTERVAL  *interval;
  STATE     *new_state;

  new_state = copy_state(s);
  interval  = make_interval(k, l, flag);
  push(new_state->Intervals, interval);

  make_pair(i, j, new_state);
  new_state->partial_energy += e;

  push(env->Stack, new_state);
  env->nopush = false;
}


PRIVATE void
fork_two_states_pair_ms(int         i,
                        int         j,
                        int         sn1,
                        int         sn2,
                        STATE       *s,
                        int         e,
                        subopt_env  *env)
{
  INTERVAL  *interval1, *interval2;
  STATE     *new_state;

  new_state = copy_state(s);
  interval1 = make_interval(i + 1, sn1, VRNA_MX_FLAG_MS5);
  interval2 = make_interval(j - 1, sn2, VRNA_MX_FLAG_MS3);
  push(new_state->Intervals, interval1);
  push(new_state->Intervals, interval2);

  make_pair(i, j, new_state);
  new_state->partial_energy += e;

  push(env->Stack, new_state);
  env->nopush = false;
}


PRIVATE void
fork_two_states(int           i,
                int           j,
                int           p,
                int           q,
                STATE         *s,
                int           e,
                unsigned int  flag1,
                unsigned int  flag2,
                subopt_env    *env)
{
  INTERVAL  *interval1, *interval2;
  STATE     *new_state;

  new_state = copy_state(s);
  interval1 = make_interval(i, j, flag1);
  interval2 = make_interval(p, q, flag2);

  if ((j - i) < (q - p)) {
    push(new_state->Intervals, interval1);
    push(new_state->Intervals, interval2);
  } else {
    push(new_state->Intervals, interval2);
    push(new_state->Intervals, interval1);
  }

  new_state->partial_energy += e;

  push(env->Stack, new_state);
  env->nopush = false;
}


PRIVATE void
scan_interval(vrna_fold_compound_t  *fc,
              int                   i,
              int                   j,
              unsigned int          array_flag,
              int                   threshold,
              STATE                 *state,
              subopt_env            *env,
              constraint_helpers    *constraints_dat)
{
  /* real backtrack routine */

  env->nopush = true;

  switch (array_flag) {
    case VRNA_MX_FLAG_F5:
      scan_ext(fc, i, j, threshold, state, env, constraints_dat);
      break;

    case VRNA_MX_FLAG_M:
      scan_mb(fc, i, j, threshold, state, env, constraints_dat);
    /* fall through */

    case VRNA_MX_FLAG_M1:
      scan_m1(fc, i, j, array_flag, threshold, state, env, constraints_dat);
      break;

    case VRNA_MX_FLAG_M2:
      scan_m2(fc, i, j, threshold, state, env, constraints_dat);
      break;

    case VRNA_MX_FLAG_C:
      scan_pair(fc, i, j, threshold, state, env, constraints_dat);
      return;

    case VRNA_MX_FLAG_MS5:
      scan_fms5(fc, i, j, threshold, state, env, constraints_dat);
      break;

    case VRNA_MX_FLAG_MS3:
      scan_fms3(fc, i, j, threshold, state, env, constraints_dat);
      break;

    case VRNA_MX_FLAG_G:
      scan_gquad(fc, i, j, threshold, state, env, constraints_dat);
      return;
  }

  if (env->nopush) {
    push_back(env->Stack, state);
    env->nopush = false;
  }
}


PRIVATE INLINE void
scan_mb(vrna_fold_compound_t  *fc,
        int                   i,
        int                   j,
        int                   threshold,
        STATE                 *state,
        subopt_env            *env,
        constraint_helpers    *constraints_dat)
{
  char                      *ptype;
  short                     *S1, s5, s3;
  unsigned int              *sn, *so;
  int                       k, type, dangle_model, element_energy, best_energy, e_gq, *c, *fML,
                            *indx, with_gquad, stopp, k1j;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  struct hc_mb_def_dat      *hc_dat;
  vrna_hc_eval_f evaluate;
  struct sc_mb_dat          *sc_dat;
  sc_mb_red_cb              sc_red_stem;
  sc_mb_red_cb              sc_decomp_ml;
  vrna_smx_csr(int)         *c_gq;

  STATE                     *temp_state;

  sn    = fc->strand_number;
  so    = fc->strand_order;
  indx  = fc->jindx;
  ptype = fc->ptype;
  S1    = fc->sequence_encoding;
  P     = fc->params;
  md    = &(P->model_details);

  dangle_model  = md->dangles;
  with_gquad    = md->gquad;

  c     = fc->matrices->c;
  fML   = fc->matrices->fML;
  c_gq  = fc->matrices->c_gq;

  hc_dat    = &(constraints_dat->hc_dat_mb);
  evaluate  = constraints_dat->hc_eval_mb;

  sc_dat        = &(constraints_dat->sc_dat_mb);
  sc_red_stem   = constraints_dat->sc_dat_mb.red_stem;
  sc_decomp_ml  = constraints_dat->sc_dat_mb.decomp_ml;

  best_energy = best_attainable_energy(fc, state);  /* .. on remaining intervals */

  if ((j < i + 1) &&
      (sn[i] == so[j])) {
    if (env->nopush) {
      push_back(env->Stack, state);
      env->nopush = false;
    }

    return;
  }

  if ((sn[i - 1] == sn[i]) && (sn[j] == sn[j + 1])) {
    /*backtrack in FML only if multiloop is possible*/
    for (k = i + 1; k <= j - 1; k++) {
      /* Multiloop decomposition if i,j contains more than 1 stack */

      if (with_gquad) {
#ifndef VRNA_DISABLE_C11_FEATURES
        int e_gq = vrna_smx_csr_get(c_gq, k + 1, j, INF);
#else
        int e_gq = vrna_smx_csr_int_get(c_gq, k + 1, j, INF);
#endif
        if ((sn[k] == sn[k + 1]) &&
            (fML[indx[k] + i] != INF) &&
            (e_gq != INF)) {
          element_energy = vrna_E_multibranch_stem(0, -1, -1, P);

          if (fML[indx[k] + i] + e_gq + element_energy + best_energy <= threshold) {
            temp_state  = derive_new_state(i, k, state, 0, VRNA_MX_FLAG_M);
            env->nopush = false;
            repeat_gquad(fc,
                         k + 1,
                         j,
                         temp_state,
                         element_energy,
                         fML[indx[k] + i],
                         best_energy,
                         threshold,
                         env,
                         constraints_dat);
            free_state_node(temp_state);
          }
        }
      }

      k1j = indx[j] + k + 1;

      if ((evaluate(i, j, k, k + 1, VRNA_DECOMP_ML_ML_STEM, hc_dat)) &&
          (fML[indx[k] + i] != INF) &&
          (c[k1j] != INF)) {
        type = vrna_get_ptype(k1j, ptype);

        switch (dangle_model) {
          case 0:
            s5 = s3 = -1;
            break;

          default:
            s5  = (sn[i - 1] == sn[i]) ? S1[k] : -1;
            s3  = (sn[j] == sn[j + 1]) ? S1[j + 1] : -1;
            break;
        }

        element_energy = vrna_E_multibranch_stem(type, s5, s3, P);

        if (sc_decomp_ml)
          element_energy += sc_decomp_ml(i, j, k, k + 1, sc_dat);

        if (sc_red_stem)
          element_energy += sc_red_stem(k + 1, j, k + 1, j, sc_dat);

        if (fML[indx[k] + i] + c[k1j] + element_energy + best_energy <= threshold) {
          temp_state  = derive_new_state(i, k, state, 0, VRNA_MX_FLAG_M);
          env->nopush = false;
          repeat(fc,
                 k + 1,
                 j,
                 temp_state,
                 element_energy,
                 fML[indx[k] + i],
                 best_energy,
                 threshold,
                 env,
                 constraints_dat);
          free_state_node(temp_state);
        }
      }
    }
  }

  stopp = j - 1;

  int up = 1;

  for (k = i; k <= stopp; k++, up++) {
    k1j = indx[j] + k + 1;

    /* Multiloop decomposition if i,j contains only 1 stack */
    if (with_gquad) {
#ifndef VRNA_DISABLE_C11_FEATURES
      e_gq = vrna_smx_csr_get(c_gq, k + 1, j, INF);
#else
      e_gq = vrna_smx_csr_int_get(c_gq, k + 1, j, INF);
#endif
      if ((e_gq != INF) &&
          (sn[i] == sn[j])) {
        element_energy = vrna_E_multibranch_stem(0, -1, -1, P) + P->MLbase * up;

        if (sc_red_stem)
          element_energy += sc_red_stem(i, j, k + 1, j, sc_dat);

        if (e_gq + element_energy + best_energy <= threshold) {
          repeat_gquad(fc,
                       k + 1,
                       j,
                       state,
                       element_energy,
                       0,
                       best_energy,
                       threshold,
                       env,
                       constraints_dat);
        }
      }
    }

    if (evaluate(i, j, k + 1, j, VRNA_DECOMP_ML_STEM, hc_dat)) {
      if (c[k1j] != INF) {
        type = vrna_get_ptype(k1j, ptype);

        switch (dangle_model) {
          case 0:
            s5 = s3 = -1;
            break;
          default:
            s5  = (sn[k - 1] == sn[k]) ? S1[k] : -1;
            s3  = (sn[j] == sn[j + 1]) ? S1[j + 1] : -1;
            break;
        }

        element_energy = vrna_E_multibranch_stem(type, s5, s3, P);

        element_energy += P->MLbase * up;

        if (sc_red_stem)
          element_energy += sc_red_stem(i, j, k + 1, j, sc_dat);

        if (c[k1j] + element_energy + best_energy <= threshold) {
          repeat(fc,
                 k + 1,
                 j,
                 state,
                 element_energy,
                 0,
                 best_energy,
                 threshold,
                 env,
                 constraints_dat);
        }
      }
    }
  }
}


PRIVATE INLINE void
scan_m1(vrna_fold_compound_t  *fc,
        int                   i,
        int                   j,
        unsigned int          array_flag,
        int                   threshold,
        STATE                 *state,
        subopt_env            *env,
        constraint_helpers    *constraints_dat)
{
  char                      *ptype;
  short                     *S1;
  unsigned int              *sn, *so;
  int                       fi, cij, e_gq, ij, type, dangle_model, element_energy, best_energy,
                            *c, *fML, *fM1, length, *indx, circular, with_gquad;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  struct hc_mb_def_dat      *hc_dat;
  vrna_hc_eval_f evaluate;
  struct sc_mb_dat          *sc_dat;
  sc_mb_red_cb              sc_red_stem;
  sc_mb_red_cb              sc_red_ml;
  vrna_smx_csr(int)         *c_gq;

  length  = fc->length;
  sn      = fc->strand_number;
  so      = fc->strand_order;
  indx    = fc->jindx;
  ptype   = fc->ptype;
  S1      = fc->sequence_encoding;
  P       = fc->params;
  md      = &(P->model_details);

  dangle_model  = md->dangles;
  circular      = md->circ;
  with_gquad    = md->gquad;

  c     = fc->matrices->c;
  fML   = fc->matrices->fML;
  fM1   = fc->matrices->fM1;
  c_gq  = fc->matrices->c_gq;

  hc_dat    = &(constraints_dat->hc_dat_mb);
  evaluate  = constraints_dat->hc_eval_mb;

  sc_dat      = &(constraints_dat->sc_dat_mb);
  sc_red_stem = constraints_dat->sc_dat_mb.red_stem;
  sc_red_ml   = constraints_dat->sc_dat_mb.red_ml;

  best_energy = best_attainable_energy(fc, state);  /* .. on remaining intervals */

  if ((j < i + 1) &&
      (sn[i] == so[j])) {
    if (env->nopush) {
      push_back(env->Stack, state);
      env->nopush = false;
    }

    return;
  }

  ij = indx[j] + i;

  if ((evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_ML, hc_dat)) &&
      (((array_flag == VRNA_MX_FLAG_M1) && (fM1[indx[j - 1] + i] != INF)) ||
       (fML[indx[j - 1] + i] != INF))) {
    element_energy = P->MLbase;

    if (sc_red_ml)
      element_energy += sc_red_ml(i, j, i, j - 1, sc_dat);

    if (array_flag == VRNA_MX_FLAG_M1)
      fi = element_energy +
           fM1[indx[j - 1] + i];
    else
      fi = element_energy +
           fML[indx[j - 1] + i];

    if (fi + best_energy <= threshold)
      fork_state(i, j - 1, state, element_energy, array_flag, env);
  }

  if (evaluate(i, j, i, j, VRNA_DECOMP_ML_STEM, hc_dat)) {
    /* i,j may pair */
    cij = c[ij];
    if (cij != INF) {
      type = vrna_get_ptype(ij, ptype);

      switch (dangle_model) {
        case 0:
          element_energy = vrna_E_multibranch_stem(type, -1, -1, P);
          break;
        default:
          element_energy = vrna_E_multibranch_stem(type,
                                    (((i > 1) && (sn[i - 1] == sn[i])) || circular) ? S1[i - 1] : -1,
                                    (((j < length) && (sn[j] == sn[j + 1])) || circular)  ? S1[j + 1] : -1,
                                    P);
          break;
      }

      if (sc_red_stem)
        element_energy += sc_red_stem(i, j, i, j, sc_dat);

      cij += element_energy;

      if (cij + best_energy <= threshold) {
        repeat(fc,
               i,
               j,
               state,
               element_energy,
               0,
               best_energy,
               threshold,
               env,
               constraints_dat);
      }
    }
  } else if (with_gquad) {
#ifndef VRNA_DISABLE_C11_FEATURES
    e_gq = vrna_smx_csr_get(c_gq, i, j, INF);
#else
    e_gq = vrna_smx_csr_int_get(c_gq, i, j, INF);
#endif
    if (e_gq != INF) {
      element_energy = vrna_E_multibranch_stem(0, -1, -1, P);

      if (sc_red_stem)
        element_energy += sc_red_stem(i, j, i, j, sc_dat);

      if (e_gq + element_energy + best_energy <= threshold) {
        repeat_gquad(fc,
                     i,
                     j,
                     state,
                     element_energy,
                     0,
                     best_energy,
                     threshold,
                     env,
                     constraints_dat);
      }
    }
  }

  return;
}


PRIVATE INLINE void
scan_m2(vrna_fold_compound_t  *fc,
        unsigned int          i,
        unsigned int          j,
        int                   threshold,
        STATE                 *state,
        subopt_env            *env,
        constraint_helpers    *constraints_dat)
{
  char                      *ptype;
  short                     *S1;
  unsigned int              *sn, *so, k, turn, type, dangle_model, length, circular, with_gquad;
  int                       fi, cij, ckj, e_gq, ij, kj, element_energy, best_energy,
                            *c, *fML, *fM2_real, *indx;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  struct hc_mb_def_dat      *hc_dat;
  vrna_hc_eval_f evaluate;
  struct sc_mb_dat          *sc_dat;
  sc_mb_red_cb              sc_red_stem;
  sc_mb_red_cb              sc_red_ml;
  sc_mb_red_cb              sc_decomp_ml;
  vrna_smx_csr(int)         *c_gq;
  STATE                     *temp_state;

  length  = fc->length;
  sn      = fc->strand_number;
  so      = fc->strand_order;
  indx    = fc->jindx;
  ptype   = fc->ptype;
  S1      = fc->sequence_encoding;
  P       = fc->params;
  md      = &(P->model_details);
  turn    = md->min_loop_size;

  dangle_model  = md->dangles;
  circular      = md->circ;
  with_gquad    = md->gquad;

  c         = fc->matrices->c;
  fML       = fc->matrices->fML;
  fM2_real  = fc->matrices->fM2_real;
  c_gq      = fc->matrices->c_gq;

  hc_dat    = &(constraints_dat->hc_dat_mb);
  evaluate  = constraints_dat->hc_eval_mb;

  sc_dat      = &(constraints_dat->sc_dat_mb);
  sc_decomp_ml = constraints_dat->sc_dat_mb.decomp_ml;
  sc_red_stem = constraints_dat->sc_dat_mb.red_stem;
  sc_red_ml   = constraints_dat->sc_dat_mb.red_ml;

  best_energy = best_attainable_energy(fc, state);  /* .. on remaining intervals */

  if ((j < i + 1) &&
      (sn[i] == so[j])) {
    if (env->nopush) {
      push_back(env->Stack, state);
      env->nopush = false;
    }

    return;
  }

  ij = indx[j] + i;

  /* case 1: j is unpaired */
  if ((evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_ML, hc_dat)) &&
      (fM2_real[indx[j - 1] + i] != INF)) {
    element_energy = P->MLbase;

    if (sc_red_ml)
      element_energy += sc_red_ml(i, j, i, j - 1, sc_dat);

    fi = element_energy +
         fM2_real[indx[j - 1] + i];

    if (fi + best_energy <= threshold)
      fork_state(i, j - 1, state, element_energy, VRNA_MX_FLAG_M2, env);
  }

  /* case 2: j pairs with some k */
  for (k = i + turn + 2; k + turn < j; k++) {
    if ((evaluate(i, j, k - 1, k, VRNA_DECOMP_ML_ML_STEM, hc_dat)) &&
        (fML[indx[k - 1] + i] != INF) &&
        (c[indx[j] + k] != INF)) {
      kj    = indx[j] + k;
      ckj   = c[kj];
      type  = vrna_get_ptype(kj, ptype);

      switch (dangle_model) {
        case 0:
          element_energy = vrna_E_multibranch_stem(type, -1, -1, P);
          break;
        default:
          element_energy = vrna_E_multibranch_stem(type,
                                    (((k > 1) && (sn[k - 1] == sn[k])) || circular) ? S1[k - 1] : -1,
                                    (((j < length) && (sn[j] == sn[j + 1])) || circular)  ? S1[j + 1] : -1,
                                    P);
          break;
      }

      if (sc_decomp_ml)
        element_energy += sc_decomp_ml(i, j, k - 1, k, sc_dat);
      if (sc_red_stem)
        element_energy += sc_red_stem(k, j, k, j, sc_dat);

      if (fML[indx[k - 1] + i] + ckj + element_energy + best_energy <= threshold) {
        temp_state  = derive_new_state(i, k - 1, state, 0, VRNA_MX_FLAG_M);
        env->nopush = false;
        repeat(fc,
               k,
               j,
               temp_state,
               element_energy,
               fML[indx[k - 1] + i],
               best_energy,
               threshold,
               env,
               constraints_dat);
        free_state_node(temp_state);
      }
    }
  }

  if (with_gquad) {
    /* case 3: gquad from k to j */
    for (k = i + turn + 2; k + VRNA_GQUAD_MIN_BOX_SIZE - 1 <= j; k++) {
#ifndef VRNA_DISABLE_C11_FEATURES
      e_gq = vrna_smx_csr_get(c_gq, k, j, INF);
#else
      e_gq = vrna_smx_csr_int_get(c_gq, k, j, INF);
#endif

      if ((evaluate(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc_dat)) &&
          (fML[indx[k - 1] + i] != INF) &&
          (e_gq != INF)) {
        element_energy = vrna_E_multibranch_stem(0, -1, -1, P);

        if (sc_decomp_ml)
          element_energy += sc_decomp_ml(i, j, k - 1, k, sc_dat);
        if (sc_red_stem)
          element_energy += sc_red_stem(k, j, k, j, sc_dat);

        if (fML[indx[k - 1] + i] + e_gq + element_energy + best_energy <= threshold) {
          temp_state  = derive_new_state(i, k - 1, state, 0, VRNA_MX_FLAG_M);
          env->nopush = false;
          repeat_gquad(fc,
                 k,
                 j,
                 temp_state,
                 element_energy,
                 fML[indx[k - 1] + i],
                 best_energy,
                 threshold,
                 env,
                 constraints_dat);
          free_state_node(temp_state);
        }
      }
    }
  }

  return;
}


PRIVATE INLINE void
scan_pair(vrna_fold_compound_t  *fc,
          int                   i,
          int                   j,
          int                   threshold,
          STATE                 *state,
          subopt_env            *env,
          constraint_helpers    *constraints_dat)
{
  unsigned int  *sn, noLP;
  int           best_energy;

  sn          = fc->strand_number;
  noLP        = fc->params->model_details.noLP;
  best_energy = best_attainable_energy(fc, state);  /* .. on remaining intervals */

  if ((j < i + 1) &&
      (sn[i] == sn[j])) {
    if (env->nopush) {
      push_back(env->Stack, state);
      env->nopush = false;
    }

    return;
  }

  repeat(fc,
         i,
         j,
         state,
         0,
         0,
         best_energy,
         threshold,
         env,
         constraints_dat);

  if (env->nopush)
    if (!noLP)
      vrna_log_warning("%d,%d\nOops, no solution in repeat!", i, j);
}


PRIVATE INLINE void
scan_ext(vrna_fold_compound_t *fc,
         int                  i,
         int                  j,
         int                  threshold,
         STATE                *state,
         subopt_env           *env,
         constraint_helpers   *constraints_dat)
{
  char                      *ptype;
  short                     *S1, s5, s3;
  unsigned int              *sn, *so;
  int                       k, type, dangle_model, element_energy, best_energy, e_gq, *f5, *c,
                            length, *indx, circular, with_gquad, kj, tmp_en;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_sc_t                 *sc;
  struct hc_ext_def_dat     *hc_dat;
  vrna_hc_eval_f evaluate;
  struct sc_f5_dat          *sc_dat;
  sc_f5_cb                  sc_red_ext;
  sc_f5_cb                  sc_red_stem;
  sc_f5_cb                  sc_decomp_stem;

  STATE                     *temp_state;
  vrna_smx_csr(int)         *c_gq;

  length  = fc->length;
  sn      = fc->strand_number;
  so      = fc->strand_order;
  indx    = fc->jindx;
  ptype   = fc->ptype;
  S1      = fc->sequence_encoding;
  P       = fc->params;
  md      = &(P->model_details);
  sc      = fc->sc;

  dangle_model  = md->dangles;
  circular      = md->circ;
  with_gquad    = md->gquad;

  f5    = fc->matrices->f5;
  c     = fc->matrices->c;
  c_gq  = fc->matrices->c_gq;

  if (circular) {
    scan_circular(fc, i, j, threshold, state, env, constraints_dat);
    return;
  }

  hc_dat    = &(constraints_dat->hc_dat_ext);
  evaluate  = constraints_dat->hc_eval_ext;

  sc_dat          = &(constraints_dat->sc_dat_ext);
  sc_red_ext      = sc_dat->red_ext5;
  sc_red_stem     = sc_dat->red_stem5;
  sc_decomp_stem  = sc_dat->decomp_stem5;

  best_energy = best_attainable_energy(fc, state);  /* .. on remaining intervals */

  if (i > 1) {
    vrna_log_error("Error while backtracking!");
    return;
  }

  if ((j < i + 1) &&
      (sn[i] == so[j])) {
    /*
     * minimal structure element
     * do not forget to add f5[j], since it may contain pseudo energies from soft constraining
     */
    state->partial_energy += f5[j];

    if (env->nopush) {
      push_back(env->Stack, state);
      env->nopush = false;
    }

    return;
  }

  if ((evaluate(1, j, 1, j - 1, VRNA_DECOMP_EXT_EXT, hc_dat)) &&
      (f5[j - 1] != INF)) {
    tmp_en = 0;

    if (sc_red_ext)
      tmp_en += sc_red_ext(j, 1, j - 1, sc_dat);

    if (f5[j - 1] + tmp_en + best_energy <= threshold)
      /* no basepair, nibbling of 3'-end */
      fork_state(i, j - 1, state, tmp_en, VRNA_MX_FLAG_F5, env);
  }

  for (k = j - 1; k > 1; k--) {
    kj = indx[j] + k;

    if (with_gquad) {
#ifndef VRNA_DISABLE_C11_FEATURES
      e_gq = vrna_smx_csr_get(c_gq, k, j, INF);
#else
      e_gq = vrna_smx_csr_int_get(c_gq, k, j, INF);
#endif
      if ((sn[k - 1] == sn[j]) &&
          (f5[k - 1] != INF) &&
          (e_gq != INF)) {
        element_energy = 0;

        if (sc_decomp_stem)
          element_energy += sc_decomp_stem(j, k - 1, k, sc_dat);

        if (f5[k - 1] + e_gq + element_energy + best_energy <= threshold) {
          temp_state  = derive_new_state(1, k - 1, state, 0, VRNA_MX_FLAG_F5);
          env->nopush = false;
          /* backtrace the quadruplex */
          repeat_gquad(fc,
                       k,
                       j,
                       temp_state,
                       element_energy,
                       f5[k - 1],
                       best_energy,
                       threshold,
                       env,
                       constraints_dat);
          free_state_node(temp_state);
        }
      }
    }

    if ((evaluate(1, j, k - 1, k, VRNA_DECOMP_EXT_EXT_STEM, hc_dat)) &&
        (f5[k - 1] != INF) &&
        (c[kj] != INF)) {
      type = vrna_get_ptype(kj, ptype);

      /* k and j pair */
      switch (dangle_model) {
        case 0:
          s5 = s3 = -1;
          break;
        default:
          s5  = (sn[k - 1] == sn[k]) ? S1[k - 1] : -1;
          s3  = ((j < length) && (sn[j] == sn[j + 1])) ? S1[j + 1] : -1;
          break;
      }

      element_energy = vrna_E_exterior_stem(type, s5, s3, P);

      if (sc_decomp_stem)
        element_energy += sc_decomp_stem(j, k - 1, k, sc_dat);

      if (f5[k - 1] + c[kj] + element_energy + best_energy <= threshold) {
        temp_state  = derive_new_state(1, k - 1, state, 0, VRNA_MX_FLAG_F5);
        env->nopush = false;
        repeat(fc,
               k,
               j,
               temp_state,
               element_energy,
               f5[k - 1],
               best_energy,
               threshold,
               env,
               constraints_dat);
        free_state_node(temp_state);
      }
    }
  }

  kj = indx[j] + 1;

  if (with_gquad) {
#ifndef VRNA_DISABLE_C11_FEATURES
    e_gq = vrna_smx_csr_get(c_gq, 1, j, INF);
#else
    e_gq = vrna_smx_csr_int_get(c_gq, 1, j, INF);
#endif
    if ((sn[1] == sn[j]) &&
        (e_gq != INF)) {
      element_energy = 0;

      if (sc_red_stem)
        element_energy += sc_red_stem(j, 1, j, sc_dat);

      if (e_gq + element_energy + best_energy <= threshold) {
        /* backtrace the quadruplex */
        repeat_gquad(fc,
                     1,
                     j,
                     state,
                     element_energy,
                     0,
                     best_energy,
                     threshold,
                     env,
                     constraints_dat);
      }
    }
  }

  if ((evaluate(1, j, 1, j, VRNA_DECOMP_EXT_STEM, hc_dat)) &&
      (c[kj] != INF)) {
    type  = vrna_get_ptype(kj, ptype);
    s5    = -1;

    switch (dangle_model) {
      case 0:
        s3 = -1;
        break;
      default:
        s3 = (j < length) && (sn[j] == sn[j + 1]) ? S1[j + 1] : -1;
        break;
    }

    element_energy = vrna_E_exterior_stem(type, s5, s3, P);

    if (sc)
      if (sc->f)
        element_energy += sc->f(1, j, 1, j, VRNA_DECOMP_EXT_STEM, sc->data);

    if (c[kj] + element_energy + best_energy <= threshold) {
      repeat(fc,
             1,
             j,
             state,
             element_energy,
             0,
             best_energy,
             threshold,
             env,
             constraints_dat);
    }
  }
}


PRIVATE INLINE void
scan_circular(vrna_fold_compound_t  *fc,
              unsigned int          i,
              unsigned int          j,
              int                   threshold,
              STATE                 *state,
              subopt_env            *env,
              constraint_helpers    *constraints_dat)
{
  unsigned char             *hard_constraints, eval;
  char                      *ptype;
  short                     *S, *S1;
  unsigned int              k, l, p, q, u1, qmin, u2, type_2, n_seq, length, turn, type;
  int                       e, e_part, tmp_en, best_energy, *c, *fML, *fM1, Fc, FcH,
                            FcI, FcM, *fM2, *fM2_real, *fM1_new, *indx, *rtype, kl, tmpE,
                            tmpE2;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  struct hc_ext_def_dat     *hc_dat_ext;
  struct hc_int_def_dat     *hc_dat_int;
  struct hc_mb_def_dat      *hc_dat_mb;
  vrna_hc_eval_f evaluate_ext;
  eval_hc                   evaluate_int;
  vrna_hc_eval_f evaluate_mb;

  struct sc_int_dat         *sc_dat_int;
  struct sc_mb_dat          *sc_dat_mb;
  sc_int_cb                 sc_int_pair_ext;
  sc_mb_red_cb              sc_mb_decomp_ml;
  sc_mb_red_cb              sc_mb_red_ml;

  STATE                     *new_state;
  INTERVAL                  *new_interval;

  length  = fc->length;
  n_seq   = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  indx    = fc->jindx;
  ptype   = fc->ptype;
  S       = fc->sequence_encoding2;
  S1      = fc->sequence_encoding;
  P       = fc->params;
  md      = &(P->model_details);
  rtype   = &(md->rtype[0]);

  turn = md->min_loop_size;

  c   = fc->matrices->c;
  fML = fc->matrices->fML;
  fM1 = fc->matrices->fM1;
  Fc  = fc->matrices->Fc;
  FcH = fc->matrices->FcH;
  FcI = fc->matrices->FcI;
  FcM = fc->matrices->FcM;
  fM2 = fc->matrices->fM2;
  fM2_real = fc->matrices->fM2_real;
  fM1_new   = fc->matrices->fM1_new;

  hc                = fc->hc;
  hard_constraints  = hc->mx;

  sc = fc->sc;

  hc_dat_ext    = &(constraints_dat->hc_dat_ext);
  hc_dat_int    = &(constraints_dat->hc_dat_int);
  hc_dat_mb     = &(constraints_dat->hc_dat_mb);
  evaluate_ext  = constraints_dat->hc_eval_ext;
  evaluate_int  = constraints_dat->hc_eval_int;
  evaluate_mb   = constraints_dat->hc_eval_mb;

  sc_dat_int      = &(constraints_dat->sc_dat_int);
  sc_dat_mb       = &(constraints_dat->sc_dat_mb);
  sc_int_pair_ext = constraints_dat->sc_dat_int.pair_ext;
  sc_mb_decomp_ml = constraints_dat->sc_dat_mb.decomp_ml;
  sc_mb_red_ml    = constraints_dat->sc_dat_mb.red_ml;

  best_energy = best_attainable_energy(fc, state);  /* .. on remaining intervals */

  if ((i != 1) ||
      (j != length)) {
    vrna_log_error("Error while backtracking!");
    return;
  }

  if (j < i + turn + 1) {
    /* minimal structure element */
    state->partial_energy += Fc;

    if (env->nopush) {
      push_back(env->Stack, state);
      env->nopush = false;
    }

    return;
  }

  /*
   * if we've done everything right, we will never reach this point more than once
   * right after the initilization of the stack with ([1,n], empty, 0)
   * lets check, if we can have an open chain without breaking the threshold
   * this is an ugly work-arround cause in case of an open chain we do not have to
   * backtrack anything further...
   */
  if (evaluate_ext(1, length, 1, length, VRNA_DECOMP_EXT_UP, hc_dat_ext)) {
    tmp_en = vrna_E_exterior_loop(length, md) * (int)n_seq;

    if (sc) {
      if (sc->energy_up)
        tmp_en += sc->energy_up[1][length];

      if (sc->f)
        tmp_en += sc->f(1, j, 1, j, VRNA_DECOMP_EXT_UP, sc->data);
    }

    if (tmp_en <= threshold) {
      new_state                 = derive_new_state(1, 2, state, 0, VRNA_MX_FLAG_F5);
      new_state->partial_energy = 0;
      push(env->Stack, new_state);
      env->nopush = false;
    }
  }

  /*
   * ok, lets check if we can do an exterior hairpin without breaking the threshold
   * best energy should be 0 if we are here
   */
  if (FcH + best_energy <= threshold) {
    /*
     * lets search for all exterior hairpin cases, that fit into our threshold barrier
     * we use index k,l to avoid confusion with i,j index of our state...
     * if we reach here, i should be 1 and j should be n respectively
     */
    for (k = i; k < j; k++) {
      if (hc->up_hp[1] < k)
        break;

      for (l = j; l >= k + turn + 1; l--) {
        kl = indx[l] + k;

        if (c[kl] != INF) {
          tmpE = vrna_eval_hairpin(fc, l, k, VRNA_EVAL_LOOP_DEFAULT);

          if (c[kl] + tmpE + best_energy <= threshold)
            /*
             * what we really have to do is something like this, isn't it?
             * we have to create a new state, with interval [k,l], then we
             * add our loop energy as initial energy of this state and put
             * the state onto the stack R... for further refinement...
             * we also denote this new interval to be scanned in C
             */
            fork_state(k, l, state, tmpE, VRNA_MX_FLAG_C, env);
        }
      }
    }
  }

  /* now lets see, if we can do an exterior internal loop without breaking the threshold */
  if (FcI + best_energy <= threshold) {
    /* now we search for our exterior internal loop possibilities */
    for (k = i; k < j; k++) {
      for (l = j; l >= k + turn + 1; l--) {
        kl = indx[l] + k;         /* just confusing these indices ;-) */

        if ((hard_constraints[length * k + l] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) &&
            (c[kl] != INF)) {
          type = rtype[vrna_get_ptype(kl, ptype)];

          for (p = l + 1; p < j; p++) {
            u1 = p - l - 1;
            if (u1 + k - 1 > MAXLOOP)
              break;

            if (hc->up_int[l + 1] < u1)
              break;

            qmin = p + turn + 1;

            if (k + j > l + MAXLOOP + turn + 2)
              qmin = u1 + k + j - MAXLOOP - 1;

            for (q = j; q >= qmin; q--) {
              u2 = j + k - q - 1;

              if (hc->up_int[q + 1] < u2)
                break;

              if ((evaluate_int(k, l, p, q, hc_dat_int)) &&
                  (c[indx[q] + p] != INF)) {
                type_2 = rtype[vrna_get_ptype(indx[q] + p, ptype)];

                if (u1 + u2 > MAXLOOP)
                  continue;

                tmpE = vrna_E_internal(u1,
                                 u2,
                                 type,
                                 type_2,
                                 S1[l + 1],
                                 S1[k - 1],
                                 S1[p - 1],
                                 S1[q + 1],
                                 P);

                if (sc_int_pair_ext)
                  tmpE += sc_int_pair_ext(k, l, p, q, sc_dat_int);

                if (c[kl] + c[indx[q] + p] + tmpE + best_energy <= threshold)
                  /*
                   * ok, similar to the hairpin stuff, we add new states onto the stack R
                   * but in contrast to the hairpin decomposition, we have to add two new
                   * intervals, enclosed by k,l and p,q respectively and we also have to
                   * add the partial energy, that comes from the exterior internal loop
                   */
                  fork_two_states(k, l, p, q, state, tmpE, VRNA_MX_FLAG_C, VRNA_MX_FLAG_C, env);
              }
            }
          }
        }
      }
    }
  }

  /* and last but not least, we have a look, if we can do an exterior multiloop within the energy threshold */
  if (FcM <= threshold) {
    /*
     * this decomposition will be somehow more complicated...so lets see what we do here...
     * first we want to find out which split inidices we can use without exceeding the threshold
     */
    for (k = turn + 1; k + 2 * turn < j; k++) {
      if ((evaluate_mb(1, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc_dat_mb)) &&
          (fM1_new[k] != INF) &&
          (fM2_real[indx[length] + k + 1] != INF)) {
        tmpE2 = fM1_new[k] +
                fM2_real[indx[length] + k + 1] +
                P->MLclosing;

        if (sc_mb_decomp_ml) {
          tmpE2 += sc_mb_decomp_ml(1,
                                   j,
                                   k,
                                   k + 1,
                                   sc_dat_mb);
        }

        if (tmpE2 + best_energy <= threshold) {
          /*
           *  1. Find stem(s) in fM1_new
           */
          for (l = 1; l + turn < k; l++) {
            eval = (hc->up_ml[1] >= (l - 1)) ?
                   (hc->mx[length * k + l] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) :
                   0;
            if ((hc->f) &&
                (!hc->f(1, k, l, k, VRNA_DECOMP_ML_ML, hc->data)))
              eval = 0;

            if (eval) {
              e = c[indx[k] + l];
              e_part = n_seq *
                       (l - 1) *
                       P->MLbase;

              switch (fc->type) {
/*
                case VRNA_FC_TYPE_COMPARATIVE:
                  for (s = 0; s < n_seq; s++) {
                    type  = vrna_get_ptype_md(SS[s][l], SS[s][k], md);
                    e     += vrna_E_multibranch_stem(type, S5[s][l], S3[s][k], P);
                  }
                  break;
*/
                default:
                  type    = vrna_get_ptype_md(S[l], S[k], md);
                  e_part  += vrna_E_multibranch_stem(type, S1[l - 1], S1[k + 1], P);
                  break;
              }

              if (sc_mb_red_ml)
                e_part += sc_mb_red_ml(1, k, l, k, sc_dat_mb);

              if (fM2_real[indx[length] + k + 1] + e + e_part + P->MLclosing <= threshold) {
                /*
                 * we've (hopefully) found a valid decomposition of fM2 and therefor we have all
                 * three intervals for our new state to be pushed on stack R
                 */
                new_state = copy_state(state);

                /* first interval leads is the single branch from fM1_new */
                new_interval = make_interval(l, k, VRNA_MX_FLAG_C);
                push(new_state->Intervals, new_interval);
                env->nopush = false;

                /* next, we put fM2 part as interval to scan further */
                new_interval = make_interval(k + 1, length, VRNA_MX_FLAG_M2);
                push(new_state->Intervals, new_interval);
                env->nopush = false;

                /* mmh, we add the energy for closing the multiloop now... */
                new_state->partial_energy += e_part +
                                             P->MLclosing;
                /* next we push our state onto the R stack */
                push(env->Stack, new_state);
                env->nopush = false;
              }
            }
          }
        }
      }
    }
  }
}


PRIVATE INLINE void
scan_fms5(vrna_fold_compound_t  *fc,
          unsigned int          i,
          unsigned int          strand,
          int                   threshold,
          STATE                 *state,
          subopt_env            *env,
          constraint_helpers    *constraints_dat)
{
  char                      *ptype;
  short                     *S1, s5, s3;
  unsigned int              k, type, *sn, *se, end;
  int                       dangle_model, element_energy, best_energy, e_gq, *c, **fms5,
                            *indx, with_gquad;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  struct hc_ext_def_dat     *hc_dat;
  vrna_hc_eval_f evaluate;
  struct sc_f5_dat          *sc_dat;
  sc_ext_red_cb             sc_red_ext;
  sc_ext_red_cb             sc_red_stem;
  sc_ext_red_cb             sc_decomp;

  STATE                     *temp_state;
  vrna_smx_csr(int)         *c_gq;

  sn    = fc->strand_number;
  se    = fc->strand_end;
  indx  = fc->jindx;
  ptype = fc->ptype;
  S1    = fc->sequence_encoding;
  P     = fc->params;
  md    = &(P->model_details);

  dangle_model  = md->dangles;
  with_gquad    = md->gquad;

  c     = fc->matrices->c;
  c_gq  = fc->matrices->c_gq;
  fms5  = fc->matrices->fms5;

  hc_dat      = &(constraints_dat->hc_dat_ext);
  evaluate    = constraints_dat->hc_eval_ext;
  sc_dat      = &(constraints_dat->sc_dat_ext);
  sc_red_ext  = sc_dat->red_ext;
  sc_red_stem = sc_dat->red_stem;
  sc_decomp   = sc_dat->decomp;

  best_energy = best_attainable_energy(fc, state);  /* .. on remaining intervals */

  end = se[strand];

  /* no more pairs if too close to strand boundary ? */
  if (i + 1 > se[strand]) {
    state->partial_energy += fms5[strand][i];

    if (env->nopush) {
      push_back(env->Stack, state);
      env->nopush = false;
    }

    return;
  }

  /* find split in fms5 */

  if ((evaluate(i, end, i + 1, end, VRNA_DECOMP_EXT_EXT, hc_dat)) &&
      (fms5[strand][i] != INF)) {
    element_energy = 0;

    if (sc_red_ext)
      element_energy += sc_red_ext(i, end, i + 1, end, sc_dat);

    if (fms5[strand][i + 1] + element_energy + best_energy <= threshold)
      /* no basepair, nibbling of 5'-end */
      fork_state(i + 1, strand, state, element_energy, VRNA_MX_FLAG_MS5, env);
  }

  if (evaluate(i, end, i, end, VRNA_DECOMP_EXT_STEM, hc_dat)) {
    type = vrna_get_ptype(indx[end] + i, ptype);

    switch (dangle_model) {
      case 2:
        s5  = ((i > 1) && (sn[i - 1] == sn[i])) ? S1[i - 1] : -1;
        s3  = -1;
        break;

      default:
        s5  = -1;
        s3  = -1;
        break;
    }

    element_energy = vrna_E_exterior_stem(type, s5, s3, P);

    if (sc_red_stem)
      element_energy += sc_red_stem(i, end, i, end, sc_dat);

    if (c[indx[end] + i] + element_energy + best_energy <= threshold) {
      repeat(fc,
             i,
             end,
             state,
             element_energy,
             0,
             best_energy,
             threshold,
             env,
             constraints_dat);
    }
  }

  if (with_gquad) {
#ifndef VRNA_DISABLE_C11_FEATURES
    e_gq = vrna_smx_csr_get(c_gq, i, end, INF);
#else
    e_gq = vrna_smx_csr_int_get(c_gq, i, end, INF);
#endif
    if(e_gq != INF) {
      element_energy = 0;

      if (sc_red_stem)
        element_energy += sc_red_stem(i, end, i, end, sc_dat);

      if (e_gq + element_energy + best_energy <= threshold) {
        repeat_gquad(fc,
                     i,
                     end,
                     state,
                     element_energy,
                     0,
                     best_energy,
                     threshold,
                     env,
                     constraints_dat);
      }
    }
  }

  for (k = i + 1; k < end; k++) {
    if (with_gquad) {
#ifndef VRNA_DISABLE_C11_FEATURES
      e_gq = vrna_smx_csr_get(c_gq, i, k, INF);
#else
      e_gq = vrna_smx_csr_int_get(c_gq, i, k, INF);
#endif
      if ((fms5[strand][k + 1] != INF) &&
          (e_gq != INF)) {
        element_energy = 0;

        if (sc_decomp)
          element_energy += sc_decomp(i, end, k, k + 1, sc_dat);

        if (sc_red_stem)
          element_energy += sc_red_stem(i, k, i, k, sc_dat);

        if (fms5[strand][k + 1] + e_gq + element_energy + best_energy <= threshold) {
          temp_state  = derive_new_state(k + 1, strand, state, 0, VRNA_MX_FLAG_MS5);
          env->nopush = false;
          repeat_gquad(fc,
                       i,
                       k,
                       temp_state,
                       element_energy,
                       fms5[strand][k + 1],
                       best_energy,
                       threshold,
                       env,
                       constraints_dat);
          free_state_node(temp_state);
        }
      }
    }

    if (evaluate(i, end, k, k + 1, VRNA_DECOMP_EXT_STEM_EXT, hc_dat)) {
      type = vrna_get_ptype(indx[k] + i, ptype);

      switch (dangle_model) {
        case 2:
          s5  = ((i > 1) && (sn[i - 1] == sn[i])) ? S1[i - 1] : -1;
          s3  = (sn[k] == sn[k + 1]) ? S1[k + 1] : -1;
          break;

        default:
          s5  = -1;
          s3  = -1;
          break;
      }

      element_energy = vrna_E_exterior_stem(type, s5, s3, P);

      if (sc_decomp)
        element_energy += sc_decomp(i, end, k, k + 1, sc_dat);

      if (sc_red_stem)
        element_energy += sc_red_stem(i, k, i, k, sc_dat);

      if (fms5[strand][k + 1] + c[indx[k] + i] + element_energy + best_energy <= threshold) {
        temp_state  = derive_new_state(k + 1, strand, state, 0, VRNA_MX_FLAG_MS5);
        env->nopush = false;
        repeat(fc,
               i,
               k,
               temp_state,
               element_energy,
               fms5[strand][k + 1],
               best_energy,
               threshold,
               env,
               constraints_dat);
        free_state_node(temp_state);
      }
    }
  }
}


PRIVATE INLINE void
scan_fms3(vrna_fold_compound_t  *fc,
          unsigned int          i,
          unsigned int          strand,
          int                   threshold,
          STATE                 *state,
          subopt_env            *env,
          constraint_helpers    *constraints_dat)
{
  char                      *ptype;
  short                     *S1, s5, s3;
  unsigned int              dangle_model, *sn, *ss, start, k, type, length, with_gquad;
  int                       element_energy, best_energy, e_gq, *c, **fms3,
                            *indx;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  struct hc_ext_def_dat     *hc_dat;
  vrna_hc_eval_f evaluate;
  struct sc_f5_dat          *sc_dat;
  sc_ext_red_cb             sc_red_ext;
  sc_ext_red_cb             sc_red_stem;
  sc_ext_red_cb             sc_decomp;
  STATE                     *temp_state;
  vrna_smx_csr(int)         *c_gq;

  length  = fc->length;
  sn      = fc->strand_number;
  ss      = fc->strand_start;
  indx    = fc->jindx;
  ptype   = fc->ptype;
  S1      = fc->sequence_encoding;
  P       = fc->params;
  md      = &(P->model_details);

  dangle_model  = md->dangles;
  with_gquad    = md->gquad;

  c     = fc->matrices->c;
  c_gq  = fc->matrices->c_gq;
  fms3  = fc->matrices->fms3;

  start = ss[strand];

  hc_dat      = &(constraints_dat->hc_dat_ext);
  evaluate    = constraints_dat->hc_eval_ext;
  sc_dat      = &(constraints_dat->sc_dat_ext);
  sc_red_ext  = sc_dat->red_ext;
  sc_red_stem = sc_dat->red_stem;
  sc_decomp   = sc_dat->decomp;

  best_energy = best_attainable_energy(fc, state);  /* .. on remaining intervals */

  /* no more pairs if too close to strand boundary ? */
  if (i < ss[strand] + 1) {
    state->partial_energy += fms3[strand][i];

    if (env->nopush) {
      push_back(env->Stack, state);
      env->nopush = false;
    }

    return;
  }

  if ((evaluate(start, i, start, i - 1, VRNA_DECOMP_EXT_EXT, hc_dat)) &&
      (fms3[strand][i - 1] != INF)) {
    element_energy = 0;

    if (sc_red_ext)
      element_energy += sc_red_ext(start, i, start, i - 1, sc_dat);

    if (fms3[strand][i - 1] + element_energy + best_energy <= threshold)
      /* no basepair, nibbling of 5'-end */
      fork_state(i - 1, strand, state, element_energy, VRNA_MX_FLAG_MS3, env);
  }

  if (evaluate(start, i, start, i, VRNA_DECOMP_EXT_STEM, hc_dat)) {
    type = vrna_get_ptype(indx[i] + start, ptype);

    switch (dangle_model) {
      case 2:
        s5  = -1;
        s3  = ((i < length) && (sn[i] == sn[i + 1])) ? S1[i + 1] : -1;
        break;

      default:
        s5 = s3 = -1;
        break;
    }

    element_energy = vrna_E_exterior_stem(type, s5, s3, P);

    if (sc_red_stem)
      element_energy += sc_red_stem(start, i, start, i, sc_dat);

    if (c[indx[i] + start] + element_energy + best_energy <= threshold) {
      repeat(fc,
             start,
             i,
             state,
             element_energy,
             0,
             best_energy,
             threshold,
             env,
             constraints_dat);
    }
  }

  if (with_gquad) {
#ifndef VRNA_DISABLE_C11_FEATURES
    e_gq = vrna_smx_csr_get(c_gq, start, i, INF);
#else
    e_gq = vrna_smx_csr_int_get(c_gq, start, i, INF);
#endif
    if (e_gq != INF) {
      element_energy = 0;

      if (sc_red_stem)
        element_energy += sc_red_stem(start, i, start, i, sc_dat);

      if (e_gq + element_energy + best_energy <= threshold) {
        repeat_gquad(fc,
                     start,
                     i,
                     state,
                     element_energy,
                     0,
                     best_energy,
                     threshold,
                     env,
                     constraints_dat);
      }
    }
  }

  for (k = start; k < i; k++) {
    if (with_gquad) {
#ifndef VRNA_DISABLE_C11_FEATURES
      e_gq = vrna_smx_csr_get(c_gq, k + 1, i, INF);
#else
      e_gq = vrna_smx_csr_int_get(c_gq, k + 1, i, INF);
#endif
      if ((fms3[strand][k] != INF) &&
          (e_gq != INF)) {
        element_energy = 0;

        if (sc_decomp)
          element_energy += sc_decomp(start, i, k, k + 1, sc_dat);

        if (sc_red_stem)
          element_energy += sc_red_stem(k + 1, i, k + 1, i, sc_dat);

        if (fms3[strand][k] + e_gq + element_energy + best_energy <= threshold) {
          temp_state  = derive_new_state(k, strand, state, 0, VRNA_MX_FLAG_MS3);
          env->nopush = false;
          repeat_gquad(fc,
                       k + 1,
                       i,
                       temp_state,
                       element_energy,
                       fms3[strand][k],
                       best_energy,
                       threshold,
                       env,
                       constraints_dat);
          free_state_node(temp_state);
        }
      }
    }

    if (evaluate(start, i, k, k + 1, VRNA_DECOMP_EXT_EXT_STEM, hc_dat)) {
      type = vrna_get_ptype(indx[i] + k + 1, ptype);

      switch (dangle_model) {
        case 2:
          s5  = (sn[k] == sn[k + 1]) ? S1[k] : -1;
          s3  = ((i < length) && (sn[i] == sn[i + 1])) ? S1[i + 1] : -1;
          break;

        default:
          s5 = s3 = -1;
          break;
      }

      element_energy = vrna_E_exterior_stem(type, s5, s3, P);

      if (sc_decomp)
        element_energy += sc_decomp(start, i, k, k + 1, sc_dat);

      if (sc_red_stem)
        element_energy += sc_red_stem(k + 1, i, k + 1, i, sc_dat);

      if (fms3[strand][k] + c[indx[i] + k + 1] + element_energy + best_energy <= threshold) {
        temp_state  = derive_new_state(k, strand, state, 0, VRNA_MX_FLAG_MS3);
        env->nopush = false;
        repeat(fc,
               k + 1,
               i,
               temp_state,
               element_energy,
               fms3[strand][k],
               best_energy,
               threshold,
               env,
               constraints_dat);
        free_state_node(temp_state);
      }
    }
  }
}


PRIVATE INLINE void
scan_gquad(vrna_fold_compound_t *fc,
           int                  i,
           int                  j,
           int                  threshold,
           STATE                *state,
           subopt_env           *env,
           constraint_helpers   *constraints_dat)
{
  int best_energy;

  best_energy = best_attainable_energy(fc, state);  /* .. on remaining intervals */

  /* we have a gquad */
  repeat_gquad(fc,
               i,
               j,
               state,
               0,
               0,
               best_energy,
               threshold,
               env,
               constraints_dat);

  if (env->nopush)
    vrna_log_warning("%d,%d\nOops, no solution in gquad-repeat!", i, j);
}


PRIVATE void
repeat_gquad(vrna_fold_compound_t *fc,
             int                  i,
             int                  j,
             STATE                *state,
             int                  part_energy,
             int                  temp_energy,
             int                  best_energy,
             int                  threshold,
             subopt_env           *env,
             constraint_helpers   *constraints_dat VRNA_UNUSED)
{
  short             *S1;
  unsigned int      *sn, cnt, *L, *l, num_gquads;
  int               element_energy;
  vrna_param_t      *P;
  vrna_smx_csr(int) *c_gq;

  sn    = fc->strand_number;
  c_gq  = fc->matrices->c_gq;
  S1    = fc->sequence_encoding;
  P     = fc->params;


  /* find all gquads that fit into the energy range and the interval [i,j] */
  STATE *new_state;

  best_energy += part_energy; /* energy of current structural element */
  best_energy += temp_energy; /* energy from unpushed interval */

  if (sn[i] == sn[j]) {
#ifndef VRNA_DISABLE_C11_FEATURES
    element_energy = vrna_smx_csr_get(c_gq, i, j, INF);
#else
    element_energy = vrna_smx_csr_int_get(c_gq, i, j, INF);
#endif
    if ((element_energy != INF) &&
        (element_energy + best_energy <= threshold)) {
      /* find out how many gquads we might expect in the interval [i,j] */
      num_gquads = get_gquad_count(S1, i, j);
      num_gquads++;
      L     = (unsigned int *)vrna_alloc(sizeof(unsigned int) * num_gquads);
      l     = (unsigned int *)vrna_alloc(sizeof(unsigned int) * num_gquads * 3);
      L[0]  = 0;

      get_gquad_pattern_exhaustive(S1, i, j, P, L, l, threshold - best_energy);

      for (cnt = 0; L[cnt] != 0; cnt++) {
        new_state = copy_state(state);
        make_gquad(i, fc->length, L[cnt], &(l[3 * cnt]), new_state);
        new_state->partial_energy += part_energy;
        /* re-compute energies */
        new_state->partial_energy += vrna_E_gquad(L[cnt], &(l[3*cnt]), P);
        /* new_state->best_energy =
         * hairpin[unpaired] + element_energy + best_energy; */
        push(env->Stack, new_state);
        env->nopush = false;
      }
      free(L);
      free(l);
    }
  }

  best_energy -= part_energy;
  best_energy -= temp_energy;
  return;
}


PRIVATE void
repeat(vrna_fold_compound_t *fc,
       unsigned int         i,
       unsigned int         j,
       STATE                *state,
       int                  part_energy,
       int                  temp_energy,
       int                  best_energy,
       int                  threshold,
       subopt_env           *env,
       constraint_helpers   *constraints_dat)
{
  /*
   * routine to find stacks, bulges, internal loops and  multiloops
   * within interval closed by basepair i,j
   */

  char                      *ptype;
  short                     *S1;
  unsigned int              n, *sn, *se, nick, k, p, q, no_close, type, type_2, rt, noGUclosure,
                            noLP, with_gquad, dangle_model, minq, cnt, u1, u2;
  int                       ij, energy, new, mm, element_energy,  *c, *fML, *fM1, **fms5, **fms3,
                            *indx, *rtype, eee, aux_eee, tmp_en;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  struct hc_int_def_dat     *hc_dat_int;
  struct hc_ext_def_dat     *hc_dat_ext;
  struct hc_mb_def_dat      *hc_dat_mb;
  eval_hc                   evaluate_int;
  vrna_hc_eval_f evaluate_ext;
  vrna_hc_eval_f evaluate_mb;

  struct sc_int_dat         *sc_dat_int;
  struct sc_mb_dat          *sc_dat_mb;
  sc_int_cb                 sc_int_pair;
  sc_mb_pair_cb             sc_mb_pair;
  sc_mb_red_cb              sc_mb_decomp_ml;
  STATE                     *new_state;

  n     = fc->length;
  S1    = fc->sequence_encoding;
  ptype = fc->ptype;
  indx  = fc->jindx;
  sn    = fc->strand_number;
  se    = fc->strand_end;
  P     = fc->params;
  md    = &(P->model_details);
  rtype = &(md->rtype[0]);

  noGUclosure   = md->noGUclosure;
  noLP          = md->noLP;
  with_gquad    = md->gquad;
  dangle_model  = md->dangles;

  c     = fc->matrices->c;
  fML   = fc->matrices->fML;
  fM1   = fc->matrices->fM1;
  fms5  = fc->matrices->fms5;
  fms3  = fc->matrices->fms3;

  hc = fc->hc;

  hc_dat_ext    = &(constraints_dat->hc_dat_ext);
  hc_dat_int    = &(constraints_dat->hc_dat_int);
  hc_dat_mb     = &(constraints_dat->hc_dat_mb);
  evaluate_ext  = constraints_dat->hc_eval_ext;
  evaluate_int  = constraints_dat->hc_eval_int;
  evaluate_mb   = constraints_dat->hc_eval_mb;

  sc              = fc->sc;
  sc_dat_int      = &(constraints_dat->sc_dat_int);
  sc_dat_mb       = &(constraints_dat->sc_dat_mb);
  sc_int_pair     = constraints_dat->sc_dat_int.pair;
  sc_mb_pair      = constraints_dat->sc_dat_mb.pair;
  sc_mb_decomp_ml = constraints_dat->sc_dat_mb.decomp_ml;

  ij = indx[j] + i;

  type = vrna_get_ptype(ij, ptype);
  /*
   * if (type==0) fprintf(stderr, "repeat: Warning: %d %d can't pair\n", i,j);
   */

  no_close = ((noGUclosure) &&
              ((type == 3) || (type == 4)));

  if ((noLP) &&
      (i + 2 < j)) {
    /* always consider the structure with additional stack */
    if (evaluate_int(i, j, i + 1, j - 1, hc_dat_int)) {
      type_2  = rtype[vrna_get_ptype(indx[j - 1] + i + 1, ptype)];
      energy  = 0;

      energy = vrna_E_internal(0, 0, type, type_2, S1[i + 1], S1[j - 1], S1[i + 1], S1[j - 1], P);

      if (sc_int_pair) {
        energy += sc_int_pair(i,
                              j,
                              i + 1,
                              j - 1,
                              sc_dat_int);
      }

      new_state = derive_new_state(i + 1, j - 1, state, part_energy + energy, VRNA_MX_FLAG_C);
      make_pair(i, j, new_state);
      make_pair(i + 1, j - 1, new_state);

      /* new_state->best_energy = new + best_energy; */
      push(env->Stack, new_state);
      env->nopush = false;
      if (i == 1 || state->structure[i - 2] != '(' || state->structure[j] != ')')
        /* adding a stack is the only possible structure */
        return;
    }
  }

  best_energy += part_energy; /* energy of current structural element */
  best_energy += temp_energy; /* energy from unpushed interval */

  if (hc->mx[n * i + j] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    for (p = i + 1; p <= MIN2(j - 2, i + MAXLOOP + 1); p++) {
      u1 = p - i - 1;

      minq = p + 1;

      if (i + MAXLOOP + 2 < j)
        minq = j + u1 - MAXLOOP - 1;

      if (hc->up_int[i + 1] < u1)
        break;

      for (q = j - 1; q >= minq; q--) {
        u2 = j - q - 1;

        if (hc->up_int[q + 1] < u2)
          break;

        /* skip stack if noLP, since we've already processed it above */
        if ((noLP) && (p == i + 1) && (q == j - 1))
          continue;

        if (!(hc->mx[n * p + q] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC))
          continue;

        if (c[indx[q] + p] == INF)
          continue;

        type_2 = vrna_get_ptype(indx[q] + p, ptype);

        if (noGUclosure)
          if (no_close || (type_2 == 3) || (type_2 == 4))
            if ((p > i + 1) || (q < j - 1))
              continue;

        if (evaluate_int(i, j, p, q, hc_dat_int)) {
          energy = vrna_E_internal(u1,
                                   u2,
                                   type,
                                   rtype[type_2],
                                   S1[i + 1],
                                   S1[j - 1],
                                   S1[p - 1],
                                   S1[q + 1],
                                   P);

          new = energy + c[indx[q] + p];

          if (sc_int_pair)
            energy += sc_int_pair(i, j, p, q, sc_dat_int);

          new = energy + c[indx[q] + p];

          if (new + best_energy <= threshold)
            /* stack, bulge, or internal loop */
            fork_int_state(i, j, p, q, state, part_energy + energy, env);
        } /*end of if block */
      }   /* end of q-loop */
    }     /* end of p-loop */
  }

  /* base pair (i,j) encloses a loop with strand nick? */
  if ((sn[i] != sn[j]) &&
      (evaluate_ext(i, j, i, j, VRNA_DECOMP_EXT_STEM, hc_dat_ext))) {
    rt = rtype[type];

    element_energy = P->DuplexInit;

    switch (dangle_model) {
      case 0:
        element_energy += vrna_E_exterior_stem(rt, -1, -1, P);
        break;
      default:
        element_energy +=
          vrna_E_exterior_stem(rt,
                          (sn[j - 1] == sn[j]) ?
                          S1[j - 1] :
                          -1,
                          (sn[i] == sn[i + 1]) ?
                          S1[i + 1] :
                          -1,
                          P);
        break;
    }

    if ((sc) && (sc->f))
      element_energy += sc->f(i, j, i, j, VRNA_DECOMP_EXT_STEM, sc->data);

    if (sn[i] != sn[i + 1]) {
      if ((sn[j - 1] != sn[j]) &&
          (i + 1 == j)) {
        if (element_energy + best_energy <= threshold) {
          fork_state_pair(i,
                          j,
                          state,
                          part_energy + element_energy,
                          env);
        }
      } else if (sn[j - 1] == sn[j]) {
        if (fms3[sn[i + 1]][j - 1] + element_energy + best_energy <= threshold) {
          /* continue backtracking in fms3[sn[i + 1]][j - 1] */
          fork_state_pair_interval(i,
                                   j,
                                   j - 1,
                                   sn[i + 1],
                                   state,
                                   part_energy + element_energy,
                                   VRNA_MX_FLAG_MS3,
                                   env);
        }
      }
    } else if (sn[j - 1] != sn[j]) {
      if (fms5[sn[j - 1]][i + 1] + element_energy + best_energy <= threshold) {
        /* continue backtracking in fms5[sn[j - 1]][i + 1] */
        fork_state_pair_interval(i,
                                 j,
                                 i + 1,
                                 sn[j - 1],
                                 state,
                                 part_energy + element_energy,
                                 VRNA_MX_FLAG_MS5,
                                 env);
      }
    } else {
      energy = 0;

      if (se[sn[i]] > i)
        energy += fms5[sn[i]][i + 1];

      if (j - 1 > se[sn[i]])
        energy += fms3[sn[se[sn[i]] + 1]][j - 1];

      if (energy + element_energy + best_energy <= threshold) {
        fork_two_states_pair_ms(i,
                                j,
                                sn[i],
                                sn[se[sn[i]] + 1],
                                state,
                                part_energy + element_energy,
                                env);
      }

      nick = se[sn[i]] + 1;

      while (sn[nick] != sn[j]) {
        energy = 0;
        if (i + 1 <= se[sn[nick]])
          energy += fms5[sn[nick]][i + 1];

        if (se[sn[nick]] + 1 <= j - 1)
          energy += fms3[sn[se[sn[nick]] + 1]][j - 1];

        if (energy + element_energy + best_energy <= threshold) {
          fork_two_states_pair_ms(i,
                                  j,
                                  sn[nick],
                                  sn[se[sn[nick]] + 1],
                                  state,
                                  part_energy + element_energy,
                                  env);
        }

        nick = se[sn[nick]] + 1;
      }
    }
  }

  mm  = P->MLclosing;
  rt  = rtype[type];

  if (evaluate_mb(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, hc_dat_mb)) {
    element_energy = mm;
    switch (dangle_model) {
      case 0:
        element_energy = vrna_E_multibranch_stem(rt, -1, -1, P) + mm;
        break;
      default:
        element_energy = vrna_E_multibranch_stem(rt, S1[j - 1], S1[i + 1], P) + mm;
        break;
    }

    if (sc_mb_pair)
      element_energy += sc_mb_pair(i, j, sc_dat_mb);

    /* multiloop decomposition */
    for (k = i + 2; k <= j - 2; k++) {
      if (evaluate_mb(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc_dat_mb)) {
        eee = fML[indx[k - 1] + i + 1];

        if ((eee != INF) &&
            (fM1[indx[j - 1] + k] != INF)) {
          eee += fM1[indx[j - 1] + k] +
                 best_energy;

          aux_eee = element_energy;

          if (sc_mb_decomp_ml)
            aux_eee += sc_mb_decomp_ml(i + 1, j - 1, k - 1, k, sc_dat_mb);

          if ((eee + aux_eee) <= threshold)
            fork_two_states_pair(i, j, k, state, part_energy + aux_eee, VRNA_MX_FLAG_M, VRNA_MX_FLAG_M1, env);
        }
      }
    }
  }

  if (sn[i] == sn[j]) {
    if (!no_close) {
      element_energy = vrna_eval_hairpin(fc, i, j, VRNA_EVAL_LOOP_DEFAULT);

      if (element_energy != INF) {
        if (element_energy + best_energy <= threshold)
          /* hairpin structure */
          fork_state_pair(i, j, state, part_energy + element_energy, env);
      }
    }

    if (with_gquad) {
      /* now we have to find all loops where (i,j) encloses a gquad in an internal loops style */
      vrna_array(int) ge = NULL;
      vrna_array(unsigned int) ps = NULL;
      vrna_array(unsigned int) qs = NULL;
      ge = vrna_gq_int_loop_subopt(fc, i, j, &ps, &qs, threshold - best_energy);
      if (ge) {
        for (cnt = 0; cnt < vrna_array_size(ge); cnt++) {
          if ((hc->up_int[i + 1] >= ps[cnt] - i - 1) &&
              (hc->up_int[qs[cnt] + 1] >= j - qs[cnt] - 1)) {
            tmp_en = ge[cnt];

            if (sc_int_pair)
              tmp_en += sc_int_pair(i, j, ps[cnt], qs[cnt], sc_dat_int);

            new_state = derive_new_state(ps[cnt], qs[cnt], state, tmp_en + part_energy, VRNA_MX_FLAG_G);

            make_pair(i, j, new_state);

            /* new_state->best_energy = new + best_energy; */
            push(env->Stack, new_state);
            env->nopush = false;
          }
        }
      }
      vrna_array_free(ge);
      vrna_array_free(ps);
      vrna_array_free(qs);
    }
  }

  best_energy -= part_energy;
  best_energy -= temp_energy;
  return;
}


PRIVATE void
old_subopt_print(const char *structure,
                 float      energy,
                 void       *data)
{
  struct old_subopt_dat *d = (struct old_subopt_dat *)data;

  if (structure && d->fp) {
    char *e_string = vrna_strdup_printf(" %6.2f", energy);
    print_structure(d->fp, structure, e_string);
    free(e_string);
  }
}


PRIVATE void
old_subopt_store(const char *structure,
                 float      energy,
                 void       *data)
{
  struct old_subopt_dat *d = (struct old_subopt_dat *)data;

  /* store solution */
  if (d->n_sol + 1 == d->max_sol) {
    d->max_sol      *= 2;
    d->SolutionList =
      (vrna_subopt_solution_t *)vrna_realloc(d->SolutionList,
                                             d->max_sol * sizeof(vrna_subopt_solution_t));
  }

  if (structure) {
    d->SolutionList[d->n_sol].energy      = energy;
    d->SolutionList[d->n_sol++].structure = strdup(structure);
  } else {
    d->SolutionList[d->n_sol].energy      = 0;
    d->SolutionList[d->n_sol++].structure = NULL;
  }
}


PRIVATE void
old_subopt_store_compressed(const char  *structure,
                            float       energy,
                            void        *data)
{
  struct old_subopt_dat *d = (struct old_subopt_dat *)data;

  /* store solution */
  if (d->n_sol + 1 == d->max_sol) {
    d->max_sol      *= 2;
    d->SolutionList =
      (vrna_subopt_solution_t *)vrna_realloc(d->SolutionList,
                                             d->max_sol * sizeof(vrna_subopt_solution_t));
  }

  if (structure) {
    d->SolutionList[d->n_sol].energy = energy;
    if (d->strands > 1) {
      char **tok  = vrna_strsplit(structure, NULL);
      char *s     = vrna_strjoin((const char **)tok, NULL);

      for (char **ptr = tok; *ptr!= NULL; ptr++)
        free(*ptr);
      free(tok);

      d->SolutionList[d->n_sol++].structure = vrna_db_pack(s);
      free(s);
    } else {
      d->SolutionList[d->n_sol++].structure = vrna_db_pack(structure);
    }
  } else {
    d->SolutionList[d->n_sol].energy      = 0;
    d->SolutionList[d->n_sol++].structure = NULL;
  }
}


/*
 * ###########################################
 * # deprecated functions below              #
 *###########################################
 */
#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

PUBLIC SOLUTION *
subopt(char *seq,
       char *structure,
       int  delta,
       FILE *fp)
{
  return wrap_subopt(seq, structure, NULL, delta, fold_constrained, 0, fp);
}


PUBLIC SOLUTION *
subopt_circ(char  *seq,
            char  *structure,
            int   delta,
            FILE  *fp)
{
  return wrap_subopt(seq, structure, NULL, delta, fold_constrained, 1, fp);
}


PUBLIC SOLUTION *
subopt_par(char         *seq,
           char         *structure,
           vrna_param_t *parameters,
           int          delta,
           int          is_constrained,
           int          is_circular,
           FILE         *fp)
{
  return wrap_subopt(seq, structure, parameters, delta, is_constrained, is_circular, fp);
}


PRIVATE SOLUTION *
wrap_subopt(char          *string,
            char          *structure,
            vrna_param_t  *parameters,
            int           delta,
            int           is_constrained,
            int           is_circular,
            FILE          *fp)
{
  vrna_fold_compound_t  *fc;
  vrna_param_t          *P;
  char                  *seq;

#ifdef _OPENMP
  /* Explicitly turn off dynamic threads */
  omp_set_dynamic(0);
#endif

  /* we need the parameter structure for hard constraints */
  if (parameters) {
    P = vrna_params_copy(parameters);
  } else {
    vrna_md_t md;
    set_model_details(&md);
    md.temperature  = temperature;
    P               = vrna_params(&md);
  }

  P->model_details.circ     = is_circular;
  P->model_details.uniq_ML  = uniq_ML = 1;

  /*
   * what about cofold sequences here? Is it safe to call the below cut_point_insert() ?
   * dirty hack to reinsert the '&' according to the global variable 'cut_point'
   */
  seq = vrna_cut_point_insert(string, cut_point);

  fc =
    vrna_fold_compound(seq,
                       &(P->model_details),
                       ((is_circular == 0) ? VRNA_OPTION_HYBRID : VRNA_OPTION_DEFAULT));

  if (parameters) {
    /* replace params if necessary */
    free(fc->params);
    fc->params = P;
  } else {
    free(P);
  }

  /* handle hard constraints in pseudo dot-bracket format if passed via simple interface */
  if (is_constrained && structure) {
    unsigned int constraint_options = 0;
    constraint_options |= VRNA_CONSTRAINT_DB
                          | VRNA_CONSTRAINT_DB_PIPE
                          | VRNA_CONSTRAINT_DB_DOT
                          | VRNA_CONSTRAINT_DB_X
                          | VRNA_CONSTRAINT_DB_ANG_BRACK
                          | VRNA_CONSTRAINT_DB_RND_BRACK
                          | VRNA_CONSTRAINT_DB_INTRAMOL
                          | VRNA_CONSTRAINT_DB_INTERMOL;

    vrna_constraints_add(fc, (const char *)structure, constraint_options);
  }

  if (backward_compat_compound && backward_compat)
    vrna_fold_compound_free(backward_compat_compound);

  backward_compat_compound  = fc;
  backward_compat           = 1;

  /* cleanup */
  free(seq);

  return vrna_subopt(fc, delta, subopt_sorted, fp);
}


#endif

/*
 * ---------------------------------------------------------------------------
 * Well, that is the end!----------------------------------------------------
 *---------------------------------------------------------------------------
 */
