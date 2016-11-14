/*
   suboptimal folding - Stefan Wuchty, Walter Fontana & Ivo Hofacker

                       Vienna RNA package
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "ViennaRNA/fold.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/fold_vars.h"
#include "list.h"
#include "ViennaRNA/eval.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/loop_energies.h"
#include "ViennaRNA/cofold.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/subopt.h"

/* hack */
#include "ViennaRNA/color_output.inc"

#ifdef _OPENMP
#include <omp.h>
#endif

#define true              1
#define false             0
#define ON_SAME_STRAND(I,J,C)  (((I)>=(C))||((J)<(C)))

/**
 *  @brief  Sequence interval stack element used in subopt.c
 */
typedef struct INTERVAL {
    int i;
    int j;
    int array_flag;
} INTERVAL;

typedef struct {
    char *structure;
    LIST *Intervals;
    int partial_energy;
    int is_duplex;
    /* int best_energy;   */ /* best attainable energy */
} STATE;

typedef struct {
  LIST          *Intervals;
  LIST          *Stack;
  int           nopush;
} subopt_env;


struct old_subopt_dat {

  unsigned long max_sol;
  unsigned long n_sol;
  SOLUTION      *SolutionList;
  FILE          *fp;
};

/*
#################################
# GLOBAL VARIABLES              #
#################################
*/
PUBLIC  int     subopt_sorted=0;                           /* output sorted by energy */
PUBLIC  int     density_of_states[MAXDOS+1];
PUBLIC  double  print_energy = 9999; /* printing threshold for use with logML */

/*
#################################
# PRIVATE VARIABLES             #
#################################
*/

/* some backward compatibility stuff */
PRIVATE int                 backward_compat           = 0;
PRIVATE vrna_fold_compound_t  *backward_compat_compound = NULL;

#ifdef _OPENMP

#pragma omp threadprivate(backward_compat_compound, backward_compat)

#endif

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/

#ifdef VRNA_BACKWARD_COMPAT

PRIVATE SOLUTION *
wrap_subopt(char *seq,
            char *structure,
            vrna_param_t *parameters,
            int delta,
            int is_constrained,
            int is_circular,
            FILE *fp);
#endif

PRIVATE void      make_pair(int i, int j, STATE *state);

/* mark a gquadruplex in the resulting dot-bracket structure */
PRIVATE void      make_gquad(int i, int L, int l[3], STATE *state);

PRIVATE INTERVAL  *make_interval (int i, int j, int ml);
/*@out@*/ PRIVATE STATE *make_state(/*@only@*/LIST *Intervals,
                                    /*@only@*/ /*@null@*/ char *structure,
                                    int partial_energy, int is_duplex, int length);
PRIVATE STATE     *copy_state(STATE * state);
PRIVATE void      print_state(STATE * state);
PRIVATE void      UNUSED print_stack(LIST * list);
/*@only@*/ PRIVATE LIST *make_list(void);
PRIVATE void      push(LIST * list, /*@only@*/ void *data);
PRIVATE void      *pop(LIST * list);
PRIVATE int       best_attainable_energy(vrna_fold_compound_t *vc, STATE * state);
PRIVATE void      scan_interval(vrna_fold_compound_t *vc, int i, int j, int array_flag, int threshold, STATE * state, subopt_env *env);
PRIVATE void      free_interval_node(/*@only@*/ INTERVAL * node);
PRIVATE void      free_state_node(/*@only@*/ STATE * node);
PRIVATE void      push_back(LIST *Stack, STATE * state);
PRIVATE char*     get_structure(STATE * state);
PRIVATE int       compare(const void *solution1, const void *solution2);
PRIVATE void      make_output(SOLUTION *SL, int cp, FILE *fp);
PRIVATE void      repeat(vrna_fold_compound_t *vc, int i, int j, STATE * state, int part_energy, int temp_energy, int best_energy, int threshold, subopt_env *env);
PRIVATE void      repeat_gquad(vrna_fold_compound_t *vc, int i, int j, STATE *state, int part_energy, int temp_energy, int best_energy, int threshold, subopt_env *env);


PRIVATE void      old_subopt_print( const char *structure, float energy, void *data);
PRIVATE void      old_subopt_store( const char *structure, float energy, void *data);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/



/*---------------------------------------------------------------------------*/
/*List routines--------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

PRIVATE void
make_pair(int i, int j, STATE *state)
{
  state->structure[i-1] = '(';
  state->structure[j-1] = ')';
}

PRIVATE void
make_gquad(int i, int L, int l[3], STATE *state)
{
  int x;
  for(x = 0; x < L; x++){
    state->structure[i - 1 + x] = '+';
    state->structure[i - 1 + x + L + l[0]] = '+';
    state->structure[i - 1 + x + 2*L + l[0] + l[1]] = '+';
    state->structure[i - 1 + x + 3*L + l[0] + l[1] + l[2]] = '+';
  }
}

/*---------------------------------------------------------------------------*/

PRIVATE INTERVAL *
make_interval(int i, int j, int array_flag)
{
  INTERVAL *interval;

  interval = lst_newnode(sizeof(INTERVAL));
  interval->i = i;
  interval->j = j;
  interval->array_flag = array_flag;
  return interval;
}

/*---------------------------------------------------------------------------*/

PRIVATE void
free_interval_node(INTERVAL * node)
{
  lst_freenode(node);
}

/*---------------------------------------------------------------------------*/

PRIVATE void
free_state_node(STATE * node)
{
  free(node->structure);
  if (node->Intervals)
    lst_kill(node->Intervals, lst_freenode);
  lst_freenode(node);
}

/*---------------------------------------------------------------------------*/

PRIVATE STATE *
make_state(LIST * Intervals,
           char *structure,
           int partial_energy,
           int is_duplex,
           int length)
{
  STATE *state;

  state = lst_newnode(sizeof(STATE));

  if (Intervals)
    state->Intervals = Intervals;
  else
    state->Intervals = lst_init();
  if (structure)
    state->structure = structure;
  else {
    int i;
    state->structure = (char *) vrna_alloc(length+1);
    for (i=0; i<length; i++)
      state->structure[i] = '.';
  }

  state->partial_energy = partial_energy;

  return state;
}

/*---------------------------------------------------------------------------*/

PRIVATE STATE *
copy_state(STATE * state)
{
  STATE *new_state;
  void *after;
  INTERVAL *new_interval, *next;

  new_state = lst_newnode(sizeof(STATE));
  new_state->Intervals = lst_init();
  new_state->partial_energy = state->partial_energy;
  /* new_state->best_energy = state->best_energy; */

  if (state->Intervals->count) {
    after = LST_HEAD(new_state->Intervals);
    for ( next = lst_first(state->Intervals); next; next = lst_next(next))
      {
        new_interval = lst_newnode(sizeof(INTERVAL));
        *new_interval = *next;
        lst_insertafter(new_state->Intervals, new_interval, after);
        after = new_interval;
      }
  }
  new_state->structure = strdup(state->structure);
  if (!new_state->structure) vrna_message_error("out of memory");
  return new_state;
}

/*---------------------------------------------------------------------------*/

/*@unused @*/ PRIVATE void
print_state(STATE * state)
{
  INTERVAL *next;

  if (state->Intervals->count)
    {
      printf("%d intervals:\n", state->Intervals->count);
      for (next = lst_first(state->Intervals); next; next = lst_next(next))
        {
          printf("[%d,%d],%d ", next->i, next->j, next->array_flag);
        }
      printf("\n");
    }
  printf("partial structure: %s\n", state->structure);
  printf("\n");
  printf(" partial_energy: %d\n", state->partial_energy);
  /* printf(" best_energy: %d\n", state->best_energy); */
  (void) fflush(stdout);
}

/*---------------------------------------------------------------------------*/

/*@unused @*/ PRIVATE void
print_stack(LIST * list)
{
  void *rec;

  printf("================\n");
  printf("%d states\n", list->count);
  for (rec = lst_first(list); rec; rec = lst_next(rec))
    {
      printf("state-----------\n");
      print_state(rec);
    }
  printf("================\n");
}

/*---------------------------------------------------------------------------*/

PRIVATE LIST *
make_list(void)
{
  return lst_init();
}

/*---------------------------------------------------------------------------*/

PRIVATE void
push(LIST * list, void *data)
{
  lst_insertafter(list, data, LST_HEAD(list));
}

/* PRIVATE void */
/* push_stack(STATE *state) { */ /* keep the stack sorted by energy */
/*   STATE *after, *next; */
/*   nopush = false; */
/*   next = after = LST_HEAD(Stack); */
/*   while ( next = lst_next(next)) { */
/*     if ( next->best_energy >= state->best_energy ) break; */
/*     after = next; */
/*   } */
/*   lst_insertafter(Stack, state, after); */
/* } */

/*---------------------------------------------------------------------------*/

PRIVATE void *
pop(LIST * list)
{
  void *data;

  data = lst_deletenext(list, LST_HEAD(list));
  return data;
}

/*---------------------------------------------------------------------------*/
/*auxiliary routines---------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

PRIVATE int
best_attainable_energy( vrna_fold_compound_t *vc,
                        STATE *state){

  /* evaluation of best possible energy attainable within remaining intervals */

  register int sum;
  INTERVAL        *next;
  vrna_md_t       *md;
  vrna_mx_mfe_t   *matrices;
  int             *indx;

  md        = &(vc->params->model_details);
  matrices  = vc->matrices;
  indx      = vc->jindx;

  sum = state->partial_energy;  /* energy of already found elements */

  for (next = lst_first(state->Intervals); next; next = lst_next(next))
    {
      if (next->array_flag == 0)
        sum += (md->circ) ? matrices->Fc : matrices->f5[next->j];
      else if (next->array_flag == 1)
        sum += matrices->fML[indx[next->j] + next->i];
      else if (next->array_flag == 2)
        sum += matrices->c[indx[next->j] + next->i];
      else if (next->array_flag == 3)
        sum += matrices->fM1[indx[next->j] + next->i];
      else if (next->array_flag == 4)
        sum += matrices->fc[next->i];
      else if (next->array_flag == 5)
        sum += matrices->fc[next->j];
      else if (next->array_flag == 6)
        sum += matrices->ggg[indx[next->j] + next->i];
    }

  return sum;
}

/*---------------------------------------------------------------------------*/

PRIVATE void
push_back(LIST *Stack, STATE * state)
{
  push(Stack, copy_state(state));
  return;
}

/*---------------------------------------------------------------------------*/

PRIVATE char*
get_structure(STATE * state)
{
  char* structure;

  structure = strdup(state->structure);
  return structure;
}

/*---------------------------------------------------------------------------*/
PRIVATE int
compare(const void *solution1, const void *solution2)
{
  if (((SOLUTION *) solution1)->energy > ((SOLUTION *) solution2)->energy)
    return 1;
  if (((SOLUTION *) solution1)->energy < ((SOLUTION *) solution2)->energy)
    return -1;
  return strcmp(((SOLUTION *) solution1)->structure,
                ((SOLUTION *) solution2)->structure);
}

/*---------------------------------------------------------------------------*/

PRIVATE void make_output(SOLUTION *SL, int cp, FILE *fp)  /* prints stuff */
{
  SOLUTION *sol;

  for (sol = SL; sol->structure!=NULL; sol++){
    char *e_string = vrna_strdup_printf(" %6.2f", sol->energy);
    print_structure(fp, sol->structure, e_string);
    free(e_string);
  }
}

PRIVATE STATE *
derive_new_state( int i,
                  int j,
                  STATE *s,
                  int e,
                  int flag){

  STATE     *s_new  = copy_state(s);
  INTERVAL  *ival   = make_interval(i, j, flag);
  push(s_new->Intervals, ival);

  s_new->partial_energy += e;

  return s_new;
}

PRIVATE void
fork_state( int i,
            int j,
            STATE *s,
            int e,
            int flag,
            subopt_env *env){

  STATE *s_new = derive_new_state(i, j, s, e, flag);
  push(env->Stack, s_new);
  env->nopush = false;
}

PRIVATE void
fork_int_state( int i, int j,
                int p, int q,
                STATE *s,
                int e,
                subopt_env *env){

  STATE *s_new = derive_new_state(p, q, s, e, 2);
  make_pair(i, j, s_new);
  make_pair(p, q, s_new);
  push(env->Stack, s_new);
  env->nopush = false;
}

PRIVATE void
fork_state_pair(int i,
                int j,
                STATE *s,
                int e,
                subopt_env *env){

  STATE     *new_state;

  new_state = copy_state(s);
  make_pair(i, j, new_state);
  new_state->partial_energy += e;
  push(env->Stack, new_state);
  env->nopush = false;
}

PRIVATE void
fork_two_states_pair( int i,
                      int j,
                      int k,
                      STATE *s,
                      int e,
                      int flag1,
                      int flag2,
                      subopt_env *env){

  INTERVAL  *interval1, *interval2;
  STATE     *new_state;

  new_state = copy_state(s);
  interval1 = make_interval(i+1, k-1, flag1);
  interval2 = make_interval(k, j-1, flag2);
  if (k-i < j-k) { /* push larger interval first */
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
fork_two_states(int i,
                int j,
                int p,
                int q,
                STATE *s,
                int e,
                int flag1,
                int flag2,
                subopt_env *env){

  INTERVAL  *interval1, *interval2;
  STATE     *new_state;

  new_state = copy_state(s);
  interval1 = make_interval(i, j, flag1);
  interval2 = make_interval(p, q, flag2);

  if((j - i) < (q - p)){
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

/*---------------------------------------------------------------------------*/
/* start of subopt backtracking ---------------------------------------------*/
/*---------------------------------------------------------------------------*/

PUBLIC SOLUTION *
vrna_subopt(vrna_fold_compound_t *vc,
            int delta,
            int sorted,
            FILE *fp){

  SOLUTION      *SolutionList;
  unsigned long max_sol, n_sol;
  struct old_subopt_dat data;

  data.SolutionList = NULL;
  data.max_sol      = 128;
  data.n_sol        = 0;
  data.fp           = fp;

  if(vc){
    /* SolutionList stores the suboptimal structures found */

    data.SolutionList = (SOLUTION *) vrna_alloc(data.max_sol*sizeof(SOLUTION));

    /* end initialize ------------------------------------------------------- */

    if (fp) {
      float min_en;
      char  *SeQ, *energies = NULL;
      if(vc->cutpoint > 0)
        min_en = vrna_mfe_dimer(vc, NULL);
      else
        min_en = vrna_mfe(vc, NULL);

      SeQ = vrna_cut_point_insert(vc->sequence, vc->cutpoint);
      energies = vrna_strdup_printf(" %6.2f %6.2f", min_en, (float)delta/100.);
      print_structure(fp, SeQ, energies);
      free(SeQ);
      free(energies);

      vrna_mx_mfe_free(vc);
    }
    /* call subopt() */
    vrna_subopt_cb(vc, delta, (!sorted && fp) ? &old_subopt_print : &old_subopt_store, (void *)&data);

    if(sorted){
      /* sort structures by energy */
      if(data.n_sol > 0)
        qsort(data.SolutionList, data.n_sol - 1, sizeof(SOLUTION), compare);
      if(fp)
        make_output(data.SolutionList, vc->cutpoint, fp);
    }

    if(fp){ /* we've printed everything -- free solutions */
      SOLUTION *sol;
      for(sol = data.SolutionList; sol->structure != NULL; sol++)
        free(sol->structure);
      free(data.SolutionList);
      data.SolutionList = NULL;
    }
  }

  return data.SolutionList;
}

PUBLIC void
vrna_subopt_cb( vrna_fold_compound_t *vc,
                int delta,
                vrna_subopt_callback *cb,
                void *data){

  subopt_env    *env;
  STATE         *state;
  INTERVAL      *interval;
  int           maxlevel, count, partial_energy, old_dangles, logML, dangle_model, length, circular, threshold, cp;
  double        structure_energy, min_en, eprint;
  char          *struc, *structure, *sequence;
  float         correction;
  vrna_param_t  *P;
  vrna_md_t     *md;
  int           minimal_energy;
  int           Fc;
  int           *f5;

  vrna_fold_compound_prepare(vc, VRNA_OPTION_MFE | VRNA_OPTION_HYBRID);

  sequence            = vc->sequence;
  length              = vc->length;
  cp                  = vc->cutpoint;
  P                   = vc->params;
  md                  = &(P->model_details);

  /* do mfe folding to get fill arrays and get ground state energy  */
  /* in case dangles is neither 0 or 2, set dangles=2 while folding */

  circular    = md->circ;
  logML       = md->logML;
  old_dangles = dangle_model = md->dangles;

  if(md->uniq_ML != 1)  /* failsafe mechanism to enforce valid fM1 array */
    md->uniq_ML = 1;

  /* temporarily set dangles to 2 if necessary */
  if((md->dangles != 0) && (md->dangles != 2))
    md->dangles = 2;

  struc = (char *)vrna_alloc(sizeof(char) * (length + 1));

  if(circular){
    min_en = vrna_mfe(vc, struc);
    Fc  = vc->matrices->Fc;
    f5    = vc->matrices->f5;

    /* restore dangle model */
    md->dangles = old_dangles;
    /* re-evaluate in case we're using logML etc */
    min_en = vrna_eval_structure(vc, struc);
  } else {
    min_en = vrna_mfe_dimer(vc, struc);

    f5    = vc->matrices->f5;

    /* restore dangle model */
    md->dangles = old_dangles;
    /* re-evaluate in case we're using logML etc */
    min_en = vrna_eval_structure(vc, struc);
  }

  free(struc);
  eprint = print_energy + min_en;

  correction = (min_en < 0) ? -0.1 : 0.1;

  /* Initialize ------------------------------------------------------------ */

  maxlevel = 0;
  count = 0;
  partial_energy = 0;

  /* Initialize the stack ------------------------------------------------- */

  minimal_energy = (circular) ? Fc : f5[length];
  threshold = minimal_energy + delta;
  if(threshold > INF){
    vrna_message_warning("Energy range too high, limiting to reasonable value");
    threshold = INF-EMAX;
  }

  /* init env data structure */
  env = (subopt_env *)vrna_alloc(sizeof(subopt_env));
  env->Stack      = NULL;
  env->nopush     = true;
  env->Stack      = make_list();                                                   /* anchor */
  env->Intervals  = make_list();                                   /* initial state: */
  interval   = make_interval(1, length, 0);          /* interval [1,length,0] */
  push(env->Intervals, interval);
  env->nopush = false;
  state  = make_state(env->Intervals, NULL, partial_energy,0, length);
  /* state->best_energy = minimal_energy; */
  push(env->Stack, state);
  env->nopush = false;

  /* end initialize ------------------------------------------------------- */


  while (1) {                    /* forever, til nothing remains on stack */

    maxlevel = (env->Stack->count > maxlevel ? env->Stack->count : maxlevel);

    if (LST_EMPTY (env->Stack))                   /* we are done! clean up and quit */
      {
        /* fprintf(stderr, "maxlevel: %d\n", maxlevel); */

        lst_kill(env->Stack, free_state_node);

        cb(NULL, 0, data); /* NULL (last time to call callback function */

        break;
      }

    /* pop the last element ---------------------------------------------- */

    state = pop(env->Stack);                       /* current state to work with */

    if (LST_EMPTY(state->Intervals))
      {
        int e;
        /* state has no intervals left: we got a solution */

        count++;
        structure = get_structure(state);
        structure_energy = state->partial_energy / 100.;

#ifdef CHECK_ENERGY
        structure_energy = vrna_eval_structure(vc, structure);

        if (!logML)
          if ((double) (state->partial_energy / 100.) != structure_energy) {
            vrna_message_error("%s %6.2f %6.2f",
                                      structure,
                                      state->partial_energy / 100.,
                                      structure_energy );
            exit(1);
          }
#endif
        if (logML || (dangle_model==1) || (dangle_model==3)) { /* recalc energy */
          structure_energy = vrna_eval_structure(vc, structure);
        }

        e = (int) ((structure_energy-min_en)*10. - correction); /* avoid rounding errors */
        if (e>MAXDOS) e=MAXDOS;
        density_of_states[e]++;
        if(structure_energy <= eprint){
          char *outstruct = vrna_cut_point_insert(structure, cp);
          cb((const char *)outstruct, structure_energy, data);
          free(outstruct);
        }
        free(structure);
      }
    else {
      /* get (and remove) next interval of state to analyze */

      interval = pop(state->Intervals);
      scan_interval(vc, interval->i, interval->j, interval->array_flag, threshold, state, env);

      free_interval_node(interval);        /* free the current interval */
    }

    free_state_node(state);                     /* free the current state */
  } /* end of while (1) */

  /* cleanup memory */
  free(env);
}


PRIVATE void
scan_interval(vrna_fold_compound_t *vc,
              int i,
              int j,
              int array_flag,
              int threshold,
              STATE * state,
              subopt_env *env){

  /* real backtrack routine */

  /* array_flag = 0:  trace back in f5-array  */
  /* array_flag = 1:  trace back in fML-array */
  /* array_flag = 2:  trace back in repeat()  */
  /* array_flag = 3:  trace back in fM1-array */

  STATE           *new_state, *temp_state;
  INTERVAL        *new_interval;
  vrna_param_t    *P;
  vrna_md_t       *md;
  register int    k, fi, cij, ij;
  register int    type;
  register int    dangle_model;
  register int    noLP;
  int             element_energy, best_energy;
  int             *fc, *f5, *c, *fML, *fM1, *ggg;
  int             FcH, FcI, FcM, *fM2;
  int             length, *indx, *rtype, circular, with_gquad, turn, cp;
  char            *ptype;
  short           *S1;
  char            *hard_constraints, hc_decompose;
  vrna_hc_t       *hc;
  vrna_sc_t       *sc;

  length    = vc->length;
  cp        = vc->cutpoint;
  indx      = vc->jindx;
  ptype     = vc->ptype;
  S1        = vc->sequence_encoding;
  P         = vc->params;
  md        = &(P->model_details);
  rtype     = &(md->rtype[0]);

  dangle_model  = md->dangles;
  noLP          = md->noLP;
  circular      = md->circ;
  with_gquad    = md->gquad;
  turn          = md->min_loop_size;

  fc  = vc->matrices->fc;
  f5  = vc->matrices->f5;
  c   = vc->matrices->c;
  fML = vc->matrices->fML;
  fM1 = vc->matrices->fM1;
  ggg = vc->matrices->ggg;
  FcH = vc->matrices->FcH;
  FcI = vc->matrices->FcI;
  FcM = vc->matrices->FcM;
  fM2 = vc->matrices->fM2;

  hc                = vc->hc;
  hard_constraints  = hc->matrix;

  sc                = vc->sc;

  best_energy = best_attainable_energy(vc, state);  /* .. on remaining intervals */
  env->nopush = true;

  if ((i > 1) && (!array_flag))
    vrna_message_error ("Error while backtracking!");

  if (j < i + turn + 1 && ON_SAME_STRAND(i,j,cp)) { /* minimal structure element */
    if(array_flag == 0){
      /* do not forget to add f5[j], since it may contain pseudo energies from soft constraining */
      state->partial_energy += f5[j];
    }
    if (env->nopush){
      push_back(env->Stack, state);
      env->nopush = false;
    }
    return;
  }

  ij = indx[j] + i;

  /* 13131313131313131313131313131313131313131313131313131313131313131313131 */
  if (array_flag == 3 || array_flag == 1) {
    /* array_flag = 3: interval i,j was generated during */
    /*                 a multiloop decomposition using array fM1 in repeat() */
    /*                 or in this block */

    /* array_flag = 1: interval i,j was generated from a */
    /*                 stack, bulge, or internal loop in repeat() */
    /*                 or in this block */

    if(hc->up_ml[j]){
      if (array_flag == 3)
        fi = fM1[indx[j-1] + i] + P->MLbase;
      else
        fi = fML[indx[j-1] + i] + P->MLbase;

      if(sc){
        if(sc->energy_up)
          fi += sc->energy_up[j][1];
        if(sc->f)
          fi += sc->f(i, j, i, j - 1, VRNA_DECOMP_ML_ML, sc->data);
      }

      if ((fi + best_energy <= threshold)&&(ON_SAME_STRAND(j-1,j, cp))) {
        /* no basepair, nibbling of 3'-end */
        fork_state(i, j-1, state, P->MLbase, array_flag, env);
      }
    }

    hc_decompose  = hard_constraints[ij];

    if (hc_decompose & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) { /* i,j may pair */
      cij = c[ij];

      type = ptype[ij];

      if(type == 0)
        type = 7;

      switch(dangle_model){
        case 0:   element_energy = E_MLstem(type, -1, -1, P);
                  break;
        default:  element_energy = E_MLstem(type,
                                  (((i > 1)&&(ON_SAME_STRAND(i-1,i,cp))) || circular)       ? S1[i-1] : -1,
                                  (((j < length)&&(ON_SAME_STRAND(j,j+1,cp))) || circular)  ? S1[j+1] : -1,
                                  P);
                  break;
      }

      if(sc){
/* should be unnecessary
        if(sc->energy_bp)
          element_energy += sc->energy_bp[ij];
*/
        if(sc->f)
          element_energy += sc->f(i, j, i, j, VRNA_DECOMP_ML_STEM, sc->data);
      }
      cij += element_energy;

      if (cij + best_energy <= threshold)
        repeat(vc, i, j, state, element_energy, 0, best_energy, threshold, env);
    } else if (with_gquad){
      element_energy = E_MLstem(0, -1, -1, P);
      cij = ggg[ij] + element_energy;
      if(cij + best_energy <= threshold)
        repeat_gquad(vc, i, j, state, element_energy, 0, best_energy, threshold, env);
    }
  }                                   /* array_flag == 3 || array_flag == 1 */

  /* 11111111111111111111111111111111111111111111111111111111111111111111111 */

  if (array_flag == 1) {
    /* array_flag = 1:                   interval i,j was generated from a */
    /*                          stack, bulge, or internal loop in repeat() */
    /*                          or in this block */

    int stopp, k1j;
    if ((ON_SAME_STRAND(i-1,i,cp))&&(ON_SAME_STRAND(j,j+1,cp))) { /*backtrack in FML only if multiloop is possible*/
      for ( k = i+turn+1 ; k <= j-1-turn ; k++) {
        /* Multiloop decomposition if i,j contains more than 1 stack */

        if(with_gquad){
          if(ON_SAME_STRAND(k, k+1, cp)){
            element_energy = E_MLstem(0, -1, -1, P);
            if(fML[indx[k]+i] + ggg[indx[j] + k + 1] + element_energy + best_energy <= threshold){
              
              temp_state = derive_new_state(i, k, state, 0, array_flag);
              env->nopush = false;
              repeat_gquad(vc, k+1, j, temp_state, element_energy, fML[indx[k]+i], best_energy, threshold, env);
              free_state_node(temp_state);
            }
          }
        }

        k1j = indx[j] + k + 1;

        if(hard_constraints[k1j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
          short s5, s3;
          type = ptype[k1j];

          if(type == 0)
            type = 7;

          switch(dangle_model){
            case 0:   s5 = s3 = -1;
                      break;
            default:  s5 = (ON_SAME_STRAND(i-1,i,cp)) ? S1[k] : -1;
                      s3 = (ON_SAME_STRAND(j,j+1,cp)) ? S1[j+1] : -1;
                      break;
          }

          element_energy = E_MLstem(type, s5, s3, P);

          if(sc){
/* should be unnecessary
            if(sc->energy_bp)
              element_energy += sc->energy_bp[k1j];
*/
            if(sc->f)
              element_energy += sc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_STEM, sc->data);
          }

          if(ON_SAME_STRAND(k, k+1, cp)){
            if(fML[indx[k]+i] + c[k1j] + element_energy + best_energy <= threshold){
              temp_state  = derive_new_state(i, k, state, 0, array_flag);
              env->nopush = false;
              repeat(vc, k+1, j, temp_state, element_energy, fML[indx[k]+i], best_energy, threshold, env);
              free_state_node(temp_state);
            }
          }
        }
      }
    }

    stopp=(cp>0)? (cp-2):(length); /*if cp -1: k on cut, => no ml*/
    stopp=MIN2(stopp, j-1-turn);
    if (i>cp) stopp=j-1-turn;
    else if (i==cp) stopp=0;   /*not a multi loop*/

    int up = 1;
    for(k = i; k <= stopp; k++, up++){

      if(hc->up_ml[i] >= up){
        k1j = indx[j] + k + 1;

        /* Multiloop decomposition if i,j contains only 1 stack */
        if(with_gquad){
          element_energy = E_MLstem(0, -1, -1, P) + P->MLbase * up;

          if(sc){
            if(sc->energy_up)
              element_energy += sc->energy_up[i][up];
          }

          if(ggg[k1j] + element_energy + best_energy <= threshold)
            repeat_gquad(vc, k+1, j, state, element_energy, 0, best_energy, threshold, env);
        }

        if(hard_constraints[k1j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
          int s5, s3;
          type = ptype[k1j];

          if(type == 0)
            type = 7;

          switch(dangle_model){
            case 0:   s5 = s3 = -1;
                      break;
            default:  s5 = (ON_SAME_STRAND(k-1,k,cp)) ? S1[k] : -1;
                      s3 = (ON_SAME_STRAND(j,j+1,cp)) ? S1[j+1] : -1;
                      break;
          }

          element_energy = E_MLstem(type, s5, s3, P);

          element_energy += P->MLbase * up;

          if(sc){
            if(sc->energy_up)
              element_energy += sc->energy_up[i][up];

/* should be unnecessary
            if(sc->energy_bp)
              element_energy += sc->energy_bp[k1j];
*/
          }

          if (c[k1j] + element_energy + best_energy <= threshold)
            repeat(vc, k+1, j, state, element_energy, 0, best_energy, threshold, env);
        }
      }
    }
  } /* array_flag == 1 */

  /* 22222222222222222222222222222222222222222222222222 */
  /*                                                    */
  /* array_flag = 2:  interval i,j was generated from a */
  /* stack, bulge, or internal loop in repeat()         */
  /*                                                    */
  /* 22222222222222222222222222222222222222222222222222 */
  if(array_flag == 2){
    repeat(vc, i, j, state, 0, 0, best_energy, threshold, env);

    if (env->nopush){
      if (!noLP){
        vrna_message_warning("%d,%d\nOops, no solution in repeat!", i, j);
      }
    }
    return;
  }

  /* 00000000000000000000000000000000000000000000000000 */
  /*                                                    */
  /* array_flag = 0:  interval i,j was found while      */
  /* tracing back through f5-array and c-array          */
  /* or within this block                               */
  /*                                                    */
  /* 00000000000000000000000000000000000000000000000000 */
  if((array_flag == 0) && !circular){
    int s5, s3, kj, tmp_en;

    if(hc->up_ext[j]){
      tmp_en = 0;
      if(sc){
        if(sc->energy_up)
          tmp_en += sc->energy_up[j][1];
      }
      if (f5[j-1] + tmp_en + best_energy <= threshold) {
        /* no basepair, nibbling of 3'-end */
        fork_state(i, j-1, state, tmp_en, 0, env);
      }
    }

    for (k = j-turn-1; k > 1; k--) {
      kj = indx[j] + k;

      if(with_gquad){
        if(ON_SAME_STRAND(k,j,cp)){
          element_energy = 0;
          if(f5[k-1] + ggg[kj] + element_energy + best_energy <= threshold){
            temp_state = derive_new_state(1, k-1, state, 0, 0);
            env->nopush = false;
            /* backtrace the quadruplex */
            repeat_gquad(vc, k, j, temp_state, element_energy, f5[k-1], best_energy, threshold, env);
            free_state_node(temp_state);
          }
        }
      }

      if(hard_constraints[kj] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP){
        type = ptype[kj];

        if(type == 0)
          type = 7;

        /* k and j pair */
        switch(dangle_model){
          case 0:   s5 = s3 = -1;
                    break;
          default:  s5 = (ON_SAME_STRAND(k-1,k,cp)) ? S1[k-1] : -1;
                    s3 = ((j < length)&&(ON_SAME_STRAND(j,j+1,cp))) ? S1[j+1] : -1;
                    break;
        }
        
        element_energy = E_ExtLoop(type, s5, s3, P);

        if (!(ON_SAME_STRAND(k,j,cp)))/*&&(state->is_duplex==0))*/ {
          element_energy+=P->DuplexInit;
          /*state->is_duplex=1;*/
        }

        if (f5[k-1] + c[kj] + element_energy + best_energy <= threshold){
          temp_state = derive_new_state(1, k-1, state, 0, 0);
          env->nopush = false;
          repeat(vc, k, j, temp_state, element_energy, f5[k-1], best_energy, threshold, env);
          free_state_node(temp_state);
        }
      }
    }

    kj = indx[j] + 1;

    if(with_gquad){
      if(ON_SAME_STRAND(k,j,cp)){
        element_energy = 0;
        if(ggg[kj] + element_energy + best_energy <= threshold){
          /* backtrace the quadruplex */
          repeat_gquad(vc, 1, j, state, element_energy, 0, best_energy, threshold, env);
        }
      }
    }

    if(hard_constraints[kj] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP){
      type  = ptype[kj];
      s5    = -1;

      if(type == 0)
        type = 7;

      switch(dangle_model){
        case 0:   s3 = -1;
                  break;
        default:  s3 = (j < length) && (ON_SAME_STRAND(j,j+1,cp)) ? S1[j+1] : -1;
                  break;
      }

      element_energy = E_ExtLoop(type, s5, s3, P);

      if(!(ON_SAME_STRAND(1,j,cp)))
        element_energy += P->DuplexInit;
      if (c[kj] + element_energy + best_energy <= threshold)
        repeat(vc, 1, j, state, element_energy, 0, best_energy, threshold, env);
    }
  } /* end array_flag == 0 && !circular*/
  /* or do we subopt circular? */
  else if(array_flag == 0){
    int k, l, p, q, tmp_en;
    /* if we've done everything right, we will never reach this case more than once   */
    /* right after the initilization of the stack with ([1,n], empty, 0)              */
    /* lets check, if we can have an open chain without breaking the threshold        */
    /* this is an ugly work-arround cause in case of an open chain we do not have to  */
    /* backtrack anything further...                                                  */
    if(hc->up_ext[1] >= length){
      tmp_en = 0;

      if(sc){
        if(sc->energy_up)
          tmp_en += sc->energy_up[1][length];
      }

      if(tmp_en <= threshold){
        new_state = derive_new_state(1,2,state,0,0);
        new_state->partial_energy = 0;
        push(env->Stack, new_state);
        env->nopush = false;
      }
    }

    /* ok, lets check if we can do an exterior hairpin without breaking the threshold */
    /* best energy should be 0 if we are here                                         */
    if(FcH + best_energy <= threshold){
      /* lets search for all exterior hairpin cases, that fit into our threshold barrier  */
      /* we use index k,l to avoid confusion with i,j index of our state...               */
      /* if we reach here, i should be 1 and j should be n respectively                   */
      for(k=i; k<j; k++){
        if(hc->up_hp[1] < k)
          break;

        for (l=j; l >= k + turn + 1; l--){
          int kl, tmpE;

          kl    = indx[l] + k;
          tmpE  = vrna_E_hp_loop(vc, l, k);

          if(c[kl] + tmpE + best_energy <= threshold){
            /* what we really have to do is something like this, isn't it?  */
            /* we have to create a new state, with interval [k,l], then we  */
            /* add our loop energy as initial energy of this state and put  */
            /* the state onto the stack R... for further refinement...      */
            /* we also denote this new interval to be scanned in C          */
            fork_state(k, l, state, tmpE, 2, env);
          }
        }
      }
    }

    /* now lets see, if we can do an exterior interior loop without breaking the threshold */
    if(FcI + best_energy <= threshold){
      /* now we search for our exterior interior loop possibilities */
      for(k=i; k<j; k++){
        for (l=j; l >= k + turn + 1; l--){
          int kl, type, tmpE;

          kl    = indx[l]+k;        /* just confusing these indices ;-) */

          if(hard_constraints[kl] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP){
            type  = ptype[kl];
            type  = rtype[type];

            if(type == 0)
              type = 7;

            for (p = l+1; p < j ; p++){
              int u1, qmin;
              u1 = p-l-1;
              if (u1+k-1>MAXLOOP) break;
              if (hc->up_int[l+1] < u1) break;

              qmin = u1+k-1+j-MAXLOOP;
              if(qmin<p+turn+1) qmin = p+turn+1;

              for(q = j; q >= qmin; q--){
                int u2, type_2;

                if(hc->up_int[q+1] < (j - q + k - 1))
                  break;

                if(hard_constraints[indx[q]+p] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP){
                  type_2 = rtype[ptype[indx[q]+p]];

                  if(type_2 == 0)
                    type_2 = 7;

                  u2 = k-1 + j-q;
                  if(u1+u2>MAXLOOP) continue;
                  tmpE = E_IntLoop(u1, u2, type, type_2, S1[l+1], S1[k-1], S1[p-1], S1[q+1], P);
                  
                  if(sc){
                    if(sc->energy_up)
                      tmpE += sc->energy_up[l+1][p-l-1]
                              + sc->energy_up[q+1][j-q]
                              + sc->energy_up[1][k-1];
                    
                    if(sc->energy_stack)
                      if(u1 + u2 == 0)
                        tmpE += sc->energy_stack[k]
                                + sc->energy_stack[l]
                                + sc->energy_stack[p]
                                + sc->energy_stack[q];
                  }

                  if(c[kl] + c[indx[q]+p] + tmpE + best_energy <= threshold){
                    /* ok, similar to the hairpin stuff, we add new states onto the stack R */
                    /* but in contrast to the hairpin decomposition, we have to add two new */
                    /* intervals, enclosed by k,l and p,q respectively and we also have to  */
                    /* add the partial energy, that comes from the exterior interior loop   */
                    fork_two_states(k, l, p, q, state, tmpE, 2, 2, env);
                  }
                }
              }
            }
          }
        }
      }
    }

    /* and last but not least, we have a look, if we can do an exterior multiloop within the energy threshold */
    if(FcM <= threshold){
      /* this decomposition will be somehow more complicated...so lets see what we do here...           */
      /* first we want to find out which split inidices we can use without exceeding the threshold */
      int tmpE2;
      for (k=turn+1; k<j-2*turn; k++){
        tmpE2 = fML[indx[k]+1]+fM2[k+1]+P->MLclosing;
        if(tmpE2 + best_energy <= threshold){
          /* grmpfh, we have found a possible split index k so we have to split fM2 and fML now */
          /* lets do it first in fM2 anyway */
          for(l=k+turn+2; l<j-turn-1; l++){
            tmpE2 = fM1[indx[l]+k+1] + fM1[indx[j]+l+1];
            if(tmpE2 + fML[indx[k]+1] + P->MLclosing <= threshold){
              /* we've (hopefully) found a valid decomposition of fM2 and therefor we have all */
              /* three intervals for our new state to be pushed on stack R */
              new_state = copy_state(state);

              /* first interval leads for search in fML array */
              new_interval = make_interval(1, k, 1);
              push(new_state->Intervals, new_interval);
              env->nopush = false;

              /* next, we have the first interval that has to be traced in fM1 */
              new_interval = make_interval(k+1, l, 3);
              push(new_state->Intervals, new_interval);
              env->nopush = false;

              /* and the last of our three intervals is also one to be traced within fM1 array... */
              new_interval = make_interval(l+1, j, 3);
              push(new_state->Intervals, new_interval);
              env->nopush = false;

              /* mmh, we add the energy for closing the multiloop now... */
              new_state->partial_energy += P->MLclosing;
              /* next we push our state onto the R stack */
              push(env->Stack, new_state);
              env->nopush = false;

            }
            /* else we search further... */
          }

          /* ok, we have to decompose fML now... */
        }
      }
    }
  }        /* thats all folks for the circular case... */

  /* 44444444444444444444444444444444444444444444444444 */
  /*                                                    */
  /* array_flag = 4:  interval i,j was found while      */
  /* tracing back through fc-array smaller than than cp */
  /* or within this block                               */
  /*                                                    */
  /* 44444444444444444444444444444444444444444444444444 */
  if (array_flag == 4) {
    int ik, s5, s3, tmp_en;

    if(hc->up_ext[i]){
      tmp_en = 0;

      if(sc){
        if(sc->energy_up)
          tmp_en += sc->energy_up[i][1];
      }

      if (fc[i+1] + tmp_en + best_energy <= threshold) {
        /* no basepair, nibbling of 5'-end */
        fork_state(i+1, j, state, tmp_en, 4, env);
      }
    }

    for (k = i+TURN+1; k < j; k++) {

      ik = indx[k] + i;

      if(with_gquad){
        if(fc[k+1] + ggg[ik] + best_energy <= threshold){
          temp_state = derive_new_state(k+1, j, state, 0, 4);
          env->nopush = false;
          repeat_gquad(vc, i, k, temp_state, 0, fc[k+1], best_energy, threshold, env);
          free_state_node(temp_state);
        }
      }

      if(hard_constraints[ik] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP){
        type = ptype[ik];

        if(type == 0)
          type = 7;

        switch(dangle_model){
          case 0:   s5 = s3 = -1;
                    break;
          default:  s5 = (i > 1) ? S1[i-1]: -1;
                    s3 = S1[k+1];
                    break;
        }

        element_energy = E_ExtLoop(type, s5, s3, P);

/* should be unnecessary
        if(sc){
          if(sc->energy_bp)
            element_energy += sc->energy_bp[ik];
        }
*/

        if (fc[k+1] + c[ik] + element_energy + best_energy <= threshold){
          temp_state = derive_new_state(k+1, j, state, 0, 4);
          env->nopush = false;
          repeat(vc, i, k, temp_state, element_energy, fc[k+1], best_energy, threshold, env);
          free_state_node(temp_state);
        }
      }
    }

    ik = indx[cp -1] + i; /* indx[j] + i; */

    if(with_gquad){
      if(ggg[ik] + best_energy <= threshold)
        repeat_gquad(vc, i, cp - 1, state, 0, 0, best_energy, threshold, env);
    }

    if(hard_constraints[ik] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP){
      type  = ptype[ik];
      s3    = -1;

      if(type == 0)
        type = 7;

      switch(dangle_model){
        case 0:   s5 = -1;
                  break;
        default:  s5 = (i>1) ? S1[i-1] : -1;
                  break;
      }

      element_energy = E_ExtLoop(type, s5, s3, P);

/* should be unnecessary
      if(sc){
        if(sc->energy_bp)
          element_energy += sc->energy_bp[ik];
      }
*/

      if(c[ik] + element_energy + best_energy <= threshold)
        repeat(vc, i, cp-1, state, element_energy, 0, best_energy, threshold, env);
    }
  } /* array_flag == 4 */

  /* 55555555555555555555555555555555555555555555555555 */
  /*                                                    */
  /* array_flag = 5:  interval cp=i,j was found while   */
  /* tracing back through fc-array greater than cp      */
  /* or within this block                               */
  /*                                                    */
  /* 55555555555555555555555555555555555555555555555555 */
  if (array_flag == 5) {
    int kj, s5, s3, tmp_en;

    if(hc->up_ext[j]){
      tmp_en = 0;

      if(sc){
        if(sc->energy_up)
          tmp_en += sc->energy_up[j][1];
      }

      if (fc[j-1] + tmp_en + best_energy <= threshold) {
        /* no basepair, nibbling of 3'-end */
        fork_state(i, j-1, state, tmp_en, 5, env);
      }
    }


    for (k = j-TURN-1; k > i; k--) {
      kj = indx[j] + k;

      if(with_gquad){
        if(fc[k-1] + ggg[kj] + best_energy <= threshold){
          temp_state = derive_new_state(i, k-1, state, 0, 5);
          env->nopush = false;
          repeat_gquad(vc, k, j, temp_state, 0, fc[k-1], best_energy, threshold, env);
          free_state_node(temp_state);
        }
      }

      if(hard_constraints[kj] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP){
        type            = ptype[kj];
        element_energy  = 0;

        if(type == 0)
          type = 7;

        switch(dangle_model){
          case 0:   s3 = s5 = -1;
                    break;
          default:  s5 = S1[k-1];
                    s3 = (j < length) ? S1[j+1] : -1;
                    break;
        }

        element_energy = E_ExtLoop(type, s5, s3, P);

/*  should be unnecessary
        if(sc){
          if(sc->energy_bp)
            element_energy += sc->energy_bp[kj];
        }
*/

        if (fc[k-1] + c[kj] + element_energy + best_energy <= threshold) {
          temp_state = derive_new_state(i, k-1, state, 0, 5);
          env->nopush = false;
          repeat(vc, k, j, temp_state, element_energy, fc[k-1], best_energy, threshold, env);
          free_state_node(temp_state);
        }
      }
    }

    kj = indx[j] + cp; /* indx[j] + i; */

    if(with_gquad){
      if(ggg[kj] + best_energy <= threshold)
        repeat_gquad(vc, cp, j, state, 0, 0, best_energy, threshold, env);
    }

    if(hard_constraints[kj] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP){
      type  = ptype[kj];
      s5    = -1;

      if(type == 0)
        type = 7;

      switch(dangle_model){
        case 0:   s3 = -1;
                  break;
        default:  s3 = (j<length) ? S1[j+1] : -1;
                  break;
      }

      element_energy = E_ExtLoop(type, s5, s3, P);

      if (c[kj] + element_energy + best_energy <= threshold)
        repeat(vc, cp, j, state, element_energy, 0, best_energy, threshold, env);
    }
  } /* array_flag == 5 */

  if (array_flag == 6) { /* we have a gquad */
    repeat_gquad(vc, i, j, state, 0, 0, best_energy, threshold, env);
    if (env->nopush){
      vrna_message_warning("%d,%d\nOops, no solution in gquad-repeat!", i, j);
    }
    return;
  }

  if (env->nopush){
    push_back(env->Stack, state);
    env->nopush = false;
  }
  return;
}

/*---------------------------------------------------------------------------*/
PRIVATE void
repeat_gquad( vrna_fold_compound_t *vc,
              int i,
              int j,
              STATE *state,
              int part_energy,
              int temp_energy,
              int best_energy,
              int threshold,
              subopt_env *env){

  int           *ggg, *indx, element_energy, cp;
  short         *S1;
  vrna_param_t  *P;

  indx  = vc->jindx;
  cp    = vc->cutpoint;
  ggg   = vc->matrices->ggg;
  S1    = vc->sequence_encoding;
  P     = vc->params;


  /* find all gquads that fit into the energy range and the interval [i,j] */
  STATE *new_state;
  best_energy += part_energy; /* energy of current structural element */
  best_energy += temp_energy; /* energy from unpushed interval */

  if(ON_SAME_STRAND(i,j,cp)){
    element_energy = ggg[indx[j] + i];
    if(element_energy + best_energy <= threshold){
      int cnt;
      int *L;
      int *l;
      /* find out how many gquads we might expect in the interval [i,j] */
      int num_gquads = get_gquad_count(S1, i, j);
      num_gquads++;
      L = (int *)vrna_alloc(sizeof(int) * num_gquads);
      l = (int *)vrna_alloc(sizeof(int) * num_gquads * 3);
      L[0] = -1;

      get_gquad_pattern_exhaustive(S1, i, j, P, L, l, threshold - best_energy);

      for(cnt = 0; L[cnt] != -1; cnt++){
        new_state = copy_state(state);

        make_gquad(i, L[cnt], &(l[3*cnt]), new_state);
        new_state->partial_energy += part_energy;
        new_state->partial_energy += element_energy;
        /* new_state->best_energy =
           hairpin[unpaired] + element_energy + best_energy; */
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
repeat( vrna_fold_compound_t *vc,
        int i,
        int j,
        STATE * state,
        int part_energy,
        int temp_energy,
        int best_energy,
        int threshold,
        subopt_env *env){

  /* routine to find stacks, bulges, internal loops and  multiloops */
  /* within interval closed by basepair i,j */

  STATE           *new_state;
  vrna_param_t    *P;
  vrna_md_t       *md;

  register int  ij, k, p, q, energy, new;
  register int  mm;
  register int  no_close, type, type_2;
  char          *ptype;
  int           element_energy;
  int             *fc, *c, *fML, *fM1, *ggg;
  int           rt, *indx, *rtype, noGUclosure, noLP, with_gquad, dangle_model, turn, cp;
  short         *S1;
  vrna_hc_t     *hc;
  vrna_sc_t     *sc;

  S1        = vc->sequence_encoding;
  ptype     = vc->ptype;
  indx      = vc->jindx;
  cp        = vc->cutpoint;
  P         = vc->params;
  md        = &(P->model_details);
  rtype     = &(md->rtype[0]);

  noGUclosure   = md->noGUclosure;
  noLP          = md->noLP;
  with_gquad    = md->gquad;
  dangle_model  = md->dangles;
  turn          = md->min_loop_size;

  fc  = vc->matrices->fc;
  c   = vc->matrices->c;
  fML = vc->matrices->fML;
  fM1 = vc->matrices->fM1;
  ggg = vc->matrices->ggg;

  hc    = vc->hc;
  sc    = vc->sc;

  ij = indx[j]+i;

  type = ptype[ij];
/*
  if (type==0) fprintf(stderr, "repeat: Warning: %d %d can't pair\n", i,j);
*/

  if(type == 0)
    type = 7;

  no_close = (((type == 3) || (type == 4)) && noGUclosure);

  if(hc->matrix[ij] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP){
    if (noLP) /* always consider the structure with additional stack */
      if(i + turn + 2 < j){
        if(hc->matrix[indx[j-1]+i+1] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC){
          type_2 = rtype[ptype[indx[j-1]+i+1]];

          if(type_2 == 0)
            type_2 = 7;

          energy = 0;

          if(ON_SAME_STRAND(i,i+1,cp) && ON_SAME_STRAND(j-1,j, cp)){
            energy = E_IntLoop(0, 0, type, type_2,S1[i+1],S1[j-1],S1[i+1],S1[j-1], P);

            if(sc){
              if(sc->energy_bp)
                energy += sc->energy_bp[ij];

              if(sc->energy_stack)
                energy += sc->energy_stack[i]
                          + sc->energy_stack[i+1]
                          + sc->energy_stack[j-1]
                          + sc->energy_stack[j];

              if(sc->f)
                energy += sc->f(i, j, i+1, j-1, VRNA_DECOMP_PAIR_IL, sc->data);
            }

            new_state = derive_new_state(i+1, j-1, state, part_energy + energy, 2);
            make_pair(i, j, new_state);
            make_pair(i+1, j-1, new_state);

            /* new_state->best_energy = new + best_energy; */
            push(env->Stack, new_state);
            env->nopush = false;
            if (i==1 || state->structure[i-2]!='('  || state->structure[j]!=')')
              /* adding a stack is the only possible structure */
              return;
          }
        }
      }
  }

  best_energy += part_energy; /* energy of current structural element */
  best_energy += temp_energy; /* energy from unpushed interval */

  if(hc->matrix[ij] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP){
    for (p = i + 1; p <= MIN2 (j-2-turn,  i+MAXLOOP+1); p++) {
      int minq = j-i+p-MAXLOOP-2;
      if (minq<p+1+turn) minq = p+1+turn;

      if(hc->up_int[i+1] < (p - i - 1))
        break;

      for (q = j - 1; q >= minq; q--) {
        if(hc->up_int[q+1] < (j - q - 1))
          break;

        /* skip stack if noLP, since we've already processed it above */
        if((noLP) && (p==i+1) && (q==j-1))
          continue;

        if(!(hc->matrix[indx[q]+p] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC))
          continue;

        type_2 = ptype[indx[q]+p];

        if(type_2 == 0)
          type_2 = 7;

        if (noGUclosure)
          if (no_close||(type_2==3)||(type_2==4))
            if ((p>i+1)||(q<j-1)) continue;  /* continue unless stack */

        if (ON_SAME_STRAND(i,p,cp) && ON_SAME_STRAND(q,j,cp)) {
          energy = E_IntLoop(p-i-1, j-q-1, type, rtype[type_2],
                              S1[i+1],S1[j-1],S1[p-1],S1[q+1], P);

          new = energy + c[indx[q]+p];

          if(sc){
            if(sc->energy_up)
              energy  +=  sc->energy_up[i+1][p-i-1]
                          + sc->energy_up[q+1][j-q-1];

            if(sc->energy_bp)
              energy  += sc->energy_bp[ij];

            if(sc->energy_stack)
              if((p == i+1) && (q == j-1))
                energy  +=  sc->energy_stack[i]
                            + sc->energy_stack[p]
                            + sc->energy_stack[q]
                            + sc->energy_stack[j];

            if(sc->f)
              energy  += sc->f(i, j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
          }

          new = energy + c[indx[q]+p];

          if (new + best_energy <= threshold) {
            /* stack, bulge, or interior loop */
            fork_int_state(i, j, p, q, state, part_energy + energy, env);
          }
        }/*end of if block */
      } /* end of q-loop */
    } /* end of p-loop */
  }

  if (!ON_SAME_STRAND(i,j,cp)) { /*look in fc*/
    if(hc->matrix[ij] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP){
      rt = rtype[type];

      if(rt == 0)
        rt = 7;

      element_energy=0;
      switch(dangle_model){
        case 0:   element_energy = E_ExtLoop(rt, -1, -1, P);
                  break;
        default:  element_energy = E_ExtLoop(rt, (ON_SAME_STRAND(j-1,j,cp)) ? S1[j-1] : -1, (ON_SAME_STRAND(i,i+1,cp)) ? S1[i+1] : -1, P);
                  break;
      }

      if (fc[i+1] + fc[j-1] +element_energy + best_energy  <= threshold)
        {
          fork_two_states_pair(i, j, cp, state, part_energy + element_energy, 4, 5, env);
        }
    }
  }

  mm = P->MLclosing;
  rt = rtype[type];

  if((hc->matrix[ij] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) && (i != cp-1) && (j != cp)){

    if(rt == 0)
      rt = 7;

    element_energy = mm;
    switch(dangle_model){
      case 0:   element_energy = E_MLstem(rt, -1, -1, P) + mm;
                break;
      default:  element_energy = E_MLstem(rt, S1[j-1], S1[i+1], P) + mm;
                break;
    }

    if(sc){
      if(sc->energy_bp)
        element_energy += sc->energy_bp[ij];
    }

    for (k = i + turn + 2; k <= j - turn - 2; k++)  {
      /* multiloop decomposition */
      if ((fML[indx[k-1] + i+1] + fM1[indx[j-1] + k] +
          element_energy + best_energy)  <= threshold)
        {
          fork_two_states_pair(i, j, k, state, part_energy + element_energy, 1, 3, env);
        }
    }
  }

  if (ON_SAME_STRAND(i,j,cp)) {
    if(hc->matrix[ij] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP){

      if(no_close)
        element_energy = FORBIDDEN;
      else
          element_energy = vrna_E_hp_loop(vc, i, j);

      if (element_energy + best_energy <= threshold) {
        /* hairpin structure */
        fork_state_pair(i, j, state, part_energy + element_energy, env);
      }
    }

    if(with_gquad){
      /* now we have to find all loops where (i,j) encloses a gquad in an interior loops style */
      int cnt, *p, *q, *en, tmp_en;
      p = q = en = NULL;
      en = E_GQuad_IntLoop_exhaustive(i, j, &p, &q, type, S1, ggg, threshold - best_energy, indx, P);
      for(cnt = 0; p[cnt] != -1; cnt++){
        if((hc->up_int[i+1] >= p[cnt] - i - 1) && (hc->up_int[q[cnt]+1] >= j - q[cnt] - 1)){
          tmp_en = en[cnt];

          if(sc){
            if(sc->energy_bp)
              tmp_en += sc->energy_bp[ij];

            if(sc->energy_up)
              tmp_en += sc->energy_up[i+1][p[cnt] - i - 1]
                        + sc->energy_up[q[cnt]+1][j - q[cnt] - 1];
          }

          new_state = derive_new_state(p[cnt], q[cnt], state, tmp_en + part_energy, 6);

          make_pair(i, j, new_state);

          /* new_state->best_energy = new + best_energy; */
          push(env->Stack, new_state);
          env->nopush = false;
        }
      }
      free(en);
      free(p);
      free(q);
    }
  }

  best_energy -= part_energy;
  best_energy -= temp_energy;
  return;
}

PRIVATE void
old_subopt_print( const char *structure,
                  float energy,
                  void *data){

  struct old_subopt_dat *d = (struct old_subopt_dat *)data;

  if(structure && d->fp){
    char *e_string = vrna_strdup_printf(" %6.2f", energy);
    print_structure(d->fp, structure, e_string);
    free(e_string);
  }
}


PRIVATE void
old_subopt_store( const char *structure,
                  float energy,
                  void *data){

  struct old_subopt_dat *d = (struct old_subopt_dat *)data;

  /* store solution */
  if(d->n_sol + 1 == d->max_sol){
    d->max_sol     *= 2;
    d->SolutionList = (SOLUTION *)vrna_realloc(d->SolutionList, d->max_sol*sizeof(SOLUTION));
  }

  if(structure){
    d->SolutionList[d->n_sol].energy      = energy;
    d->SolutionList[d->n_sol++].structure = strdup(structure);
  } else {
    d->SolutionList[d->n_sol].energy      = 0;
    d->SolutionList[d->n_sol++].structure = NULL;
  }
}

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/
#ifdef VRNA_BACKWARD_COMPAT

PUBLIC SOLUTION *
subopt( char *seq,
        char *structure,
        int delta,
        FILE *fp){

  return wrap_subopt(seq, structure, NULL, delta, fold_constrained, 0, fp);
}

PUBLIC SOLUTION *
subopt_circ(char *seq,
            char *structure,
            int delta,
            FILE *fp){

  return wrap_subopt(seq, structure, NULL, delta, fold_constrained, 1, fp);
}

PUBLIC SOLUTION *subopt_par(char *seq,
                            char *structure,
                            vrna_param_t *parameters,
                            int delta,
                            int is_constrained,
                            int is_circular,
                            FILE *fp){

  return wrap_subopt(seq, structure, parameters, delta, is_constrained, is_circular, fp);
}

PRIVATE SOLUTION *
wrap_subopt(char *string,
            char *structure,
            vrna_param_t *parameters,
            int delta,
            int is_constrained,
            int is_circular,
            FILE *fp){

  vrna_fold_compound_t  *vc;
  vrna_param_t        *P;
  char                *seq;

#ifdef _OPENMP
/* Explicitly turn off dynamic threads */
  omp_set_dynamic(0);
#endif

  /* we need the parameter structure for hard constraints */
  if(parameters){
    P = vrna_params_copy(parameters);
  } else {
    vrna_md_t md;
    set_model_details(&md);
    md.temperature = temperature;
    P = vrna_params(&md);
  }
  P->model_details.circ     = is_circular;
  P->model_details.uniq_ML  = uniq_ML = 1;

  /* what about cofold sequences here? Is it safe to call the below cut_point_insert() ? */
  /* dirty hack to reinsert the '&' according to the global variable 'cut_point' */
  seq = vrna_cut_point_insert(string, cut_point);

  vc = vrna_fold_compound(seq, &(P->model_details), ((is_circular == 0) ? VRNA_OPTION_HYBRID : VRNA_OPTION_DEFAULT));

  if(parameters){ /* replace params if necessary */
    free(vc->params);
    vc->params = P;
  } else {
    free(P);
  }

  /* handle hard constraints in pseudo dot-bracket format if passed via simple interface */
  if(is_constrained && structure){
    unsigned int constraint_options = 0;
    constraint_options |= VRNA_CONSTRAINT_DB
                          | VRNA_CONSTRAINT_DB_PIPE
                          | VRNA_CONSTRAINT_DB_DOT
                          | VRNA_CONSTRAINT_DB_X
                          | VRNA_CONSTRAINT_DB_ANG_BRACK
                          | VRNA_CONSTRAINT_DB_RND_BRACK
                          | VRNA_CONSTRAINT_DB_INTRAMOL
                          | VRNA_CONSTRAINT_DB_INTERMOL;

    vrna_constraints_add(vc, (const char *)structure, constraint_options);
  }

  if(backward_compat_compound && backward_compat)
    vrna_fold_compound_free(backward_compat_compound);

  backward_compat_compound  = vc;
  backward_compat           = 1;

  /* cleanup */
  free(seq);

  return vrna_subopt(vc, delta, subopt_sorted, fp);
}

#endif

/*---------------------------------------------------------------------------*/
/* Well, that is the end!----------------------------------------------------*/
/*---------------------------------------------------------------------------*/
