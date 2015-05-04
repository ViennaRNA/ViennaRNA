/*
  $Log: subopt.c,v $
  Revision 2.0  2010/12/06 20:04:20  ronny
  repaired subopt for cofolding

  Revision 1.24  2008/11/01 21:10:20  ivo
  avoid rounding errors when computing DoS

  Revision 1.23  2008/03/31 15:06:49  ivo
  Add cofolding support in subopt

  Revision 1.22  2008/02/23 09:42:35  ivo
  fix circular folding bugs with dangles that cross the origin

  Revision 1.21  2008/01/08 15:08:51  ivo
  circular fold would fail for open chain

  Revision 1.20  2008/01/08 14:08:20  ivo
  add an option to compute the density of state

  Revision 1.19  2007/12/05 13:04:04  ivo
  add various circfold variants from Ronny

  Revision 1.18  2003/10/06 08:56:45  ivo
  use P->TerminalAU

  Revision 1.17  2003/08/26 09:26:08  ivo
  don't modify print_energy in subopt(); use doubles instead of floats

  Revision 1.16  2001/10/01 13:50:00  ivo
  sorted -> subopt_sorted

  Revision 1.15  2001/09/17 10:30:42  ivo
  move scale_parameters() into params.c
  returns pointer to paramT structure

  Revision 1.14  2001/08/31 15:02:19  ivo
  Let subopt either write to file pointer or return a list of structures,
  so we can nicely integrate it into the library

  Revision 1.13  2001/04/05 07:35:08  ivo
  remove uneeded declaration of TETRA_ENERGY

  Revision 1.12  2000/10/10 08:53:20  ivo
  adapted for new Turner energy parameters
  supports all constraints that forbid pairs

  Revision 1.11  2000/04/08 15:56:18  ivo
  with noLonelyPairs=1 will produce no structures with isolated base pairs
  (Giegerich's canonical structures)

  Revision 1.10  1999/05/06 10:13:35  ivo
  recalculte energies before printing if logML is set
  + cosmetic changes

  Revision 1.9  1998/05/19 16:31:52  ivo
  added support for constrained folding

  Revision 1.8  1998/03/30 14:44:54  ivo
  cleanup of make_printout etc.

  Revision 1.7  1998/03/30 14:39:31  ivo
  replaced BasePairs list with structure string in STATE
  save memory by not storing (and sorting) structures
  modified for use with ViennaRNA-1.2.1

  Revision 1.6  1997/10/21 11:34:09  walter
  steve update

  Revision 1.1  1997/08/04 21:05:32  walter
  Initial revision

*/
/*
   suboptimal folding - Stefan Wuchty, Walter Fontana & Ivo Hofacker

                       Vienna RNA package
*/
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "fold.h"
#include "utils.h"
#include "energy_par.h"
#include "fold_vars.h"
#include "pair_mat.h"
#include "list.h"
#include "params.h"
#include "loop_energies.h"
#include "cofold.h"
#include "gquad.h"
#include "subopt.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define true              1
#define false             0
#define SAME_STRAND(I,J)  (((I)>=cut_point)||((J)<cut_point))
#define NEW_NINIO         1         /* use new asymetry penalty */
#define STACK_BULGE1      1         /* stacking energies for bulges of size 1 */
#define MAXALPHA          20        /* maximal length of alphabet */

/*
#################################
# GLOBAL VARIABLES              #
#################################
*/
PUBLIC  int     subopt_sorted=0;                           /* output sorted by energy */
PUBLIC  int     density_of_states[MAXDOS+1];
PUBLIC  double  print_energy = 9999; /* printing threshold for use with logML */

typedef struct {
    char *structure;
    LIST *Intervals;
    int partial_energy;
    int is_duplex;
    /* int best_energy;   */ /* best attainable energy */
} STATE;

/*
#################################
# PRIVATE VARIABLES             #
#################################
*/
PRIVATE int     turn;
PRIVATE LIST    *Stack = NULL;
PRIVATE int     nopush;
PRIVATE int     best_energy;          /* best_energy = remaining energy */
PRIVATE int     *f5 = NULL;           /* energy of 5 end */
PRIVATE int     *c = NULL;            /* energy array, given that i-j pair */
PRIVATE int     *fML = NULL;          /* multi-loop auxiliary energy array */
PRIVATE int     *fM1 = NULL;          /* another multi-loop auxiliary energy array */
PRIVATE int     *fc = NULL;           /*energy array, from i (j)  to cut*/
PRIVATE int     *indx = NULL;         /* index for moving in the triangle matrices c[] and f[] */
PRIVATE short   *S=NULL, *S1=NULL;
PRIVATE char    *ptype=NULL;
PRIVATE paramT  *P = NULL;
PRIVATE int     length;
PRIVATE int     minimal_energy;       /* minimum free energy */
PRIVATE int     element_energy;       /* internal energy of a structural element */
PRIVATE int     threshold;            /* minimal_energy + delta */
PRIVATE char    *sequence = NULL;
PRIVATE int     circular            = 0;
PRIVATE int     struct_constrained  = 0;
PRIVATE int     *fM2 = NULL;                 /* energies of M2 */
PRIVATE int     Fc, FcH, FcI, FcM;    /* parts of the exterior loop energies */
PRIVATE int     with_gquad          = 0;


PRIVATE int     *ggg = NULL;

#ifdef _OPENMP

#pragma omp threadprivate(turn, Stack, nopush, best_energy, f5, c, fML, fM1, fc, indx, S, S1,\
                          ptype, P, length, minimal_energy, element_energy, threshold, sequence,\
                          fM2, Fc, FcH, FcI, FcM, circular, struct_constrained,\
                          ggg, with_gquad)

#endif

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/
PRIVATE void      make_pair(int i, int j, STATE *state);

/* mark a gquadruplex in the resulting dot-bracket structure */
PRIVATE void      make_gquad(int i, int L, int l[3], STATE *state);

PRIVATE INTERVAL  *make_interval (int i, int j, int ml);
/*@out@*/ PRIVATE STATE *make_state(/*@only@*/LIST *Intervals,
                                    /*@only@*/ /*@null@*/ char *structure,
                                    int partial_energy, int is_duplex);
PRIVATE STATE     *copy_state(STATE * state);
PRIVATE void      print_state(STATE * state);
PRIVATE void      UNUSED print_stack(LIST * list);
/*@only@*/ PRIVATE LIST *make_list(void);
PRIVATE void      push(LIST * list, /*@only@*/ void *data);
PRIVATE void      *pop(LIST * list);
PRIVATE int       best_attainable_energy(STATE * state);
PRIVATE void      scan_interval(int i, int j, int array_flag, STATE * state);
PRIVATE void      free_interval_node(/*@only@*/ INTERVAL * node);
PRIVATE void      free_state_node(/*@only@*/ STATE * node);
PRIVATE void      push_back(STATE * state);
PRIVATE char*     get_structure(STATE * state);
PRIVATE int       compare(const void *solution1, const void *solution2);
PRIVATE void      make_output(SOLUTION *SL, FILE *fp);
PRIVATE char      *costring(char *string);
PRIVATE void      repeat(int i, int j, STATE * state, int part_energy, int temp_energy);
PRIVATE void      repeat_gquad( int i, int j, STATE *state, int part_energy, int temp_energy);

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
           int is_duplex)
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
    state->structure = (char *) space(length+1);
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
  if (!new_state->structure) nrerror("out of memory");
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
  nopush = false;
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
best_attainable_energy(STATE * state)
{
  /* evaluation of best possible energy attainable within remaining intervals */

  register int sum;
  INTERVAL *next;

  sum = state->partial_energy;  /* energy of already found elements */

  for (next = lst_first(state->Intervals); next; next = lst_next(next))
    {
      if (next->array_flag == 0)
        sum += (circular) ? Fc : f5[next->j];
      else if (next->array_flag == 1)
        sum += fML[indx[next->j] + next->i];
      else if (next->array_flag == 2)
        sum += c[indx[next->j] + next->i];
      else if (next->array_flag == 3)
        sum += fM1[indx[next->j] + next->i];
      else if (next->array_flag == 4)
        sum += fc[next->i];
      else if (next->array_flag == 5)
        sum += fc[next->j];
      else if (next->array_flag == 6)
        sum += ggg[indx[next->j] + next->i];
    }

  return sum;
}

/*---------------------------------------------------------------------------*/

PRIVATE void
push_back(STATE * state)
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

PRIVATE void make_output(SOLUTION *SL, FILE *fp)  /* prints stuff */
{
  SOLUTION *sol;

  for (sol = SL; sol->structure!=NULL; sol++)
    if (cut_point<0) fprintf(fp, "%s %6.2f\n", sol->structure, sol->energy);
    else {
      char *tStruc;
      tStruc=costring(sol->structure);
      fprintf(fp, "%s %6.2f\n", tStruc, sol->energy);
      free(tStruc);
    }
}

/*---------------------------------------------------------------------------*/
/* start of subopt backtracking ---------------------------------------------*/
/*---------------------------------------------------------------------------*/

PUBLIC SOLUTION *subopt(char *seq, char *structure, int delta, FILE *fp){
  return subopt_par(seq, structure, NULL, delta, fold_constrained, 0, fp);
}

PUBLIC SOLUTION *subopt_circ(char *seq, char *structure, int delta, FILE *fp){
  return subopt_par(seq, structure, NULL, delta, fold_constrained, 1, fp);
}

PUBLIC SOLUTION *subopt_par(char *seq,
                            char *structure,
                            paramT *parameters,
                            int delta,
                            int is_constrained,
                            int is_circular,
                            FILE *fp){

  STATE         *state;
  LIST          *Intervals;
  INTERVAL      *interval;
  SOLUTION      *SolutionList;
  unsigned long max_sol, n_sol;
  int           maxlevel, count, partial_energy, old_dangles, logML, dangle_model;
  double        structure_energy, min_en, eprint;
  char          *struc;

  max_sol             = 128;
  n_sol               = 0;
  sequence            = seq;
  length              = strlen(sequence);
  circular            = is_circular;
  struct_constrained  = is_constrained;

  struc = (char *) space(sizeof(char)*(length+1));
  if (struct_constrained) strncpy(struc, structure, length);

  /* do mfe folding to get fill arrays and get ground state energy  */
  /* in case dangles is neither 0 or 2, set dangles=2 while folding */

  if(P) free(P);
  if(parameters){
    P = get_parameter_copy(parameters);
  } else {
    model_detailsT md;
    set_model_details(&md);
    P = get_scaled_parameters(temperature, md);
  }

  logML       = P->model_details.logML;
  old_dangles = dangle_model = P->model_details.dangles;
  with_gquad  = P->model_details.gquad;

  /* temporarily set dangles to 2 if necessary */
  if((P->model_details.dangles != 0) && (P->model_details.dangles != 2))
    P->model_details.dangles = 2;


  turn = (cut_point<0) ? 3 : 0;
  uniq_ML = 1;
  if(circular){
    min_en = fold_par(sequence, struc, P, struct_constrained, circular);
    export_circfold_arrays(&Fc, &FcH, &FcI, &FcM, &fM2, &f5, &c, &fML, &fM1, &indx, &ptype);
    /* restore dangle model */
    P->model_details.dangles = old_dangles;
    /* re-evaluate in case we're using logML etc */
    min_en = energy_of_circ_struct_par(sequence, struc, P, 0);
  } else {
    min_en = cofold_par(sequence, struc, P, struct_constrained);

    if(with_gquad){
      export_cofold_arrays_gq(&f5, &c, &fML, &fM1, &fc, &ggg, &indx, &ptype);
    } else {
      export_cofold_arrays(&f5, &c, &fML, &fM1, &fc, &indx, &ptype);
    }

    /* restore dangle model */
    P->model_details.dangles = old_dangles;
    /* re-evaluate in case we're using logML etc */
    min_en = energy_of_struct_par(sequence, struc, P, 0);
  }

  free(struc);
  eprint = print_energy + min_en;
  if (fp) {
    char *SeQ;
    SeQ=costring(sequence);
    fprintf(fp, "%s %6d %6d\n", SeQ, (int) (-0.1+100*min_en), delta);
    free(SeQ);
  }
  make_pair_matrix();
  S   = encode_sequence(sequence, 0);
  S1  = encode_sequence(sequence, 1);

  /* Initialize ------------------------------------------------------------ */

  maxlevel = 0;
  count = 0;
  partial_energy = 0;

  /* Initialize the stack ------------------------------------------------- */

  minimal_energy = (circular) ? Fc : f5[length];
  threshold = minimal_energy + delta;
  if(threshold > INF){
    warn_user("energy range too high, limiting to reasonable value");
    threshold = INF-EMAX;
  }
  Stack = make_list();                                                   /* anchor */
  Intervals = make_list();                                   /* initial state: */
  interval = make_interval(1, length, 0);          /* interval [1,length,0] */
  push(Intervals, interval);
  state = make_state(Intervals, NULL, partial_energy,0);
  /* state->best_energy = minimal_energy; */
  push(Stack, state);

  /* SolutionList stores the suboptimal structures found */

  SolutionList = (SOLUTION *) space(max_sol*sizeof(SOLUTION));

  /* end initialize ------------------------------------------------------- */


  while (1) {                    /* forever, til nothing remains on stack */

    maxlevel = (Stack->count > maxlevel ? Stack->count : maxlevel);

    if (LST_EMPTY (Stack))                   /* we are done! clean up and quit */
      {
        /* fprintf(stderr, "maxlevel: %d\n", maxlevel); */

        lst_kill(Stack, free_state_node);

        SolutionList[n_sol].structure = NULL; /* NULL terminate list */

        if (subopt_sorted) {
          /* sort structures by energy */
          qsort(SolutionList, n_sol, sizeof(SOLUTION), compare);

          if (fp) make_output(SolutionList, fp);
        }

        break;
      }

    /* pop the last element ---------------------------------------------- */

    state = pop(Stack);                       /* current state to work with */

    if (LST_EMPTY(state->Intervals))
      {
        int e;
        /* state has no intervals left: we got a solution */

        count++;
        structure = get_structure(state);
        structure_energy = state->partial_energy / 100.;

#ifdef CHECK_ENERGY
        structure_energy = (circular) ? energy_of_circ_struct_par(sequence, structure, P, 0) : (with_gquad) ? energy_of_gquad_struct_par(sequence, structure, P, 0) : energy_of_struct_par(sequence, structure, P, 0);

        if (!logML)
          if ((double) (state->partial_energy / 100.) != structure_energy) {
            fprintf(stderr, "%s %6.2f %6.2f\n", structure,
                    state->partial_energy / 100., structure_energy );
            exit(1);
          }
#endif
        if (logML || (dangle_model==1) || (dangle_model==3)) { /* recalc energy */
          structure_energy = (circular) ? energy_of_circ_struct_par(sequence, structure, P, 0) : (with_gquad) ? energy_of_gquad_struct_par(sequence, structure, P, 0) : energy_of_struct_par(sequence, structure, P, 0);
        }

        e = (int) ((structure_energy-min_en)*10. + 0.1); /* avoid rounding errors */
        if (e>MAXDOS) e=MAXDOS;
        density_of_states[e]++;
        if (structure_energy>eprint) {
          free(structure);
        } else {
          if (!subopt_sorted && fp) {
            /* print and forget */
            if (cut_point<0)
              fprintf(fp, "%s %6.2f\n", structure, structure_energy);
            else {
              char * outstruc;
              /*make ampersand seperated output if 2 sequences*/
              outstruc=costring(structure);
              fprintf(fp, "%s %6.2f\n", outstruc, structure_energy);
              free(outstruc);
            }
            free(structure);
          }
          else {
            /* store solution */
            if (n_sol+1 == max_sol) {
              max_sol *= 2;
              SolutionList = (SOLUTION *)
                xrealloc(SolutionList, max_sol*sizeof(SOLUTION));
            }
            SolutionList[n_sol].energy =  structure_energy;
            SolutionList[n_sol++].structure = structure;
          }
        }
      }
    else {
      /* get (and remove) next interval of state to analyze */

      interval = pop(state->Intervals);

      scan_interval(interval->i, interval->j, interval->array_flag, state);

      free_interval_node(interval);        /* free the current interval */
    }

    free_state_node(state);                     /* free the current state */
  } /* end of while (1) */

  /* free arrays left over from cofold() */
  free(S); free(S1);
  (circular) ? free_arrays():free_co_arrays();
  if (fp) { /* we've printed everything -- free solutions */
    SOLUTION *sol;
    for (sol=SolutionList; sol->structure != NULL; sol++)
      free(sol->structure);
    free(SolutionList);
    SolutionList = NULL;
  }

  return SolutionList;
}


PRIVATE void
scan_interval(int i, int j, int array_flag, STATE * state)
{
  /* real backtrack routine */

  /* array_flag = 0:  trace back in f5-array  */
  /* array_flag = 1:  trace back in fML-array */
  /* array_flag = 2:  trace back in repeat()  */
  /* array_flag = 3:  trace back in fM1-array */

  STATE *new_state, *temp_state;
  INTERVAL *new_interval;
  register int k, fi, cij;
  register int type;
  register int dangle_model = P->model_details.dangles;
  register int noGUclosure  = P->model_details.noGUclosure;
  register int noLP         = P->model_details.noLP;

  best_energy = best_attainable_energy(state);  /* .. on remaining intervals */
  nopush = true;

  if ((i > 1) && (!array_flag))
    nrerror ("Error while backtracking!");

  if (j < i + turn + 1 && SAME_STRAND(i,j)) { /* minimal structure element */
    if (nopush)
      push_back(state);
    return;
  }

  /* 13131313131313131313131313131313131313131313131313131313131313131313131 */
  if (array_flag == 3 || array_flag == 1) {
    /* array_flag = 3: interval i,j was generated during */
    /*                 a multiloop decomposition using array fM1 in repeat() */
    /*                 or in this block */

    /* array_flag = 1: interval i,j was generated from a */
    /*                 stack, bulge, or internal loop in repeat() */
    /*                 or in this block */

    if (array_flag == 3)
      fi = fM1[indx[j-1] + i] + P->MLbase;
    else
      fi = fML[indx[j-1] + i] + P->MLbase;

    if ((fi + best_energy <= threshold)&&(SAME_STRAND(j-1,j))) {
      /* no basepair, nibbling of 3'-end */

      new_state = copy_state(state);
      new_interval = make_interval(i, j-1, array_flag);
      push(new_state->Intervals, new_interval);
      new_state->partial_energy += P->MLbase;
      /* new_state->best_energy = fi + best_energy; */
      push(Stack, new_state);
    }

    type = ptype[indx[j]+i];

    if (type) { /* i,j may pair */

      if(dangle_model)
        element_energy = E_MLstem(type,
                                  (((i > 1)&&(SAME_STRAND(i-1,i))) || circular)       ? S1[i-1] : -1,
                                  (((j < length)&&(SAME_STRAND(j,j+1))) || circular)  ? S1[j+1] : -1,
                                  P);
      else
        element_energy = E_MLstem(type, -1, -1, P);

      cij = c[indx[j] + i] + element_energy;
      if (cij + best_energy <= threshold)
        repeat(i, j, state, element_energy, 0);
    } else if (with_gquad){
      element_energy = E_MLstem(0, -1, -1, P);
      cij = ggg[indx[j] + i] + element_energy;
      if(cij + best_energy <= threshold)
        repeat_gquad(i, j, state, element_energy, 0);
    }
  }                                   /* array_flag == 3 || array_flag == 1 */

  /* 11111111111111111111111111111111111111111111111111111111111111111111111 */

  if (array_flag == 1) {
    /* array_flag = 1:                   interval i,j was generated from a */
    /*                          stack, bulge, or internal loop in repeat() */
    /*                          or in this block */

    int stopp;
    if ((SAME_STRAND(i-1,i))&&(SAME_STRAND(j,j+1))) { /*backtrack in FML only if multiloop is possible*/
      for ( k = i+turn+1 ; k <= j-1-turn ; k++) {
        /* Multiloop decomposition if i,j contains more than 1 stack */
        
        if(with_gquad){
          if(SAME_STRAND(k, k+1)){
            element_energy = E_MLstem(0, -1, -1, P);
            if(fML[indx[k]+i] + ggg[indx[j] + k + 1] + element_energy + best_energy <= threshold){
              temp_state = copy_state (state);
              new_interval = make_interval (i, k, 1);
              push (temp_state->Intervals, new_interval);
              repeat_gquad(k+1, j, temp_state, element_energy, fML[indx[k]+i]);
              free_state_node(temp_state);
            }
          }
        }

        type = ptype[indx[j]+k+1];
        if (type==0) continue;

        if(dangle_model)
          element_energy = E_MLstem(type,
                                    (SAME_STRAND(i-1,i)) ? S1[k] : -1,
                                    (SAME_STRAND(j,j+1)) ? S1[j+1] : -1,
                                    P);
        else
          element_energy = E_MLstem(type, -1, -1, P);

        if (SAME_STRAND(k,k+1)) {
          if (fML[indx[k]+i] + c[indx[j] + k+1] +
              element_energy + best_energy <= threshold) {
            temp_state = copy_state (state);
            new_interval = make_interval (i, k, 1);
            push (temp_state->Intervals, new_interval);
            repeat(k+1, j, temp_state, element_energy, fML[indx[k]+i]);
            free_state_node(temp_state);
          }
        }
      }
    }

    stopp=(cut_point>0)? (cut_point-2):(length); /*if cut_point -1: k on cut, => no ml*/
    stopp=MIN2(stopp, j-1-turn);
    if (i>cut_point) stopp=j-1-turn;
    else if (i==cut_point) stopp=0;   /*not a multi loop*/
    for (k = i ; k <= stopp; k++) {
      /* Multiloop decomposition if i,j contains only 1 stack */
      if(with_gquad){
        element_energy = E_MLstem(0, -1, -1, P) + P->MLbase*(k-i+1);
        if(ggg[indx[j] + k + 1] + element_energy + best_energy <= threshold)
          repeat_gquad(k+1, j, state, element_energy, 0);
      }

      type = ptype[indx[j]+k+1];
      if (type==0) continue;

      if(dangle_model)
        element_energy = E_MLstem(type,
                                  (SAME_STRAND(k-1,k)) ? S1[k] : -1,
                                  (SAME_STRAND(j,j+1)) ? S1[j+1] : -1,
                                  P);
      else
        element_energy = E_MLstem(type, -1, -1, P);

      element_energy += P->MLbase*(k-i+1);

      if (c[indx[j]+k+1] + element_energy + best_energy <= threshold)
        repeat(k+1, j, state, element_energy, 0);
    }
  }                                                    /* array_flag == 1 */

  /* 2222222222222222222222222222222222222222222222222222222222222222222222 */

  if (array_flag == 2)
    {
      /* array_flag = 2:                  interval i,j was generated from a */
      /* stack, bulge, or internal loop in repeat() */

      repeat(i, j, state, 0, 0);

      if (nopush){
        if (!noLP){
          fprintf(stderr, "%d,%d", i, j);
          fprintf(stderr, "Oops, no solution in repeat!\n");
        }
      }
      return;
    }

  /* 0000000000000000000000000000000000000000000000000000000000000000000000 */

  if ((array_flag == 0) && !circular)
    {
      /* array_flag = 0:                        interval i,j was found while */
      /* tracing back through f5-array and c-array */
      /* or within this block */

      if (f5[j-1] + best_energy <= threshold) {
        /* no basepair, nibbling of 3'-end */

        new_state = copy_state(state);
        new_interval = make_interval(i, j-1 , 0);
        push(new_state->Intervals, new_interval);
        /* new_state->best_energy = f5[j-1] + best_energy; */
        push(Stack, new_state);
      }

      for (k = j-turn-1; k > 1; k--) {

        if(with_gquad){
          if(SAME_STRAND(k,j)){
            element_energy = 0;
            if(f5[k-1] + ggg[indx[j]+k] + element_energy + best_energy <= threshold){
              temp_state = copy_state(state);
              new_interval = make_interval(1,k-1,0);
              push(temp_state->Intervals, new_interval);
              /* backtrace the quadruplex */
              repeat_gquad(k, j, temp_state, element_energy, f5[k-1]);
              free_state_node(temp_state);
            }
          }
        }

        type = ptype[indx[j]+k];
        if (type==0)   continue;

        /* k and j pair */
        if(dangle_model)
          element_energy = E_ExtLoop(type,
                                    (SAME_STRAND(k-1,k)) ? S1[k-1] : -1,
                                    ((j < length)&&(SAME_STRAND(j,j+1))) ? S1[j+1] : -1,
                                    P);
        else
          element_energy = E_ExtLoop(type, -1, -1, P);

        if (!(SAME_STRAND(k,j)))/*&&(state->is_duplex==0))*/ {
          element_energy+=P->DuplexInit;
          /*state->is_duplex=1;*/
        }

        if (f5[k-1] + c[indx[j]+k] + element_energy + best_energy <= threshold)
          {
            temp_state = copy_state(state);
            new_interval = make_interval(1, k-1, 0);
            push(temp_state->Intervals, new_interval);
            repeat(k, j, temp_state, element_energy, f5[k-1]);
            free_state_node(temp_state);
          }
      }
      type = ptype[indx[j]+1];
      if (type) {
        if (dangle_model && (j < length)&&(SAME_STRAND(j,j+1)))
          element_energy = E_ExtLoop(type, -1, S1[j+1], P);
        else
          element_energy = E_ExtLoop(type, -1, -1, P);

        if (!(SAME_STRAND(1,j))) element_energy+=P->DuplexInit;

        if (c[indx[j]+1] + element_energy + best_energy <= threshold)
          repeat(1, j, state, element_energy, 0);
      } else if (with_gquad){
        if(SAME_STRAND(k,j)){
          element_energy = 0;
          if(ggg[indx[j]+1] + element_energy + best_energy <= threshold){
            /* backtrace the quadruplex */
            repeat_gquad(1, j, state, element_energy, 0);
          }
        }
      }
    }/* end array_flag == 0 && !circular*/
  /* or do we subopt circular? */
  else if(array_flag == 0){
    int k, l, p, q;
    /* if we've done everything right, we will never reach this case more than once   */
    /* right after the initilization of the stack with ([1,n], empty, 0)              */
    /* lets check, if we can have an open chain without breaking the threshold        */
    /* this is an ugly work-arround cause in case of an open chain we do not have to  */
    /* backtrack anything further...                                                  */
    if(0 <= threshold){
      new_state = copy_state(state);
      new_interval = make_interval(1,2,0);
      push(new_state->Intervals, new_interval);
      new_state->partial_energy = 0;
      push(Stack, new_state);
    }
    /* ok, lets check if we can do an exterior hairpin without breaking the threshold */
    /* best energy should be 0 if we are here                                         */
    if(FcH + best_energy <= threshold){
      /* lets search for all exterior hairpin cases, that fit into our threshold barrier  */
      /* we use index k,l to avoid confusion with i,j index of our state...               */
      /* if we reach here, i should be 1 and j should be n respectively                   */
      for(k=i; k<j; k++)
        for (l=k+turn+1; l <= j; l++){
          int kl, type, u, tmpE, no_close;
          u = j-l + k-1;        /* get the hairpin loop length */
          if(u<turn) continue;

          kl = indx[l]+k;        /* just confusing these indices ;-) */
          type = ptype[kl];
          no_close = ((type==3)||(type==4))&&noGUclosure;
          type=rtype[type];
          if (!type) continue;
          if (!no_close){
            /* now lets have a look at the hairpin energy */
            char loopseq[10];
            if (u<7){
              strcpy(loopseq , sequence+l-1);
              strncat(loopseq, sequence, k);
            }
            tmpE = E_Hairpin(u, type, S1[l+1], S1[k-1], loopseq, P);
          }
          if(c[kl] + tmpE + best_energy <= threshold){
            /* what we really have to do is something like this, isn't it?  */
            /* we have to create a new state, with interval [k,l], then we  */
            /* add our loop energy as initial energy of this state and put  */
            /* the state onto the stack R... for further refinement...      */
            /* we also denote this new interval to be scanned in C          */
            new_state = copy_state(state);
            new_interval = make_interval(k,l,2);
            push(new_state->Intervals, new_interval);
            /* hopefully we add this energy in the right way... */
            new_state->partial_energy += tmpE;
            push(Stack, new_state);
          }
        }
    }

    /* now lets see, if we can do an exterior interior loop without breaking the threshold */
    if(FcI + best_energy <= threshold){
      /* now we search for our exterior interior loop possibilities */
      for(k=i; k<j; k++)
        for (l=k+turn+1; l <= j; l++){
          int kl, type, tmpE;

          kl = indx[l]+k;        /* just confusing these indices ;-) */
          type = ptype[kl];
          type=rtype[type];
          if (!type) continue;

          for (p = l+1; p < j ; p++){
            int u1, qmin;
            u1 = p-l-1;
            if (u1+k-1>MAXLOOP) break;
            qmin = u1+k-1+j-MAXLOOP;
            if(qmin<p+turn+1) qmin = p+turn+1;
            for(q = qmin; q <=j; q++){
              int u2, type_2;
              type_2 = rtype[ptype[indx[q]+p]];
              if(!type_2) continue;
              u2 = k-1 + j-q;
              if(u1+u2>MAXLOOP) continue;
              tmpE = E_IntLoop(u1, u2, type, type_2, S1[l+1], S1[k-1], S1[p-1], S1[q+1], P);
              if(c[kl] + c[indx[q]+p] + tmpE + best_energy <= threshold){
                /* ok, similar to the hairpin stuff, we add new states onto the stack R */
                /* but in contrast to the hairpin decomposition, we have to add two new */
                /* intervals, enclosed by k,l and p,q respectively and we also have to  */
                /* add the partial energy, that comes from the exterior interior loop   */
                new_state = copy_state(state);
                new_interval = make_interval(k, l, 2);
                push(new_state->Intervals, new_interval);
                new_interval = make_interval(p,q,2);
                push(new_state->Intervals, new_interval);
                new_state->partial_energy += tmpE;
                push(Stack, new_state);
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

              /* next, we have the first interval that has to be traced in fM1 */
              new_interval = make_interval(k+1, l, 3);
              push(new_state->Intervals, new_interval);

              /* and the last of our three intervals is also one to be traced within fM1 array... */
              new_interval = make_interval(l+1, j, 3);
              push(new_state->Intervals, new_interval);

              /* mmh, we add the energy for closing the multiloop now... */
              new_state->partial_energy += P->MLclosing;
              /* next we push our state onto the R stack */
              push(Stack, new_state);

            }
            /* else we search further... */
          }

          /* ok, we have to decompose fML now... */
        }
      }
    }
  }        /* thats all folks for the circular case... */

  /* 44444444444444444444444444444444444444444444444444444444444444 */
  if (array_flag == 4) {
    /* array_flag = 4:                        interval i,j was found while */
    /* tracing back through fc-array smaller than than cut_point*/
    /* or within this block */

    if (fc[i+1] + best_energy <= threshold) {
      /* no basepair, nibbling of 5'-end */
      new_state = copy_state(state);
      new_interval = make_interval(i+1, j , 4);
      push(new_state->Intervals, new_interval);
      push(Stack, new_state);
    }

    for (k = i+TURN+1; k < j; k++) {

      if(with_gquad){
        if(fc[k+1] + ggg[indx[k]+i] + best_energy <= threshold){
          temp_state = copy_state(state);
          new_interval = make_interval(k+1,j, 4);
          push(temp_state->Intervals, new_interval);
          repeat_gquad(i, k, temp_state, 0, fc[k+1]);
          free_state_node(temp_state);
        }
      }

      type = ptype[indx[k]+i];
      if (type==0)   continue;

      /* k and j pair */
      if (dangle_model)
        element_energy = E_ExtLoop(type, (i > 1) ? S1[i-1]: -1, S1[k+1], P);
      else  /* no dangles */
        element_energy = E_ExtLoop(type, -1, -1, P);

      if (fc[k+1] + c[indx[k]+i] + element_energy + best_energy <= threshold) {
        temp_state = copy_state(state);
        new_interval = make_interval(k+1,j, 4);
        push(temp_state->Intervals, new_interval);
        repeat(i, k, temp_state, element_energy, fc[k+1]);
        free_state_node(temp_state);
      }
    }
    type = ptype[indx[j]+i];
    if (type) {
      if (dangle_model)
        element_energy = E_ExtLoop(type, (i>1) ? S1[i-1] : -1, -1, P);
      else
        element_energy = E_ExtLoop(type, -1, -1, P);

      if (c[indx[cut_point-1]+i] + element_energy + best_energy <= threshold)
        repeat(i, cut_point-1, state, element_energy, 0);
    } else if(with_gquad){
      if(ggg[indx[cut_point -1] + i] + best_energy <= threshold)
        repeat_gquad(i, cut_point - 1, state, 0, 0);
    }
  } /* array_flag == 4 */

  /*55555555555555555555555555555555555555555555555555555555555555555555555*/
  if (array_flag == 5) {
    /* array_flag = 5:  interval cut_point=i,j was found while  */
    /* tracing back through fc-array greater than cut_point     */
    /* or within this block                                     */

    if (fc[j-1] + best_energy <= threshold) {
      /* no basepair, nibbling of 3'-end */
      new_state = copy_state(state);
      new_interval = make_interval(i, j-1 , 5);
      push(new_state->Intervals, new_interval);
      push(Stack, new_state);
    }

    for (k = j-TURN-1; k > i; k--) {

      if(with_gquad){
        if(fc[k-1] + ggg[indx[j] + k] + best_energy <= threshold){
          temp_state = copy_state(state);
          new_interval = make_interval(i, k-1, 5);
          push(temp_state->Intervals, new_interval);
          repeat_gquad(k, j, temp_state, 0, fc[k-1]);
          free_state_node(temp_state);
        }
      }

      type = ptype[indx[j]+k];
      if (type==0)   continue;
      element_energy = 0;

      /* k and j pair */
      if (dangle_model)
        element_energy = E_ExtLoop(type, S1[k-1], (j < length) ? S1[j+1] : -1, P);
      else
        element_energy = E_ExtLoop(type, -1, -1, P);

      if (fc[k-1] + c[indx[j]+k] + element_energy + best_energy <= threshold) {
        temp_state = copy_state(state);
        new_interval = make_interval(i, k-1, 5);
        push(temp_state->Intervals, new_interval);
        repeat(k, j, temp_state, element_energy, fc[k-1]);
        free_state_node(temp_state);
      }
    }
    type = ptype[indx[j]+i];
    if (type) {
      if(dangle_model)
        element_energy = E_ExtLoop(type, -1, (j<length) ? S1[j+1] : -1, P);

      if (c[indx[j]+cut_point] + element_energy + best_energy <= threshold)
        repeat(cut_point, j, state, element_energy, 0);
    } else if (with_gquad){
      if(ggg[indx[j] + cut_point] + best_energy <= threshold)
        repeat_gquad(cut_point, j, state, 0, 0);
    }
  } /* array_flag == 5 */

  if (array_flag == 6) { /* we have a gquad */
      repeat_gquad(i, j, state, 0, 0);
      if (nopush){
        fprintf(stderr, "%d,%d", i, j);
        fprintf(stderr, "Oops, no solution in gquad-repeat!\n");
      }
      return;
  
  
  }

  if (nopush)
    push_back(state);
  return;
}

/*---------------------------------------------------------------------------*/
PRIVATE void
repeat_gquad( int i,
              int j,
              STATE *state,
              int part_energy,
              int temp_energy){

  /* find all gquads that fit into the energy range and the interval [i,j] */
  STATE *new_state;
  best_energy += part_energy; /* energy of current structural element */
  best_energy += temp_energy; /* energy from unpushed interval */

  if(SAME_STRAND(i,j)){
    element_energy = ggg[indx[j] + i];
    if(element_energy + best_energy <= threshold){
      int cnt;
      int *L;
      int *l;
      /* find out how many gquads we might expect in the interval [i,j] */
      int num_gquads = get_gquad_count(S1, i, j);
      num_gquads++;
      L = (int *)space(sizeof(int) * num_gquads);
      l = (int *)space(sizeof(int) * num_gquads * 3);
      L[0] = -1;

      get_gquad_pattern_exhaustive(S1, i, j, P, L, l, threshold - best_energy);

      for(cnt = 0; L[cnt] != -1; cnt++){
        new_state = copy_state(state);

        make_gquad(i, L[cnt], &(l[3*cnt]), new_state);
        new_state->partial_energy += part_energy;
        new_state->partial_energy += element_energy;
        /* new_state->best_energy =
           hairpin[unpaired] + element_energy + best_energy; */
        push(Stack, new_state);
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
repeat(int i, int j, STATE * state, int part_energy, int temp_energy)
{
  /* routine to find stacks, bulges, internal loops and  multiloops */
  /* within interval closed by basepair i,j */

  STATE *new_state;
  INTERVAL *new_interval;

  register int  k, p, q, energy, new;
  register int  mm;
  register int  no_close, type, type_2;
  int           rt;
  int           dangle_model  = P->model_details.dangles;
  int           noLP          = P->model_details.noLP;
  int           noGUclosure   = P->model_details.noGUclosure;

  type = ptype[indx[j]+i];
  if (type==0) fprintf(stderr, "repeat: Warning: %d %d can't pair\n", i,j);

  no_close = (((type == 3) || (type == 4)) && noGUclosure);

  if (noLP) /* always consider the structure with additional stack */
    if ((i+turn+2<j) && ((type_2 = ptype[indx[j-1]+i+1]))) {
      new_state = copy_state(state);
      make_pair(i, j, new_state);
      make_pair(i+1, j-1, new_state);
      new_interval = make_interval(i+1, j-1, 2);
      push(new_state->Intervals, new_interval);
      if(SAME_STRAND(i,i+1) && SAME_STRAND(j-1,j))
        energy = E_IntLoop(0, 0, type, rtype[type_2],S1[i+1],S1[j-1],S1[i+1],S1[j-1], P);

      new_state->partial_energy += part_energy;
      new_state->partial_energy += energy;
      /* new_state->best_energy = new + best_energy; */
      push(Stack, new_state);
      if (i==1 || state->structure[i-2]!='('  || state->structure[j]!=')')
        /* adding a stack is the only possible structure */
        return;
    }

  best_energy += part_energy; /* energy of current structural element */
  best_energy += temp_energy; /* energy from unpushed interval */

  for (p = i + 1; p <= MIN2 (j-2-turn,  i+MAXLOOP+1); p++) {
    int minq = j-i+p-MAXLOOP-2;
    if (minq<p+1+turn) minq = p+1+turn;
    for (q = j - 1; q >= minq; q--) {
      if ((noLP) && (p==i+1) && (q==j-1)) continue;

      type_2 = ptype[indx[q]+p];
      if (type_2==0) continue;

      if (noGUclosure)
        if (no_close||(type_2==3)||(type_2==4))
          if ((p>i+1)||(q<j-1)) continue;  /* continue unless stack */

      if (SAME_STRAND(i,p) && SAME_STRAND(q,j)) {
        energy = E_IntLoop(p-i-1, j-q-1, type, rtype[type_2],
                            S1[i+1],S1[j-1],S1[p-1],S1[q+1], P);

        new = energy + c[indx[q]+p];

        if (new + best_energy <= threshold) {
          /* stack, bulge, or interior loop */

          new_state = copy_state(state);
          make_pair(i, j, new_state);
          make_pair(p, q, new_state);

          new_interval = make_interval(p, q, 2);
          push(new_state->Intervals, new_interval);
          new_state->partial_energy += part_energy;
          new_state->partial_energy += energy;
          /* new_state->best_energy = new + best_energy; */
          push(Stack, new_state);
        }
      }/*end of if block */
    } /* end of q-loop */
  } /* end of p-loop */

  if (!SAME_STRAND(i,j)) { /*look in fc*/
    rt = rtype[type];
    element_energy=0;
    if (dangle_model)
      element_energy = E_ExtLoop(rt, (SAME_STRAND(j-1,j)) ? S1[j-1] : -1, (SAME_STRAND(i,i+1)) ? S1[i+1] : -1, P);
    else
      element_energy = E_ExtLoop(rt, -1, -1, P);

    if (fc[i+1] + fc[j-1] +element_energy + best_energy  <= threshold)
      {
        INTERVAL *interval1, *interval2;

        new_state = copy_state(state);
        interval1 = make_interval(i+1, cut_point-1, 4);
        interval2 = make_interval(cut_point, j-1, 5);
        if (cut_point-i < j-cut_point) { /* push larger interval first */
          push(new_state->Intervals, interval1);
          push(new_state->Intervals, interval2);
        } else {
          push(new_state->Intervals, interval2);
          push(new_state->Intervals, interval1);
        }
        make_pair(i, j, new_state);
        new_state->partial_energy += part_energy;
        new_state->partial_energy += element_energy;
        push(Stack, new_state);
      }
  }

  mm = P->MLclosing;
  rt = rtype[type];

  for (k = i + 1 + turn; k <= j - 2 - turn; k++)  {
    /* multiloop decomposition */

    element_energy = mm;
    if (dangle_model)
      element_energy = E_MLstem(rt, S1[j-1], S1[i+1], P) + mm;
    else
      element_energy = E_MLstem(rt, -1, -1, P) + mm;

    if ((fML[indx[k] + i+1] + fM1[indx[j-1] + k+1] +
        element_energy + best_energy)  <= threshold)
      {
        INTERVAL *interval1, *interval2;
        new_state = copy_state(state);
        interval1 = make_interval(i+1, k, 1);
        interval2 = make_interval(k+1, j-1, 3);
        if (k-i+1 < j-k-2) { /* push larger interval first */
          push(new_state->Intervals, interval1);
          push(new_state->Intervals, interval2);
        } else {
          push(new_state->Intervals, interval2);
          push(new_state->Intervals, interval1);
        }
        make_pair(i, j, new_state);
        new_state->partial_energy += part_energy;
        new_state->partial_energy += element_energy;
        /* new_state->best_energy = fML[indx[k] + i+1] + fM1[indx[j-1] + k+1]
           + element_energy + best_energy; */
        push(Stack, new_state);
      }
  } /* end of k-loop */


  if (SAME_STRAND(i,j)) {
    if (no_close) element_energy = FORBIDDEN;
    else
      element_energy = E_Hairpin(j-i-1, type, S1[i+1], S1[j-1], sequence+i-1, P);
    if (element_energy + best_energy <= threshold) {
      /* hairpin structure */

      new_state = copy_state(state);
      make_pair(i, j, new_state);
      new_state->partial_energy += part_energy;
      new_state->partial_energy += element_energy;
      /* new_state->best_energy =
         hairpin[unpaired] + element_energy + best_energy; */
      push(Stack, new_state);
    }

    if(with_gquad){
      /* now we have to find all loops where (i,j) encloses a gquad in an interior loops style */
      int cnt, *p, *q, *en;
      p = q = en = NULL;
      en = E_GQuad_IntLoop_exhaustive(i, j, &p, &q, type, S1, ggg, threshold - best_energy, indx, P);
      for(cnt = 0; p[cnt] != -1; cnt++){
          new_state = copy_state(state);
          make_pair(i, j, new_state);

          new_interval = make_interval(p[cnt], q[cnt], 6);
          push(new_state->Intervals, new_interval);
          new_state->partial_energy += part_energy;
          new_state->partial_energy += en[cnt];
          /* new_state->best_energy = new + best_energy; */
          push(Stack, new_state);
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

PRIVATE char *costring(char *string)
{
  char *ctmp;
  int len;
  len = strlen(string);
  ctmp = (char *)space((len+2) * sizeof(char));
  /* first sequence */
  if (cut_point<=0) {
    (void) strncpy(ctmp, string, len);
    return ctmp;
  }
  (void) strncpy(ctmp, string, cut_point-1);
  /* spacer */
  ctmp[cut_point-1] = '&';
  /* second sequence */
  (void) strcat(ctmp, string+cut_point-1);

  return ctmp;
}

/*---------------------------------------------------------------------------*/
/* Well, that is the end!----------------------------------------------------*/
/*---------------------------------------------------------------------------*/
