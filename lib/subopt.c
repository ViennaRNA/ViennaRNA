/*
  $Log: subopt.c,v $
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
#include "subopt.h"
#include "params.h"

#define true	  1
#define false	  0

#define PUBLIC
#define PRIVATE	  static
#define SAME_STRAND(I,J) (((I)>=cut_point)||((J)<cut_point))

/*@unused@*/
PRIVATE char UNUSED rcsid[] = "$Id: subopt.c,v 1.23 2008/03/31 15:06:49 ivo Exp $";

/*Typedefinitions ---------------------------------------------------------- */

typedef struct {
  char *structure;
  LIST *Intervals;
  int partial_energy;
  int is_duplex;
  /* int best_energy;   */ /* best attainable energy */
} STATE;

typedef struct {
  int i;
  int j;
} PAIR;

typedef struct {
  int i;
  int j;
  int array_flag;
} INTERVAL;

#define MAXDOS 1000
PUBLIC   int density_of_states[MAXDOS+1];
PRIVATE  void make_pair(int i, int j, STATE *state);
PRIVATE  INTERVAL *make_interval (int i, int j, int ml);
/*@out@*/ PRIVATE STATE *make_state(/*@only@*/LIST *Intervals,
				    /*@only@*/ /*@null@*/ char *structure,
				    int partial_energy, int is_duplex);
PRIVATE  STATE *copy_state(STATE * state);
PRIVATE  void print_state(STATE * state);
PRIVATE  void UNUSED print_stack(LIST * list);
/*@only@*/ PRIVATE LIST *make_list(void);
PRIVATE  void push(LIST * list, /*@only@*/ void *data);
PRIVATE  void *pop(LIST * list);
PUBLIC   SOLUTION *subopt (char *seq, char *sequence, int delta, FILE *fp);
PRIVATE  int  best_attainable_energy(STATE * state);
PRIVATE  void scan_interval(int i, int j, int array_flag, STATE * state);
PRIVATE  void free_interval_node(/*@only@*/ INTERVAL * node);
PRIVATE  void free_state_node(/*@only@*/ STATE * node);
PRIVATE  void push_back(STATE * state);
PRIVATE  char* get_structure(STATE * state);
PRIVATE  int compare(const void *solution1, const void *solution2);
PRIVATE  void make_output(SOLUTION *SL, FILE *fp);
PRIVATE  char *costring(char *string);
PRIVATE  void repeat(int i, int j, STATE * state,
		     int part_energy, int temp_energy);

/*Globals ------------------------------------------------------------------ */
/* options that may be modified by RNAsubopt.c */
int subopt_sorted=0;                           /* output sorted by energy */

#define MAXALPHA 20	                /* maximal length of alphabet */

PRIVATE int turn;
PRIVATE LIST *Stack;
PRIVATE int nopush;
PRIVATE int best_energy;                /* best_energy = remaining energy */

PRIVATE int *f5;                         /* energy of 5 end */
PRIVATE int *c;		                /* energy array, given that i-j pair */
PRIVATE int *fML;		        /* multi-loop auxiliary energy array */
PRIVATE int *fM1;                /* another multi-loop auxiliary energy array */
PRIVATE int *fc;                     /*energy array, from i (j)  to cut*/
PRIVATE int *indx;   /* index for moving in the triangle matrices c[] and f[] */
PRIVATE short *S, *S1;

PRIVATE char *ptype;

PRIVATE const paramT *P;
/* extern float energy_of_struct(char *, char *); */
extern int  LoopEnergy(int n1, int n2, int type, int type_2,
		       int si1, int sj1, int sp1, int sq1);
extern int  HairpinE(int size, int type, int si1, int sj1, const char *string);
extern void export_fold_arrays(int **f5_p, int **c_p, int **fML_p,
			       int **fM1_p, int **indx_p, char **ptype_p);
extern void export_cofold_arrays(int **f5_p, int **c_p, int **fML_p,
				 int **fM1_p, int **fc_p, int **indx_p,
				 char **ptype_p);
extern int  uniq_ML;
extern  int cut_point;   /* set to first pos of second seq for cofolding */

PRIVATE int length;
PRIVATE int minimal_energy;                           /* minimum free energy */
PRIVATE int element_energy;       /* internal energy of a structural element */
PRIVATE int threshold;                             /* minimal_energy + delta */
PRIVATE char *sequence;
PUBLIC  double print_energy = 9999; /* printing threshold for use with logML */

/* some needful things for subopt_circ */
/* take int circ from fold.c	 */
extern	int circ;
PUBLIC	SOLUTION *subopt_circ(char *seq, char *sequence, int delta, FILE *fp);
PRIVATE int *fM2;	 /* energies of M2 */
PUBLIC	int	Fc, FcH, FcI, FcM;		/* parts of the exterior loop energies */

PRIVATE void encode_seq(char *sequence) {
  unsigned int i,l;

  l = strlen(sequence);
  S = (short *) space(sizeof(short)*(l+2));
  S1= (short *) space(sizeof(short)*(l+2));
  /* S1 exists only for the special X K and I bases and energy_set!=0 */
  S[0] = (short) l;

  for (i=1; i<=l; i++) { /* make numerical encoding of sequence */
    S[i]= (short) encode_char(toupper(sequence[i-1]));
    S1[i] = alias[S[i]];   /* for mismatches of nostandard bases */
  }
  /* for circular folding add first base at position n+1 and last base at position 0 in S1*/
  S[l+1] = S[1]; S1[l+1]=S1[1]; S1[0] = S1[l];
}

/*---------------------------------------------------------------------------*/
/*List routines--------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

PRIVATE void
make_pair(int i, int j, STATE *state)
{
  state->structure[i-1] = '(';
  state->structure[j-1] = ')';
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
	sum += (circ) ? Fc : f5[next->j];
      else if (next->array_flag == 1)
	sum += fML[indx[next->j] + next->i];
      else if (next->array_flag == 2)
	sum += c[indx[next->j] + next->i];
      else if (next->array_flag == 3)
	sum += fM1[indx[next->j] + next->i];
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

PUBLIC SOLUTION *subopt_circ(char *seq, char *structure, int delta, FILE *fp){
  circ = 1;
  return subopt(seq, structure, delta, fp);
}

PUBLIC SOLUTION *subopt(char *seq, char *structure, int delta, FILE *fp)
{
  STATE *state;
  LIST *Intervals;
  INTERVAL *interval;
  SOLUTION *SolutionList;

  unsigned long max_sol = 128, n_sol = 0;

  int maxlevel, count, partial_energy, old_dangles;
  double structure_energy, min_en, eprint;
  char* struc;

  sequence = seq;
  length = strlen(sequence);

  struc = (char *) space(sizeof(char)*(length+1));
  if (fold_constrained) strncpy(struc, structure, length);
  /* do mfe folding to get fill arrays and get ground state energy  */
  /* in case dangles is neither 0 or 2, set dangles=2 while folding */
  old_dangles = dangles;
  if ((dangles!=0) && (dangles != 2)) dangles = 2;

  turn = (cut_point<0) ? 3 : 0;
  uniq_ML = 1;
  (circ)? initialize_fold(length) : initialize_cofold(length);
  min_en = (circ) ? circfold(sequence, struc) : cofold(sequence, struc);
  (circ) ? export_circfold_arrays(&Fc, &FcH, &FcI, &FcM, &fM2, &f5, &c, &fML, &fM1, &indx, &ptype) : export_cofold_arrays(&f5, &c, &fML, &fM1, &fc, &indx, &ptype);

  dangles = old_dangles;
  /* re-evaluate in case we're using logML etc */
  min_en = (circ) ? energy_of_circ_struct(sequence, struc) : energy_of_struct(sequence, struc);
  free(struc);
  eprint = print_energy + min_en;
  if (fp) {
    char *SeQ;
    SeQ=costring(sequence);
    fprintf(fp, "%s %6d %6d\n", SeQ, (int) (-0.1+100*min_en), delta);
    free(SeQ);
  }
  make_pair_matrix();
  encode_seq(sequence);
  P = scale_parameters();

  /* Initialize ------------------------------------------------------------ */

  maxlevel = 0;
  count = 0;
  partial_energy = 0;

  /* Initialize the stack ------------------------------------------------- */

  minimal_energy = (circ) ? Fc : f5[length];
  threshold = minimal_energy + delta;

  Stack = make_list();		                                   /* anchor */
  Intervals = make_list();	                           /* initial state: */
  interval = make_interval(1, length, 0);          /* interval [1,length,0] */
  push(Intervals, interval);
  state = make_state(Intervals, NULL, partial_energy,0);
  /* state->best_energy = minimal_energy; */
  push(Stack, state);

  /* SolutionList stores the suboptimal structures found */

  SolutionList = (SOLUTION *) space(max_sol*sizeof(SOLUTION));

  /* end initialize ------------------------------------------------------- */


  while (1) {		    /* forever, til nothing remains on stack */

    maxlevel = (Stack->count > maxlevel ? Stack->count : maxlevel);

    if (LST_EMPTY (Stack))	           /* we are done! clean up and quit */
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

    state = pop(Stack);	               /* current state to work with */

    if (LST_EMPTY(state->Intervals))
      {
	int e;
	/* state has no intervals left: we got a solution */

	count++;
	structure = get_structure(state);
	structure_energy = state->partial_energy / 100.;

#ifdef CHECK_ENERGY
	structure_energy = (circ) ? energy_of_circ_struct(sequence, structure) : energy_of_struct(sequence, structure);

	if (!logML)
	  if ((double) (state->partial_energy / 100.) != structure_energy) {
	    fprintf(stderr, "%s %6.2f %6.2f\n", structure,
		    state->partial_energy / 100., structure_energy );
	    exit(1);
	  }
#endif
	if (logML || (dangles==1) || (dangles==3)) { /* recalc energy */
	  structure_energy = (circ) ? energy_of_circ_struct(sequence, structure) : energy_of_struct(sequence, structure);
	}

	e = (int) ((structure_energy-min_en)*10.);
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

      free_interval_node(interval);	/* free the current interval */
    }

    free_state_node(state);	             /* free the current state */
  } /* end of while (1) */

  /* free arrays left over from cofold() */
  free(S); free(S1);
  (circ) ? free_arrays():free_co_arrays();
  if (fp) { /* we've printed everything -- free solutions */
    SOLUTION *sol;
    for (sol=SolutionList; sol->structure != NULL; sol++)
      free(sol->structure);
    free(SolutionList);
    SolutionList = NULL;
  }

  return SolutionList;
}

/*---------------------------------------------------------------------------*/
/* Definitions---------------------------------------------------------------*/
/* these have to identical to defines in fold.c */
#define NEW_NINIO         1        /* use new asymetry penalty */
#define STACK_BULGE1      1	   /* stacking energies for bulges of size 1 */
#define MIN2(A, B)        ((A) < (B) ? (A) : (B))

/*---------------------------------------------------------------------------*/

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

      element_energy = P->MLintern[type];

      if ( type && dangles ) {                        /* dangling ends */
	if ((i > 1)&&(SAME_STRAND(i-1,i)) || circ)
	  element_energy +=  P->dangle5[type][S1[i-1]];
	if ((j < length)&&(SAME_STRAND(j,j+1)) || circ)
	  element_energy += P->dangle3[type][S1[j+1]];
      }

      cij = c[indx[j] + i] + element_energy;
      if (cij + best_energy <= threshold)
	repeat(i, j, state, element_energy, 0);
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

	type = ptype[indx[j]+k+1];
	if (type==0) continue;

	element_energy = P->MLintern[type];
	if (dangles) {
	  if (SAME_STRAND(j,j+1)) element_energy += P->dangle3[type][S1[j+1]];
	  if (SAME_STRAND(i-1,i)) element_energy +=P->dangle5[type][S1[k]];
	}
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
      type = ptype[indx[j]+k+1];
      if (type==0) continue;

      element_energy = P->MLintern[type] + P->MLbase*(k-i+1);
      if (dangles) {
	if (SAME_STRAND(j,j+1)) element_energy += P->dangle3[type][S1[j+1]];
	if (SAME_STRAND(k-1,k)) element_energy += P->dangle5[type][S1[k]];
      }

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

      if (nopush)
	if (!noLonelyPairs)
	  fprintf(stderr, "Oops, no solution in repeat!\n");
      return;
    }

  /* 0000000000000000000000000000000000000000000000000000000000000000000000 */

  if ((array_flag == 0) && !circ)
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

	type = ptype[indx[j]+k];
	if (type==0)   continue;

	/* k and j pair */
	if (dangles) {
	  element_energy =(SAME_STRAND(k-1,k))? P->dangle5[type][S1[k - 1]]:0;
	  if ((j < length)&&(SAME_STRAND(j,j+1)))
	    element_energy += P->dangle3[type][S1[j+1]];
	}
	else                                                 /* no dangles */
	  element_energy = 0;

	if (type>2)
	  element_energy += P->TerminalAU;


	if (!(SAME_STRAND(k,j)))/*&&(state->is_duplex==0))*/ {
	  element_energy+=DuplexInit;
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
	if (dangles && (j < length)&&(SAME_STRAND(j,j+1))) {
	  element_energy = P->dangle3[type][S1[j+1]];
	}
	else
	  element_energy = 0;

	if (!(SAME_STRAND(1,j))) element_energy+=DuplexInit;
	if (type>2)
	  element_energy += P->TerminalAU;

	if (c[indx[j]+1] + element_energy + best_energy <= threshold)
	  repeat(1, j, state, element_energy, 0);
      }
    }/* end array_flag == 0 && !circ*/
  /* or do we subopt circ? */
  else if(array_flag == 0){
    int k, l, p, q;
    /* if we've done everything right, we will never reach this case more than once	 */
    /* right after the initilization of the stack with ([1,n], empty, 0)						 */
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
    /* best energy should be 0 if we are here																				 */
    if(FcH + best_energy <= threshold){
      /* lets search for all exterior hairpin cases, that fit into our threshold barrier */
      /* we use index k,l to avoid confusion with i,j index of our state...							 */
      /* if we reach here, i should be 1 and j should be n respectively									 */
      for(k=i; k<j; k++)
	for (l=k+turn+1; l <= j; l++){
	  int kl, type, u, new_c, tmpE, no_close;
	  u = j-l + k-1;	/* get the hairpin loop length */
	  if(u<turn) continue;

	  kl = indx[l]+k;	/* just confusing these indices ;-) */
	  type = ptype[kl];
	  no_close = ((type==3)||(type==4))&&no_closingGU;
	  type=rtype[type];
	  if (!type) continue;
	  if (no_close) new_c = FORBIDDEN;
	  else{
	    /* now lets have a look at the hairpin energy */
	    char loopseq[10];
	    if (u<7){
	      strcpy(loopseq , sequence+l-1);
	      strncat(loopseq, sequence, k);
	    }
	    tmpE = HairpinE(u, type, S1[l+1], S1[k-1], loopseq);
	  }
	  if(c[kl] + tmpE + best_energy <= threshold){
	    /* what we really have to do is something like this, isn't it? */
	    /* we have to create a new state, with interval [k,l], then we */
	    /* add our loop energy as initial energy of this state and put */
	    /* the state onto the stack R... for further refinement...	 */
	    /* we also denote this new interval to be scanned in C	 */
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
	  int kl, type, u, new_c, tmpE, no_close;

	  kl = indx[l]+k;	/* just confusing these indices ;-) */
	  type = ptype[kl];
	  no_close = ((type==3)||(type==4))&&no_closingGU;
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
	      tmpE = LoopEnergy(u1, u2, type, type_2, S1[l+1], S1[k-1], S1[p-1], S1[q+1]);
	      if(c[kl] + c[indx[q]+p] + tmpE + best_energy <= threshold){
		/* ok, similar to the hairpin stuff, we add new states onto the stack R */
		/* but in contrast to the hairpin decomposition, we have to add two new */
		/* intervals, enclosed by k,l and p,q respectively and we also have to */
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
      /* this decomposition will be somehow more complicated...so lets see what we do here...	   */
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
  }	/* thats all folks for the circular case... */

  if (nopush)
    push_back(state);
  return;
}

/*---------------------------------------------------------------------------*/

PRIVATE void
repeat(int i, int j, STATE * state, int part_energy, int temp_energy)
{
  /* routine to find stacks, bulges, internal loops and  multiloops */
  /* within interval closed by basepair i,j */

  STATE *new_state;
  INTERVAL *new_interval;

  register int k, p, q, energy, new;
  register int mm;
  register int no_close, no_close_2, type, type_2;
  int  rt;

  no_close_2 = 0;

  type = ptype[indx[j]+i];
  if (type==0) fprintf(stderr, "repeat: Warning: %d %d can't pair\n", i,j);

  no_close = (((type == 3) || (type == 4)) && no_closingGU);

  if (noLonelyPairs) /* always consider the structure with additional stack */
    if ((i+turn+2<j) && ((type_2 = ptype[indx[j-1]+i+1]))) {
      new_state = copy_state(state);
      make_pair(i, j, new_state);
      make_pair(i+1, j-1, new_state);
      new_interval = make_interval(i+1, j-1, 2);
      push(new_state->Intervals, new_interval);
      if (SAME_STRAND(i,i+1) && SAME_STRAND(j-1,j)) energy = LoopEnergy(0, 0, type, rtype[type_2],S1[i+1],S1[j-1],S1[i+1],S1[j-1]);

      new_state->partial_energy += part_energy;
      new_state->partial_energy += energy;
      /* new_state->best_energy = new + best_energy; */
      push(Stack, new_state);
      if (i==1 || state->structure[i-2]!='('  || state->structure[j]!=')')
	/* adding a stack is the only possible structure */
	return;
    }

  best_energy += part_energy;        /* energy of current structural element */
  best_energy += temp_energy;               /* energy from unpushed interval */

  for (p = i + 1; p <= MIN2 (j-2-turn,  i+MAXLOOP+1); p++) {
    int minq = j-i+p-MAXLOOP-2;
    if (minq<p+1+turn) minq = p+1+turn;
    for (q = j - 1; q >= minq; q--) {
      if ((noLonelyPairs) && (p==i+1) && (q==j-1)) continue;

      type_2 = ptype[indx[q]+p];
      if (type_2==0) continue;

      if (no_closingGU)
	if (no_close||(type_2==3)||(type_2==4))
	  if ((p>i+1)||(q<j-1)) continue;  /* continue unless stack */

      if (SAME_STRAND(i,p) && SAME_STRAND(q,j)) {
	energy = LoopEnergy(p-i-1, j-q-1, type, rtype[type_2],
			    S1[i+1],S1[j-1],S1[p-1],S1[q+1]);

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
    }                                                  /* end of q-loop */
  }                                                    /* end of p-loop */

  if (!SAME_STRAND(i,j)) { /*look in fc*/
    int fci, fcj;
    rt = rtype[type];
    element_energy=0;
    if (dangles) {
      if (SAME_STRAND(i,i+1))
	element_energy += P->dangle3[rt][S1[i+1]];
      if (SAME_STRAND(j-1,j))
	element_energy+=P->dangle5[rt][S1[j-1]];
    }
    if (type>2) element_energy += P->TerminalAU;
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

  mm = P->MLclosing + P->MLintern[type];
  rt = rtype[type];

  for (k = i + 1 + turn; k <= j - 2 - turn; k++)  {
    /* multiloop decomposition */

    element_energy = mm;
    if (dangles)
      element_energy += P->dangle3[rt][S1[i+1]] + P->dangle5[rt][S1[j-1]];

    if (fML[indx[k] + i+1] + fM1[indx[j-1] + k+1] +
	element_energy + best_energy  <= threshold)
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
  }                                                     /* end of k-loop */


  if (SAME_STRAND(i,j)) {
    if (no_close) energy = FORBIDDEN;
    else
      energy = HairpinE(j-i-1, type, S1[i+1], S1[j-1], sequence+i-1);

    if (energy + best_energy <= threshold) {
      /* hairpin structure */

      new_state = copy_state(state);
      make_pair(i, j, new_state);
      new_state->partial_energy += part_energy;
      new_state->partial_energy += energy;
      /* new_state->best_energy =
	 hairpin[unpaired] + element_energy + best_energy; */
      push(Stack, new_state);
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
