/*
  $Log: subopt.c,v $
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

#define true	  1
#define false	  0

#define PUBLIC
#define PRIVATE	  static

/*@unused@*/
PRIVATE char rcsid[] = "$Id: subopt.c,v 1.13 2001/04/05 07:35:08 ivo Exp $";

/*Typedefinitions ---------------------------------------------------------- */

typedef struct {
    char *structure;
    LIST *Intervals;
    int partial_energy;
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

typedef struct {
  float e_o_s;                            /* e_o_s ... energy of structure */
  char *structure;
} SOLUTION;


PRIVATE void make_pair(int i, int j, STATE *state);
PRIVATE INTERVAL *make_interval (int i, int j, int ml);
/*@out@*/ PRIVATE STATE *make_state(/*@only@*/LIST *Intervals,
			  /*@only@*/ /*@null@*/ char *structure,
			  int partial_energy);
PRIVATE SOLUTION *make_solution (float e, /*@only@*/ char* structure);
PRIVATE STATE *copy_state(STATE * state);
PRIVATE void print_state(STATE * state);
PRIVATE void print_stack(LIST * list);
/*@only@*/ PRIVATE LIST *make_list(void);
PRIVATE  void push(LIST * list, /*@only@*/ void *data);
PRIVATE  void *pop(LIST * list);
PUBLIC  void subopt (char *seq, char *sequence, int delta);
PRIVATE int  best_attainable_energy(STATE * state);
PRIVATE void scan_interval(int i, int j, int array_flag, STATE * state);
PRIVATE  void free_pair_node(/*@only@*/ PAIR * node);
PRIVATE  void free_interval_node(/*@only@*/ INTERVAL * node);
PRIVATE  void free_state_node(/*@only@*/ STATE * node);
PRIVATE  void free_solution_node(/*@only@*/ SOLUTION* node);
PRIVATE void push_back(STATE * state);
PRIVATE char* get_structure(STATE * state);
PRIVATE int compare(SOLUTION * solution1, SOLUTION* solution2);
PRIVATE  void make_output(void);
PRIVATE void repeat(int i, int j, STATE * state,
		     int part_energy, int temp_energy);

/*Globals ------------------------------------------------------------------ */
/* options that may be modified by RNAsubopt.c */
int sorted=0;                           /* output sorted by energy */
int LODOS_ONLY=0; 

#define MAXALPHA 20	                /* maximal length of alphabet */

PRIVATE LIST *Stack;
PRIVATE LIST *SolutionList;
PRIVATE int nopush;
PRIVATE int best_energy;                /* best_energy = remaining energy */

extern int *f5;                         /* energy of 5 end */
extern int *c;		                /* energy array, given that i-j pair */
extern int *fML;		        /* multi-loop auxiliary energy array */
extern int *fM1;                /* another multi-loop auxiliary energy array */
extern int *indx;   /* index for moving in the triangle matrices c[] and f[] */
extern short *S, *S1;

extern int stack[NBPAIRS+1][NBPAIRS+1];
extern int hairpin[31];
extern int bulge[MAXLOOP+1];
extern int internal_loop[MAXLOOP+1];
extern int mismatchI[NBPAIRS+1][5][5];
extern int mismatchH[NBPAIRS+1][5][5];
extern int dangle5[NBPAIRS+1][5];
extern int dangle3[NBPAIRS+1][5];
extern int F_ninio[5];
extern double lxc;
extern int MLbase;
extern int MLintern[];
extern int MLclosing;
extern char *ptype;

extern float energy_of_struct(char *, char *);
extern int   LoopEnergy(int n1, int n2, int type, int type_2,
                         int si1, int sj1, int sp1, int sq1);
extern int   HairpinE(int i, int j, int type, const char *string);

PRIVATE int length;
PRIVATE int minimal_energy;                           /* minimum free energy */
PRIVATE int element_energy;       /* internal energy of a structural element */
PRIVATE int threshold;                             /* minimal_energy + delta */
PRIVATE char *sequence;
PUBLIC  float print_energy = 9999;  /* printing threshold for use with logML */

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

PRIVATE SOLUTION *
make_solution(float e_o_s, char* structure)
{
  SOLUTION* solution;
  
  solution = lst_newnode(sizeof(SOLUTION));
  solution->e_o_s = e_o_s;
  solution->structure = structure;

  return solution;
}

/*---------------------------------------------------------------------------*/
/*@unused@*/
PRIVATE void
free_pair_node(PAIR * node)
{
  lst_freenode(node);
}

/*---------------------------------------------------------------------------*/

PRIVATE void
free_interval_node(INTERVAL * node)
{
  lst_freenode(node);
}

/*---------------------------------------------------------------------------*/

PRIVATE void
free_solution_node(SOLUTION * node)
{
  free(node->structure);
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
	    int partial_energy)
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
	  sum += f5[next->j];
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
compare(SOLUTION * solution1, SOLUTION * solution2)
{
  if (solution1->e_o_s > solution2->e_o_s)
    return 1;
  if (solution1->e_o_s < solution2->e_o_s)
    return -1;
  return strcmp(solution1->structure,solution2->structure);
}

/*---------------------------------------------------------------------------*/

PRIVATE void 
make_output(void)  /* prints stuff */
{
  SOLUTION *next;
  int degeneracy, state;
  float this_energy, old_energy;

#ifdef BLAH
  if (!LODOS_ONLY && !GAPSTAT) {
    if (GAPS >= 0) {
      printf("\n\n the %d best structures between %4.2f and %4.2f kcal",
	      GAPS+1, minimal_energy/100., threshold/100.);
    }
    else {
      printf("\n\n all best structures between %4.2f and %4.2f kcal", 
	      minimal_energy/100.,
	      threshold/100.);
    }
    printf("\n\n    %s\n", sequence);
  }
#endif
  old_energy = minimal_energy/100. - 1.;
  state = -1;
  degeneracy = 1;
  for (next = lst_first(SolutionList); next; next = lst_next(next)) {
    state++;
    this_energy = next->e_o_s;
    if (LODOS_ONLY)
      if (this_energy == old_energy)
	degeneracy++;
      else { 
	printf("%4.2f %d\n", this_energy, degeneracy);
	degeneracy = 1;
      }
    else 
      printf("%s %6.2f\n", next->structure, this_energy);
    old_energy = this_energy;
  }
}

/*---------------------------------------------------------------------------*/
/* start of subopt backtracking ---------------------------------------------*/
/*---------------------------------------------------------------------------*/

PUBLIC void subopt(char *seq, char *structure, int delta)
{
  STATE *state;
  LIST *Intervals;
  INTERVAL *interval;
  SOLUTION *new_solution;

  int maxlevel, count, partial_energy, old_dangles;
  float structure_energy, min_en;
  char* struc;

  sequence = seq;
  length = strlen(sequence);

  struc = (char *) space(sizeof(char)*(length+1));
  if (fold_constrained) strncpy(struc, structure, length);
  /* do mfe folding to get fill arrays and get ground state energy  */
  /* in case dangles is neither 0 or 2, set dangles=2 while folding */
  old_dangles = dangles;
  if ((dangles!=0) && (dangles != 2)) dangles = 2;
  
  initialize_fold(length);
  min_en = fold(sequence, struc);
  
  dangles = old_dangles;
  /* re-evaluate in case we're using logML etc */
  min_en = energy_of_struct(sequence, struc);
  free(struc);
  print_energy += min_en;

  printf("%s %6d %6d\n", sequence, (int) (-0.1+100*min_en), delta); 

  /* Initialize ------------------------------------------------------------ */
  
  maxlevel = 0;
  count = 0;
  partial_energy = 0;

  make_pair_matrix();
  
  /* Initialize the stack ------------------------------------------------- */

  minimal_energy = f5[length];
  threshold = minimal_energy + delta;
  
  Stack = make_list();		                                   /* anchor */
  Intervals = make_list();	                           /* initial state: */
  interval = make_interval(1, length, 0);          /* interval [1,length,0] */
  push(Intervals, interval);
  state = make_state(Intervals, NULL, partial_energy);
  /* state->best_energy = minimal_energy; */
  push(Stack, state);

  /* SolutionList stores the suboptimal structures found */
                            
  SolutionList = make_list(); 
  
  /* end initialize ------------------------------------------------------- */

  
  while (1)			    /* forever, til nothing remains on stack */
    {
      maxlevel = (Stack->count > maxlevel ? Stack->count : maxlevel);
      
      if (LST_EMPTY (Stack))	           /* we are done! clean up and quit */
	{
	  fprintf(stderr, "maxlevel: %d\n", maxlevel);

	  lst_kill(Stack, free_state_node);

	  if (sorted) { 
	    /* sort structures by energy */
	    lst_mergesort(SolutionList, compare);    

	    make_output(); 
	  }

	  /* listkillroutines  -------------------------------------------- */ 
	     
	  lst_kill(SolutionList, free_solution_node); 

	  break;
	}
	    
      /* pop the last element ---------------------------------------------- */
	  
      state = pop(Stack);	               /* current state to work with */
      
      if (LST_EMPTY(state->Intervals))   
	{
	  /* state has no intervals left: we got a solution */

	  count++;		                          
	  structure = get_structure(state);
	  structure_energy = state->partial_energy / 100.;

#ifdef CHECK_ENERGY
	  structure_energy = energy_of_struct(sequence, structure);

	  if (!logML)
	    if ((float) (state->partial_energy / 100.) != structure_energy) {
	      fprintf(stderr, "%s %6.2f %6.2f\n", structure,
		      state->partial_energy / 100., structure_energy );
	      exit(1);
	    }
#endif
	  if (logML || (dangles==1) || (dangles==3)) { /* recalc energy */
	    structure_energy = energy_of_struct(sequence, structure);
	  }
	  
	  if (structure_energy>print_energy) {
	    free(structure);
	  } else {
	    if (!sorted) { 
	      /* print and forget */
	      printf("%s %6.2f\n", structure, structure_energy);
	      free(structure);
	    }
	    else {
	      /* store solution */
	      new_solution = make_solution(structure_energy, structure);
	      push(SolutionList, new_solution);
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

  /* free arrays left over from fold() */
  free(S); free(S1);
  free_arrays();
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
  
  if (j < i + TURN + 1)	{                       /* minimal structure element */
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
      fi = fM1[indx[j-1] + i] + MLbase;
    else
      fi = fML[indx[j-1] + i] + MLbase;
    
    if (fi + best_energy <= threshold) {
      /* no basepair, nibbling of 3'-end */
      
      new_state = copy_state(state);
      new_interval = make_interval(i, j-1, array_flag);
      push(new_state->Intervals, new_interval);
      new_state->partial_energy += MLbase;
      /* new_state->best_energy = fi + best_energy; */
      push(Stack, new_state);
    }

    type = ptype[indx[j]+i];
    
    if (type) { /* i,j may pair */
      
      element_energy = MLintern[type];
      
      if ( type && dangles ) {                        /* dangling ends */
	if (i > 1)
	  element_energy +=  dangle5[type][S1[i-1]];
	if (j < length)
	  element_energy += dangle3[type][S1[j+1]];
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

    for ( k = i+TURN+1 ; k <= j-1-TURN ; k++) {
      /* Multiloop decomposition if i,j contains more than 1 stack */
	   
      type = ptype[indx[j]+k+1];
      if (type==0) continue;
      
      element_energy = MLintern[type];
      if (dangles)
	element_energy += dangle3[type][S1[j+1]] + dangle5[type][S1[k]];
      
      
      if (fML[indx[k]+i] + c[indx[j] + k+1] +
	  element_energy + best_energy <= threshold) {
	temp_state = copy_state (state);
	new_interval = make_interval (i, k, 1);
	push (temp_state->Intervals, new_interval);
	repeat(k+1, j, temp_state, element_energy, fML[indx[k]+i]);
	free_state_node(temp_state);
      }
    }
    
    for (k = i ; k <= j-1-TURN; k++) {
      /* Multiloop decomposition if i,j contains only 1 stack */
      type = ptype[indx[j]+k+1];
      if (type==0) continue;
      
      element_energy = MLintern[type] + MLbase*(k-i+1);
      if (dangles)
	element_energy += dangle3[type][S1[j+1]] + dangle5[type][S1[k]];
      
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
  
  if (array_flag == 0)
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
      
      for (k = j-TURN-1; k > 1; k--) {
	
	type = ptype[indx[j]+k];
	if (type==0)   continue;
	
	/* k and j pair */
	if (dangles) {
	  element_energy =  dangle5[type][S1[k - 1]];
	  if (j < length)
	    element_energy += dangle3[type][S1[j+1]];
	}
	else                                                 /* no dangles */
	  element_energy = 0;
	
	if (type>2)
	  element_energy += TerminalAU;

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
	if (dangles && (j < length)) {
	  element_energy = dangle3[type][S1[j+1]];
	}
	else
	  element_energy = 0;
	
	if (type>2)
	  element_energy += TerminalAU;
	
	if (c[indx[j]+1] + element_energy + best_energy <= threshold)
	  repeat(1, j, state, element_energy, 0);
      }
    }                                                     /* array_flag == 0 */
  
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
    if ((i+TURN+2<j) && ((type_2 = ptype[indx[j-1]+i+1]))) { 
      new_state = copy_state(state);
      make_pair(i, j, new_state);
      make_pair(i+1, j-1, new_state);
      new_interval = make_interval(i+1, j-1, 2);
      push(new_state->Intervals, new_interval);
      energy = LoopEnergy(0, 0, type, rtype[type_2], 
			  S1[i+1],S1[j-1],S1[i+1],S1[j-1]);
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
  
  for (p = i + 1; p <= MIN2 (j-2-TURN,  i+MAXLOOP+1); p++) {
    int minq = j-i+p-MAXLOOP-2;
    if (minq<p+1+TURN) minq = p+1+TURN;
    for (q = j - 1; q >= minq; q--) {
      if ((noLonelyPairs) && (p==i+1) && (q==j-1)) continue;
      
      type_2 = ptype[indx[q]+p];
      if (type_2==0) continue;

      if (no_closingGU) 
	if (no_close||(type_2==3)||(type_2==4))
	  if ((p>i+1)||(q<j-1)) continue;  /* continue unless stack */
      
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
    }                                                  /* end of q-loop */ 
  }                                                    /* end of p-loop */
  
  mm = MLclosing + MLintern[type];
  rt = rtype[type];
  
  for (k = i + 1 + TURN; k <= j - 1 - TURN; k++)  {
    /* multiloop decomposition */
    
    element_energy = mm;
    if (dangles)
      element_energy += dangle3[rt][S1[i+1]] + dangle5[rt][S1[j-1]];
    
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
  
  
  if (no_close) energy = FORBIDDEN;
  else
    energy = HairpinE(i, j, type, sequence);
  
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
  
  best_energy -= part_energy;
  best_energy -= temp_energy;
  return;
}

/*---------------------------------------------------------------------------*/
/* Well, that is the end!----------------------------------------------------*/
/*---------------------------------------------------------------------------*/
