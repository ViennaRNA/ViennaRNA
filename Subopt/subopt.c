/*
  $Log: subopt.c,v $
  Revision 1.7  1998/03/30 14:39:31  ivo
  replaced BasePairs list with structure string in STATE
  save memory by not storing (and sorting) structures
  modified for use with ViennaRNA-1.2.1

  
  Revision 1.5  1997/08/26 12:57:55  walter
  *** empty log message ***

  Revision 1.4  1997/08/20 15:06:02  walter
  *** empty log message ***

  Revision 1.4  1997/08/19 14:09:52  walter
  *** empty log message ***

  Revision 1.3  1997/08/18 19:09:32  walter
  *** empty log message ***

  Revision 1.2  1997/08/05 00:02:44  walter
  *** empty log message ***

  Revision 1.1  1997/08/04 21:05:32  walter
  Initial revision

*/
/*
   suboptimal folding - Stefan Wuchty, Walter Fontana & Ivo Hofacker

                      subopt.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "utilities.h"
#include "energy_par.h"
#include "fold_vars.h"
#include "pair_mat.h"
#include "list.h"

PRIVATE char rcsid[] = "$Id: subopt.c,v 1.7 1998/03/30 14:39:31 ivo Exp $";

/*Typedefinitions ----------------------------------------------------------- */

#define NEW_NINIO   1  /* new asymetry penalty, must be the same as in fold.c */
#define BPLIST 0
#define NOSORT

typedef struct
  {
#if BPLIST
    LIST *BasePairs;
#else
    char *structure;
#endif
    LIST *Intervals;
    int partial_energy;
  }
STATE;

typedef struct
  {
    int i;
    int j;
  }
PAIR;

typedef struct
  {
    int i;
    int j;
    int array_flag;
  }
INTERVAL;

typedef struct
  {
    float e_o_s;                              /* e_o_s ... energy of structure */
    char *structure;
  }
SOLUTION;


PUBLIC void make_pair (int i, int j, STATE *state);
PUBLIC INTERVAL *make_interval (int i, int j, int ml);
PUBLIC STATE *make_state (LIST * Intervals,
#if BPLIST
			  LIST * BasePairs,
#else
			  char *structure,
#endif
			  int partial_energy);
PUBLIC SOLUTION *make_solution (float e, char* structure);
PUBLIC STATE *copy_state (STATE * state);
PUBLIC void print_state (STATE * state);
PUBLIC void print_stack (LIST * list);
PUBLIC LIST *make_list (void);
PUBLIC void push (LIST * list, void *data);
PUBLIC void *pop (LIST * list);
PUBLIC void subopt (void);
PUBLIC int best_attainable_energy (STATE * state);
PUBLIC void scan_interval (int i, int j, int array_flag, STATE * state);
PUBLIC void free_pair_node (PAIR * node);
PUBLIC void free_interval_node (INTERVAL * node);
PUBLIC void free_state_node (STATE * node);
PUBLIC void free_solution_node (SOLUTION* node);
PUBLIC void push_back (STATE * state);
PUBLIC char* get_structure (STATE * state);
PRIVATE int compare (SOLUTION * solution1, SOLUTION* solution2);
PUBLIC void make_output (void);
PRIVATE void repeat (int i, int j, STATE * state, int part_energy, int temp_energy);

/*Globals -------------------------------------------------------------------- */

#define MAXALPHA 20		                 /* maximal length of alphabet */

PRIVATE LIST *Stack;
PUBLIC LIST *SolutionList;
PRIVATE int nopush;
PRIVATE int best_energy;		     /* best_energy = remaining energy */

extern int MLbase;
extern int length;
extern int delta; 
extern int *f5;                                             /* energy of 5 end */
extern int *c;		                  /* energy array, given that i-j pair */
extern int *fML;		          /* multi-loop auxiliary energy array */
extern int *fM1;                  /* another multi-loop auxiliary energy array */
extern int *indx;     /* index for moving in the triangle matrices c[] and f[] */
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
extern int MLintern;
extern int MLclosing;
extern int TETRA_ENERGY[40];
extern int *BP;

extern char* sequence;
extern int GAPS; 
extern int GAPSTAT; 
extern int LODOS_ONLY; 
extern FILE *struc[3];
extern int **Gap;

extern float energy_of_struct (char *, char *);

PRIVATE int minimal_energy;                             /* minimum free energy */
PRIVATE int element_energy;         /* internal energy of a structural element */
PRIVATE int threshold;                               /* minimal_energy + delta */

/*-----------------------------------------------------------------------------*/
/*List routines----------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/

PUBLIC void
make_pair (int i, int j, STATE *state)
{
#if BPLIST
  PAIR *pair;
  
  pair = lst_newnode (sizeof (PAIR));
  pair->i = i;
  pair->j = j;
  push (state->BasePairs, pair);
#else
  state->structure[i-1] = '(';
  state->structure[j-1] = ')';
#endif
}

/*-----------------------------------------------------------------------------*/

PUBLIC INTERVAL *
make_interval (int i, int j, int array_flag)
{
  INTERVAL *interval;

  interval = lst_newnode (sizeof (INTERVAL));
  interval->i = i;
  interval->j = j;
  interval->array_flag = array_flag;
  return interval;
}

/*-----------------------------------------------------------------------------*/

PUBLIC SOLUTION *
make_solution (float e_o_s, char* structure)
{
  SOLUTION* solution;
  
  solution = lst_newnode (sizeof (SOLUTION));
  solution->e_o_s = e_o_s;
  solution->structure = structure;

  return solution;
}

/*-----------------------------------------------------------------------------*/

PUBLIC void
free_pair_node (PAIR * node)
{
  lst_freenode (node);
}

/*-----------------------------------------------------------------------------*/

PUBLIC void
free_interval_node (INTERVAL * node)
{
  lst_freenode (node);
}

/*-----------------------------------------------------------------------------*/

PUBLIC void
free_solution_node (SOLUTION * node)
{
  free (node->structure);
  lst_freenode (node);
}

/*-----------------------------------------------------------------------------*/

PUBLIC void
free_state_node (STATE * node)
{
#if BPLIST
  if (node->BasePairs)
    lst_kill (node->BasePairs, lst_freenode);
#endif
  free(node->structure);
  if (node->Intervals)
    lst_kill (node->Intervals, lst_freenode);
  lst_freenode (node);
}

/*-----------------------------------------------------------------------------*/

PUBLIC STATE *
make_state (LIST * Intervals,
#if BPLIST
	    LIST * BasePairs,
#else
	    char *structure,
#endif
	    int partial_energy)
{
  STATE *state;

  state = lst_newnode (sizeof (STATE));

  if (Intervals)
    state->Intervals = Intervals;
  else
    state->Intervals = lst_init ();
#if BPLIST
  if (BasePairs)
     state->BasePairs = BasePairs;
  else
     state->BasePairs = lst_init ();
#else
  if (structure)
    state->structure = structure;
  else {
    int i;
    state->structure = (char *) space(length+1);
    for (i=0; i<length; i++)
      state->structure[i] = '.';
  }
#endif
  
  state->partial_energy = partial_energy;

  return state;
}

/*-----------------------------------------------------------------------------*/

PUBLIC STATE *
copy_state (STATE * state)
{
  STATE *new_state;
  void *after;
#if BPLIST
  PAIR *new_basepair, *rec;
#endif
  INTERVAL *new_interval, *next;
  
  new_state = lst_newnode (sizeof (STATE));
  new_state->Intervals = lst_init ();
  new_state->partial_energy = state->partial_energy;

  if (state->Intervals->count)
    {
      after = LST_HEAD (new_state->Intervals);
      for ( next = lst_first (state->Intervals); next; next = lst_next (next))
	{
	  new_interval = lst_newnode (sizeof (INTERVAL));
	  /* new_interval->i = next->i; */
	  /* new_interval->j = next->j; */
	  /* new_interval->array_flag = next->array_flag; */
	  *new_interval = *next;
	  lst_insertafter (new_state->Intervals, new_interval, after);
	  after = new_interval;
	}
    }
#if BPLIST
  new_state->BasePairs = lst_init ();
  if (state->BasePairs->count)
    {
      after = LST_HEAD (new_state->BasePairs);
      for (rec = lst_first (state->BasePairs); rec; rec = lst_next (rec))
	{
	  new_basepair = lst_newnode (sizeof (PAIR));
	  /* new_basepair->i = rec->i; */
	  /* new_basepair->j = rec->j; */
	  *new_basepair = *rec;
	  lst_insertafter (new_state->BasePairs, new_basepair, after);
	  after = new_basepair;
	}
    }
#else
  new_state->structure = strdup(state->structure);
  if (!new_state->structure) nrerror("out of memory");
#endif
  return new_state;
}

/*-----------------------------------------------------------------------------*/

PUBLIC void
print_state (STATE * state)
{
  PAIR *rec;
  INTERVAL *next;

  if (state->Intervals->count)
    {
      printf ("%d intervals:\n", state->Intervals->count);
      for (next = lst_first (state->Intervals); next; next = lst_next (next))
	{
	  printf ("[%d,%d],%d ", next->i, next->j, next->array_flag);
	}
      printf ("\n");
    }
#if BPLIST
  if (state->BasePairs->count)
     {
       printf ("%d base pairs:\n", state->BasePairs->count);
       for (rec = lst_first (state->BasePairs); rec; rec = lst_next (rec))
	  {
	    printf ("(%d,%d) ", rec->i, rec->j);
	  }
       printf ("\n");
     }
#else
  printf ("partial structure: %s\n", state->structure);
#endif
  printf ("\n");
  printf (" partial_energy: %d\n", state->partial_energy);
  fflush (stdout);
}

/*-----------------------------------------------------------------------------*/

PUBLIC void
print_stack (LIST * list)
{
  void *rec;

  printf ("================\n");
  printf ("%d states\n", list->count);
  for (rec = lst_first (list); rec; rec = lst_next (rec))
    {
      printf ("state-----------\n");
      print_state (rec);
    }
  printf ("================\n");
}

/*-----------------------------------------------------------------------------*/

PUBLIC LIST *
make_list (void)
{
  return lst_init ();
}

/*-----------------------------------------------------------------------------*/

PUBLIC void
push (LIST * list, void *data)
{
  nopush = false;
  lst_insertafter (list, data, LST_HEAD (list));
}

/*-----------------------------------------------------------------------------*/

PUBLIC void *
pop (LIST * list)
{
  void *data;

  data = lst_first (list);
  lst_deletenext (list, LST_HEAD (list));
  return data;
}

/*-----------------------------------------------------------------------------*/
/*auxiliary routines-----------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/

PUBLIC int
best_attainable_energy (STATE * state)
{
   /* evaluation of best possible energy attainable within remaining intervals */
  
  register int sum;
  INTERVAL *next;

  sum = state->partial_energy;/* energy of already found elements of structure */

  for (next = lst_first (state->Intervals); next; next = lst_next (next))
    {
      if (next->array_flag == 0)
	  sum += f5[next->j];
      else if (next->array_flag == 1)
	  sum += fML[indx[next->j] + next->i];
      else if (next->array_flag == 2)
	  sum += c[indx[next->j] + next->i];
    }
    
  return sum;
}

/*-----------------------------------------------------------------------------*/

PUBLIC void
push_back (STATE * state)
{
  push (Stack, copy_state (state));
  return;
}

/*-----------------------------------------------------------------------------*/

PUBLIC char*
get_structure (STATE * state)
{
  PAIR *rec;
  char* structure;
  int x;

#if BPLIST
  structure = (char *) space (sizeof (char) * (length + 1)); 
  
  for (x = 0; x < length; x++)
     structure[x] = '.';
  
  for (rec = lst_first (state->BasePairs); rec; rec = lst_next (rec))
    {
      structure[rec->i - 1] = '(';
      structure[rec->j - 1] = ')';
    }
#else
  structure = strdup(state->structure);
#endif
  return structure;
}

/*-----------------------------------------------------------------------------*/
PUBLIC int 
compare (SOLUTION * solution1, SOLUTION * solution2)
{
  if (solution1->e_o_s > solution2->e_o_s)
    return 1;
  if (solution1->e_o_s < solution2->e_o_s)
    return -1;
  return strcmp(solution1->structure,solution2->structure);
}

/*-----------------------------------------------------------------------------*/

PUBLIC void 
make_output (void)  /* prints stuff */
{
  SOLUTION *next, *nnext, *old_ptr;
  int degeneracy, state;
  float gap, gap_to_ground, this_energy;

  if (!LODOS_ONLY && !GAPSTAT)
    {
      if (GAPS >= 0)
	{
	  printf ("\n\n the %d best structures between %4.2f and %4.2f kcal", GAPS+1,
		  minimal_energy/100.,
		  threshold/100.);
	}
      else
	{
	  printf ("\n\n all best structures between %4.2f and %4.2f kcal", 
		  minimal_energy/100.,
		  threshold/100.);
	}
      printf ("\n\n    %s\n", sequence);
    }
  
  gap_to_ground = 0.;
  gap = 0.;
  state = -1;
  for (next = lst_first(SolutionList); next; next = lst_next (next))
    {
      state++;
      degeneracy = 1;
      this_energy = next->e_o_s;
      old_ptr = next;

      for (nnext = lst_next(next); nnext; nnext = lst_next (nnext))
	{
	  if (nnext->e_o_s == this_energy)
	    {
	      degeneracy++;

	      if (GAPSTAT)
		{
		  if (state)
		    {
		      /* Gap[state][(int)(gap_to_ground*100+0.5)]++; */
		      fprintf (struc[state], "%s %4.2f %4.2f %4.2f %4.2f\n",
			       nnext->structure,
			       minimal_energy/100., this_energy, gap, gap_to_ground);
		      fflush (struc[state]);
		    }
		}
	      else if (!LODOS_ONLY)
		{
		  printf("%3d %s\n", state, nnext->structure);
		}
	      old_ptr = nnext;
	    }
	  else
	    break;
	}

      if (GAPSTAT)
	{
	  if (state)
	    {
	      Gap[state][(int)(gap_to_ground*100+0.5)]++;
	      fprintf (struc[state], "%s %4.2f %4.2f %4.2f %4.2f\n",
		       next->structure,
		       minimal_energy/100., this_energy, gap, gap_to_ground);
	      fflush (struc[state]);
	    }
	}

      if (!GAPSTAT)
	{
	  if (LODOS_ONLY)
	    {
	      printf("%4.2f %d", this_energy, degeneracy);
	    }
	  else
	    {
	      printf("%3d %s %4.2f %d", state, next->structure, this_energy, degeneracy);
	    }
	  printf(" %4.2f %4.2f", gap, gap_to_ground);
	}
      
      if (GAPS)
	{
	  if (state >= GAPS && GAPS != -1)
	    goto quit;
	  if (nnext)
	    {
	      gap = nnext->e_o_s - this_energy;
	      gap_to_ground += gap;
	    }
	}
      
      if (!GAPSTAT)
	  printf("\n");
      next = old_ptr;
    }
 quit:
  if (!GAPSTAT)
      printf ("\n\n");
  fflush (stdout);
}

/*-----------------------------------------------------------------------------*/
/* start of subopt backtracking -----------------------------------------------*/
/*-----------------------------------------------------------------------------*/

PUBLIC void
subopt (void)
{
  STATE *state;
  LIST *Intervals;
  INTERVAL *interval;
  SOLUTION *new_solution;

  int partial_energy;
  float structure_energy;
  register int maxlevel, count;
  register int scan_level;
  char* structure;

  /* Initialize -------------------------------------------------------------- */
  
  maxlevel = 0;
  count = 0;
  partial_energy = 0;

  make_pair_matrix();
  
  /* Initialize the stack ---------------------------------------------------- */

  minimal_energy = f5[length];
  threshold = minimal_energy + delta;
  
  Stack = make_list ();		                                     /* anchor */
  Intervals = make_list ();	                             /* initial state: */
  interval = make_interval (1, length, 0);            /* interval [1,length,0] */
  push (Intervals, interval);
  state = make_state (Intervals, NULL, partial_energy);
  push (Stack, state);

  /* SolutionList stores the suboptimal structures found */
                            
  SolutionList = make_list (); 
  
  /* end initialize ---------------------------------------------------------- */

#ifdef NOSORT
  printf("%s\n", sequence);   /* first line of output */
#endif
  
  while (1)			      /* forever, til nothing remains on stack */
    {
      maxlevel = (Stack->count > maxlevel ? Stack->count : maxlevel);
      
      if (LST_EMPTY (Stack))	     /* we are done! clean up and quit */
	{
	  fprintf(stderr, "maxlevel: %d\n", maxlevel);

	  lst_kill (Stack, free_state_node);
#ifndef NOSORT	  
	  lst_mergesort (SolutionList, compare);    /* sorting of structures */
	                                            /* by their energies */
	    
	  /* outputroutines ----------------------------------------------- */ 

	  make_output (); 
#endif
	  /* listkillroutines  -------------------------------------------- */ 
	     
	  lst_kill (SolutionList, free_solution_node); 

	  return;
	}
	    
      /* pop the last element -------------------------------------------- */
	  
      state = pop (Stack);	               /* current state to work with */
      
      if (LST_EMPTY (state->Intervals))   
	{
	  /* state has no intervals left: we got a solution */
	  short *ss, *ss1;

	  count++;		                          
	  structure = get_structure (state);
	  structure_energy = state->partial_energy / 100.;

	  ss = S; ss1=S1;  /* save the pointers */
	  structure_energy = energy_of_struct (sequence, structure);
	  S = ss; S1=ss1;
	  
	  if ((float) (state->partial_energy / 100.) != structure_energy)
	    {
	      printf ("%s %6.2f %6.2f\n", structure,
		      state->partial_energy / 100., structure_energy );
	      exit (1);
	    }

#ifdef NOSORT
	  printf("%s %6.2f\n", structure, structure_energy);
	  free(structure);
#else
	  new_solution = make_solution (structure_energy, structure);
	  
	  push (SolutionList, new_solution);
#endif
	}
      else
	{
	  /* get (and remove) next interval of state to analyze */

	  interval = pop (state->Intervals);  
	  
	  scan_interval (interval->i, interval->j, interval->array_flag, state);
	  
	  free_interval_node (interval);	/* free the current interval */
	}
      
      free_state_node (state);	             /* free the current state */
    } /* end of while (1) */
}
 
/*---------------------------------------------------------------------------*/
/* Definitions---------------------------------------------------------------*/

#define STACK_BULGE1      1	     /* stacking energies for bulges of size 1 */
#define MIN2(A, B)        ((A) < (B) ? (A) : (B))

/*-----------------------------------------------------------------------------*/

PUBLIC void
scan_interval (int i, int j, int array_flag, STATE * state)
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
  
  best_energy = best_attainable_energy (state);  /* ... on remaining intervals */
  nopush = true;
   
  if ((i > 1) && (!array_flag))
      nrerror ("Error while backtracking!");
  
  if (j < i + TURN + 1)	                         /* minimal structure element */
    {
      if (nopush)
	push_back (state);
      return;
    }

   /* 13131313131313131313131313131313131313131313131313131313131313131313131 */

  if (array_flag == 3 || array_flag == 1)
    {
      /* array_flag = 3:                    interval i,j was generated during */
                     /* a multiloop decomposition using array fM1 in repeat() */
                                                          /* or in this block */
      
      /* array_flag = 1:                    interval i,j was generated from a */
                                /* stack, bulge, or internal loop in repeat() */
                                                          /* or in this block */

      if (array_flag == 3)
	{
	  fi = fM1[indx[j-1] + i] + MLbase;
	}
      else
	{
	  fi = fML[indx[j-1] + i] + MLbase;
	}
      
      if (fi + best_energy <= threshold)                                      
	{
	                                   /* no basepair, nibbling of 3'-end */
	  
	  new_state = copy_state (state);
	  new_interval = make_interval (i, j-1, array_flag);
	  push (new_state->Intervals, new_interval);
	  new_state->partial_energy += MLbase;
	  push (Stack, new_state);
	}

      type = pair [S[i]][S[j]];
      element_energy = MLintern;

      if (( type) && (dangles) ) {                        /* dangling ends */
	 if (i > 1)
	    element_energy +=  dangle5[type][S1[i-1]];
	 if (j < length)
	    element_energy += dangle3[type][S1[j+1]];
      }

                                                              /* i,j may pair */
      
      cij = c[indx[j] + i] + element_energy;           
      if (cij + best_energy <= threshold)                   
	  repeat (i, j, state, element_energy, 0);

    }                                   /* array_flag == 3 || array_flag == 1 */

   /* 11111111111111111111111111111111111111111111111111111111111111111111111 */

   if (array_flag == 1) {
       /* array_flag = 1:                   interval i,j was generated from a */
                                /* stack, bulge, or internal loop in repeat() */
                                                          /* or in this block */

      for ( k = i+TURN+1 ; k <= j-1-TURN ; k++) {
	 /* Multiloop decomposition if i,j contains more than 1 stack */
	   
	 type = pair [S[k+1]][S[j]];
	 if (type==0) continue;
	 
	 element_energy = MLintern;
	 if (dangles)
	    element_energy += dangle3[type][S1[j+1]] + dangle5[type][S1[k]];
	 
	 
	 if (fML[indx[k]+i] + c[indx[j] + k+1] + element_energy + best_energy <=
	     threshold) 
	    {
	       temp_state = copy_state (state);
	       new_interval = make_interval (i, k, 1);
	       push (temp_state->Intervals, new_interval);
	       repeat (k+1, j, temp_state, element_energy, fML[indx[k]+i]);
	       free_state_node (temp_state);
	    }
      }
      
       for (k = i ; k <= j-1-TURN; k++) {
	  /* Multiloop decomposition if i,j contains only 1 stack */
	   
	  type = pair [S[k+1]][S[j]];
	  if (type==0) continue;

	  element_energy = MLintern + MLbase*(k-i+1);
	  if (dangles)
	     element_energy += dangle3[type][S1[j+1]] + dangle5[type][S1[k]];
	  
	  if (c[indx[j]+k+1] + element_energy + best_energy <= threshold)
	     repeat (k+1, j, state, element_energy, 0);
	 }
     }                                                    /* array_flag == 1 */

   /* 2222222222222222222222222222222222222222222222222222222222222222222222 */
   
   if (array_flag == 2)   	  
     {
       /* array_flag = 2:                  interval i,j was generated from a */
                               /* stack, bulge, or internal loop in repeat() */            
       
       repeat (i, j, state, 0, 0);

       if (nopush)
	  push_back (state);
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
	  
	 new_state = copy_state (state);
	 new_interval = make_interval (i, j-1 , 0);
	 push (new_state->Intervals, new_interval);
	 push (Stack, new_state);
      }
      
      for (k = j-TURN-1; k > 1; k--) {

	 type = pair[S[k]][S[j]]; if ((type==0)&&(BP[k]==j)) type=7;
	 if (type==0)
	    continue;

	                                                     /* k and j pair */
	 if (dangles) {
	    element_energy =  dangle5[type][S1[k - 1]];
	    if (j < length)
	       element_energy += dangle3[type][S1[j+1]];
	 }
	 else                                                 /* no dangles */
	    element_energy = 0;
	  
	 if (f5[k-1] + c[indx[j]+k] + element_energy + best_energy <= threshold)
	    {
	       temp_state = copy_state (state);
	       new_interval = make_interval (1, k-1, 0);
	       push (temp_state->Intervals, new_interval);
	       repeat (k, j, temp_state, element_energy, f5[k-1]);
	       free_state_node (temp_state);
	    }
      }
      
      if (dangles && (j < length)) {
	 type = pair[S[1]][S[j]]; if ((type==0)&&(BP[1]==j)) type=7;
	 element_energy = dangle3[type][S1[j+1]];
      }
      else
	 element_energy = 0;

      if (c[indx[j]+1] + element_energy + best_energy <= threshold)
	 repeat (1, j, state, element_energy, 0);
    }                                                       /* array_flag == 0 */
   
  if (nopush)
    push_back (state);
  return;
}

/*-----------------------------------------------------------------------------*/

PRIVATE void
repeat (int i, int j, STATE * state, int part_energy, int temp_energy)
{
             /* routine to find stacks, bulges, internal loops and  multiloops */
                                     /* within interval closed by basepair i,j */
  
  STATE *new_state;
  INTERVAL *new_interval;

  register int k, p, q, m, energy, new;
  register int unpaired, sizecorr, mm;
  register int no_close, no_close_2, type, type_2;
  register int n1, n2, tetracorr;
  int bonus=0, rt;
    
  no_close_2 = 0;
  
  best_energy += part_energy;          /* energy of current structural element */ 
  best_energy += temp_energy;                 /* energy from unpushed interval */
  
  type = pair[S[i]][S[j]]; if ((type==0)&&(BP[i]==j)) type=7;
  no_close = (((type == 3) || (type == 4)) && no_closingGU);


  for (p = i + 1; p <= MIN2 (j - 2 - TURN, i + MAXLOOP + 1); p++)
    {
      for (q = j - 1; q >= p + 1 + TURN; q--)
	{
	  /* stack, bulge, interior loop */
	  
	  n1 = p - i - 1;
	  n2 = j - q - 1;
	  
	  if (n1 + n2 > MAXLOOP) continue;
	  
	  if (n1 > n2) { m = n1; n1 = n2; n2 = m; }

	  type_2 = pair[S[p]][S[q]];
	  if ((BP[p]==q) && (type_2==0)) type_2=7; /* nonstandard */
	  
	  if (type_2 == 0) continue;

	  if (no_closingGU)
	    no_close_2 = (no_close || ((type_2 == 3) || (type_2 == 4)));

	  if (n2 == 0)
	    energy = stack[type][type_2];
	  
	  else if (n1 == 0) {  	     /* bulge */
	      
	     if (!no_close_2)
		energy = bulge[n2];
	     else
		energy = FORBIDDEN;
#if STACK_BULGE1
	     if (n2 == 1)
		energy += stack[type][type_2];
#endif
	  }
	  else {  	     /* interior loop */
	     if (!no_close_2) {
		if ((n1+n2==2)&&(james_rule)) {
		   /* special case for size 2 loop */
		   energy = internal2_energy;
		 } else {
		    register int rt;
		    energy = internal_loop[n1+n2];
		    
#if NEW_NINIO
		    energy += MIN2(MAX_NINIO,(abs(n1-n2)*F_ninio[2]));
#else
		    m       = MIN2(4, n1);
		    energy += MIN2(MAX_NINIO,(abs(n1-n2)*F_ninio[m]));
#endif	  
		    rt  = rtype[type_2];
		    energy += mismatchI[type][S1[i+1]][S1[j-1]]+
		       mismatchI[rt][S1[q+1]][S1[p-1]];
		 }
	     }
	     else
		energy = FORBIDDEN;
	  }

	  new = energy + c[indx[q]+p] + bonus;

	  if (new + best_energy <= threshold) {
	     /* stack, bulge, or interior loop */
	     
	     new_state = copy_state (state);
	     make_pair (i, j, new_state);
	     make_pair (p, q, new_state);
	     
	     new_interval = make_interval (p, q, 2);
	     push (new_state->Intervals, new_interval);
	     new_state->partial_energy += part_energy;
	     new_state->partial_energy += energy;
	     push (Stack, new_state);
	  }
       }                                                  /* end of q-loop */ 
   }                                                      /* end of p-loop */
  
  mm = MLclosing + MLintern;
  rt = rtype[type];
  
  for (k = i + 1 + TURN; k <= j - 1 - TURN; k++)  {
     /* multiloop decomposition */

     element_energy = mm;
     if (dangles)
	element_energy += dangle3[rt][S1[i+1]] + dangle5[rt][S1[j-1]];
      
     if (fML[indx[k] + i+1] + fM1[indx[j-1] + k+1] + element_energy + best_energy
	 <= threshold)
	{
	   new_state = copy_state (state);
	   new_interval = make_interval (i+1, k, 1);
	   push (new_state->Intervals, new_interval);
	   new_interval = make_interval (k+1, j-1, 3);
	   push (new_state->Intervals, new_interval);
	   make_pair (i, j, new_state);
	   new_state->partial_energy += part_energy;
	   new_state->partial_energy += element_energy;
	   push (Stack, new_state);
	}
  }                                                     /* end of k-loop */

  unpaired = j - i - 1;
  sizecorr = 0;
  tetracorr = 0;
  
  if (unpaired > 30) {
     unpaired = 30;
     sizecorr = (int) (lxc * log ((double) (j - i - 1) / 30.));
  }
  else if (tetra_loop)
     if (unpaired == 4) {
	char tl[5] = {0,0,0,0,0}, *ts;
	strncpy(tl, sequence+i, 4);
	if (ts=strstr(Tetraloops, tl))
	  tetracorr = TETRA_ENERGY[(ts-Tetraloops)/5];
     }
  
  mm = (unpaired > 3) ? mismatchH[type][S1[i+1]][S1[j-1]] : 0;
  
  if (no_close) {
     if (c[indx[j]+i] == FORBIDDEN)  {
	fprintf (stderr, "A BUG!\n\n");
	exit (1);
     }
  }
  else {
     element_energy = sizecorr + tetracorr + mm;
     if (hairpin[unpaired] + element_energy + best_energy <= threshold) {
	/* hairpin structure */
	
	new_state = copy_state (state);
	make_pair (i, j, new_state);
	new_state->partial_energy += part_energy;
	new_state->partial_energy += hairpin[unpaired] + element_energy;
	push (Stack, new_state);
     }
  }
  
  best_energy -= part_energy;
  best_energy -= temp_energy;
  return;
}

/*-----------------------------------------------------------------------------*/
/* Well, that is the end!------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/
