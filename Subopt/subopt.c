/*
  $Log: subopt.c,v $
  Revision 1.2  1997/08/05 00:02:44  walter
  *** empty log message ***

  Revision 1.1  1997/08/04 21:05:32  walter
  Initial revision

*/
/*
   suboptimal folding - Stefan Wuchty & Walter Fontana & Ivo Hofacker

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
#include "list.h"

PRIVATE char rcsid[] = "$Id: subopt.c,v 1.2 1997/08/05 00:02:44 walter Exp $";

/*Typedefinitions ----------------------------------------------------------- */

typedef struct state
  {
    LIST *BasePairs;
    LIST *Intervals;
    int partial_energy;
  }
STATE;

typedef struct _pair
  {
    int i;
    int j;
  }
PAIR;

typedef struct interval
  {
    int i;
    int j;
    int array_flag;
  }
INTERVAL;


PUBLIC PAIR *make_pair (int i, int j);
PUBLIC INTERVAL *make_interval (int i, int j, int ml);
PUBLIC STATE *make_state (LIST * Intervals, LIST * BasePairs, int partial_energy);
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
PUBLIC void push_back (STATE * state);
PUBLIC void get_structure (STATE * state, char *solution);
PRIVATE void repeat (int i, int j, STATE * state, int part_energy, int temp_energy);
PRIVATE int check_ML (int i, int j, int k);

/*Globals -------------------------------------------------------------------- */

#define MAXALPHA 20		                 /* maximal length of alphabet */

PRIVATE LIST *Stack;
PRIVATE int nopush;
PRIVATE int best_energy;		     /* best_energy = remaining energy */

extern int pair[MAXALPHA + 1][MAXALPHA + 1];
extern int MLbase;
extern int length;
extern int delta; 
extern int *f3;                                             /* energy of 3 end */
extern int *f5;                                             /* energy of 5 end */
extern int *c;		                  /* energy array, given that i-j pair */
extern int *fML;		          /* multi-loop auxiliary energy array */
extern int *fM1;                  /* another multi-loop auxiliary energy array */
extern int *indx;     /* index for moving in the triangle matrices c[] and f[] */
extern short *S, *S1;

extern float energy_of_struct (char *, char *);

int minimal_energy;                                     /* minimum free energy */
int element_energy;                /* internal energy of a structural  element */


/*-----------------------------------------------------------------------------*/
/*List routines----------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/

PUBLIC PAIR *
make_pair (int i, int j)
{
  PAIR *pair;

  pair = lst_newnode (sizeof (PAIR));
  pair->i = i;
  pair->j = j;
  return pair;
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
free_state_node (STATE * node)
{
  if (node->BasePairs)
    lst_kill (node->BasePairs, lst_freenode);
  if (node->Intervals)
    lst_kill (node->Intervals, lst_freenode);
  lst_freenode (node);
}

/*-----------------------------------------------------------------------------*/

PUBLIC STATE *
make_state (LIST * Intervals, LIST * BasePairs, int partial_energy)
{
  STATE *state;

  state = lst_newnode (sizeof (STATE));

  if (Intervals)
    state->Intervals = Intervals;
  else
    state->Intervals = lst_init ();

  if (BasePairs)
    state->BasePairs = BasePairs;
  else
    state->BasePairs = lst_init ();

  state->partial_energy = partial_energy;

  return state;
}

/*-----------------------------------------------------------------------------*/

PUBLIC STATE *
copy_state (STATE * state)
{
  STATE *new_state;
  void *after;
  PAIR *new_basepair, *rec;
  INTERVAL *new_interval, *next;
  
  new_state = lst_newnode (sizeof (STATE));
  new_state->Intervals = lst_init ();
  new_state->BasePairs = lst_init ();
  new_state->partial_energy = state->partial_energy;

  if (state->Intervals->count)
    {
      after = LST_HEAD (new_state->Intervals);
      for ( next = lst_first (state->Intervals); next; next = lst_next (next))
	{
	  new_interval = lst_newnode (sizeof (INTERVAL));
	  new_interval->i = next->i;
	  new_interval->j = next->j;
	  new_interval->array_flag = next->array_flag;
	  lst_insertafter (new_state->Intervals, new_interval, after);
	  after = new_interval;
	}
    }
  if (state->BasePairs->count)
    {
      after = LST_HEAD (new_state->BasePairs);
      for (rec = lst_first (state->BasePairs); rec; rec = lst_next (rec))
	{
	  new_basepair = lst_newnode (sizeof (PAIR));
	  new_basepair->i = rec->i;
	  new_basepair->j = rec->j;
	  lst_insertafter (new_state->BasePairs, new_basepair, after);
	  after = new_basepair;
	}
    }
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
      fflush (stdout);
    }
  if (state->BasePairs->count)
    {
      printf ("%d base pairs:\n", state->BasePairs->count);
      for (rec = lst_first (state->BasePairs); rec; rec = lst_next (rec))
	{
	  printf ("(%d,%d) ", rec->i, rec->j);
	}
      printf ("\n");
      fflush (stdout);
    }
  printf ("\n");
  fflush (stdout);
  printf (" partial_energy: %d\n", state->partial_energy);
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

  sum = state->partial_energy; /* energy of already found elements of structure */

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

PUBLIC void
get_structure (STATE * state, char *solution)
{
  int x;
  PAIR *rec;

  for (x = 0; x < length; x++)
    solution[x] = '.';
  for (rec = lst_first (state->BasePairs); rec; rec = lst_next (rec))
    {
      solution[rec->i - 1] = '(';
      solution[rec->j - 1] = ')';
    }
  
  return;
}

/*-----------------------------------------------------------------------------*/
/* start of subopt backtracking- ----------------------------------------------*/
/*-----------------------------------------------------------------------------*/

PUBLIC void
subopt (void)
{
  extern char *sequence;
  extern int PRINT;
  
  STATE *state;
  LIST *Intervals, *BasePairs;
  INTERVAL *interval;
  PAIR *rec;
  int partial_energy;
  char *solution, *empty;
  float structure_energy;
  register int l, maxlevel, count;

  /* Initialize -------------------------------------------------------------- */
  
  solution = (char *) space (sizeof (char) * (length + 1));
  empty = (char *) space (sizeof (char) * (length + 1));
  maxlevel = 0;
  count = 0;
  partial_energy = 0;

  for (l = 0; l < length; l++)
    empty[l] = '.';

  /* Initialize the stack ---------------------------------------------------- */

  minimal_energy = f5[length]; 
  
  Stack = make_list ();		                                     /* anchor */
  Intervals = make_list ();	                             /* initial state: */
  interval = make_interval (1, length, 0);            /* interval [1,length,0] */
  push (Intervals, interval);
  state = make_state (Intervals, NULL, partial_energy);
  push (Stack, state);

  /* end initialize ---------------------------------------------------------- */

  while (1)			      /* forever, til nothing remains on stack */
    {
      maxlevel = (Stack->count > maxlevel ? Stack->count : maxlevel);
      
      if (LST_EMPTY (Stack))	             /* we are done! clean up and quit */
	{
	  lst_kill (Stack, free_state_node);
	  free (solution);
	  free (empty);

	  printf ("\n\n%d solutions \n", count);
	  printf ("(max stack level was %d)\n", maxlevel);
	  fflush (stdout);
	  return;
	}

  /* pop the last element ---------------------------------------------------- */
      
      state = pop (Stack);	                 /* current state to work with */
      
      if (LST_EMPTY (state->Intervals))   
	{
	  /* state has no intervals left: we got a solution */
	  
	  count++;		                          
	  get_structure (state, solution);
	  structure_energy = 0.;
	  
	  structure_energy = energy_of_struct (sequence, solution);

	  if (PRINT)
	    {
	      printf ("%d %s %6.2f %6.2f\n",
		      count,
		      solution, (float) (state->partial_energy / 100.),
		      structure_energy);
	      fflush (stdout);
	    }
	  else /* print a dot for every 1000 solutions found */
	    {
	      if (count % 1000 == 0)
		{
		  printf (".");
		  fflush (stdout);
		}
	    }
	  
	  if ((float) (state->partial_energy / 100.) != structure_energy)
	    {
	      printf ("                A BUG ! ! ! ! ! !\n");
	      printf ("partial_energy and energy of structure are not equal!!!\n");
	      printf ("         Calculation terminated!\n");
	      exit (1);
	    } 
	}
      else
	{
	  /* get (and remove) next interval of state to analyze */
	  
	  interval = pop (state->Intervals);  
	  
	  scan_interval (interval->i, interval->j, interval->array_flag, state);
	         	  
	  free_interval_node (interval);	  /* free the current interval */
	}
      
      free_state_node (state);	                     /* free the current state */
    }                                                      /* end of while (1) */
}

/*-----------------------------------------------------------------------------*/
/* Definitions-----------------------------------------------------------------*/

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

  extern int MLintern[NBPAIRS + 1];

  STATE *new_state, *temp_state;
  PAIR *new_basepair;
  INTERVAL *new_interval;
  register int k, fi, cij;
  register int n1, n2, tetracorr;
  register int energy1, energy2, energy3;
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

   /* 1313131313131313131313131313131313131313131313131313131313131313131313 */

  if (array_flag == 3 || array_flag == 1)
    {
      /* array_flag = 3: interval i,j was generated during */
      /* a multiloop decomposition using array fM1 in repeat() */
      /* or in this block */
      
      /* array_flag = 1: interval i,j was generated from a */
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
      
      if (fi + best_energy <= minimal_energy + delta)                                      
	{
	  /* no basepair, nibbling of 3'-end */
	  
	  new_state = copy_state (state);
	  new_interval = make_interval (i, j-1, array_flag);
	  push (new_state->Intervals, new_interval);
	  new_state->partial_energy += MLbase;
	  push (Stack, new_state);
	}

      type = pair [S[i]][S[j]];
      element_energy = MLintern[type];

      if (type)
	{
	  if (dangles)  /* dangling ends */
	    {
	      if ((i > 1) && (j < length))
		{
		  element_energy +=  dangle5[type][S1[i-1]] + dangle3[type][S1[j+1]];
		}
	      else
		{
		  if (i > 1)
		    {
		      element_energy +=  dangle5[type][S1[i-1]];
		    }
		  else if (j < length)
		    {
		      element_energy += dangle3[type][S1[j+1]];
		    }
		}
	    }
	}

      /* i,j may pair */
      
      cij = c[indx[j] + i] + element_energy;           
      if (cij + best_energy <= minimal_energy + delta)                   
	{
	  repeat (i, j, state, element_energy, 0);
	}
    }  /* array_flag == 3 || array_flag == 1 */

   /* 1111111111111111111111111111111111111111111111111111111111111111111111 */

   if (array_flag == 1)
     {
       /* array_flag = 1: interval i,j was generated from a */
       /* stack, bulge, or internal loop in repeat() */
       /* or in this block */

       for ( k = i+TURN+1 ; k <= j-1-TURN ; k++)    
	 {
	   /* Multiloop decomposition if i,j contains more than 1 stack */
	   
	   type = pair [S[k+1]][S[j]];
	   if (type==0) continue;

	   if (dangles)
	     {
	       element_energy = MLintern[type] + dangle3[type][S1[j+1]] + dangle5[type][S1[k]];
	     }
	   else
	     {
	       element_energy = MLintern[type];
	     }
	   
	   if (fML[indx[k]+i] + c[indx[j] + k+1] + element_energy + best_energy <=
	       minimal_energy + delta)
	     {
	       temp_state = copy_state (state);
	       new_interval = make_interval (i, k, 1);
	       push (temp_state->Intervals, new_interval);
	       repeat (k+1, j, temp_state, element_energy, fML[indx[k]+i]);
	       free_state_node (temp_state);
	     }
	 }
      
       for (k = i ; k <= j-1-TURN; k++)       
	 {
	   /* Multiloop decomposition if i,j contains only 1 stack */
	   
	   type = pair [S[k+1]][S[j]];
	   if (type==0) continue;
	  
	   if (dangles)
	     {
	       element_energy = MLintern[1] + MLbase*(k-i+1) + dangle3[type][S1[j+1]] +
		 dangle5[type][S1[k]];
	     }
	   else
	     {
	       element_energy = MLintern[1] + MLbase*(k-i+1);
	     }
	  
	   if (c[indx[j]+k+1] + element_energy + best_energy <= minimal_energy + delta)
	     {
	       repeat (k+1, j, state, element_energy, 0);
	     }
	 }
     } /* array_flag == 1 */

   /* 2222222222222222222222222222222222222222222222222222222222222222222222 */
   
   if (array_flag == 2)   	  
     {
       /* array_flag = 2: interval i,j was generated from a */
       /* stack, bulge, or internal loop in repeat() */            
       
       repeat (i, j, state, 0, 0);

       if (nopush)
	 push_back (state);
       return;
     }
  
   /* 0000000000000000000000000000000000000000000000000000000000000000000000 */
   
   if (array_flag == 0)
    {
      /* array_flag = 0: interval i,j was found while */
      /* tracing back through f5-array and c-array */
      /* or within this block */
      
      if (f5[j-1] + best_energy <= minimal_energy + delta)
        {
	  /* no basepair, nibbling of 3'-end */
	  
	  new_state = copy_state (state);
	  new_interval = make_interval (i, j-1 , 0);
	  push (new_state->Intervals, new_interval);
	  push (Stack, new_state);
	}
      
      for (k = j-TURN-1; k > 1; k--)
	{
	  type = pair[S[k]][S[j]];
	  if (type==0)
	      continue;

	  /* k and j pair */
	  
	  if (dangles)
	    {                   
	      if (j < length)
		{
		  element_energy = dangle3[type][S1[j+1]] + dangle5[type][S1[k-1]];
		}
	      else  /* j == length */
		{
		  element_energy =  dangle5[type][S1[k - 1]];
		}
	    }
	  else  /* no dangles */
	    {
	      element_energy = 0;
	    }
	  
	  if (f5[k-1] + c[indx[j]+k] + element_energy + best_energy <= minimal_energy + delta)
	    {
	      temp_state = copy_state (state);
	      new_interval = make_interval (1, k-1, 0);
	      push (temp_state->Intervals, new_interval);
	      repeat (k, j, temp_state, element_energy, f5[k-1]);
	      free_state_node (temp_state);
	    }
	}
      
      if (dangles && (j < length))
	{
	  type = pair[S[1]][S[j]];
	  element_energy = dangle3[type][S1[j+1]];
	}
      else
	{
	  element_energy = 0;
	}

      if (c[indx[j]+1] + element_energy + best_energy <= minimal_energy + delta)
	{
	  repeat (1, j, state, element_energy, 0);
	}
    } /* array_flag == 0 */
   
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
  
  extern int stack[NBPAIRS + 1][NBPAIRS + 1];
  extern int hairpin[31];
  extern int bulge[MAXLOOP + 1];
  extern int internal_loop[MAXLOOP + 1];
  extern int F_ninio[5];
  extern double lxc;
  extern int MLclosing[NBPAIRS + 1];
  extern int TETRA_ENERGY;
  extern char *sequence;

  STATE *new_state;
  PAIR *new_basepair;
  INTERVAL *new_interval;

  register int k, p, q, m, energy, new;
  register int unpaired, sizecorr, mm;
  register int no_close, no_close_2, type, type_2;
  register int n1, n2, tetracorr;
    
  no_close_2 = 0;
  
  best_energy += part_energy;             /* energy of current structural element */ 
  best_energy += temp_energy;             /* energy from unpushed interval */
  
  type = pair[S[i]][S[j]];
  no_close = (((type == 3) || (type == 4)) && no_closingGU);

  unpaired = j - i - 1;
  sizecorr = 0;
  tetracorr = 0;
  
  if (unpaired > 30)
    {
      unpaired = 30;
      sizecorr = (int) (lxc * log ((double) (j - i - 1) / 30.));
    }
  else if (tetra_loop)
    {
      if (unpaired == 4)
	{
	  for (k = 0; k < N_TETRALOOPS; k++)
	    {
	      if (strncmp (sequence+i, Tetraloops[k], 4) == 0)
		tetracorr = TETRA_ENERGY;
	    }
	}
    }

  /* these are the wrong mismatches !! */
  mm = (unpaired > 3) ? mismatch[S1[i]][S1[j]][S1[i + 1]][S1[j - 1]] : 0;

  if (no_close)
    {
      if (c[indx[j]+i] == FORBIDDEN)
	{
	  printf ("A BUG!\n\n");
	  fflush (stdout);
	  exit (1);
	}
    }
  else
    {
      element_energy = sizecorr + tetracorr + mm;
      if (hairpin[unpaired] + element_energy + best_energy <= minimal_energy + delta)
	{
	  /* hairpin structure */
	  
	  new_state = copy_state (state);
	  new_basepair = make_pair (i, j);
	  push (new_state->BasePairs, new_basepair);
	  new_state->partial_energy += part_energy;
	  new_state->partial_energy += hairpin[unpaired] + element_energy;
	  push (Stack, new_state);
	}
    }

  for (p = i + 1; p <= MIN2 (j - 2 - TURN, i + MAXLOOP + 1); p++)
    {
      for (q = j - 1; q >= p + 1 + TURN; q--)
	{
	  /* stack, bulge, interior loop */
	  
	  n1 = p - i - 1;
	  n2 = j - q - 1;
	  
	  if (n1 + n2 > MAXLOOP)
	    continue;
	  
	  if (n1 > n2)
	    {
	      m = n1;
	      n1 = n2;
	      n2 = m;
	    }

	  type_2 = pair[S[p]][S[q]];
	  if (type_2 == 0)
	    continue;

	  if (no_closingGU)
	    no_close_2 = (no_close || ((type_2 == 3) || (type_2 == 4)));

	  if (n2 == 0)
	    energy = stack[type][type_2];
	  
	  else if (n1 == 0)
	    {
	      /* bulge */
	      
	      if (!no_close_2)
		energy = bulge[n2];
	      else
		energy = FORBIDDEN;
#if STACK_BULGE1
	      if (n2 == 1)
		energy += stack[type][type_2];
#endif
	    }
	  else
	    {
	      /* interior loop */
	      
	      if (!no_close_2)
		{
		  energy = internal_loop[n1+n2];

		  m = MIN2 (4, n1);
		  energy += MIN2 (MAX_NINIO, (abs (n1-n2) * F_ninio[m]));
		  /* wrong mismatches !! */
		  energy += mismatch[S1[i]][S1[j]][S1[i+1]][S1[j-1]] +
		    mismatch[S1[p-1]][S1[q+1]][S1[p]][S1[q]];
		}
	      else
		energy = FORBIDDEN;
	    }

	  new = energy + c[indx[q]+p];

	  if (new + best_energy <= minimal_energy + delta)
	    {
	      /* stack, bulge, or interior loop */
	      
	      new_state = copy_state (state);
	      new_basepair = make_pair (i, j);
	      push (new_state->BasePairs, new_basepair);
	      new_basepair = make_pair (p, q);
	      push (new_state->BasePairs, new_basepair);
	      new_interval = make_interval (p, q, 2);
	      push (new_state->Intervals, new_interval);
	      new_state->partial_energy += part_energy;
	      new_state->partial_energy += energy;
	      push (Stack, new_state);
	    }
	} /* end of q-loop */ 
    } /* end of p-loop */
 
  mm = MLclosing[type];
  type = pair[S[j]][S[i]];

  for (k = i + 1 + TURN; k <= j - 1 - TURN; k++)
    {
      /* multiloop decomposition */
      
      if (dangles)
	{
	  element_energy = dangle3[type][S1[i+1]] + dangle5[type][S1[j-1]] + mm;
	}
      else
	{
	  element_energy = mm;
	}
      
      if (fML[indx[k] + i+1] + fM1[indx[j-1] + k+1] + element_energy + best_energy <=
	  minimal_energy + delta)
	{
	  new_state = copy_state (state);
	  new_interval = make_interval (i+1, k, 1);
	  push (new_state->Intervals, new_interval);
	  new_interval = make_interval (k+1, j-1, 3);
	  push (new_state->Intervals, new_interval);
	  new_basepair = make_pair (i, j);
	  push (new_state->BasePairs, new_basepair);
	  new_state->partial_energy += part_energy;
	  new_state->partial_energy += element_energy;
	  push (Stack, new_state);
	}
    }                                                         /* end of k-loop */
  
  best_energy -= part_energy;
  best_energy -= temp_energy;
  return;
}

/*-----------------------------------------------------------------------------*/
/* Well, that is the end!------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/
