/*
  $Log: subopt.c,v $
  Revision 1.1  1997/08/04 21:05:32  walter
  Initial revision

*/
/*
   suboptimal folding - Stefan Wuchty & Walter Fontana

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

#define PUBLIC
#define PRIVATE static

PRIVATE char rcsid[] = "$Id: subopt.c,v 1.1 1997/08/04 21:05:32 walter Exp $";

/*Typedefinitions ----------------------------------------------------------- */

typedef struct state
  {
    LIST *BasePairs;
    LIST *Intervals;
    int cur_en;
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
    int ml;
  }
INTERVAL;


PUBLIC PAIR *make_pair (int i, int j);
PUBLIC INTERVAL *make_interval (int i, int j, int ml);
PUBLIC STATE *make_state (LIST * Intervals, LIST * BasePairs, int cur_en);
PUBLIC STATE *copy_state (STATE * state);
PUBLIC void print_state (STATE * state);
PUBLIC void print_stack (LIST * list);
PUBLIC LIST *make_list (void);
PUBLIC void push (LIST * list, void *data);
PUBLIC void *pop (LIST * list);
PUBLIC void subopt (void);
PUBLIC int remaining_energy (STATE * state);
PUBLIC void scan_interval (int i, int j, int ml, STATE * state);
PUBLIC void free_pair_node (PAIR * node);
PUBLIC void free_interval_node (INTERVAL * node);
PUBLIC void free_state_node (STATE * node);
PUBLIC void no_push (int nopush, STATE * state);
PUBLIC char *print_structure (STATE * state, char *solution);
PRIVATE void repeat (int i, int j, STATE * state, int part_energy, int temp_energy);
PRIVATE int check_ML (int i, int j, int k);

/*Globals -------------------------------------------------------------------- */

#define MAXALPHA 20		                 /* maximal length of alphabet */

PRIVATE LIST *Stack;
PRIVATE int nopush;
PRIVATE int rem_en;		                  /* rem_en = remaining energy */
PRIVATE int cij, ci1j, cij1, ci1j1;


extern int pair[MAXALPHA + 1][MAXALPHA + 1];
extern int MLbase;
extern int length;
extern int *BP;
extern int delta; 
extern int *f3;                                             /* energy of 3 end */
extern int *f5;                                             /* energy of 5 end */
extern int *c;		                  /* energy array, given that i-j pair */
extern int *fML;		          /* multi-loop auxiliary energy array */
extern int *fM1;                  /* another multi-loop auxiliary energy array */
                                             /* used in double_dangle_subopt.c */
extern int *indx;     /* index for moving in the triangle matrices c[] and f[] */
extern short *S, *S1;
extern float energy_of_struct (char *, char *);
int min_en;                                    /* min_en = minimal free energy */

int en_intern;          /* internal energy of an element of structure */


/*-----------------------------------------------------------------------------*/
/*Subroutines------------------------------------------------------------------*/
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
make_interval (int i, int j, int ml)
{
  INTERVAL *interval;

  interval = lst_newnode (sizeof (INTERVAL));
  interval->i = i;
  interval->j = j;
  interval->ml = ml;
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
make_state (LIST * Intervals, LIST * BasePairs, int cur_en)
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

  state->cur_en = cur_en;

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
  new_state->cur_en = state->cur_en;

  if (state->Intervals->count)
    {
      after = LST_HEAD (new_state->Intervals);
      for ( next = lst_first (state->Intervals); next; next = lst_next (next))
	{
	  new_interval = lst_newnode (sizeof (INTERVAL));
	  new_interval->i = next->i;
	  new_interval->j = next->j;
	  new_interval->ml = next->ml;
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
	  printf ("[%d,%d],%d ", next->i, next->j, next->ml);
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
  printf (" cur_en: %d\n", state->cur_en);
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

PUBLIC int
remaining_energy (STATE * state)
{			    /* evaluation of energy within remaining intervals */
  register int sum;
  INTERVAL *next;

  sum = state->cur_en;        /* energy of already found elements of structure */

  for (next = lst_first (state->Intervals); next; next = lst_next (next))
    {
      if (next->ml == 0)
	  sum += f5[next->j];
      else if (next->ml == 1)
	  sum += fML[indx[next->j] + next->i];
      else if (next->ml == 2)
	  sum += c[indx[next->j] + next->i];
    }
    
  return sum;
}

/*-----------------------------------------------------------------------------*/
/* start of new backtrackroutine ----------------------------------------------*/
/*-----------------------------------------------------------------------------*/

PUBLIC void
subopt (void)
{
  STATE *state;
  LIST *Intervals, *BasePairs;
  INTERVAL *interval;
  PAIR *rec;
  
  int cur_en = 0;
  extern char sequence[1000];
  char *solution, *empty;

  register int x, maxlevel, count;
  register float struc_en;
  double time;
    
  time = cpu_time (0.);

  solution = (char *) space (sizeof (char) * (length + 1));
  empty = (char *) space (sizeof (char) * (length + 1));

  maxlevel = 0;
  count = 0;

  for (x = 0; x < length; x++)
    empty[x] = '.';


  /* Initialize the stack ---------------------------------------------------- */

  ci1j = cij1 = ci1j1 = INF;
  min_en = f5[length]; 
  
  Stack = make_list ();		                                     /* anchor */
  Intervals = make_list ();	                             /* initial state: */
  interval = make_interval (1, length, 0);           /* interval [1,length,ml] */
                                           /* ml = 0:  trace back in f5-array  */
                                           /* ml = 1:  trace back in fML-array */
                                           /* ml = 2:  trace back in repeat()  */
                                           /* ml = 3:  trace back in fM1-array */
  push (Intervals, interval);
  state = make_state (Intervals, NULL, cur_en);
  push (Stack, state);

  /* end initialize ---------------------------------------------------------- */

  while (1)			      /* forever, til nothing remains on stack */
    {
      maxlevel = (Stack->count > maxlevel ? Stack->count : maxlevel);
      
      if (LST_EMPTY (Stack))	                               /* we are done! */
	{
	                                                  /* clean up and quit */
	  lst_kill (Stack, free_state_node);
	  free (solution);
	  free (empty);

	  printf ("%d solutions \n", count);
	  printf ("(max stack level was %d)\n", maxlevel);
	  printf ("%.1f seconds cpu time\n\n", cpu_time (time));
	  fflush (stdout);
	  return;
	}

  /* pop the last element ---------------------------------------------------- */

      state = pop (Stack);	                 /* current state to work with */
      
      if (LST_EMPTY (state->Intervals))                   /* we got a solution */
	{
	  count++;		                          
	  solution = print_structure (state, solution);
	  struc_en = 0.;
	  
	  struc_en = energy_of_struct(sequence, solution);
	  printf ("%d %s %6.2f %6.2f\n", count, solution, (float) (state->cur_en / 100.),
		  struc_en);
	  
	  if ((float) (state->cur_en/ 100.) != struc_en)
	    {
	      printf ("      M A L F U N C T I O N ! ! ! ! ! !\n");
	      printf ("cur_en and energy of structure are not equal!!!\n");
	      printf ("         Calculation terminated!\n");
	      printf ("(max stack level was %d)\n", maxlevel);
	      break;
	    } 
	  
	  fflush (stdout);
	}
      else
	{
	  interval = pop (state->Intervals);  /* current interval to work with */
	  
	  scan_interval (interval->i, interval->j, interval->ml, state);
	         	  
	  free_interval_node (interval);	  /* free the current interval */
	}
      free_state_node (state);	                     /* free the current state */
      
    }                                                      /* end of while (1) */
}

/* Definitions-----------------------------------------------------------------*/

#define STACK_BULGE1      1	     /* stacking energies for bulges of size 1 */
#define MAXSECTORS      500	            /* dimension for a backtrack array */
#define LOCALITY         0.	          /* locality parameter for base-pairs */
#define MIN2(A, B)      ((A) < (B) ? (A) : (B))

/*-----------------------------------------------------------------------------*/

PUBLIC void
scan_interval (int i, int j, int ml, STATE * state)  /* real backtrackroutine: */
                                         /* trace back within i,j with flag ml */
{
  STATE *new_state, *temp_state;
  PAIR *new_basepair;
  INTERVAL *new_interval;

  extern int MLintern[NBPAIRS + 1];
  
  register int k, fi;
  register int n1, n2, tetracorr;
  register int energy1, energy2, energy3;
  register int type;
  
  /* start the routine ------------------------------------------------------- */

  rem_en = remaining_energy (state);     /* calculates best possible energy of */
                                           /* remaining intervals of the state */
  nopush = 1;
   
  if ((i > 1) && (!ml))
      nrerror ("Error while backtracking!");
  
  if (((j - i + 1) == length) && (backtrack_type == 'C'))
    {
      repeat (i, j, state, 0., 0.);
    }

  if (j < i + TURN + 1)	         /* minimal demand for an element of structure */
    {
      no_push (nopush, state);
      return;
    }
  
  
   if (ml == 3)                       /* ml = 3: interval i,j was found during */
                       /* a multiloopdecomposition using fM1-array in repeat() */
                                     /* or within the following interrogations */
    {
      fi = fM1[indx[j - 1] + i] + MLbase;
      if (fi + rem_en <= min_en + delta)                                      
	{                                    /* no basepair, nibbling of 3-end */
	  new_state = copy_state (state);
	  new_interval = make_interval (i, j - 1, 3);
	  push (new_state->Intervals, new_interval);
	  new_state->cur_en += MLbase;
	  push (Stack, new_state);
	  nopush = 0;
	}

      type = pair [S[i]][S[j]];
      if (( type == 0 ) && (BP[i] == j))
	  type = 7;
      
      en_intern = MLintern[type];

      if (type)
	  if (dangles)                           /* dangling end contributions */
	    {
	      if ((i>1)&&(j<length))
		  en_intern +=  dangle5[type][S1[i-1]] + dangle3[type][S1[j+1]];
	      else                                    /* just for completeness */
		  if (i>1)
		      en_intern +=  dangle5[type][S1[i-1]];
		  else if (j<length)
		      en_intern += dangle3[type][S1[j+1]];
	    }

      cij = c[indx[j] + i] + en_intern;           
            
      if (cij + rem_en <= min_en + delta)                   
	{                                                /* basepair i,j found */      
	  repeat (i, j, state, en_intern, 0.);
	}
    }                                                   /* end of if (ml == 3) */


  if (ml == 2)                        /* ml = 2: interval i,j was found during */
                    /* stack-, bulge, -internal loop decomposition in repeat() */            
	  
    {
      repeat (i, j, state, 0., 0.);
      no_push (nopush, state);
      return;
    }
  
  if (ml == 1)                        /* ml = 1: interval i,j was found during */
                    /* stack-, bulge, -internal loop decomposition in repeat() */
                                     /* or within the following interrogations */
    {
      fi = fML[indx[j - 1] + i] + MLbase;
      if (fi + rem_en <= min_en + delta)     /* no basepair, nibbling of 3-end */
	{
	  new_state = copy_state (state);
	  new_interval = make_interval (i, j - 1, 1);
	  push (new_state->Intervals, new_interval);
	  new_state->cur_en += MLbase;
	  push (Stack, new_state);
	  nopush = 0;
	}

      type = pair [S[i]][S[j]];
      if (( type == 0 ) && (BP[i] == j))
	  type = 7;
      
      en_intern = MLintern[type];

      if (type)
	  if (dangles)                           /* dangling end contributions */
	    {
	      if ((i>1)&&(j<length))
		  en_intern +=  dangle5[type][S1[i-1]] + dangle3[type][S1[j+1]];
	      else                                     /* just for completness */
		  if (i>1)
		      en_intern +=  dangle5[type][S1[i-1]];
		  else if (j<length)
		      en_intern += dangle3[type][S1[j+1]];
	    }

      cij = c[indx[j] + i] + en_intern;                  
            
      if (cij + rem_en <= min_en + delta)                /* basepair i,j found */              
	{                                                
	  repeat (i, j, state, en_intern, 0.);
	}
     
      for ( k = i+TURN+1 ; k <= j - 1 - TURN ; k++)    
	{
	  type = pair [S[k+1]][S[j]];
	  if (( type == 0 ) && (BP[k+1] == j))
	      type = 7;
	  if (type==0) continue;

	                             /* Multiloopdecomposition if i,j contains */
	                                                  /* more than 1 stack */
	  if (dangles)
	    {
	      en_intern = MLintern[type] + dangle3[type][S1[j+1]] + dangle5[type][S1[k]];
	      if (fML[indx[k] + i] + c[indx[j] + k + 1] + en_intern + rem_en <= min_en + delta)
		{
		  temp_state = copy_state (state);
		  new_interval = make_interval (i, k, 1);
		  push (temp_state->Intervals, new_interval);
		  repeat (k + 1, j, temp_state, en_intern, fML[indx[k] + i]);
		  free_state_node (temp_state);
		}
	    }
	  else
	      if (fML[indx[k] + i] + c[indx[j] + k + 1] + MLintern[type]
		  + rem_en <= min_en + delta)
		{                                   
		  temp_state = copy_state (state);
		  new_interval = make_interval (i, k, 1);
		  push (temp_state->Intervals, new_interval);
		  repeat (k + 1, j, temp_state, MLintern[1], fML[indx[k] + i]);
		  free_state_node (temp_state);
		}
	  
	}
      
      for ( k = i ; k <= j - 1 - TURN; k++)       
	                              
	{
	  type = pair [S[k+1]][S[j]];
	  if (( type == 0 ) && (BP[k+1] == j))
	      type = 7;
	  if (type==0) continue;
                                     /* Multiloopdecomposition if i,j contains */
	                                                       /* only 1 stack */
	  
	  if (dangles)
	    {
	      en_intern = MLintern[1] + MLbase*(k-i+1) + dangle3[type][S1[j+1]] +
		  dangle5[type][S1[k]];
	      if (c[indx[j] + k + 1] + en_intern + rem_en <= min_en + delta)
		{
		  repeat (k+1, j, state, en_intern, 0.);
		}
	    }
	  else
	    {
	      en_intern = MLintern[1] + MLbase*(k-i+1);
	      if (c[indx[j] + k + 1] + en_intern + rem_en <= min_en + delta)
		{                                         /*  cik + k-i+1 unpaired bases */
		  repeat (k+1, j, state, en_intern, 0.); 
		}
	    }
	  
	}                                                     
      
    }                                                   /* end of if (ml == 1) */
  
  
  if (ml == 0)                        /* ml = 0: interval i,j was found during */
                                  /* tracing back through f5-array and c-array */
                                        /* within the following interrogations */
    {
      if (f5[j-1] + rem_en <= min_en + delta)/* no basepair, nibbling of 3-end */     
        {
	  new_state = copy_state (state);
	  new_interval = make_interval (i, j - 1 , 0);
	  push (new_state->Intervals, new_interval);
	  push (Stack, new_state);
	  nopush = 0;
	}
      
      for (k=j-TURN-1; k>1; k--)
	{
	  type = pair[S[k]][S[j]];
	  if ((type==0)&&(BP[k]==j))
	      type=7;
	  if (type==0)
	      continue;
	                            
	  if (dangles)                                   /* basepair k,j found */
	    {                   
	      if (j<length)
		{
		  en_intern = dangle3[type][S1[j + 1]] + dangle5[type][S1[k - 1]];
		  if (f5[k - 1] + c[indx[j] + k] + en_intern + rem_en <= min_en + delta)
		    {
		      temp_state = copy_state (state);
		      new_interval = make_interval (1, k - 1, 0);
		      push (temp_state->Intervals, new_interval);
		      repeat (k, j, temp_state, en_intern, f5[k - 1]);
		      free_state_node (temp_state);
		    }
		}
	      else
		  {
		    en_intern =  dangle5[type][S1[k - 1]];        /* j==length */
		    if (f5[k - 1]+c[indx[j]+k] + en_intern + rem_en <= min_en +
			delta)
		    {
		      temp_state = copy_state (state);
		      new_interval = make_interval (1, k - 1, 0);
		      push (temp_state->Intervals, new_interval);
		      repeat (k, j, temp_state, en_intern, f5[k - 1] );
		      free_state_node (temp_state);
		    }
		  }
	      
	    }
	  else                                                   /* no dangles */
	    {
	      if (f5[k - 1] + c[indx[j] + k] + rem_en <= min_en + delta)
		{
		  temp_state = copy_state (state);
		  new_interval = make_interval (1, k - 1, 0);
		  push (temp_state->Intervals, new_interval);
		  repeat (k, j, temp_state, 0., f5[k - 1]);
		  free_state_node (temp_state);
		}
	    }
	}
      
      if ((dangles)&&(j<length))
	{
	  type = pair[S[1]][S[j]];
	  if ((type==0)&&(BP[1]==j))
	      type=7;
	  
	  en_intern = dangle3[type][S1[j+1]];
	  if (c[indx[j]+1] + en_intern <= min_en + delta)
	    {
	      repeat (1, j, state, en_intern, 0.);
	    }
	}
      else
	  if (c[indx[j] + 1] + rem_en <= min_en + delta)
	    {
	      repeat (1, j, state, 0., 0.);
	    }
    }                                                   /* end of if (ml == 0) */
  no_push (nopush, state);
  return;
}                                                   /* end of scan_interval () */

/*-----------------------------------------------------------------------------*/

PRIVATE void
repeat (int i, int j, STATE * state, int part_energy, int temp_energy)
             /* routine to find stacks, bulges, internal loops and  multiloops */
                                     /* within interval closed by basepair i,j */
{
  STATE *new_state;
  PAIR *new_basepair;
  INTERVAL *new_interval;

  extern int stack[NBPAIRS + 1][NBPAIRS + 1];
  extern int hairpin[31];
  extern int bulge[MAXLOOP + 1];
  extern int internal_loop[MAXLOOP + 1];
  extern int F_ninio[5];
  extern double lxc;
  extern int MLclosing[NBPAIRS + 1];
  extern int TETRA_ENERGY;
  extern char sequence[1000];

  register int k, p, q, m, energy, new;
  register int unpaired, sizecorr, mm;
  register int no_close, no_close_2, type, type_2;
  register int n1, n2, tetracorr;
  register int bonus = 0;
    
  no_close_2 = 0;
  
  rem_en += part_energy;                        /* energy basepair i,j carries */ 
  rem_en += temp_energy; /* best possible energycontribution unpushed interval */
  
  type = pair[S[i]][S[j]];
  if ((BP[i] == j) && (type == 0))
    type = 7;			                                /* nonstandard */

  no_close = (((type == 3) || (type == 4)) && no_closingGU && (bonus == 0));

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
	      if (strncmp (sequence + i, Tetraloops[k], 4) == 0)
		tetracorr = TETRA_ENERGY;
	    }
	}
    }
  mm = (unpaired > 3) ? mismatch[S1[i]][S1[j]][S1[i + 1]][S1[j - 1]] : 0;

  bonus = 0;
  
  if ((BP[i] == j) || (BP[i] == -1) || (BP[j] == -1) || (BP[i] == -2) || (BP[j] == -3))
    bonus -= BONUS;

  if (no_close)
    {
      if (c[indx[j] + i] == FORBIDDEN)
	{
	  no_push (nopush, state);
	  printf ("bullshit!\n\n");
	  fflush (stdout);
	  rem_en -= part_energy;
	  rem_en -= temp_energy;
	  return;
	}
    }
  else
    {
      en_intern = sizecorr + tetracorr + bonus + mm;
      if (hairpin[unpaired] + en_intern + rem_en <= min_en + delta)
	                                             /* hairpinstructure found */
	{
	  new_state = copy_state (state);
	  new_basepair = make_pair (i, j);
	  push (new_state->BasePairs, new_basepair);
	  new_state->cur_en += part_energy;
	  new_state->cur_en += hairpin[unpaired] + en_intern;
	  push (Stack, new_state);
	  nopush = 0;
	}
    }

  for (p = i + 1; p <= MIN2 (j - 2 - TURN, i + MAXLOOP + 1); p++)
                                  /* stack-, bulge-, interiorloopdecomposition */
    {
      for (q = j - 1; q >= p + 1 + TURN; q--)
	{
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
	  
	  if ((BP[p] == q) && (type_2 == 0))
	    type_2 = 7;		                                /* nonstandard */

	  if (type_2 == 0)
	    continue;

	  if (no_closingGU)
	    no_close_2 = (no_close || ((type_2 == 3) || (type_2 == 4)));

	  if (n2 == 0)
	    energy = stack[type][type_2];
	  
	  else if (n1 == 0)
	    {			/* bulge */
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
	    {			                              /* interior loop */
	      if (!no_close_2)
		{
		  energy = internal_loop[n1 + n2];

		  m = MIN2 (4, n1);
		  energy += MIN2 (MAX_NINIO, (abs (n1 - n2) * F_ninio[m]));

		  energy += mismatch[S1[i]][S1[j]][S1[i + 1]][S1[j - 1]] +
		    mismatch[S1[p - 1]][S1[q + 1]][S1[p]][S1[q]];
		}
	      else
		energy = FORBIDDEN;
	    }

	  new = energy + c[indx[q] + p] + bonus;

	  if (new + rem_en <= min_en + delta)
	                                   /* stack, bulge, interiorloop found */
	    {
	      new_state = copy_state (state);
	      new_basepair = make_pair (i, j);
	      push (new_state->BasePairs, new_basepair);
	      new_basepair = make_pair (p, q);
	      push (new_state->BasePairs, new_basepair);
	      new_interval = make_interval (p, q, 2);
	      push (new_state->Intervals, new_interval);
	      new_state->cur_en += part_energy;
	      new_state->cur_en += energy + bonus;
	      push (Stack, new_state);
	      nopush = 0;
	    }
	}                                                     /* end of q-loop */ 
    }                                                         /* end of p-loop */

 
  mm = bonus + MLclosing[type];
  type = pair[S[j]][S[i]];
  if ((type == 0) && (BP[i] == j))
    type = 7;

  for (k = i + 1 + TURN; k <= j - 1 - TURN; k++)
                                                     /* multiloopdecomposition */
    {
      if (dangles)
	{
	  en_intern = dangle3[type][S1[i + 1]] + dangle5[type][S1[j - 1]] + mm;
	  if (fML[indx[k] + i + 1] + fM1[indx[j - 1] + k + 1] + en_intern + rem_en <= min_en +
	      delta)
	    {
	      new_state = copy_state (state);
	      new_interval = make_interval (i + 1, k, 1);
	      push (new_state->Intervals, new_interval);
	      new_interval = make_interval (k + 1, j - 1, 3);
	      push (new_state->Intervals, new_interval);
	      new_basepair = make_pair (i, j);
	      push (new_state->BasePairs, new_basepair);
	      new_state->cur_en += part_energy;
	      new_state->cur_en += en_intern;
	      push (Stack, new_state);
	      nopush = 0; 
	    }
	}
      else
	{
	  if (fML[indx[k] + i + 1] + fM1[indx[j - 1] + k + 1] + mm + rem_en <= min_en + delta)
	    {
	      new_state = copy_state (state);
	      new_interval = make_interval (i + 1, k, 1);
	      push (new_state->Intervals, new_interval);
	      new_interval = make_interval (k + 1, j - 1, 3);
	      push (new_state->Intervals, new_interval);
	      new_basepair = make_pair (i, j);
	      push (new_state->BasePairs, new_basepair);
	      new_state->cur_en += part_energy;
	      new_state->cur_en += mm;
	      push (Stack, new_state);
	      nopush = 0;
	    }
	}
    }                                                         /* end of k-loop */
  
  rem_en -= part_energy;
  rem_en -= temp_energy;
  return;
}				                           /* end of repeat () */

/*-----------------------------------------------------------------------------*/

PUBLIC void
no_push (int nopush, STATE * state)

{
  if (nopush)			                             /* nothing pushed */
    {
      push (Stack, copy_state (state));
    } 
  return;
}

/*-----------------------------------------------------------------------------*/

PUBLIC char *
print_structure (STATE * state, char *solution)
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
  
  return (solution);
}

/*-----------------------------------------------------------------------------*/
/* Well, that is the end!------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/
