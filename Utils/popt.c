/*
   extract Zuker's p-optimal folds from subopt output.
   prints p-optimal structures followed by the list of pairs
   with respect to which they are optimal

   input must be sorted by energy! 

   Ivo Hofacker, Walter Fontana 
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "utils.h"
#define PRIVATE static
#define PUBLIC


PRIVATE void
help ()
{
  printf ("p-optimal filter to subopt output. usage:\n"
	  "RNAsubopt -s < seq | popt \n");
  exit(EXIT_FAILURE);
}

/*-------------------------------------------------------------------------*/

int main (int argc, char *argv[])
{
  char *line, *sequence, *structure;
  int i, j, length;
  float energy;

  short *ptable;
  char *pair_seen;
  int  *iindx;

  if (argc>1) help();
  
  line = get_line (stdin);
  if (*line=='>') {
    free (line);
    line = get_line (stdin);
  }

  /* read the sequence */
  sequence = (char *) space (sizeof (char) * (strlen(line)+1));
  (void) sscanf(line, "%s", sequence);
  free (line);
  length = (int) strlen(sequence);

  iindx = (int *) space(sizeof(int)*(length+1));
  for (i=1; i<=length; i++) 
    iindx[i] = ((length+1-i)*(length-i))/2 +length+1;
  pair_seen = (char *) space(((length+1)*(length+2))/2 * sizeof(char));

  structure = (char *) space (sizeof (char) * (length+1));
  
  /* get of suboptimal structures */

  while ((line = get_line (stdin))) {
    int r, popt;
    
    r = sscanf(line, "%s %f", structure, &energy);
    free(line);
    if (r==0) continue;
    ptable = make_pair_table(structure);

    popt = 0;
    for (i = 0; i < length; i++) {
      if ((j=ptable[i])>i)
	if (pair_seen[iindx[i]-j]==0) {
	  pair_seen[iindx[i]-j]=(char) 1;
	  if (popt==0) {
	    printf("%s %6.2f", structure, energy);
	    popt = 1;
	  }
	  printf(" (%d,%d)", i,j);
	}
    }
    if (popt==1) printf("\n");
    free(ptable);
  }
  free(sequence); free(structure); free(iindx); free(pair_seen);
  
  return 0;
}
