/* Last changed Time-stamp: <1998-04-08 22:44:40 ivo> */
/*                
		Ineractive Access to suboptimal folding

			   c Ivo L Hofacker
			  Vienna RNA package
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "fold.h"
#include "fold_vars.h"
#include "utils.h"
static char rcsid[] = "$Id: RNAsubopt.c,v 1.2 1998/05/19 16:31:17 ivo rel $";

#define PRIVATE static

static char  scale[] = "....,....1....,....2....,....3....,....4"
                       "....,....5....,....6....,....7....,....8";

PRIVATE void usage(void);

int sorted=0;
int delta=100;
char *sequence;

int LODOS_ONLY = 0;

extern void subopt(int delta);

/*--------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
   char  *line;
   char *structure = NULL, *cstruc = NULL;
   char  fname[21];
   char  ParamFile[256] = "";
   char  ns_bases[33] = "", *c;
   int   i, length, l, sym, r;
   float energy, min_en;
   int   istty;
   float deltaf;
   
   do_backtrack = 1;
   dangles = 2;
   sequence=NULL;
   for (i=1; i<argc; i++) {
      if (argv[i][0]=='-') 
	switch ( argv[i][1] )
	  {
	  case 'T':  if (argv[i][2]!='\0') usage();
	    if(i==argc-1) usage();
	    r=sscanf(argv[++i], "%f", &temperature);
	    if (!r) usage();
	    break;
	  case 'n':
	    if ( strcmp(argv[i], "-noGU" )==0) noGU=1;
	    if ( strcmp(argv[i], "-noCloseGU" ) ==0) no_closingGU=1;
	    if ( strcmp(argv[i], "-nsp") ==0) {
	      if (i==argc-1) usage();
	      r=sscanf(argv[++i], "%32s", ns_bases);
	      if (!r) usage();
	    }
	    break;
	  case '4':
	    tetra_loop=0;
	    break;
 	  case 'C': 
 	    fold_constrained=1; 
 	    break; 
	  case 'd': dangles=0;
	    if (strcmp(argv[i],"-d2")==0) dangles=2;
	    /* danlges == 1 not allowed in subopt() */
	    if (strcmp(argv[i],"-d1")==0) usage();
	    break;
	  case 'P':
	    if (i==argc-1) usage();
	    r=sscanf(argv[++i], "%255s", ParamFile);
	    if (!r) usage();
	    break;
	  case 's':
	    sorted=1;
	    break;
	  case 'l':
	    if (strcmp(argv[i],"-lodos")==0) {
	      LODOS_ONLY=1;
	      sorted = 1;
	    }
	    else usage();
	    break;
	  case 'e':
	    r=sscanf(argv[++i], "%f", &deltaf);
	    if (!r) usage();
	    delta = (int) (0.1+deltaf*100);
	    break;
	  default: usage();
	  } 
   }

   if (ParamFile[0])
     read_parameter_file(ParamFile);
   
   if (ns_bases[0]) {
      nonstandards = space(33);
      c=ns_bases;
      i=sym=0;
      if (*c=='-') {
	 sym=1; c++;
      }
      while (*c) {
	 if (*c!=',') {
	    nonstandards[i++]=*c++;
	    nonstandards[i++]=*c;
	    if ((sym)&&(*c!=*(c-1))) {
	       nonstandards[i++]=*c;
	       nonstandards[i++]=*(c-1);
	    }
	 }
	 c++;
      }
   }
   istty = isatty(fileno(stdout))&&isatty(fileno(stdin));
   if ((fold_constrained)&&(istty)) {
      printf("Input constraints using the following notation:\n");
      /* printf("| : paired with another base\n"); */
      printf(". : no constraint at all\n");
      printf("x : base must not pair\n");
      /* printf("< : base i is paired with a base j<i\n"); */
      /* printf("> : base i is paired with a base j>i\n"); */
      /* printf("matching brackets ( ): base i pairs base j\n"); */
   } 
	
   do {				/* main loop: continue until end of file */
      if (istty) {
	 printf("\nInput string (upper or lower case); @ to quit\n");
	 printf("%s\n", scale);
      }
      fname[0]='\0';
      if ((line = get_line(stdin))==NULL) break;

      /* skip comment lines and get filenames */
      while ((*line=='*')||(*line=='\0')||(*line=='>')) {
	if (*line=='>')
	  sscanf(line, ">%20s", fname);
	free(line);
	if ((line = get_line(stdin))==NULL) line = "@";
      } 

      if (strcmp(line, "@") == 0) break;

      sequence = (char *) space(strlen(line)+1);
      sscanf(line,"%s",sequence);
      free(line);
      length = strlen(sequence);

      if (fold_constrained) 
	 cstruc = get_line(stdin);
      
      structure = (char *) space(length+1);
      
      for (l = 0; l < length; l++) sequence[l] = toupper(sequence[l]);
      if (istty)
	printf("length = %d\n", length);

      initialize_fold(length);
      if (fold_constrained) {
	strncpy(structure, cstruc, length+1);
	for (i=0; i<length; i++)
	  if ((structure[i]!='.')&&(structure[i]!='x'))
	    nrerror("only constraints of type 'x' allowed");
      }
      min_en = fold(sequence, structure);
      /* first lines of output (suitable  for sort +1n) */
      if (fname[0] != '\0')
	printf("> %s [%6.2f to %6.2f]\n", fname, min_en, min_en+delta/100.);
      printf("%s %6d %6d\n", sequence, (int) (-0.1+100*min_en), delta); 
       
      subopt(delta);
      
      fflush(stdout);

      free_arrays();
      if (fold_constrained)
	free(cstruc);
      free(sequence);
      free(structure); 
   } while (1);
   return 0;
}

PRIVATE void usage(void)
{
   nrerror("usage: "
	   "RNAsubopt [-e range] [-s] [-lodos]\n"
	   "          [-C] [-T temp] [-4] [-d[2]] [-noGU] [-noCloseGU]\n" 
	   "          [-P paramfile] [-nsp pairs]");
}
