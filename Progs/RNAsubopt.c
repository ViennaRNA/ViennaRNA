/* Last changed Time-stamp: <2001-09-15 10:56:22 ivo> */
/*                
		Ineractive Access to suboptimal folding

			   c Ivo L Hofacker
			  Vienna RNA package
*/
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "fold_vars.h"
#include "utils.h"
#include "subopt.h"
extern void  read_parameter_file(const char fname[]);
/*@unused@*/
static char UNUSED rcsid[] = "$Id: RNAsubopt.c,v 1.6 2001/09/17 10:30:42 ivo Exp $";

#define PRIVATE static

static char  scale[] = "....,....1....,....2....,....3....,....4"
                       "....,....5....,....6....,....7....,....8";

extern float print_energy;
PRIVATE void usage(void);

/*--------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
   char  *line;
   char *sequence;
   char *structure = NULL;
   char  fname[21];
   char  *ParamFile = NULL;
   char  *ns_bases = NULL, *c;
   int   i, length, l, sym, r;
   int   istty;
   float deltaf, deltap=0;
   int delta=100;
   
   do_backtrack = 1;
   dangles = 2;
   for (i=1; i<argc; i++) {
      if (argv[i][0]=='-') 
	switch ( argv[i][1] )
	  {
	  case 'T':  if (argv[i][2]!='\0') usage();
	    if(i==argc-1) usage();
	    r=sscanf(argv[++i], "%lf", &temperature);
	    if (r!=1) usage();
	    break;
	  case 'n':
	    if ( strcmp(argv[i], "-noGU" )==0) noGU=1;
	    if ( strcmp(argv[i], "-noCloseGU" ) ==0) no_closingGU=1;
	    if ( strcmp(argv[i], "-noLP")==0) noLonelyPairs=1;
	    if ( strcmp(argv[i], "-nsp") ==0) {
	      if (i==argc-1) usage();
	      ns_bases = argv[++i];
	    }
	    break;
	  case '4':
	    tetra_loop=0;
	    break;
 	  case 'C': 
 	    fold_constrained=1; 
 	    break; 
	  case 'd': dangles=0;
	    if (argv[i][2]!='\0') {
              r=sscanf(argv[i]+2, "%d", &dangles);
	      if (r!=1) usage();
	    }
	    break;
	  case 'P':
	    if (i==argc-1) usage();
	    ParamFile = argv[++i];
	    break;
	  case 's':
	    sorted=1;
	    break;
	  case 'l':
	    if (strcmp(argv[i],"-logML")==0) {
	      logML=1;
	      break;
	    }
	    else usage();
	    break;
	  case 'e':
	    if (i>=argc-1) usage();
	    if (strcmp(argv[i],"-ep")==0) 
	      r=sscanf(argv[++i], "%f", &deltap);
	    else {
	      r=sscanf(argv[++i], "%f", &deltaf);
	      delta = (int) (0.1+deltaf*100);
	    }
	    if (r!=1) usage();
	    break;
	  default: usage();
	  } 
   }

   if (ParamFile != NULL)
     read_parameter_file(ParamFile);
   
   if (ns_bases != NULL) {
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
	  (void) sscanf(line, ">%20s", fname);
	free(line);
	if ((line = get_line(stdin))==NULL) break;;
      } 

      if ((line==NULL) || strcmp(line, "@") == 0) break;

      sequence = (char *) space(strlen(line)+1);
      (void) sscanf(line,"%s",sequence);
      free(line);
      length = (int) strlen(sequence);

      structure = (char *) space((unsigned) length+1);
      if (fold_constrained) {
	char *cstruc;
	cstruc = get_line(stdin);
	if (cstruc!=NULL) {
	  strncpy(structure, cstruc, length);
	  for (i=0; i<length; i++)
	    if (structure[i]=='|')
	      nrerror("constraints of type '|' not allowed");
	  free(cstruc);
	}
      }      
      
      for (l = 0; l < length; l++) sequence[l] = toupper(sequence[l]);
      if (istty)
	printf("length = %d\n", length);

      if (logML!=0 || dangles==1 || dangles==3) 
	if (deltap<=0) deltap=delta/100. +0.001;
      if (deltap>0)
	print_energy = deltap;

      /* first lines of output (suitable  for sort +1n) */
      if (fname[0] != '\0')
	printf("> %s [%d]\n", fname, delta);

      subopt(sequence, structure, delta, stdout);
      
      (void)fflush(stdout);
      
      free(sequence);
      free(structure); 
   } while (1);
   return 0;
}

PRIVATE void usage(void)
{
   nrerror("usage: "
	   "RNAsubopt [-e range] [-ep prange] [-s] [-logML]\n"
	   "          [-C] [-T temp] [-4] [-d[2]] [-noGU] [-noCloseGU]\n" 
	   "          [-noLP] [-P paramfile] [-nsp pairs]");
}
