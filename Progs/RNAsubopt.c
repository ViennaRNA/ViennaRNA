/* Last changed Time-stamp: <2008-12-03 16:38:01 ivo> */
/*
		Ineractive access to suboptimal folding

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
#include "part_func.h"
#include "fold.h"
#include "cofold.h"
#include "fold_vars.h"
#include "utils.h"
#include "subopt.h"
extern void  read_parameter_file(const char fname[]);
extern int   st_back;
/*@unused@*/
static char UNUSED rcsid[] = "$Id: RNAsubopt.c,v 1.20 2008/12/03 16:55:44 ivo Exp $";

#define PRIVATE static

static char  scale[] = "....,....1....,....2....,....3....,....4"
		       "....,....5....,....6....,....7....,....8";

extern double print_energy;
#define MAXDOS 1000
extern int    density_of_states[MAXDOS+1];
PRIVATE char *tokenize(char *line);
PRIVATE void usage(void);
PRIVATE void putoutzuker(SOLUTION* zukersolution);

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
   double deltaf, deltap=0;
   int delta=100;
   int n_back = 0;
   int noconv = 0;
   int circ=0;
   int dos=0;
   int zuker=0;
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
	  case 'p':
	    if (argv[i][2]!='\0') usage();
	    if(i==argc-1) usage();
	    (void) sscanf(argv[++i], "%d", &n_back);
	    init_rand();
	    break;
	  case 'n':
	    if ( strcmp(argv[i], "-noGU" )==0) noGU=1;
	    if ( strcmp(argv[i], "-noCloseGU" ) ==0) no_closingGU=1;
	    if ( strcmp(argv[i], "-noLP")==0) noLonelyPairs=1;
	    if ( strcmp(argv[i], "-nsp") ==0) {
	      if (i==argc-1) usage();
	      ns_bases = argv[++i];
	    }
	    if ( strcmp(argv[i], "-noconv")==0) noconv=1;
	    break;
	  case '4':
	    tetra_loop=0;
	    break;
	  case 'C':
	    fold_constrained=1;
	    break;
	  case 'D':
	    dos=1;
	    print_energy = -999999;
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
	    subopt_sorted=1;
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
	      r=sscanf(argv[++i], "%lf", &deltap);
	    else {
	      r=sscanf(argv[++i], "%lf", &deltaf);
	      delta = (int) (0.1+deltaf*100);
	    }
	    if (r!=1) usage();
	    break;
	  case 'c':
	    if ( strcmp(argv[i], "-circ")==0) circ=1;
	    break;
	  case 'z':
	    zuker=1;
	    break;
	  default: usage();
	  }
   }

   if ((zuker)&&(circ)) {
     printf("Sorry, zuker subopts not yet implemented for circfold\n");
     usage();
   }
   if ((zuker)&&(n_back>0)) {
     printf("Cna't do zuker subopts and stochastic subopts at the same time\n");
     usage();
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
   }

   do {				/* main loop: continue until end of file */
     cut_point = -1;
     if (istty) {
       printf("\nInput string (upper or lower case); @ to quit\n");
       if (!zuker)printf("Use '&' to connect 2 sequences that shall form a complex.\n");
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

     if ((line==NULL)||strcmp(line,"@")==0) break;

     sequence = tokenize(line); /* frees line */
     length = (int) strlen(sequence);
     structure = (char *) space((unsigned) length+1);

     if (fold_constrained) {
       char *cstruc;
       cstruc = tokenize(get_line(stdin));
       if (cstruc!=NULL) {
	 strncpy(structure, cstruc, length);
	 for (i=0; i<length; i++)
	   if (structure[i]=='|')
	     nrerror("constraints of type '|' not allowed");
	 free(cstruc);
       }
     }

     for (l = 0; l < length; l++) {
       sequence[l] = toupper(sequence[l]);
       if (!noconv && sequence[l] == 'T') sequence[l] = 'U';
     }
     if (istty) {
       if (cut_point == -1)
	 printf("length = %d\n", length);
       else
	 printf("length1 = %d\nlength2 = %d\n",
		cut_point-1, length-cut_point+1);
     }


     if ((logML!=0 || dangles==1 || dangles==3) && dos==0)
       if (deltap<=0) deltap=delta/100. +0.001;
     if (deltap>0)
       print_energy = deltap;

     /* first lines of output (suitable  for sort +1n) */
     if (fname[0] != '\0')
       printf("> %s [%d]\n", fname, delta);

     if (n_back>0) {  /* stochastic backtrack */
       double mfe, kT;
       char *ss;
       st_back=1;
       ss = (char *) space(strlen(sequence)+1);
       strncpy(ss, structure, length);
       mfe = fold(sequence, ss);
       kT = (temperature+273.15)*1.98717/1000.; /* in Kcal */
       pf_scale = exp(-(1.03*mfe)/kT/length);
       strncpy(ss, structure, length);
       /* ignore return value, we are not interested in the free energy */
       (circ) ? (void) pf_circ_fold(sequence, ss) : (void) pf_fold(sequence, ss);
       free(ss);
       for (i=0; i<n_back; i++) {
	 char *s;
	 s =(circ) ? pbacktrack_circ(sequence) : pbacktrack(sequence);
	 printf("%s\n", s);
	 free(s);
       }
       free_pf_arrays();
     } else if (!zuker) { /* normal subopt */
       (circ) ? subopt_circ(sequence, structure, delta, stdout) : subopt(sequence, structure, delta, stdout);
       if (dos) {
	 int i;
	 for (i=0; i<= MAXDOS && i<=delta/10; i++) {
	   printf("%4d %6d\n", i, density_of_states[i]);
	 }
       }
     } else { /* Zuker suboptimals */
       SOLUTION *zr;
       int i;
       if (cut_point!=-1) {
	 printf("Sorry, zuker subopts not yet implemented for cofold\n");
	 usage();
       }
       zr = zukersubopt(sequence);
       putoutzuker(zr);
       (void)fflush(stdout);
       for (i=0; zr[i].structure; i++) {
	 free(zr[i].structure);
       }
       free(zr);
     }
     (void)fflush(stdout);
     free(sequence);
     free(structure);
   } while (1);
   return 0;
}

PRIVATE void usage(void)
{
   nrerror("usage: "
	   "RNAsubopt [-e range] [-ep prange] [-s] [-p num] [-logML]\n"
	   "          [-C] [-T temp] [-4] [-d[2]] [-noGU] [-noCloseGU]\n"
	   "          [-noLP] [-P paramfile] [-nsp pairs] [-circ] [-z]");
}
PRIVATE char *tokenize(char *line)
{
  char *pos, *copy;
  int cut = -1;

  copy = (char *) space(strlen(line)+1);
  (void) sscanf(line, "%s", copy);
  pos = strchr(copy, '&');
  if (pos) {
    cut = (int) (pos-copy)+1;
    if (cut >= strlen(copy)) cut = -1;
    if (strchr(pos+1, '&')) nrerror("more than one cut-point in input");
    for (;*pos;pos++) *pos = *(pos+1); /* splice out the & */
  }
  if (cut > -1) {
    if (cut_point==-1) cut_point = cut;
    else if (cut_point != cut) {
      fprintf(stderr,"cut_point = %d cut = %d\n", cut_point, cut);
      nrerror("Sequence and Structure have different cut points.");
    }
  }

  free(line);
  return copy;
}
PRIVATE void putoutzuker(SOLUTION* zukersolution) {
  int i;
  printf("%s [%.2f]\n",zukersolution[0].structure,zukersolution[0].energy/100.);
  for(i=1; zukersolution[i].structure; i++) {
    printf("%s [%.2f]\n", zukersolution[i].structure, zukersolution[i].energy/100.);
  }
  return;
}
