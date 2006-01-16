/* Last changed Time-stamp: <2006-01-16 10:38:56 ivo> */
/*
		  Ineractive Access to folding Routines

		  c Ivo L Hofacker
		  Vienna RNA package
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include "fold.h"
#include "part_func.h"
#include "fold_vars.h"
#include "utils.h"
#include "PS_dot.h"

extern float Lfold(char *string, char *structure, int maxdist);
extern void  read_parameter_file(const char fname[]);
extern int pfl_fold(char *sequence, int winSize, float cutoff, struct plist **pl);
extern void  init_pf_foldLP(int length);
extern void  free_pf_arraysLP(void);

/*@unused@*/
static char rcsid[] = "$Id: RNAplfold.c,v 1.2 2006/01/16 09:49:07 ivo Exp $";

#define PRIVATE static

static char  scale[] = "....,....1....,....2....,....3....,....4"
	"....,....5....,....6....,....7....,....8";

PRIVATE void usage(void);

/*--------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
  char *string, *line;
  char *structure=NULL, *cstruc=NULL;
  char  fname[30], ffname[20];
  char  *ParamFile=NULL;
  char  *ns_bases=NULL, *c;
  int   i, length, l, sym,r;
  int   istty;
  int noconv=0;
  int maxdist=70;
  int winSize;
  float cutoff=0.01;
  int hit;
  plist *pl;
  do_backtrack = 1;
  string=NULL;
  dangles=2;
  for (i=1; i<argc; i++) {
    if (argv[i][0]=='-')
      switch ( argv[i][1] )
	{
	case 'T':  if (argv[i][2]!='\0') usage();
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i], "%lf", &temperature);
	  if (!r) usage();
	  break;
	case 'n':
	  if ( strcmp(argv[i], "-noGU")==0) noGU=1;
	  if ( strcmp(argv[i], "-noCloseGU")==0) no_closingGU=1;
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
	case 'e':
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i],"%d", &energy_set);
	  if (!r) usage();
	  break;
	case 'c':
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i],"%f", &cutoff);
	  if (!r) usage();
	  break;
/*	case 'C':
	  fold_constrained=1;
	  break;
	case 'S':
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i],"%lf", &sfact);
	  if (!r) usage();
	  break;
*/
	case 'd': dangles=0;
	  if (argv[i][2]!='\0') {
	    r=sscanf(argv[i]+2, "%d", &dangles);
	    if (r!=1) usage();
	    if (!dangles%2) usage(); /*1,3 not possible*/
	  }
	  break;
	case 'P':
	  if (i==argc-1) usage();
	  ParamFile = argv[++i];
	  break;
	case 'L':
	  if (i==argc-1) usage();
	  r=sscanf(argv[++i], "%d", &maxdist);
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
    while (*c!='\0') {
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
  winSize=maxdist;

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
	(void) sscanf(line, ">%12s", fname);
      printf("%s\n", line);
      free(line);
      if ((line = get_line(stdin))==NULL) break;
    }

    if ((line ==NULL) || (strcmp(line, "@") == 0)) break;

    string = (char *) space(strlen(line)+1);
    (void) sscanf(line,"%s",string);
    free(line);
    length = (int) strlen(string);

    structure = (char *) space((unsigned) length+1);

    for (l = 0; l < length; l++) {
      string[l] = toupper(string[l]);
      if (!noconv && string[l] == 'T') string[l] = 'U';
    }
    if (istty)
      printf("length = %d\n", length);

    /* initialize_fold(length); */
    update_fold_params();
    maxdist=winSize;
    if (length<maxdist) {
      fprintf(stderr, "WARN: window size %d larger than sequence length %d\n",
	      maxdist, length);
      maxdist=length;
    }
    if (length >= 5) {
      pf_scale = -1;

      init_pf_foldLP(length);

      hit=pfl_fold(string, maxdist, cutoff, &pl);
      free_pf_arraysLP();
      if (fname[0]!='\0') {
	strcpy(ffname, fname);
	strcat(ffname, "_dp.ps");
      }
      else strcpy(ffname, "plfold_dp.ps");
      PS_dot_plot_turn(string, pl, ffname, maxdist);
      free(pl);

      if (cstruc!=NULL) free(cstruc);
    (void) fflush(stdout);
    }
    free(string);
    free(structure);
  } while (1);
  return 0;
}

PRIVATE void usage(void)
{
  nrerror("usage:\n"
	  "RNAplfold [-L span]\n"
	  "          [-T temp] [-4] [-d[0|1|2]] [-noGU] [-noCloseGU]\n"
	  "          [-noLP] [-P paramfile] [-nsp pairs] [-noconv]\n");
}
