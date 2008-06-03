/* Last changed Time-stamp: <2008-03-25 23:13:39 ivo> */
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
#include "LPfold.h"

extern float Lfold(char *string, char *structure, int winsize);
extern void  read_parameter_file(const char fname[]);

/*@unused@*/
static char rcsid[] = "$Id: RNAplfold.c,v 1.9 2008/06/03 21:19:42 ivo Exp $";

#define PRIVATE static

static char  scale[] = "....,....1....,....2....,....3....,....4"
	"....,....5....,....6....,....7....,....8";

PRIVATE void usage(void);
PRIVATE void putout_pup(double *pup,int length, int winsize);
int unpaired;
/*--------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
  char *string=NULL, *line;
  char *structure=NULL, *cstruc=NULL;
  char  fname[30], ffname[20];
  char  *ParamFile=NULL;
  char  *ns_bases=NULL, *c;
  int   i, length, l, sym,r;
  int   istty;
  int noconv=0;
  int winsave, winsize=70;
  int pairdist=0;
  float cutoff=0.01;
  double *pup=NULL; /*prob of being unpaired*/
  plist *pl;

  dangles=2;
  unpaired=0;
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
	case 'W':
	  if (i==argc-1) usage();
	  r=sscanf(argv[++i], "%d", &winsize);
	  if (r!=1) usage();
	  break;
	case 'u':
	  if (i==argc-1) usage();
	  r=sscanf(argv[++i], "%d", &unpaired);
	  if (r!=1) usage();
	  if (unpaired<0) usage();
	  break;
	case 'L':
	  if (i==argc-1) usage();
	  r=sscanf(argv[++i], "%d", &pairdist);
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
   if (pairdist==0) pairdist=winsize;
  if (pairdist>winsize) {
    fprintf(stderr, "pairdist (-L %d) should be <= winsize (-W %d);"
	    "Setting pairdist=winsize\n",pairdist, winsize);
    pairdist=winsize;
  }
  winsave = winsize;
  do {				/* main loop: continue until end of file */
    winsize = winsave;
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
    if (unpaired) {
      pup=(double *)space((length+1)*sizeof(double));
      pup[0]=unpaired;
    }

    structure = (char *) space((unsigned) length+1);

    for (l = 0; l < length; l++) {
      string[l] = toupper(string[l]);
      if (!noconv && string[l] == 'T') string[l] = 'U';
    }
    if (istty)
      printf("length = %d\n", length);

    /* initialize_fold(length); */
    update_fold_params();
    if (length<winsize) {
      fprintf(stderr, "WARN: window size %d larger than sequence length %d\n",
	      winsize, length);
      winsize=length;
    }
    if (length >= 5) {
      pf_scale = -1;

      pl=pfl_fold(string, winsize, pairdist, cutoff, pup);
      if (fname[0]!='\0') {
	strcpy(ffname, fname);
	strcat(ffname, "_dp.ps");
      }
      else strcpy(ffname, "plfold_dp.ps");
      PS_dot_plot_turn(string, pl, ffname, pairdist);
      free(pl);
      if (unpaired) {
	putout_pup(pup,length,winsize);
	free(pup);
      }
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
	  "RNAplfold [-L span] [-W winsize] [-u size of unpaired region]\n"
	  "          [-T temp] [-4] [-d[0|1|2]] [-noGU] [-noCloseGU]\n"
	  "          [-noLP] [-P paramfile] [-nsp pairs] [-noconv]\n");
}

#define MAX(A,B) (A)<(B)?(B):(A)
#define MIN(A,B) (A)>(B)?(B):(A)
PRIVATE void putout_pup(double *pup,int length, int winsize) {
  int i;
  double factor;

  printf("# prob of being unpaired between i-%d and i\n",unpaired);
  fflush(NULL);
  for (i=unpaired; i<=length; i++) {
    int leftmost, rightmost;

    leftmost = MAX(i,winsize);
    rightmost = MIN(length,i-unpaired+winsize);

    factor = 1./(rightmost-leftmost+1);

    printf("%d %.6g\n",i,pup[i]*factor);
  }
}
