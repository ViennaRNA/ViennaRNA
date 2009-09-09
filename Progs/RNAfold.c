/* Last changed Time-stamp: <2009-09-08 19:54:31 ivo> */
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
#include "PS_dot.h"
#include "utils.h"
#include "MEA.h"
extern void  read_parameter_file(const char fname[]);
extern plist * stackProb(double cutoff);

/*@unused@*/
static char UNUSED rcsid[] = "$Id: RNAfold.c,v 1.25 2009/02/24 14:22:21 ivo Exp $";

#define PRIVATE static

static char  scale1[] = "....,....1....,....2....,....3....,....4";
static char  scale2[] = "....,....5....,....6....,....7....,....8";

PRIVATE struct plist *b2plist(const char *struc);
PRIVATE struct plist *make_plist(int length, double pmin);
PRIVATE void usage(void);

/*--------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
  char *string, *line;
  char *structure=NULL, *cstruc=NULL;
  char  fname[13], ffname[20], gfname[20];
  char  *ParamFile=NULL;
  char  *ns_bases=NULL, *c;
  int   i, length, l, sym, r;
  double energy, min_en;
  double kT, sfact=1.07;
  int   pf=0, noPS=0, istty;
  int noconv=0;
  int circ=0;
  int doMEA=0;
  double MEAgamma = 1.;

  do_backtrack = 1;
  string=NULL;
  for (i=1; i<argc; i++) {
    if (argv[i][0]=='-')
      switch ( argv[i][1] )
	{
	case 'T':  if (argv[i][2]!='\0') usage();
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i], "%lf", &temperature);
	  if (!r) usage();
	  break;
	case 'p':  pf=1;
	  if (argv[i][2]!='\0')
	    (void) sscanf(argv[i]+2, "%d", &do_backtrack);
	  break;
	case 'n':
	  if ( strcmp(argv[i], "-noGU")==0) noGU=1;
	  if ( strcmp(argv[i], "-noCloseGU")==0) no_closingGU=1;
	  if ( strcmp(argv[i], "-noLP")==0) noLonelyPairs=1;
	  if ( strcmp(argv[i], "-noPS")==0) noPS=1;
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
	case 'C':
	  fold_constrained=1;
	  break;
	case 'c':
	  if ( strcmp(argv[i], "-circ")==0) circ=1;
	  break;
	case 'S':
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i],"%lf", &sfact);
	  if (!r) usage();
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
	case 'M':
	  if (strcmp(argv[i], "-MEA")==0) pf=doMEA=1;
	  if (i<argc-1)
	    if (isdigit(argv[i+1][0]) || argv[i+1][0]=='.') {
	      r=sscanf(argv[++i], "%lf", &MEAgamma);
	      if (r!=1) usage();
	    }
	  break;

	default: usage();
	}
  }

  if (circ && noLonelyPairs)
    fprintf(stderr, "warning, depending on the origin of the circular sequence, some structures may be missed when using -noLP\nTry rotating your sequence a few times\n");
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
  if ((fold_constrained)&&(istty)) {
    printf("Input constraints using the following notation:\n");
    printf("| : paired with another base\n");
    printf(". : no constraint at all\n");
    printf("x : base must not pair\n");
    printf("< : base i is paired with a base j<i\n");
    printf("> : base i is paired with a base j>i\n");
    printf("matching brackets ( ): base i pairs base j\n");
  }

  do {				/* main loop: continue until end of file */
    if (istty) {
      printf("\nInput string (upper or lower case); @ to quit\n");
      printf("%s%s\n", scale1, scale2);
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
    if (fold_constrained) {
      cstruc = get_line(stdin);
      if (cstruc!=NULL)
	strncpy(structure, cstruc, length);
      else
	fprintf(stderr, "constraints missing\n");
    }
    for (l = 0; l < length; l++) {
      string[l] = toupper(string[l]);
      if (!noconv && string[l] == 'T') string[l] = 'U';
    }
    if (istty)
      printf("length = %d\n", length);

    /* initialize_fold(length); */
    if (circ)
      min_en = circfold(string, structure);
    else
      min_en = fold(string, structure);
    printf("%s\n%s", string, structure);
    if (istty)
      printf("\n minimum free energy = %6.2f kcal/mol\n", min_en);
    else
      printf(" (%6.2f)\n", min_en);

    (void) fflush(stdout);

    if (fname[0]!='\0') {
      strcpy(ffname, fname);
      strcat(ffname, "_ss.ps");
      strcpy(gfname, fname);
      strcat(gfname, "_ss.g");
    } else {
      strcpy(ffname, "rna.ps");
      strcpy(gfname, "rna.g");
    }
    if (!noPS) (void) PS_rna_plot(string, structure, ffname);
    if (length>2000) free_arrays();
    if (pf) {
      char *pf_struc;
      pf_struc = (char *) space((unsigned) length+1);
	if (dangles==1) {
	  dangles=2;   /* recompute with dangles as in pf_fold() */
	  min_en = (circ) ? energy_of_circ_struct(string, structure) : energy_of_struct(string, structure);
	  dangles=1;
      }

      kT = (temperature+273.15)*1.98717/1000.; /* in Kcal */
      pf_scale = exp(-(sfact*min_en)/kT/length);
      if (length>2000) fprintf(stderr, "scaling factor %f\n", pf_scale);

      (circ) ? init_pf_circ_fold(length) : init_pf_fold(length);

      if (cstruc!=NULL)
	strncpy(pf_struc, cstruc, length+1);
      energy = (circ) ? pf_circ_fold(string, pf_struc) : pf_fold(string, pf_struc);

      if (do_backtrack) {
	printf("%s", pf_struc);
	if (!istty) printf(" [%6.2f]\n", energy);
	else printf("\n");
      }
      if ((istty)||(!do_backtrack))
	printf(" free energy of ensemble = %6.2f kcal/mol\n", energy);
      if (do_backtrack) {
	plist *pl1,*pl2;
	char *cent;
	double dist, cent_en;
	cent = centroid(length, &dist);
	cent_en = (circ) ? energy_of_circ_struct(string, cent) :energy_of_struct(string, cent);
	printf("%s {%6.2f d=%.2f}\n", cent, cent_en, dist);
	free(cent);
	if (fname[0]!='\0') {
	  strcpy(ffname, fname);
	  strcat(ffname, "_dp.ps");
	} else strcpy(ffname, "dot.ps");
	pl1 = make_plist(length, 1e-5);
	pl2 = b2plist(structure);
	(void) PS_dot_plot_list(string, ffname, pl1, pl2, "");
	free(pl2);
	if (do_backtrack==2) {
	  pl2 = stackProb(1e-5);
	  if (fname[0]!='\0') {
	    strcpy(ffname, fname);
	    strcat(ffname, "_dp2.ps");
	  } else strcpy(ffname, "dot2.ps");
	  PS_dot_plot_list(string, ffname, pl1, pl2,
			   "Probabilities for stacked pairs (i,j)(i+1,j-1)");
	  free(pl2);
	}
	free(pl1);
	free(pf_struc);
	if (doMEA) {
	  float mea, mea_en;
	  plist *pl;
	  pl = make_plist(length, 1e-4/(1+MEAgamma));
	  mea = MEA(pl, structure, MEAgamma);
	  mea_en = (circ) ? energy_of_circ_struct(string, structure) :
			    energy_of_struct(string, structure);
	  printf("%s {%6.2f MEA=%.2f}\n", structure, mea_en, mea);
	  free(pl);
	}
      }
      printf(" frequency of mfe structure in ensemble %g; ",
	     exp((energy-min_en)/kT));
      if (do_backtrack)
	printf("ensemble diversity %-6.2f", mean_bp_dist(length));

      printf("\n");


      free_pf_arrays();

    }
    if (cstruc!=NULL) free(cstruc);
    (void) fflush(stdout);
    free(string);
    free(structure);
  } while (1);
  if (length<=2000)
    free_arrays();
  return 0;
}

PRIVATE struct plist *b2plist(const char *struc) {
  /* convert bracket string to plist */
  short *pt;
  struct plist *pl;
  int i,k=0;
  pt = make_pair_table(struc);
  pl = (struct plist *)space(strlen(struc)/2*sizeof(struct plist));
  for (i=1; i<strlen(struc); i++) {
    if (pt[i]>i) {
      pl[k].i = i;
      pl[k].j = pt[i];
      pl[k++].p = 0.95*0.95;
    }
  }
  free(pt);
  pl[k].i=0;
  pl[k].j=0;
  pl[k++].p=0.;
  return pl;
}


PRIVATE struct plist *make_plist(int length, double pmin) {
  /* convert matrix of pair probs to plist */
  struct plist *pl;
  int i,j,k=0,maxl;
  maxl = 2*length;
  pl = (struct plist *)space(maxl*sizeof(struct plist));
  k=0;
  for (i=1; i<length; i++)
    for (j=i+1; j<=length; j++) {
      if (pr[iindx[i]-j]<pmin) continue;
      if (k>=maxl-1) {
	maxl *= 2;
	pl = (struct plist *)xrealloc(pl,maxl*sizeof(struct plist));
      }
      pl[k].i = i;
      pl[k].j = j;
      pl[k++].p = pr[iindx[i]-j];
    }
  pl[k].i=0;
  pl[k].j=0;
  pl[k++].p=0.;
  return pl;
}

PRIVATE void usage(void)
{
  nrerror("usage:\n"
	  "RNAfold [-p[0|2]] [-C] [-T temp] [-4] [-d[2|3]] [-noGU] [-noCloseGU]\n"
	  "        [-noLP] [-e e_set] [-P paramfile] [-nsp pairs] [-S scale]\n"
	  "        [-noconv] [-noPS] [-circ] [-MEA [gamma]]\n");
}
