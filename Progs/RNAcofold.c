/* Last changed Time-stamp: <2004-02-05 15:29:47 ivo> */
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
#include "cofold.h"
#include "fold.h"
#include "part_func.h"
#include "fold_vars.h"
#include "PS_dot.h"
#include "utils.h"
extern void  read_parameter_file(const char fname[]);

/*@unused@*/
static char rcsid[] = "$Id: RNAcofold.c,v 1.2 2004/02/09 16:52:51 ivo Exp $";

#define PRIVATE static

static char  scale1[] = "....,....1....,....2....,....3....,....4";
static char  scale2[] = "....,....5....,....6....,....7....,....8";

PRIVATE char *costring(char *string);
PRIVATE char *tokenize(char *line);
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
  int   pf=0, istty;
  int noconv=0;
   
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
    cut_point = -1;
    if (istty) {
      printf("\nInput string (upper or lower case); @ to quit\n");
      printf("Use '&' to connect 2 sequences that shall form a complex.\n");
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

    string = tokenize(line); /* frees line */

    length = (int) strlen(string);

    structure = (char *) space((unsigned) length+1);
    if (fold_constrained) { 
      cstruc = tokenize(get_line(stdin));
      if (cstruc!=NULL)
        strncpy(structure, cstruc, length);
      else
        fprintf(stderr, "constraints missing\n");
    }
    
    for (l = 0; l < length; l++) {
      string[l] = toupper(string[l]);
      if (!noconv && string[l] == 'T') string[l] = 'U';
    }
    if (istty) {
      if (cut_point == -1)
	printf("length = %d\n", length);
      else
	printf("length1 = %d\nlength2 = %d\n",
	       cut_point-1, length-cut_point+1);
    }

    /* initialize_fold(length); */
    min_en = cofold(string, structure);
    if (cut_point == -1)
      printf("%s\n%s", string, structure);
    else {
      char *pstring, *pstruct;
      pstring = costring(string);
      pstruct = costring(structure);
      printf("%s\n%s", pstring,  pstruct);
      free(pstring);
      free(pstruct);
    }
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
    if (length<2000)
      (void) PS_rna_plot(string, structure, ffname);
    else { 
      struct bond  *bp;
      fprintf(stderr,"INFO: structure too long, not doing xy_plot\n");

      /* free mfe arrays but preserve base_pair for PS_dot_plot */
      bp = base_pair; base_pair = space(16);
      free_co_arrays();  /* free's base_pair */
      base_pair = bp;
    } 
    if (pf) {

      if (dangles==1) {
	dangles=2;   /* recompute with dangles as in pf_fold() */
	min_en = energy_of_struct(string, structure);
	dangles=1;
      }
	 
      kT = (temperature+273.15)*1.98717/1000.; /* in Kcal */
      pf_scale = exp(-(sfact*min_en)/kT/length);
      if (length>2000) fprintf(stderr, "scaling factor %f\n", pf_scale);

      init_pf_fold(length);

      if (cstruc!=NULL)
	strncpy(structure, cstruc, length+1);
      energy = pf_fold(string, structure);
	 
      if (do_backtrack) {
	printf("%s", structure);
	if (!istty) printf(" [%6.2f]\n", energy);
	else printf("\n");
      }
      if ((istty)||(!do_backtrack)) 
	printf(" free energy of ensemble = %6.2f kcal/mol\n", energy);
      printf(" frequency of mfe structure in ensemble %g\n",
	     exp((energy-min_en)/kT));
	 
      if (do_backtrack) {
	if (fname[0]!='\0') {
	  strcpy(ffname, fname);
	  strcat(ffname, "_dp.ps");
	} else strcpy(ffname, "dot.ps");
	(void) PS_dot_plot(string, ffname);
      }
      free_pf_arrays();

    }
    if (cstruc!=NULL) free(cstruc);
    if (length>=2000) free(base_pair);
    (void) fflush(stdout);
    free(string);
    free(structure); 
  } while (1);
  return 0;
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

PRIVATE char *costring(char *string)
{
  char *ctmp;
  int len;
  
  len = strlen(string);
  ctmp = (char *)space((len+2) * sizeof(char));
  /* first sequence */
  (void) strncpy(ctmp, string, cut_point-1);
  /* spacer */
  ctmp[cut_point-1] = '&';
  /* second sequence */
  (void) strcat(ctmp, string+cut_point-1);

  return ctmp;
}

PRIVATE void usage(void)
{
  nrerror("usage:\n"
	  "RNAfold [-p[0]] [-C] [-T temp] [-4] [-d[2|3]] [-noGU] [-noCloseGU]\n" 
	  "        [-noLP] [-e e_set] [-P paramfile] [-nsp pairs] [-S scale] "
	  "[-noconv]\n");
}
