/* Last changed Time-stamp: <1999-09-17 16:17:28 ivo> */
/*                
		Ineractive Access to folding Routines

			   c Ivo L Hofacker
			  Vienna RNA package
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "fold.h"
#include "part_func.h"
#include "fold_vars.h"
#include "PS_dot.h"
#include "utils.h"
extern void  read_parameter_file(const char fname[]);

static char rcsid[] = "$Id: RNAfold.c,v 1.11 1999/11/04 12:15:35 ivo Exp $";

#define PRIVATE static

static char  scale1[] = "....,....1....,....2....,....3....,....4";
static char  scale2[] = "....,....5....,....6....,....7....,....8";

PRIVATE void usage(void);

static struct bond  *bp, *bpp;


/*--------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
   char *string, *line;
   char *structure=NULL, *cstruc=NULL;
   char  fname[13], ffname[20], gfname[20];
   char  ParamFile[256]="";
   char  ns_bases[33]="", *c;
   int   i, length, l, sym, r;
   float energy, min_en;
   float kT, sfact=1.07;
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
	    r=sscanf(argv[++i], "%f", &temperature);
	    if (!r) usage();
	    break;
	  case 'p':  pf=1;
	    if (argv[i][2]!='\0')
	      sscanf(argv[i]+2, "%d", &do_backtrack);
	    break;
	  case 'n':
	    if ( strcmp(argv[i], "-noGU")==0) noGU=1;
	    if ( strcmp(argv[i], "-noCloseGU")==0) no_closingGU=1;
	    if ( strcmp(argv[i], "-noLP")==0) noLonelyPairs=1;
	    if ( strcmp(argv[i], "-nsp") ==0) {
	      if (i==argc-1) usage();
	      r=sscanf(argv[++i], "%32s", ns_bases);
	      if (!r) usage();
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
	    r=sscanf(argv[++i],"%f", &sfact);
	    if (!r) usage();
	    break;
	  case 'd': dangles=0;
	    if (strcmp(argv[i],"-d2")==0) dangles=2;
	    if (strcmp(argv[i],"-d1")==0) dangles=1;
	    break;
	  case 'P':
	    if (i==argc-1) usage();
	    r=sscanf(argv[++i], "%255s", ParamFile);
	    if (!r) usage();
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
	    sscanf(line, ">%12s", fname);
	 printf("%s\n", line);
	 free(line);
	 if ((line = get_line(stdin))==NULL) line = "@";
      } 

      if (strcmp(line, "@") == 0) break;

      string = (char *) space(strlen(line)+1);
      sscanf(line,"%s",string);
      free(line);
      length = strlen(string);

      if (fold_constrained) 
	 cstruc = get_line(stdin);
      structure = (char *) space(length+1);
      
      for (l = 0; l < length; l++) {
	string[l] = toupper(string[l]);
	if (!noconv && string[l] == 'T') string[l] = 'U';
      }
      if (istty)
	 printf("length = %d\n", length);

      initialize_fold(length);
      if (fold_constrained) 
	strncpy(structure, cstruc, length+1);
      min_en = fold(string, structure);
      printf("%s\n%s", string, structure);
      if (istty)
	 printf("\n minimum free energy = %6.2f kcal/mol\n", min_en);
      else
	 printf(" (%6.2f)\n", min_en);
       
      fflush(stdout);
       
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
	PS_rna_plot(string, structure, ffname);
      else 
	fprintf(stderr,"INFO: structure too long, not doing xy_plot\n");

      bp = base_pair;
      bpp= space(16);
      base_pair=bpp;
      free_arrays();
      base_pair = bp;
       
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

	 if (fold_constrained)
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
	    PS_dot_plot(string, ffname);
	 }
	 free_pf_arrays();

      }
      if (fold_constrained)
	free(cstruc);
      free(base_pair);
      fflush(stdout);
      free(string);
      free(structure); 
   } while (1);
   return 0;
}

PRIVATE void usage(void)
{
   nrerror("usage:\n"
	   "RNAfold [-p[0]] [-C] [-T temp] [-4] [-d[2]] [-noGU] [-noCloseGU]\n" 
	   "        [-noLP] [-e e_set] [-P paramfile] [-nsp pairs] [-S scale] "
	   "[-noconv]\n");
}
