/*
	    Interactive access to inverse folding routines
		    c Ivo Hofacker, Peter Stadler
			  Vienna RNA Package
*/
/* Last changed Time-stamp: <1998-05-12 15:12:53 ivo> */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "inverse.h"
#include "fold_vars.h"
#include "fold.h"
#include "part_func.h"
#include "utils.h"
#ifdef dmalloc
#include  "/usr/local/include/dmalloc.h"
#define space(X) calloc(1,(X))
#endif

#define  PUBLIC
#define  PRIVATE   static
static char rcsid[] = "$Id: RNAinverse.c,v 1.9 1998/05/19 17:44:29 ivo Exp $";
static char scale[] = "....,....1....,....2....,....3....,....4"
                      "....,....5....,....6....,....7....,....8";

#define  REPEAT_DEFAULT  100
#define  INFINITY        100000

extern void  read_parameter_file(const char fname[]);

PRIVATE void usage(void);
extern int inv_verbose;

int main(int argc, char *argv[])
{
   char *string, *start, *structure, *str2, *line;
   char  ParamFile[256]="";
   int   i,j,length1, length2, l, hd;
   float energy, kT;
   int   pf, mfe, istty, rstart=0;
   int   repeat, found;

   do_backtrack = 0; pf = 0; mfe = 1;
   repeat = 0;
   init_rand();
   strcpy(symbolset,"GCAU");
   string=NULL;
   for (i=1; i<argc; i++) {
      if (argv[i][0]=='-')
	 switch ( argv[i][1] )
	   {
	   case 'a':
	     i++;
	     strncpy(symbolset,argv[i],20);
	     break;
	   case 'T':  if (argv[i][2]!='\0') usage(); 
	     if (sscanf(argv[++i], "%f", &temperature)==0)
	       usage();
	     break;
	   case 'F':
	     mfe = 0; pf = 0;
	     for(j=2;j<strlen(argv[i]);j++){
	       switch( argv[i][j] ) {
	       case 'm' :  mfe = 1;
		 break;
	       case 'p' :  pf = 1; /* old version had dangles=0 here */
		 break;
	       default : usage();
	       }
	     }
	     break;
	   case 'R': repeat = REPEAT_DEFAULT;
	     if(++i<argc) sscanf(argv[i], "%d", &repeat);
	     break;
	   case 'n':
	     if ( strcmp(argv[i], "-noGU" )==0) noGU=1;
	     break;
	   case '4':
	     tetra_loop=0;
	     break;
	   case 'e':
	     if (sscanf(argv[++i],"%d", &energy_set)==0)
	       usage();
	     break;
	   case 'd': dangles=0;
	     if (strcmp(argv[i],"-d2")==0) dangles=2;
	     if (strcmp(argv[i],"-d1")==0) dangles=1;
	     break;
	   case 'f': /* when to stop RNAfold -p */
	     if (sscanf(argv[++i],"%f", &final_cost)==0)
	       usage(); 
	     break;
	   case 'P':
	     if (sscanf(argv[++i], "%255s", ParamFile)==0)
	       usage();
	     break;
	   case 'v':
	     inv_verbose = 1;
	     break;
	   default: usage();
	   }
   }

   kT = (temperature+273.15)*1.98717/1000.0;

   istty = (isatty(fileno(stdout))&&isatty(fileno(stdin)));

   if (ParamFile[0])
     read_parameter_file(ParamFile);

   give_up = (repeat<0);

   do {
      if (istty) {
	 printf("\nInput structure & start string"
		" (lower case letters for const positions)\n"
		"    @ to quit, and 0 for random start string\n");
	 printf("%s\n", scale);
      }
      
      if ((line = get_line(stdin))==NULL) break;

      /* skip comment lines and get filenames */
      while ((*line=='*')||(*line=='\0')||(*line=='>')) {
	 printf("%s\n", line);
	 free(line);
	 if ((line = get_line(stdin))==NULL) line = "@";
      } 

      if (strcmp(line, "@") == 0) break;

      structure = (char *) space(strlen(line)+1);
      sscanf(line,"%s",structure);
      free(line);
      
      length2    = strlen(structure);
      str2 = (char *) space(length2+1);

      if ((line = get_line(stdin))==NULL) break;

      /* skip comment lines and get filenames */
      while ((*line=='*')||(*line=='\0')||(*line=='>')) {
	 printf("%s\n", line);
	 free(line);
	 if ((line = get_line(stdin))==NULL) line = "@";
      } 
      
      if (strcmp(line, "@") == 0) break;

      string = (char *) space(strlen(line)+1);
      sscanf(line,"%s",string);
      free(line);
      length1 = strlen(string);

      if (string[0]=='0') {
	 rstart = 1;
	 length1 = length2;
      }

      if(length1!=length2)
	 nrerror("Sequence and Structure have unequal length.");

      for (l = 0; l < strlen(symbolset); l++)
	 symbolset[l] = toupper(symbolset[l]);
      if (istty) printf("length = %d\n", length1);

      start = (char *) space(sizeof(char)*(length2+1));
      if (!rstart) strcpy(start, string);

      if (repeat!=0) found = (repeat>0)? repeat : (-repeat);
      else found = 1;

      initialize_fold(length2);

      while(found>0){
	 if (rstart) {
	    free(string);
	    string = random_string(length2, symbolset);
	    strcpy(start, string);
	 } else
	    strcpy(string, start);
	 
	 if (mfe) {
            energy = inverse_fold(string, structure);
            if( (repeat>=0) || (energy==0.0) ) {
	       found--;
               hd = hamming(start, string);
	       printf("%s  %3d", string, hd);
	       if (energy>0.) {
	          printf("   d= %g\n", energy);
                  if(istty) {
                     energy = fold(string,str2);
	             printf("%s\n", str2);
                  }
               } else printf("\n");
            }
	 }
	 if (pf) {
	    if( (!mfe) || (repeat >= 0) || (energy==0.) ) {
	       float prob, min_en, sfact=1.07;

	       /* get a reasonable pf_scale */
	       min_en = fold(string,str2); 
	       pf_scale = exp(-(sfact*min_en)/kT/length2);
	       init_pf_fold(length2);

	       
	       if (dangles) dangles = 2; /* for energy_of_struct */
	       energy = inverse_pf_fold(string, structure);
	       prob = exp(-energy/kT);
	       hd = hamming(start, string);
	       printf("%s  %3d  (%g)\n", string, hd, prob);
	       free_pf_arrays();
	    }
	    if (!mfe) found--;
	 }
	 fflush(stdout);
      }
      free_arrays();

      free(string);
      free(structure);
      free(str2);
      free(start);
      fflush(stdout);
   } while (1);
   return 0;
}


PRIVATE void usage(void)
{
  nrerror("usage: RNAinverse [-F[mp]] [-a ALPHABET] [-R [repeats]] [-f final]\n"
          "                  [-T temp] [-4] [-d[2]] [-noGU] [-P paramfile] [-e e_set] [-v]");
}
