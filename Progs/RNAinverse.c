/*
	    Interactive access to inverse folding routines
		    c Ivo Hofacker, Peter Stadler
			  Vienna RNA Package
*/
/* Last changed Time-stamp: <1999-11-03 19:27:01 ivo> */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
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
static char rcsid[] = "$Id: RNAinverse.c,v 1.10 1999/11/04 12:16:37 ivo Exp $";
static char scale[] = "....,....1....,....2....,....3....,....4"
                      "....,....5....,....6....,....7....,....8";

#define  REPEAT_DEFAULT  100
#define  INFINITY        100000

extern void  read_parameter_file(const char fname[]);

PRIVATE void usage(void);
extern int inv_verbose;

int main(int argc, char *argv[])
{
  char *start, *structure, *rstart, *str2, *line;
  char  ParamFile[256]="";
  int   i,j, length, l, hd;
  float energy=0., kT;
  int   pf, mfe, istty;
  int   repeat, found;
  
  do_backtrack = 0; pf = 0; mfe = 1;
  repeat = 0;
  init_rand();
  strcpy(symbolset,"GCAU");
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

      /* read structure, skipping over comment lines */
      while ((*line=='*')||(*line=='\0')||(*line=='>')) {
	printf("%s\n", line);
	free(line);
	if ((line = get_line(stdin))==NULL) line = "@";
      } 
      /* stop at eof or '@' */
      if (strcmp(line, "@") == 0) break;

      structure = (char *) space(strlen(line)+1);
      sscanf(line,"%s",structure); /* scanf gets rid of trailing junk */
      free(line);
      
      length    = strlen(structure);
      str2 = (char *) space(length+1);

      if ((line = get_line(stdin))==NULL)
	line = "@";
      if (strcmp(line, "@") == 0) break;

      start = (char *) space(strlen(line)+1);
      sscanf(line,"%s",start);
      free(line);
      
      if(strlen(start)>length)
	nrerror("Sequence is longer than Structure");

      /* symbolset should only have uppercase characters */
      for (l = 0; l < strlen(symbolset); l++)
	symbolset[l] = toupper(symbolset[l]);
      if (istty) printf("length = %d\n", length);

      if (repeat!=0) found = (repeat>0)? repeat : (-repeat);
      else found = 1;
      
      initialize_fold(length);

      rstart = (char *) space(length+1);
      while(found>0) {
	char *string;
	string = (char *) space(length+1);
	strcpy(string, start);
	for (i=0; i<length; i++) {
	  /* lower case characters are kept fixed, any other character
	     not in symbolset is replaced by a random character */
	  if (islower(string[i])) continue;

	  if (string[i]=='\0' || (strchr(symbolset,string[i])==NULL))
	    string[i]=symbolset[int_urn(0,strlen(symbolset)-1)];
	}
	strcpy(rstart, string); /* remember start string */
	
	if (mfe) {
	  energy = inverse_fold(string, structure);
	  if( (repeat>=0) || (energy==0.0) ) {
	    found--;
	    hd = hamming(rstart, string);
	    printf("%s  %3d", string, hd);
	    if (energy>0.) { /* no solution found */
	      printf("   d= %g\n", energy);
	      if(istty) {
		energy = fold(string,str2);
		printf("%s\n", str2);
	      }
	    } else printf("\n");
	  }
	}
	if (pf) {
	  if (!(mfe && give_up && (energy>0))) {
	    /* unless we gave up in the mfe part */
	    float prob, min_en, sfact=1.07;
	    
	    /* get a reasonable pf_scale */
	    min_en = fold(string,str2); 
	    pf_scale = exp(-(sfact*min_en)/kT/length);
	    init_pf_fold(length);
	    
	    energy = inverse_pf_fold(string, structure);
	    prob = exp(-energy/kT);
	    hd = hamming(rstart, string);
	    printf("%s  %3d  (%g)\n", string, hd, prob);
	    free_pf_arrays();
	  }
	  if (!mfe) found--;
	}
	fflush(stdout);
	free(string);
      }
      free(rstart);
      free_arrays();
      
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
