/* Last changed Time-stamp: <2000-09-28 11:20:30 ivo> */
/*

	  Calculate Energy of given Sequences and Structures
			   c Ivo L Hofacker
			  Vienna RNA Pckage
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include "fold_vars.h"
#include "fold.h"
#include "utils.h"

/*@unused@*/
static char rcsid[] = "$Id: RNAeval.c,v 1.7 2000/09/28 09:30:29 ivo Rel $";

#define  PUBLIC
#define  PRIVATE   static

static char  scale[] = "....,....1....,....2....,....3....,....4"
                       "....,....5....,....6....,....7....,....8";

PRIVATE void usage(void);
extern int logML;
extern void  read_parameter_file(const char fname[]);

int main(int argc, char *argv[])
{
   char *line, *string, *structure;
   char  fname[12];
   char  *ParamFile=NULL;
   int   i, l, length1, length2;
   float energy;
   int   istty;
   int   noconv=0;
      
   string=NULL;
   for (i=1; i<argc; i++) {
     if (argv[i][0]=='-')
       switch ( argv[i][1] )
	 {
	 case 'T':  if (argv[i][2]!='\0') usage();
	   if (sscanf(argv[++i], "%lf", &temperature)==0)
	     usage();
	   break;
	 case '4':
	   tetra_loop=0;
	   break;
	 case 'e':
	   if (sscanf(argv[++i],"%d", &energy_set)==0)
	     usage();
	   break;
	 case 'd': dangles=0;
	   if (argv[i][2]!='\0')
	     if (sscanf(argv[i]+2, "%d", &dangles)==0)
	       usage();
	   break;
	 case 'P':
	   if (++i <= argc) 
	     ParamFile = argv[i];
	   else usage();
	   break;
	 case 'l':
	   if (strcmp(argv[i],"-logML")==0) 
	     logML=1;
	   else usage();
	   break;
	 case 'n':
	   if ( strcmp(argv[i], "-noconv")==0) noconv=1;
	   break;
	 default: usage();
	 }
   }

   istty = isatty(fileno(stdout))&&isatty(fileno(stdin));

   if (ParamFile!=NULL)
     read_parameter_file(ParamFile);

   update_fold_params();

   do {
      if (istty) {
	 printf("\nInput sequence & structure;   @ to quit\n");
	 printf("%s\n", scale);
      }

      fname[0]='\0';
      if ((line = get_line(stdin))==NULL) break;

     /* skip comment lines and get filenames */
      while ((*line=='*')||(*line=='\0')||(*line=='>')) {
	 if (*line=='>')
	    (void) sscanf(line, ">%s", fname);
	 printf("%s\n", line);
	 free(line);
	 if ((line = get_line(stdin))==NULL) break;
      }  

      if (line==NULL) break;
      if (strcmp(line, "@") == 0) {free(line); break;}

      string = (char *) space(strlen(line)+1);
      (void) sscanf(line,"%s",string);
      free(line);
      length2 = (int) strlen(string);

      if ((line = get_line(stdin))==NULL) break;
      /* skip comment lines */
      while ((*line=='*')||(*line=='\0')||(*line=='>')) {
	 printf("%s\n", line);
	 free(line);
	 if ((line = get_line(stdin))==NULL) break;
      }  
      if (line==NULL) break;
      if (strcmp(line, "@") == 0) {free(line); break;}
      
      structure = (char *) space(strlen(line)+1);
      (void) sscanf(line, "%s", structure);
      free(line);
      length1 = (int) strlen(structure);
      
      if(length1!=length2)
	 nrerror("Sequence and Structure have unequal length.");

      for (l = 0; l < length1; l++) {
        string[l] = toupper((int)string[l]);
        if (!noconv && string[l] == 'T') string[l] = 'U';
      }

      if (istty) printf("length = %d\n", length1);
      
      energy = energy_of_struct(string, structure);
      printf("%s\n%s", string, structure);
      if (istty)
	 printf("\n energy = %6.2f\n", energy);
      else
	 printf(" (%6.2f)\n", energy);
      
      free(string);
      free(structure);
      (void) fflush(stdout);
   } while (1);
   return 0;
}


PRIVATE void usage(void)
{
  nrerror("usage: RNAeval  [-T temp] [-4] [-d[0|1|2]] [-e e_set] [-logML] [-P paramfile]");
}
