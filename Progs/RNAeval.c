/* Last changed Time-stamp: <2006-02-25 19:55:55 ivo> */
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
#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/*@unused@*/
static char UNUSED rcsid[]="$Id: RNAeval.c,v 1.9 2006/02/28 19:11:20 ivo Exp $";

#define  PUBLIC
#define  PRIVATE   static

static char  scale[] = "....,....1....,....2....,....3....,....4"
		       "....,....5....,....6....,....7....,....8";

PRIVATE char *costring(char *string);
PRIVATE char *tokenize(char *line);
PRIVATE void usage(void);
extern int logML;
extern int cut_point;
extern int eos_debug;
extern void  read_parameter_file(const char fname[]);
extern float energy_of_circ_struct(const char *seq, const char *str);

int main(int argc, char *argv[])
{
   char *line, *string, *structure;
   char  fname[12];
   char  *ParamFile=NULL;
   int   i, l, length1, length2;
   float energy;
   int   istty;
   int circ=0;
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
	 case 'c':
	   if ( strcmp(argv[i], "-circ")==0) circ=1;
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
	   else usage();
	   break;
	 case 'v':
	   if ( strcmp(argv[i], "-v")==0) eos_debug=1;
	   else usage();
	   break;
	 default: usage();
	 }
   }

   istty = isatty(fileno(stdout))&&isatty(fileno(stdin));

   if (ParamFile!=NULL)
     read_parameter_file(ParamFile);

   update_fold_params();

   do {
     cut_point = -1;
      if (istty) {
	 printf("\nInput sequence & structure;   @ to quit\n");
       printf("Use '&' to connect 2 sequences that shall form a complex.\n");
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

      string = tokenize(line);
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

      structure = tokenize(line);
      free(line);
      length1 = (int) strlen(structure);

      if(length1!=length2)
	 nrerror("Sequence and Structure have unequal length.");

      for (l = 0; l < length1; l++) {
	string[l] = toupper((int)string[l]);
	if (!noconv && string[l] == 'T') string[l] = 'U';
      }

      if (istty) {
	if (cut_point == -1)
	  printf("length = %d\n", length1);
	else
	  printf("length1 = %d\nlength2 = %d\n",
		 cut_point-1, length1-cut_point+1);
      }
      if (circ)
	energy = energy_of_circ_struct(string, structure);
      else
	energy = energy_of_struct(string, structure);
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
	printf("\n energy = %6.2f\n", energy);
      else
	printf(" (%6.2f)\n", energy);

      free(string);
      free(structure);
      (void) fflush(stdout);
   } while (1);
   return 0;
}

PRIVATE char *tokenize(char *line)
{
  char *token, *copy, *ctmp;
  int cut = -1;

  copy = (char *) space(strlen(line)+1);
  ctmp = (char *) space(strlen(line)+1);
  (void) sscanf(line, "%s", copy);
  ctmp[0] = '\0';
  token = strtok(copy, "&");
  cut = strlen(token)+1;
  while (token) {
    strcat(ctmp, token);
    token = strtok(NULL, "&");
  }
  if (cut > strlen(ctmp)) cut = -1;
  if (cut > -1) {
    if (cut_point==-1) cut_point = cut;
    else if (cut_point != cut) {
      fprintf(stderr,"cut_point = %d cut = %d\n", cut_point, cut);
      nrerror("Sequence and Structure have different cut points.");
    }
  }
  free(copy);

  return ctmp;
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
  nrerror("usage: RNAeval  [-T temp] [-4] [-d[0|1|2]] [-e e_set] [-logML] [-P paramfile]");
}
