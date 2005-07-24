/* Last changed Time-stamp: <2005-07-23 16:50:24 ivo> */
/*                
	     Compute duplex structure of two RNA strands

			   c Ivo L Hofacker
			  Vienna RNA package
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include "duplex.h"
#include "fold.h"
#include "fold_vars.h"
#include "utils.h"

extern void  read_parameter_file(const char fname[]);
extern int subopt_sorted;
static void  print_struc(duplexT const *dup);
/*@unused@*/
static char rcsid[] = "$Id: RNAduplex.c,v 1.4 2005/07/24 08:35:15 ivo Exp $";

static char  scale[] = "....,....1....,....2....,....3....,....4"
                       "....,....5....,....6....,....7....,....8";

static void usage(void);

/*--------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
  char *s1, *s2, *line;
  char  fname[13];
  char  *ParamFile=NULL;
  char  *ns_bases=NULL, *c;
  int   i, l, sym, r;
  double deltaf;
  double kT, sfact=1.07;
  int   pf=0, istty, delta=-1;
  int noconv=0;
   
  for (i=1; i<argc; i++) {
    if (argv[i][0]=='-') 
      switch ( argv[i][1] )
	{
	case 'T':  if (argv[i][2]!='\0') usage();
	  if(i==argc-1) usage();
	  r=sscanf(argv[++i], "%lf", &temperature);
	  if (!r) usage();
	  break;
	case 'p':
	  fprintf(stderr, "partition function folding not yet implemented\n");
	  usage();
	  pf=1;
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
	  if (i>=argc-1) usage();
	  r=sscanf(argv[++i], "%lf", &deltaf);
	  if (r!=1) usage();
	  delta = (int) (0.1+deltaf*100);
	  break;
	case 's': subopt_sorted=1;
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
	
  do {				/* main loop: continue until end of file */
    duplexT mfe, *subopt;
    if (istty) {
      printf("\nInput two sequences (one line each); @ to quit\n");
      printf("%s\n", scale);
    }
    fname[0]='\0';
    
    if ((line = get_line(stdin))==NULL) break;
    /* skip empty lines, comment lines, name lines */
    while (line && ((*line=='*')||(*line=='\0')||(*line=='>'))) {
      printf("%s\n", line); free(line);
      if ((line = get_line(stdin))==NULL) break;
    } 
    if ((line ==NULL) || (strcmp(line, "@") == 0)) break;
    
    s1 = (char *) space(strlen(line)+1);
    (void) sscanf(line,"%s",s1);  free(line);
    
    if ((line = get_line(stdin))==NULL) break;
    /* skip comment lines and get filenames */
    while ((*line=='*')||(*line=='\0')||(*line=='>')) {
      printf("%s\n", line); free(line);
      if ((line = get_line(stdin))==NULL) break;
    } 
    if ((line ==NULL) || (strcmp(line, "@") == 0)) break;
    
    s2 = (char *) space(strlen(line)+1);
    (void) sscanf(line,"%s",s2); free(line);

    for (l = 0; l < strlen(s1); l++) {
      s1[l] = toupper(s1[l]);
      if (!noconv && s1[l] == 'T') s1[l] = 'U';
    }
    for (l = 0; l < strlen(s2); l++) {
      s2[l] = toupper(s2[l]);
      if (!noconv && s2[l] == 'T') s2[l] = 'U';
    }
    if (istty)
      printf("lengths = %d,%d\n", strlen(s1), strlen(s2));

    /* initialize_fold(length); */
    update_fold_params();
    if (delta>=0) {
      duplexT *sub;
      subopt = duplex_subopt(s1, s2, delta, 0);
      for (sub=subopt; sub->i >0; sub++) {
	print_struc(sub);
	free(sub->structure);
      }
      free(subopt);
    }
    else {
      mfe = duplexfold(s1, s2);
      print_struc(&mfe);
      free(mfe.structure);
    }
    (void) fflush(stdout);
       
    (void) fflush(stdout);
    free(s1); free(s2);
  } while (1);
  return 0;
}

static void print_struc(duplexT const *dup) {
  int l1;
  l1 = strchr(dup->structure, '&')-dup->structure;
  printf("%s %3d,%-3d : %3d,%-3d (%5.2f)\n", dup->structure, dup->i+1-l1,
	 dup->i, dup->j, dup->j+strlen(dup->structure)-l1-2, dup->energy);
}
    
static void usage(void)
{
  nrerror("usage:\n"
	  "RNAduplex [-e range] [-s]\n"
	  "          [-T temp] [-4] [-d] [-noGU] [-noCloseGU]\n" 
	  "          [-noLP] [-P paramfile] [-nsp pairs] [-noconv]\n");
}
