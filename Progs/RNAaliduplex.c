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
#include "aln_util.h"

extern void  read_parameter_file(const char fname[]);
extern int subopt_sorted;
extern duplexT aliduplexfold(const char *s1[], const char *s2[]);
static void  print_struc(duplexT const *dup);
/*@unused@*/
static char rcsid[] = "$Id: RNAaliduplex.c,v 1.1 2007/08/26 10:08:44 ivo Exp $";

static char  scale[] = "....,....1....,....2....,....3....,....4"
                       "....,....5....,....6....,....7....,....8";

static void usage(void);

/*--------------------------------------------------------------------------*/
#define MAX_NUM_NAMES    500

int main(int argc, char *argv[])
{
  char *AS1[MAX_NUM_NAMES], *AS2[MAX_NUM_NAMES], 
    *names1[MAX_NUM_NAMES], *names2[MAX_NUM_NAMES];
  char  fname[13];
  char  *ParamFile=NULL;
  int   i, r, n_seq, n_seq2;
  double deltaf;
  double kT, sfact=1.07;
  int   pf=0, istty, delta=-1;
  int noconv=0;
  duplexT mfe, *subopt;
  FILE  *file1=NULL, *file2=NULL; /* input alignments */

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
    else { /* doesn't start with '-' should be filename */
      if (i>argc-1) usage();
      file1 = fopen(argv[i], "r");
      if (file1 == NULL) {
        fprintf(stderr, "can't open %s\n", argv[i]);
        usage();
      }
      file2 = fopen(argv[++i], "r");
      if (file1 == NULL) {
        fprintf(stderr, "can't open %s\n", argv[i]);
        usage();
      }
    }
  }

  if (!(file1 && file2)) {
    fprintf(stderr, "No input files");
    usage();
  }
  if (ParamFile != NULL)
    read_parameter_file(ParamFile);
   
  istty = isatty(fileno(stdout))&&isatty(fileno(stdin));
	
  n_seq = read_clustal(file1, AS1, names1);
  n_seq2 = read_clustal(file2, AS2, names2);
  fclose(file1); fclose(file2);

  if (n_seq != n_seq2) nrerror("unequal number of seqs in alignments");

  update_fold_params();

  if (delta>=0) {
    duplexT *sub;
    subopt = aliduplex_subopt(AS1, AS2, delta, 5);
    for (sub=subopt; sub->i >0; sub++) {
      print_struc(sub);
      free(sub->structure);
    }
    free(subopt);
  } else {
    mfe = aliduplexfold(AS1, AS2);
    print_struc(&mfe);
    free(mfe.structure);
  }
  for (i=0; AS1[i]; i++) {
    free(AS1[i]); free(names1[i]);
    free(AS2[i]); free(names2[i]);
  }
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
	  "RNAaliduplex [-e range] [-s] [standard options] file1.aln file2.aln\n");
}
