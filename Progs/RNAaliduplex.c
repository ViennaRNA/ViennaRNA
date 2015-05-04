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
#include "fold.h"
#include "fold_vars.h"
#include "utils.h"
#include "aln_util.h"
#include "read_epars.h"
#include "subopt.h"
#include "duplex.h"
#include "RNAaliduplex_cmdl.h"

/*@unused@*/
static char rcsid[] = "$Id: RNAaliduplex.c,v 1.1 2007/08/26 10:08:44 ivo Exp $";

PRIVATE void  print_struc(duplexT const *dup);

#define MAX_NUM_NAMES    500

int main(int argc, char *argv[]){
  struct        RNAaliduplex_args_info args_info;
  char          *AS1[MAX_NUM_NAMES], *AS2[MAX_NUM_NAMES], *names1[MAX_NUM_NAMES], *names2[MAX_NUM_NAMES];
  char          *ParamFile=NULL, *c, *ns_bases=NULL;
  int           i, r, n_seq, n_seq2;
  int           istty, delta=-1, sym;
  int           noconv=0;
  duplexT       mfe, *subopt;
  FILE          *file1=NULL, *file2=NULL; /* input alignments */


  /*
  #############################################
  # check the command line parameters
  #############################################
  */
  if(RNAaliduplex_cmdline_parser (argc, argv, &args_info) != 0) exit(1);
  /* temperature */
  if(args_info.temp_given)        temperature = args_info.temp_arg;
  /* do not take special tetra loop energies into account */
  if(args_info.noTetra_given)     tetra_loop=0;
  /* set dangle model */
  if(args_info.dangles_given)     dangles = args_info.dangles_arg;
  /* do not allow weak pairs */
  if(args_info.noLP_given)        noLonelyPairs = 1;
  /* do not allow wobble pairs (GU) */
  if(args_info.noGU_given)        noGU = 1;
  /* do not allow weak closing pairs (AU,GU) */
  if(args_info.noClosingGU_given) no_closingGU = 1;
  /* do not convert DNA nucleotide "T" to appropriate RNA "U" */
  if(args_info.noconv_given)      noconv = 1;
  /* take another energy parameter set */
  if(args_info.paramFile_given)   ParamFile = strdup(args_info.paramFile_arg);
  /* Allow other pairs in addition to the usual AU,GC,and GU pairs */
  if(args_info.nsp_given)         ns_bases = strdup(args_info.nsp_arg);
  /*energy range */
  if(args_info.deltaEnergy_given) delta = (int) (0.1+args_info.deltaEnergy_arg*100);
  /* sorted output */
  if(args_info.sorted_given)      subopt_sorted = 1;

  /* check unnamed options a.k.a. filenames of input alignments */
  if(args_info.inputs_num == 2){
    file1 = fopen(args_info.inputs[0], "r");
    if(file1 == NULL){
      fprintf(stderr, "can't open %s\n", args_info.inputs[0]);
    }
    file2 = fopen(args_info.inputs[1], "r");
    if(file2 == NULL) {
      fprintf(stderr, "can't open %s\n", args_info.inputs[1]);
    }
  }
  else{
    RNAaliduplex_cmdline_parser_print_help();
    exit(1);
  }

  if (!(file1 && file2)){
    fprintf(stderr, "No input files");
    RNAaliduplex_cmdline_parser_print_help();
    exit(1);
  }
  else{
    n_seq   = read_clustal(file1, AS1, names1);
    n_seq2  = read_clustal(file2, AS2, names2);
    fclose(file1);
    fclose(file2);
    if(n_seq != n_seq2) nrerror("unequal number of seqs in alignments");
  }

  /* free allocated memory of command line data structure */
  RNAaliduplex_cmdline_parser_free (&args_info);

  /*
  #############################################
  # begin initializing
  #############################################
  */


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

  /*
  #############################################
  # begin calculations
  #############################################
  */
  update_fold_params();

  if (delta>=0) {
    duplexT *sub;
    subopt = aliduplex_subopt((const char **)AS1, (const char **)AS2, delta, 5);
    for (sub=subopt; sub->i >0; sub++) {
      print_struc(sub);
      free(sub->structure);
    }
    free(subopt);
  } else {
    mfe = aliduplexfold((const char **)AS1, (const char **)AS2);
    print_struc(&mfe);
    free(mfe.structure);
  }
  for (i=0; AS1[i]; i++) {
    free(AS1[i]); free(names1[i]);
    free(AS2[i]); free(names2[i]);
  }
  return 0;
}

PRIVATE void print_struc(duplexT const *dup) {
  int l1;
  l1 = strchr(dup->structure, '&')-dup->structure;
  printf("%s %3d,%-3d : %3d,%-3d (%5.2f)\n", dup->structure, dup->i+1-l1,
         dup->i, dup->j, dup->j+strlen(dup->structure)-l1-2, dup->energy);
}
