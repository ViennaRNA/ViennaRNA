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
#include "read_epars.h"
#include "subopt.h"
#include "duplex.h"
#include "RNAduplex_cmdl.h"

PRIVATE void  print_struc(duplexT const *dup);
/*@unused@*/
static char rcsid[] = "$Id: RNAduplex.c,v 1.5 2007/08/26 09:41:12 ivo Exp $";

/*--------------------------------------------------------------------------*/

int main(int argc, char *argv[]){
  struct        RNAduplex_args_info args_info;
  char          fname[FILENAME_MAX_LENGTH], *input_string, *s1, *s2, *orig_s1, *orig_s2, *c, *ParamFile=NULL, *ns_bases=NULL;
  unsigned int  input_type;
  int           i, l, sym, r;
  int           istty, delta=-1;
  int           noconv=0;

  s1 = s2 = orig_s1 = orig_s2 = NULL;

  dangles = 2;

  /*
  #############################################
  # check the command line parameters
  #############################################
  */
  if(RNAduplex_cmdline_parser (argc, argv, &args_info) != 0) exit(1);
  /* temperature */
  if(args_info.temp_given)        temperature = args_info.temp_arg;
  /* do not take special tetra loop energies into account */
  if(args_info.noTetra_given)     tetra_loop=0;
  /* set dangle model */
  if(args_info.dangles_given){
    if((args_info.dangles_arg < 0) || (args_info.dangles_arg > 3))
      warn_user("required dangle model not implemented, falling back to default dangles=2");
    else
      dangles = args_info.dangles_arg;
  }
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

  /* free allocated memory of command line data structure */
  RNAduplex_cmdline_parser_free (&args_info);

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
  # main loop: continue until end of file
  #############################################
  */
  do{
    duplexT mfe, *subopt;
    /*
    ########################################################
    # handle user input from 'stdin'
    ########################################################
    */
    if(istty) print_tty_input_seq_str("Input two sequences (one line each)");

    /* extract filename from fasta header if available */
    fname[0] = '\0';
    while((input_type = get_input_line(&input_string, 0)) == VRNA_INPUT_FASTA_HEADER){
      printf(">%s\n", input_string);
      (void) sscanf(input_string, "%" XSTR(FILENAME_ID_LENGTH) "s", fname);
      free(input_string);
    }

    /* break on any error, EOF or quit request */
    if(input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR)){ break;}
    /* else assume a proper sequence of letters of a certain alphabet (RNA, DNA, etc.) */
    else{
      s1 = strdup(input_string);
      free(input_string);
    }

    /* get second sequence */
    while((input_type = get_input_line(&input_string, 0)) == VRNA_INPUT_FASTA_HEADER){
      printf(">%s\n", input_string);
      free(input_string);
    }
    /* break on any error, EOF or quit request */
    if(input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR)){ break;}
    /* else assume a proper sequence of letters of a certain alphabet (RNA, DNA, etc.) */
    else{
      s2 = strdup(input_string);
      free(input_string);
    }

    /* convert DNA alphabet to RNA if not explicitely switched off */
    if(!noconv){
      str_DNA2RNA(s1);
      str_DNA2RNA(s2);
    }
    /* store case-unmodified sequence */
    orig_s1 = strdup(s1);
    orig_s2 = strdup(s2);
    /* convert sequence to uppercase letters only */
    str_uppercase(s1);
    str_uppercase(s2);

    if (istty) printf("lengths = %d,%d\n", strlen(s1), strlen(s2));

    /*
    ########################################################
    # done with 'stdin' handling, now init everything properly
    ########################################################
    */
    update_fold_params();

     /*
    ########################################################
    # begin actual computations
    ########################################################
    */
    if (delta>=0) {
      duplexT *sub;
      subopt = duplex_subopt(s1, s2, delta, 5);
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
    free(s1);
    free(s2);
    free(orig_s1);
    free(orig_s2);
    s1 = s2 = orig_s1 = orig_s2 = NULL;
  } while (1);
  return 0;
}

PRIVATE void print_struc(duplexT const *dup) {
  int l1;
  l1 = strchr(dup->structure, '&')-dup->structure;
  printf("%s %3d,%-3d : %3d,%-3d (%5.2f)\n", dup->structure, dup->i+1-l1,
         dup->i, dup->j, dup->j+strlen(dup->structure)-l1-2, dup->energy);
}
