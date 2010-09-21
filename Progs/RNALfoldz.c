/* Last changed Time-stamp: <2003-04-23 11:56:44 ivo> */
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
#include "fold_vars.h"
#include "utils.h"
#include "read_epars.h"
#include "Lfold.h"
#include "RNALfoldz_cmdl.h"

int main(int argc, char *argv[]){
  struct  RNALfoldz_args_info args_info;
  char                        *input_string, *c, *string, *structure, *ParamFile, *ns_bases;
  int                         i, length, l, sym, r, istty, noconv, maxdist, zsc;
  double                      energy, min_en, min_z;
  unsigned int                input_type;

  string = structure = ParamFile = ns_bases = NULL;
  do_backtrack  = 1;
  noconv        = 0;
  maxdist       = 150;
  zsc           = 1;
  min_z         = -2.0;

  /*
  #############################################
  # check the command line parameters
  #############################################
  */
  if(RNALfoldz_cmdline_parser (argc, argv, &args_info) != 0) exit(1);
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
  /* set energy model */
  if(args_info.energyModel_given) energy_set = args_info.energyModel_arg;
  /* take another energy parameter set */
  if(args_info.paramFile_given)   ParamFile = strdup(args_info.paramFile_arg);
  /* Allow other pairs in addition to the usual AU,GC,and GU pairs */
  if(args_info.nsp_given)         ns_bases = strdup(args_info.nsp_arg);
  /* set the maximum base pair span */
  if(args_info.span_given)        maxdist = args_info.span_arg;
  if(args_info.zscore_given)      min_z = args_info.zscore_arg;
  
  /* check for errorneous parameter options */
  if(maxdist < 0){
    RNALfoldz_cmdline_parser_print_help();
    exit(EXIT_FAILURE);
  }

  /* free allocated memory of command line data structure */
  RNALfoldz_cmdline_parser_free (&args_info);

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
    /*
    ########################################################
    # handle user input from 'stdin'
    ########################################################
    */
    if(istty) print_tty_input_seq();

    /* skip fasta header and comment lines */
    while((input_type = get_input_line(&input_string, 0)) & VRNA_INPUT_FASTA_HEADER){
      printf(">%s\n", input_string);
      free(input_string);
    }

    /* break on any error, EOF or quit request */
    if(input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR)){ break;}
    /* else assume a proper sequence of letters of a certain alphabet (RNA, DNA, etc.) */
    else{
      length = (int)    strlen(input_string);
      string = strdup(input_string);
      free(input_string);
    }

    if(noconv)  str_RNA2RNA(string);
    else        str_DNA2RNA(string);

    if(istty) printf("length = %d\n", length);
    /*
    ########################################################
    # done with 'stdin' handling
    ########################################################
    */

    min_en = (zsc) ? Lfoldz((const char *)string, NULL, maxdist, zsc, min_z) : Lfold((const char *)string, NULL, maxdist);
    printf("%s\n", string);

    if (istty)
      printf("\n minimum free energy = %6.2f kcal/mol\n", min_en);
    else
      printf(" (%6.2f)\n", min_en);

    (void) fflush(stdout);
    free(string);
  } while (1);
  return 0;
}

