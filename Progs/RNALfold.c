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
#include "RNALfold_cmdl.h"

/*@unused@*/
static char rcsid[] = "$Id: RNALfold.c,v 1.2 2003/07/14 13:38:47 ivo Exp $";

int main(int argc, char *argv[]){
  struct  RNALfold_args_info  args_info;
  char                        *input_string, *c, *string, *structure, *ParamFile, *ns_bases, *rec_sequence, *rec_id, **rec_rest, *orig_sequence;
  int                         i, length, l, sym, r, istty, noconv, maxdist, zsc;
  double                      energy, min_en, min_z;
  unsigned int                input_type;
  unsigned int                rec_type, read_opt;

  string        = structure = ParamFile = ns_bases = NULL;
  do_backtrack  = 1;
  noconv        = 0;
  dangles       = 2;
  maxdist       = 150;
  zsc           = 0;
  min_z         = -2.0;
  gquad         = 0;
  rec_type      = read_opt = 0;
  rec_id        = rec_sequence = orig_sequence = NULL;
  rec_rest      = NULL;

  /*
  #############################################
  # check the command line parameters
  #############################################
  */
  if(RNALfold_cmdline_parser (argc, argv, &args_info) != 0) exit(1);
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
  /* set energy model */
  if(args_info.energyModel_given) energy_set = args_info.energyModel_arg;
  /* take another energy parameter set */
  if(args_info.paramFile_given)   ParamFile = strdup(args_info.paramFile_arg);
  /* Allow other pairs in addition to the usual AU,GC,and GU pairs */
  if(args_info.nsp_given)         ns_bases = strdup(args_info.nsp_arg);
  /* set the maximum base pair span */
  if(args_info.span_given)        maxdist = args_info.span_arg;
  if(args_info.zscore_given){
#ifdef USE_SVM
    zsc = 1;
    if(args_info.zscore_arg != -2)
      min_z = args_info.zscore_arg;
#else
  nrerror("\'z\' option is available only if compiled with SVM support!");
#endif
  }
  /* gquadruplex support */
  if(args_info.gquad_given)      gquad = 1;

  /* check for errorneous parameter options */
  if(maxdist < 0){
    RNALfold_cmdline_parser_print_help();
    exit(EXIT_FAILURE);
  }

  /* free allocated memory of command line data structure */
  RNALfold_cmdline_parser_free (&args_info);

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
  read_opt |= VRNA_INPUT_NO_REST;
  if(istty){
    print_tty_input_seq();
    read_opt |= VRNA_INPUT_NOSKIP_BLANK_LINES;
  }
  /*
  #############################################
  # main loop: continue until end of file
  #############################################
  */
  while(
    !((rec_type = read_record(&rec_id, &rec_sequence, &rec_rest, read_opt))
        & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT))){

    /*
    ########################################################
    # init everything according to the data we've read
    ########################################################
    */
    if(rec_id && !istty) printf("%s\n", rec_id);

    length = (int)strlen(rec_sequence);

    /* convert DNA alphabet to RNA if not explicitely switched off */
    if(!noconv) str_DNA2RNA(rec_sequence);
    /* store case-unmodified sequence */
    orig_sequence = strdup(rec_sequence);
    /* convert sequence to uppercase letters only */
    str_uppercase(rec_sequence);

    if(istty) printf("length = %d\n", length);
    /*
    ########################################################
    # done with 'stdin' handling
    ########################################################
    */

    min_en = (zsc) ? Lfoldz((const char *)rec_sequence, NULL, maxdist, zsc, min_z) : Lfold((const char *)rec_sequence, NULL, maxdist);
    printf("%s\n", orig_sequence);

    if (istty)
      printf("\n minimum free energy = %6.2f kcal/mol\n", min_en);
    else
      printf(" (%6.2f)\n", min_en);

    (void) fflush(stdout);

    /* clean up */
    if(rec_id) free(rec_id);
    free(rec_sequence);
    free(orig_sequence);
    rec_id = rec_sequence = orig_sequence = NULL;
    rec_rest = NULL;
    /* print user help for the next round if we get input from tty */

    if(istty) print_tty_input_seq();
  }
  return EXIT_SUCCESS;
}

