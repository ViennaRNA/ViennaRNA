/*
 ###################################
 # Access to RNAlib routines to    #
 # convert energy parameter files  #
 # from ViennaRNAPackage 1.8.4 to  #
 # 2.0 format                      #
 #                                 #
 # Ronny Lorenz                    #
 ###################################
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include "ViennaRNA/utils.h"
#include "ViennaRNA/convert_epars.h"
#include "RNAparconv_cmdl.h"

int
main(int  argc,
     char *argv[])
{
  struct  RNAparconv_args_info  args_info;
  char                          *ifileName, *ofileName;
  int                           i, silent;
  unsigned int                  options;

  ifileName = ofileName = NULL;
  silent    = 0;
  options   = VRNA_CONVERT_OUTPUT_ALL;

  /*
   #############################################
   # check the command line parameters
   #############################################
   */
  if (RNAparconv_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  /* use input file */
  if (args_info.input_given)
    ifileName = strdup(args_info.input_arg);

  /* use output file */
  if (args_info.output_given)
    ofileName = strdup(args_info.output_arg);

  /* output without (special) hairpin energies + enthalpies */
  if (args_info.without_HairpinE_given)
    options &= ~(VRNA_CONVERT_OUTPUT_HP
                 | VRNA_CONVERT_OUTPUT_SPECIAL_HP);

  /* output without stacking energies + enthalpies */
  if (args_info.without_StackE_given)
    options &= ~VRNA_CONVERT_OUTPUT_STACK;

  /* output without 5' dangle energies + enthalpies */
  if (args_info.without_Dangle5_given)
    options &= ~VRNA_CONVERT_OUTPUT_DANGLE5;

  /* output without 3' dangle energies + enthalpies */
  if (args_info.without_Dangle3_given)
    options &= ~VRNA_CONVERT_OUTPUT_DANGLE3;

  /* output without interior loop energies + enthalpies */
  if (args_info.without_IntE_given) {
    options &= ~(VRNA_CONVERT_OUTPUT_NINIO
                 | VRNA_CONVERT_OUTPUT_INT
                 | VRNA_CONVERT_OUTPUT_INT_11
                 | VRNA_CONVERT_OUTPUT_INT_21
                 | VRNA_CONVERT_OUTPUT_INT_22);
  }

  /* output without bulge loop energies + enthalpies */
  if (args_info.without_BulgeE_given)
    options &= ~VRNA_CONVERT_OUTPUT_BULGE;

  /* output without multi loop energies + enthalpies */
  if (args_info.without_MultiE_given)
    options &= ~VRNA_CONVERT_OUTPUT_ML;

  /* output without misc energies + enthalpies */
  if (args_info.without_Misc_given)
    options &= ~VRNA_CONVERT_OUTPUT_MISC;

  /* output without hairpin mismatch energies + enthalpies */
  if (args_info.without_MismatchH_given)
    options &= ~VRNA_CONVERT_OUTPUT_MM_HP;

  /* output without interior loop mismatch energies + enthalpies */
  if (args_info.without_MismatchI_given)
    options &= ~(VRNA_CONVERT_OUTPUT_MM_INT
                 | VRNA_CONVERT_OUTPUT_MM_INT_1N
                 | VRNA_CONVERT_OUTPUT_MM_INT_23);

  /* output without multi loop mismatch energies + enthalpies */
  if (args_info.without_MismatchM_given)
    options &= ~VRNA_CONVERT_OUTPUT_MM_MULTI;

  /* output without exterior loop mismatch energies + enthalpies */
  if (args_info.without_MismatchE_given)
    options &= ~VRNA_CONVERT_OUTPUT_MM_EXT;

  /* output given energy parameters only */
  if (args_info.vanilla_given)
    options |= VRNA_CONVERT_OUTPUT_VANILLA;

  /* just dump Vienna 1.8.4 parameters in new format */
  if (args_info.dump_given)
    options |= VRNA_CONVERT_OUTPUT_DUMP;

  /* file names given as unnamed option? */
  for (i = 0; i < args_info.inputs_num; i++) {
    if (ifileName == NULL) {
      ifileName = strdup(args_info.inputs[i]);
    } else if (ofileName == NULL) {
      ofileName = strdup(args_info.inputs[i]);
    } else {
      RNAparconv_cmdline_parser_print_help();
      vrna_message_error("unrecognized or too many parameter options given!");
    }
  }

  /* be silent */
  if (args_info.silent_given)
    silent = 1;

  /* free allocated memory of command line data structure */
  RNAparconv_cmdline_parser_free(&args_info);

  /*
   #############################################
   # do something
   #############################################
   */

  if (ifileName == NULL)
    if (!(silent || (options & VRNA_CONVERT_OUTPUT_DUMP)))
      fprintf(stdout, "Enter parameters:\n(Indicate last line of input parameters by typing \"# END\" or EOF\n\n");

  convert_parameter_file(ifileName, ofileName, options);

  return EXIT_SUCCESS;
}
