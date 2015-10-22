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
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/data_structures.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/read_epars.h"
#include "ViennaRNA/Lfold.h"
#include "ViennaRNA/file_formats.h"
#include "RNALfold_cmdl.h"

/*@unused@*/
static char rcsid[] = "$Id: RNALfold.c,v 1.2 2003/07/14 13:38:47 ivo Exp $";

int main(int argc, char *argv[]){
  FILE                        *input, *output;
  struct  RNALfold_args_info  args_info;
  char                        *input_string, *c, *string, *structure, *ParamFile, *ns_bases, *rec_sequence, *rec_id, **rec_rest, *orig_sequence;
  char                        fname[FILENAME_MAX_LENGTH], *infile, *outfile;
  int                         i, length, l, sym, r, istty, noconv, maxdist, zsc;
  double                      energy, min_en, min_z;
  unsigned int                input_type;
  unsigned int                rec_type, read_opt;
  vrna_md_t                   md;

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
  outfile       = NULL;
  infile        = NULL;
  input         = NULL;
  output        = NULL;

  /* apply default model details */
  vrna_md_set_default(&md);

  /*
  #############################################
  # check the command line parameters
  #############################################
  */
  if(RNALfold_cmdline_parser (argc, argv, &args_info) != 0) exit(1);
  /* temperature */
  if(args_info.temp_given)
    md.temperature = temperature = args_info.temp_arg;
  /* do not take special tetra loop energies into account */
  if(args_info.noTetra_given)
    md.special_hp = tetra_loop = 0;
  /* set dangle model */
  if(args_info.dangles_given){
    if((args_info.dangles_arg < 0) || (args_info.dangles_arg > 3))
      vrna_message_warning("required dangle model not implemented, falling back to default dangles=2");
    else
      md.dangles = dangles = args_info.dangles_arg;
  }
  /* do not allow weak pairs */
  if(args_info.noLP_given)
    md.noLP = noLonelyPairs = 1;
  /* do not allow wobble pairs (GU) */
  if(args_info.noGU_given)
    md.noGU = noGU = 1;
  /* do not allow weak closing pairs (AU,GU) */
  if(args_info.noClosingGU_given)
    md.noGUclosure = no_closingGU = 1;
  /* do not convert DNA nucleotide "T" to appropriate RNA "U" */
  if(args_info.noconv_given)
    noconv = 1;
  /* set energy model */
  if(args_info.energyModel_given)
    md.energy_set = energy_set = args_info.energyModel_arg;
  /* take another energy parameter set */
  if(args_info.paramFile_given)
    ParamFile = strdup(args_info.paramFile_arg);
  /* Allow other pairs in addition to the usual AU,GC,and GU pairs */
  if(args_info.nsp_given)
    ns_bases = strdup(args_info.nsp_arg);
  /* set the maximum base pair span */
  if(args_info.span_given)
    maxdist = args_info.span_arg;
  if(args_info.zscore_given){
#ifdef USE_SVM
    zsc = 1;
    if(args_info.zscore_arg != -2)
      min_z = args_info.zscore_arg;
#else
  vrna_message_error("\'z\' option is available only if compiled with SVM support!");
#endif
  }
  /* gquadruplex support */
  if(args_info.gquad_given)
    md.gquad = gquad = 1;

  if(args_info.outfile_given){
    outfile = strdup(args_info.outfile_arg);
  }
  if(args_info.infile_given){
    infile = strdup(args_info.infile_arg);
  }

  /* check for errorneous parameter options */
  if(maxdist <= 0){
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

  md.max_bp_span = md.window_size = maxdist;

  if(infile){
    input = fopen((const char *)infile, "r");
    if(!input)
      vrna_message_error("Could not read input file");
  }

  if (ParamFile != NULL)
    read_parameter_file(ParamFile);

  if (ns_bases != NULL) {
    /* nonstandards = vrna_alloc(33); */
    c=ns_bases;
    i=sym=0;
    if (*c=='-') {
      sym=1; c++;
    }
    while (*c!='\0') {
      if (*c!=',') {
        md.nonstandards[i++]=*c++;
        md.nonstandards[i++]=*c;
        if ((sym)&&(*c!=*(c-1))) {
          md.nonstandards[i++]=*c;
          md.nonstandards[i++]=*(c-1);
        }
      }
      c++;
    }
  }

  istty = (!infile) && isatty(fileno(stdout)) && isatty(fileno(stdin));
  read_opt |= VRNA_INPUT_NO_REST;
  if(istty){
    vrna_message_input_seq_simple();
    read_opt |= VRNA_INPUT_NOSKIP_BLANK_LINES;
  }
  /*
  #############################################
  # main loop: continue until end of file
  #############################################
  */
  while(
    !((rec_type = vrna_read_fasta_record(&rec_id, &rec_sequence, &rec_rest, NULL, read_opt))
        & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT))){

    /*
    ########################################################
    # init everything according to the data we've read
    ########################################################
    */
    char *prefix      = NULL;
    char *v_file_name = NULL;
    /*
    ########################################################
    # init everything according to the data we've read
    ########################################################
    */
    if(rec_id){
      if(!istty && !outfile) printf("%s\n", rec_id);
      (void) sscanf(rec_id, ">%" XSTR(FILENAME_ID_LENGTH) "s", fname);
    }
    else fname[0] = '\0';

    if(outfile){
      /* prepare the file prefix */
      if(fname[0] != '\0'){
        prefix = (char *)vrna_alloc(sizeof(char) * (strlen(fname) + strlen(outfile) + 1));
        strcpy(prefix, outfile);
        strcat(prefix, "_");
        strcat(prefix, fname);
      } else {
        prefix = (char *)vrna_alloc(sizeof(char) * (strlen(outfile) + 1));
        strcpy(prefix, outfile);
      }
    }

    length = (int)strlen(rec_sequence);

    /* convert DNA alphabet to RNA if not explicitely switched off */
    if(!noconv) vrna_seq_toRNA(rec_sequence);
    /* store case-unmodified sequence */
    orig_sequence = strdup(rec_sequence);
    /* convert sequence to uppercase letters only */
    vrna_seq_toupper(rec_sequence);

    if(!outfile && istty)
      printf("length = %d\n", length);

    /*
    ########################################################
    # done with 'stdin' handling
    ########################################################
    */
    vrna_fold_compound *vc = vrna_get_fold_compound((const char *)rec_sequence, &md, VRNA_OPTION_MFE | VRNA_OPTION_WINDOW);

    if(outfile){
      v_file_name = (char *)vrna_alloc(sizeof(char) * (strlen(prefix) + 8));
      strcpy(v_file_name, prefix);
      strcat(v_file_name, ".lfold");

      if(infile && !strcmp(infile, v_file_name))
        vrna_message_error("Input and output file names are identical");

      output = fopen((const char *)v_file_name, "a");
      if(!output)
        vrna_message_error("Failed to open file for writing");
    } else {
      output = stdout;
    }

    /*
    ########################################################
    # begin actual computations
    ########################################################
    */

    min_en = (zsc) ? vrna_Lfoldz(vc, min_z, output) : vrna_Lfold(vc, output);
    fprintf(output, "%s\n", orig_sequence);

    if(!outfile && istty)
      printf("\n minimum free energy = %6.2f kcal/mol\n", min_en);
    else
      fprintf(output, " (%6.2f)\n", min_en);

    if(output)
      (void) fflush(output);

    if(outfile && output){
      fclose(output);
      output = NULL;
    }

    /* clean up */
    vrna_free_fold_compound(vc);
    if(rec_id) free(rec_id);
    free(rec_sequence);
    free(orig_sequence);
    rec_id = rec_sequence = orig_sequence = NULL;
    rec_rest = NULL;
    /* print user help for the next round if we get input from tty */

    if(istty) vrna_message_input_seq_simple();
  }

  if(input)
    fclose(input);

  return EXIT_SUCCESS;
}

