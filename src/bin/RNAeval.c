/* Last changed Time-stamp: <2008-09-02 10:47:24 ivo> */
/*

          Calculate Energy of given Sequences and Structures
                           c Ivo L Hofacker
                          Vienna RNA Pckage
*/
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/data_structures.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/read_epars.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/file_formats.h"
#include "ViennaRNA/eval.h"
#include "ViennaRNA/cofold.h"
#include "RNAeval_cmdl.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

static void
add_shape_constraints(vrna_fold_compound *vc,
                      const char *shape_method,
                      const char *shape_conversion,
                      const char *shape_file,
                      int verbose,
                      unsigned int constraint_type){

  float p1, p2;
  char method;
  char *sequence;
  double *values;
  int length = vc->length;

  if(!vrna_sc_SHAPE_parse_method(shape_method, &method, &p1, &p2)){
    vrna_message_warning("Method for SHAPE reactivity data conversion not recognized!");
    return;
  }

  if(verbose){
    fprintf(stderr, "Using SHAPE method '%c'", method);
    if(method != 'W'){
      if(method == 'Z')
        fprintf(stderr, " with parameter p1=%f", p1);
      else
        fprintf(stderr, " with parameters p1=%f and p2=%f", p1, p2);
    }
    fputc('\n', stderr);
  }

  sequence = vrna_alloc(sizeof(char) * (length + 1));
  values = vrna_alloc(sizeof(double) * (length + 1));
  vrna_read_SHAPE_file(shape_file, length, method == 'W' ? 0 : -1, sequence, values);

  if(method == 'D'){
    (void)vrna_sc_add_SHAPE_deigan(vc, (const double *)values, p1, p2, constraint_type);
  }
  else if(method == 'Z'){
    (void)vrna_sc_add_SHAPE_zarringhalam(vc, (const double *)values, p1, 0.5, shape_conversion, constraint_type);
  } else {
    assert(method == 'W');
    vrna_sc_add_up(vc, values, constraint_type);
  }

  free(values);
  free(sequence);
}

int main(int argc, char *argv[]){
  struct RNAeval_args_info  args_info;
  char                      *string, *structure, *orig_sequence, *tmp;
  char                      *rec_sequence, *rec_id, **rec_rest;
  char                      *shape_file, *shape_method, *shape_conversion;
  char                      fname[FILENAME_MAX_LENGTH];
  char                      *ParamFile;
  int                       i, length1, length2, with_shapes;
  float                     energy;
  int                       istty;
  int                       circular=0;
  int                       noconv=0;
  int                       verbose = 0;
  unsigned int              rec_type, read_opt;
  vrna_md_t                 md;
  vrna_param_t              *P;

  string  = orig_sequence = ParamFile = NULL;
  gquad   = 0;
  dangles = 2;
  shape_file    = NULL;
  shape_method  = NULL;
  with_shapes   = 0;

  /* apply default model details */
  vrna_md_set_default(&md);

  /*
  #############################################
  # check the command line parameters
  #############################################
  */
  if(RNAeval_cmdline_parser (argc, argv, &args_info) != 0) exit(1);

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
  /* do not convert DNA nucleotide "T" to appropriate RNA "U" */
  if(args_info.noconv_given)      noconv = 1;
  /* set energy model */
  if(args_info.energyModel_given)
    md.energy_set = energy_set = args_info.energyModel_arg;
  /* take another energy parameter set */
  if(args_info.paramFile_given)   ParamFile = strdup(args_info.paramFile_arg);
  /* assume RNA sequence to be circular */
  if(args_info.circ_given)
    md.circ = circular = 1;
  /* logarithmic multiloop energies */
  if(args_info.logML_given)
    md.logML = logML = 1;
  /* be verbose */
  if(args_info.verbose_given)     verbose = 1;
  /* gquadruplex support */
  if(args_info.gquad_given)
    md.gquad = gquad = 1;

  if(args_info.shape_given){
    with_shapes = 1;
    shape_file = strdup(args_info.shape_arg);
  }

  shape_method = strdup(args_info.shapeMethod_arg);
  shape_conversion = strdup(args_info.shapeConversion_arg);
  if(args_info.verbose_given){
    verbose = 1;
  }

  /* free allocated memory of command line data structure */
  RNAeval_cmdline_parser_free (&args_info);

  /*
  #############################################
  # begin initializing
  #############################################
  */

  if (ParamFile!=NULL) read_parameter_file(ParamFile);

  rec_type      = read_opt = 0;
  rec_id        = rec_sequence = NULL;
  rec_rest      = NULL;
  istty         = isatty(fileno(stdout)) && isatty(fileno(stdin));

  if(circular && gquad){
    vrna_message_error("G-Quadruplex support is currently not available for circular RNA structures");
  }

  P = vrna_params_get(&md);

  /* set options we wanna pass to vrna_read_fasta_record() */

  if(istty){
    read_opt |= VRNA_INPUT_NOSKIP_BLANK_LINES;
    vrna_message_input_seq("Use '&' to connect 2 sequences that shall form a complex.\n"
                            "Input sequence (upper or lower case) followed by structure");
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
    if(rec_id){
      (void) sscanf(rec_id, ">%" XSTR(FILENAME_ID_LENGTH) "s", fname);
    }
    else fname[0] = '\0';

    /* convert DNA alphabet to RNA if not explicitely switched off */
    if(!noconv) vrna_seq_toRNA(rec_sequence);
    /* store case-unmodified sequence */
    orig_sequence = strdup(rec_sequence);
    /* convert sequence to uppercase letters only */
    vrna_seq_toupper(rec_sequence);

    vrna_fold_compound *vc = vrna_get_fold_compound(rec_sequence, &md, VRNA_OPTION_MFE | VRNA_OPTION_EVAL_ONLY);

    tmp       = vrna_extract_record_rest_structure((const char **)rec_rest, 0, (rec_id) ? VRNA_OPTION_MULTILINE : 0);

    if(!tmp)
      vrna_message_error("structure missing");

    int cp = -1;
    structure = vrna_cut_point_remove(tmp, &cp);
    if(cp != vc->cutpoint){
      fprintf(stderr,"cut_point = %d cut = %d\n", vc->cutpoint, cp);
      vrna_message_error("Sequence and Structure have different cut points.");
    }

    length1   = (int) strlen(structure);
    if(length1 != vc->length)
      vrna_message_error("structure and sequence differ in length!");

    free(tmp);

    if(with_shapes)
      add_shape_constraints(vc, shape_method, shape_conversion, shape_file, verbose, VRNA_CONSTRAINT_SOFT_MFE);

    if(istty){
      if (vc->cutpoint == -1)
        printf("length = %d\n", length1);
      else
        printf("length1 = %d\nlength2 = %d\n", vc->cutpoint-1, length1-vc->cutpoint+1);
    }

    /*
    ########################################################
    # begin actual computations
    ########################################################
    */

    if(verbose)
      energy = vrna_eval_structure_verbose(vc, structure, NULL);
    else
      energy = vrna_eval_structure(vc, structure);

    if (vc->cutpoint == -1)
      printf("%s\n%s", orig_sequence, structure);
    else {
      char *pstruct = vrna_cut_point_insert(structure, vc->cutpoint);
      printf("%s\n%s", orig_sequence,  pstruct);
      free(pstruct);
    }
    if (istty)
      printf("\n energy = %6.2f\n", energy);
    else
      printf(" (%6.2f)\n", energy);

    /* clean up */
    (void) fflush(stdout);
    if(rec_id) free(rec_id);
    free(rec_sequence);
    free(structure);
    /* free the rest of current dataset */
    if(rec_rest){
      for(i=0;rec_rest[i];i++) free(rec_rest[i]);
      free(rec_rest);
    }
    rec_id = rec_sequence = structure = NULL;
    rec_rest = NULL;

    free(string);
    free(orig_sequence);
    string = orig_sequence = NULL;

    vrna_free_fold_compound(vc);

    if(with_shapes)
      break;

    /* print user help for the next round if we get input from tty */
    if(istty){
      vrna_message_input_seq("Use '&' to connect 2 sequences that shall form a complex.\n"
                              "Input sequence (upper or lower case) followed by structure");
    }
  }
  return EXIT_SUCCESS;
}
