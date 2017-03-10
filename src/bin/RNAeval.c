/*
 *
 *        Calculate Energy of given Sequences and Structures
 *                         c Ivo L Hofacker
 *                        Vienna RNA Pckage
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <limits.h>
#include <string.h>
#include <sys/types.h>
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/data_structures.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/read_epars.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/constraints_SHAPE.h"
#include "ViennaRNA/file_formats.h"
#include "ViennaRNA/eval.h"
#include "ViennaRNA/cofold.h"
#include "RNAeval_cmdl.h"
#include "gengetopt_helper.h"
#include "input_id_helper.h"

#include "ViennaRNA/color_output.inc"

int
main(int  argc,
     char *argv[])
{
  struct RNAeval_args_info  args_info;
  char                      *string, *structure, *orig_sequence, *tmp, *rec_sequence,
                            *rec_id, **rec_rest, *shape_file, *shape_method, *id_delim,
                            *shape_conversion, *id_prefix;
  unsigned int              rec_type, read_opt;
  int                       i, length1, with_shapes, istty, noconv, verbose,
                            auto_id, id_digits, filename_full;
  long int                  seq_number;
  float                     energy;
  vrna_md_t                 md;

  string        = orig_sequence = NULL;
  noconv        = 0;
  verbose       = 0;
  gquad         = 0;
  dangles       = 2;
  seq_number    = 1;
  id_prefix     = NULL;
  auto_id       = 0;
  id_digits     = 4;
  filename_full = 0;

  /* apply default model details */
  vrna_md_set_default(&md);

  /*
   #############################################
   # check the command line parameters
   #############################################
   */
  if (RNAeval_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  /* get basic set of model details */
  ggo_get_md_eval(args_info, md);
  ggo_get_circ(args_info, md.circ);

  /* check dangle model */
  if ((md.dangles < 0) || (md.dangles > 3)) {
    vrna_message_warning("required dangle model not implemented, falling back to default dangles=2");
    md.dangles = dangles = 2;
  }

  /* SHAPE reactivity data */
  ggo_get_SHAPE(args_info, with_shapes, shape_file, shape_method, shape_conversion);

  ggo_get_ID_manipulation(args_info,
                          auto_id,
                          id_prefix, "sequence",
                          id_delim, "_",
                          id_digits, 4,
                          seq_number, 1);

  /* do not convert DNA nucleotide "T" to appropriate RNA "U" */
  if (args_info.noconv_given)
    noconv = 1;

  /* logarithmic multiloop energies */
  if (args_info.logML_given)
    md.logML = logML = 1;

  /* be verbose */
  if (args_info.verbose_given)
    verbose = 1;

  /* free allocated memory of command line data structure */
  RNAeval_cmdline_parser_free(&args_info);

  /*
   #############################################
   # begin initializing
   #############################################
   */

  rec_type  = read_opt = 0;
  rec_id    = rec_sequence = NULL;
  rec_rest  = NULL;
  istty     = isatty(fileno(stdout)) && isatty(fileno(stdin));

  if (md.circ && md.gquad)
    vrna_message_error("G-Quadruplex support is currently not available for circular RNA structures");

  /* set options we wanna pass to vrna_file_fasta_read_record() */

  if (istty) {
    read_opt |= VRNA_INPUT_NOSKIP_BLANK_LINES;
    vrna_message_input_seq("Use '&' to connect 2 sequences that shall form a complex.\n"
                           "Input sequence (upper or lower case) followed by structure");
  }

  /*
   #############################################
   # main loop: continue until end of file
   #############################################
   */
  while (
    !((rec_type = vrna_file_fasta_read_record(&rec_id, &rec_sequence, &rec_rest, NULL, read_opt))
      & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT))) {
    /*
     ########################################################
     # init everything according to the data we've read
     ########################################################
     */
    char  *SEQ_ID         = NULL, *msg = NULL;
    int   maybe_multiline = 0;

    if (rec_id) {
      maybe_multiline = 1;
      /* remove '>' from FASTA header */
      rec_id = memmove(rec_id, rec_id + 1, strlen(rec_id));
    }

    /* construct the sequence ID */
    ID_generate(SEQ_ID, rec_id, auto_id, id_prefix, id_delim, id_digits, seq_number, filename_full);

    /* convert DNA alphabet to RNA if not explicitely switched off */
    if (!noconv)
      vrna_seq_toRNA(rec_sequence);

    /* store case-unmodified sequence */
    orig_sequence = strdup(rec_sequence);
    /* convert sequence to uppercase letters only */
    vrna_seq_toupper(rec_sequence);

    vrna_fold_compound_t *vc = vrna_fold_compound(rec_sequence, &md, VRNA_OPTION_MFE | VRNA_OPTION_EVAL_ONLY);

    tmp = vrna_extract_record_rest_structure((const char **)rec_rest, 0, (maybe_multiline) ? VRNA_OPTION_MULTILINE : 0);

    if (!tmp)
      vrna_message_error("structure missing");

    int cp = -1;
    structure = vrna_cut_point_remove(tmp, &cp);
    if (cp != vc->cutpoint) {
      vrna_message_warning("cut_point = %d cut = %d", vc->cutpoint, cp);
      vrna_message_error("Sequence and Structure have different cut points.");
    }

    length1 = (int)strlen(structure);
    if (length1 != vc->length)
      vrna_message_error("structure and sequence differ in length!");

    free(tmp);

    if (with_shapes)
      vrna_constraints_add_SHAPE(vc, shape_file, shape_method, shape_conversion, verbose, VRNA_OPTION_MFE);

    if (istty) {
      if (vc->cutpoint == -1)
        vrna_message_info(stdout, "length = %d", length1);
      else
        vrna_message_info(stdout, "length1 = %d\nlength2 = %d", vc->cutpoint - 1, length1 - vc->cutpoint + 1);
    }

    /*
     ########################################################
     # begin actual computations
     ########################################################
     */

    print_fasta_header(stdout, rec_id);

    energy = vrna_eval_structure_v(vc, structure, verbose, NULL);

    fprintf(stdout, "%s\n", orig_sequence);

    if (istty)
      msg = vrna_strdup_printf("\n energy = %6.2f kcal/mol", energy);
    else
      msg = vrna_strdup_printf(" (%6.2f)", energy);

    if (vc->cutpoint == -1) {
      print_structure(stdout, structure, msg);
    } else {
      char *pstruct = vrna_cut_point_insert(structure, vc->cutpoint);
      print_structure(stdout, pstruct, msg);
      free(pstruct);
    }

    (void)fflush(stdout);

    /* clean up */
    free(rec_id);
    free(SEQ_ID);
    free(msg);
    free(rec_sequence);
    free(structure);
    /* free the rest of current dataset */
    if (rec_rest) {
      for (i = 0; rec_rest[i]; i++)
        free(rec_rest[i]);
      free(rec_rest);
    }

    rec_id    = rec_sequence = structure = NULL;
    rec_rest  = NULL;

    free(string);
    free(orig_sequence);
    string = orig_sequence = NULL;

    vrna_fold_compound_free(vc);

    if (with_shapes)
      break;

    ID_number_increase(seq_number, "Sequence");

    /* print user help for the next round if we get input from tty */
    if (istty)
      vrna_message_input_seq("Use '&' to connect 2 sequences that shall form a complex.\n"
                             "Input sequence (upper or lower case) followed by structure");
  }

  free(id_prefix);
  free(id_delim);

  return EXIT_SUCCESS;
}
