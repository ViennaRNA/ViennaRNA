/*
 * Plot RNA structures using different layout algorithms
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>
#include "ViennaRNA/utils.h"
#include "ViennaRNA/PS_dot.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/file_formats.h"
#include "ViennaRNA/file_utils.h"
#include "RNAplot_cmdl.h"
#include "gengetopt_helper.h"
#include "input_id_helpers.h"

#define PRIVATE static

int
main(int  argc,
     char *argv[])
{
  struct RNAplot_args_info  args_info;
  char                      *structure, *pre, *post, *ffname, *tmp_string,
                            *filename_delim, *rec_sequence, *rec_id, **rec_rest,
                            format[5] = "ps";
  unsigned int              rec_type, read_opt;
  int                       i, istty, filename_full;
  vrna_md_t                 md;
  dataset_id                id_control;

  structure     = pre = post = NULL;
  filename_full = 0;
  vrna_md_set_default(&md);

  /*
   #############################################
   # check the command line parameters
   #############################################
   */
  if (RNAplot_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  /* parse options for ID manipulation */
  ggo_get_id_control(args_info, id_control, "Sequence", "sequence", "_", 4, 1);

  if (args_info.layout_type_given)
    rna_plot_type = args_info.layout_type_arg;

  if (args_info.pre_given)
    pre = strdup(args_info.pre_arg);

  if (args_info.post_given)
    post = strdup(args_info.post_arg);

  if (args_info.output_format_given) {
    strncpy(format, args_info.output_format_arg, 4);
    format[4] = '\0';
  }

  /* filename sanitize delimiter */
  if (args_info.filename_delim_given)
    filename_delim = strdup(args_info.filename_delim_arg);
  else if (get_id_delim(id_control))
    filename_delim = strdup(get_id_delim(id_control));

  if ((filename_delim) && isspace(*filename_delim)) {
    free(filename_delim);
    filename_delim = NULL;
  }

  /* full filename from FASTA header support */
  if (args_info.filename_full_given)
    filename_full = 1;

  /* free allocated memory of command line data structure */
  RNAplot_cmdline_parser_free(&args_info);

  /*
   #############################################
   # begin initializing
   #############################################
   */
  rec_type  = read_opt = 0;
  rec_id    = rec_sequence = NULL;
  rec_rest  = NULL;
  istty     = isatty(fileno(stdout)) && isatty(fileno(stdin));

  /* set options we wanna pass to vrna_file_fasta_read_record() */
  if (istty) {
    read_opt |= VRNA_INPUT_NOSKIP_BLANK_LINES;
    vrna_message_input_seq("Input sequence (upper or lower case) followed by structure");
  }

  /*
   #############################################
   # main loop: continue until end of file
   #############################################
   */
  while (
    !((rec_type = vrna_file_fasta_read_record(&rec_id, &rec_sequence, &rec_rest, NULL, read_opt))
      & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT))) {
    char  *SEQ_ID         = NULL;
    int   maybe_multiline = 0;

    if (rec_id) {
      maybe_multiline = 1;
      /* remove '>' from FASTA header */
      rec_id = memmove(rec_id, rec_id + 1, strlen(rec_id));
    }

    /* construct the sequence ID */
    set_next_id(&rec_id, id_control);
    SEQ_ID = fileprefix_from_id(rec_id, id_control, filename_full);

    structure = vrna_extract_record_rest_structure((const char **)rec_rest,
                                                   0,
                                                   (maybe_multiline) ? VRNA_OPTION_MULTILINE : 0);

    if (!structure)
      vrna_message_error("structure missing");

    if (strlen(rec_sequence) != strlen(structure))
      vrna_message_error("sequence and structure have unequal length");

    if (SEQ_ID)
      ffname = vrna_strdup_printf("%s%sss", SEQ_ID, filename_delim);
    else
      ffname = vrna_strdup_printf("rna");

    switch (format[0]) {
      case 'p':
        tmp_string = vrna_strdup_printf("%s.ps", ffname);
        free(ffname);
        ffname = vrna_filename_sanitize(tmp_string, filename_delim);
        free(tmp_string);

        (void)vrna_file_PS_rnaplot_a(rec_sequence, structure, ffname, pre, post, &md);

        break;
      case 'g':
        tmp_string = vrna_strdup_printf("%s.gml", ffname);
        free(ffname);
        ffname = vrna_filename_sanitize(tmp_string, filename_delim);
        free(tmp_string);

        gmlRNA(rec_sequence, structure, ffname, 'x');
        break;
      case 'x':
        tmp_string = vrna_strdup_printf("%s.ss", ffname);
        free(ffname);
        ffname = vrna_filename_sanitize(tmp_string, filename_delim);
        free(tmp_string);

        xrna_plot(rec_sequence, structure, ffname);
        break;
      case 's':
        tmp_string = vrna_strdup_printf("%s.svg", ffname);
        free(ffname);
        ffname = vrna_filename_sanitize(tmp_string, filename_delim);
        free(tmp_string);

        svg_rna_plot(rec_sequence, structure, ffname);
        break;
      default:
        RNAplot_cmdline_parser_print_help();
        exit(EXIT_FAILURE);
    }

    fflush(stdout);

    /* clean up */
    free(rec_id);
    free(rec_sequence);
    free(structure);
    structure = NULL;
    /* free the rest of current dataset */
    if (rec_rest) {
      for (i = 0; rec_rest[i]; i++)
        free(rec_rest[i]);
      free(rec_rest);
    }

    rec_id    = rec_sequence = structure = NULL;
    rec_rest  = NULL;

    free(SEQ_ID);
    free(ffname);
    ffname = NULL;

    /* print user help for the next round if we get input from tty */
    if (istty)
      vrna_message_input_seq("Input sequence (upper or lower case) followed by structure");
  }

  free(filename_delim);

  free_id_data(id_control);

  return EXIT_SUCCESS;
}
