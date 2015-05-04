/*
  Plot RNA structures using different layout algorithms
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>
#include "utils.h"
#include "PS_dot.h"
#include "RNAplot_cmdl.h"

#define PRIVATE static

int main(int argc, char *argv[]){
  struct        RNAplot_args_info args_info;
  char          *structure, *pre, *post;
  char          fname[FILENAME_MAX_LENGTH], ffname[FILENAME_MAX_LENGTH];
  char          *rec_sequence, *rec_id, **rec_rest;
  int           i, length;
  int           istty;
  char          format[5]="ps";
  unsigned int  rec_type, read_opt;

  structure = pre = post = NULL;
  length = 0;

  /*
  #############################################
  # check the command line parameters
  #############################################
  */
  if(RNAplot_cmdline_parser (argc, argv, &args_info) != 0) exit(1);

  if(args_info.layout_type_given) rna_plot_type = args_info.layout_type_arg;
  if(args_info.pre_given)         pre           = strdup(args_info.pre_arg);
  if(args_info.post_given)        post          = strdup(args_info.post_arg);
  if(args_info.output_format_given){
    strncpy(format, args_info.output_format_arg, 4); format[4] = '\0';
  }

  /* free allocated memory of command line data structure */
  RNAplot_cmdline_parser_free (&args_info);

  /*
  #############################################
  # begin initializing
  #############################################
  */
  rec_type      = read_opt = 0;
  rec_id        = rec_sequence = NULL;
  rec_rest      = NULL;
  istty         = isatty(fileno(stdout)) && isatty(fileno(stdin));

  /* set options we wanna pass to read_record */
  if(istty){
    read_opt |= VRNA_INPUT_NOSKIP_BLANK_LINES;
    print_tty_input_seq_str("Input sequence (upper or lower case) followed by structure");
  }

  /*
  #############################################
  # main loop: continue until end of file
  #############################################
  */
  while(
    !((rec_type = read_record(&rec_id, &rec_sequence, &rec_rest, read_opt))
        & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT))){

    if(rec_id){
      (void) sscanf(rec_id, ">%" XSTR(FILENAME_ID_LENGTH) "s", fname);
    }
    else fname[0] = '\0';

    length = (int)strlen(rec_sequence);

    structure = extract_record_rest_structure((const char **)rec_rest, 0, (rec_id) ? VRNA_OPTION_MULTILINE : 0);

    if(!structure) nrerror("structure missing");
    if((int)strlen(structure) != length) nrerror("structure and sequence differ in length!");

    if (fname[0]!='\0'){
      strcpy(ffname, fname);
      strcat(ffname, "_ss");
    } else
      strcpy(ffname, "rna");

    structure = NULL;
    unsigned int struct_options = (rec_id) ? VRNA_CONSTRAINT_MULTILINE : 0;
    struct_options |= VRNA_CONSTRAINT_ALL;
    getConstraint(&structure, (const char **)rec_rest, struct_options);

    if(strlen(rec_sequence) != strlen(structure))
      nrerror("sequence and structure have unequal length");

    switch (format[0]) {
      case 'p':
        strcat(ffname, ".ps");

        (void) PS_rna_plot_a_gquad(rec_sequence, structure, ffname, pre, post);

        /* PS_rna_plot_a(rec_sequence, structure, ffname, pre, post); */

        break;
      case 'g':
        strcat(ffname, ".gml");
        gmlRNA(rec_sequence, structure, ffname, 'x');
        break;
      case 'x':
        strcat(ffname, ".ss");
        xrna_plot(rec_sequence, structure, ffname);
        break;
      case 's':
        strcat(ffname, ".svg");
        svg_rna_plot(rec_sequence, structure, ffname);
        break;
      default:
        RNAplot_cmdline_parser_print_help(); exit(EXIT_FAILURE);
    }

    fflush(stdout);

    /* clean up */
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

    /* print user help for the next round if we get input from tty */
    if(istty){
      print_tty_input_seq_str("Input sequence (upper or lower case) followed by structure");
    }
  }

  return EXIT_SUCCESS;
}
