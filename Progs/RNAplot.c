/*
  Plot RNA structures using different layout algorithms
  Last changed Time-stamp: <2003-09-10 13:55:01 ivo>
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
  char          *string, *input_string, *structure, *pre, *post;
  char          fname[FILENAME_MAX_LENGTH], ffname[FILENAME_MAX_LENGTH];
  int           i, r, fasta, length;
  float         energy;
  int           istty;
  char          format[5]="ps";
  unsigned int  input_type;

  string = input_string = structure = pre = post = NULL;
  fasta = length = 0;
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
  istty = isatty(fileno(stdin));

  /*
  #############################################
  # main loop: continue until end of file
  #############################################
  */
  do{
    if (istty) print_tty_input_seq();

    fname[0]='\0';

    input_type  = get_multi_input_line(&input_string, ((fasta) ? VRNA_INPUT_FASTA_HEADER : 0) | (istty ? VRNA_INPUT_NOSKIP_COMMENTS : 0));

    /* skip everything we are not interested in at this moment */
    while(input_type & (VRNA_INPUT_MISC | VRNA_INPUT_CONSTRAINT)){
      free(input_string); input_string  = NULL;
      /* get more input */
      input_type    = get_multi_input_line(&input_string, ((fasta) ? VRNA_INPUT_FASTA_HEADER : 0) | (istty ? VRNA_INPUT_NOSKIP_COMMENTS : 0));
    }
    if(input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR)) break;

    if(input_type & VRNA_INPUT_FASTA_HEADER){
      fasta = 1;
      (void) sscanf(input_string, ">%FILENAME_ID_LENGTHs", fname);
      if(!istty) printf("%s\n", input_string);
      free(input_string); input_string = NULL;
      input_type = get_multi_input_line(&input_string, VRNA_INPUT_FASTA_HEADER | (istty ? VRNA_INPUT_NOSKIP_COMMENTS : 0));
      if(input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR)) break;
    }
    else if(fasta) warn_user("fasta header missing");
    if(input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR)) break;

    if(input_type & VRNA_INPUT_SEQUENCE){
      length        = (int)strlen(input_string);
      string        = input_string;
      input_string  = NULL;
      input_type = get_multi_input_line(&input_string, ((fasta) ? VRNA_INPUT_FASTA_HEADER : 0) | (istty ? VRNA_INPUT_NOSKIP_COMMENTS : 0));
    }
    else nrerror("sequence missing");
    if(input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR)) break;

    if(input_type & VRNA_INPUT_CONSTRAINT){
      structure = (char *) space(strlen(input_string)+1);
      sscanf(input_string,"%s (%f)", structure, &energy);
      free(input_string); input_string  = NULL;
      if((int)strlen(structure) != length)  nrerror("structure and sequence differ in length!");
    }
    else nrerror("structure missing");

    if (fname[0]!='\0'){
      strcpy(ffname, fname);
      strcat(ffname, "_ss");
    } else
      strcpy(ffname, "rna");

    switch (format[0]) {
      case 'p':
        strcat(ffname, ".ps");
        PS_rna_plot_a(string, structure, ffname, pre, post);
        break;
      case 'g':
        strcat(ffname, ".gml");
        gmlRNA(string, structure, ffname, 'x');
        break;
      case 'x':
        strcat(ffname, ".ss");
        xrna_plot(string, structure, ffname);
        break;
      case 's':
        strcat(ffname, ".svg");
        svg_rna_plot(string, structure, ffname);
        break;
      default:
        RNAplot_cmdline_parser_print_help(); exit(EXIT_FAILURE);
     }
     fflush(stdout);
     free(string);
     free(structure);
   } while (1);
   return 0;
}
