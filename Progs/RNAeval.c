/* Last changed Time-stamp: <2008-09-02 10:47:24 ivo> */
/*

          Calculate Energy of given Sequences and Structures
                           c Ivo L Hofacker
                          Vienna RNA Pckage
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include "fold_vars.h"
#include "fold.h"
#include "utils.h"
#include "read_epars.h"
#include "RNAeval_cmdl.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/*@unused@*/
static char UNUSED rcsid[]="$Id: RNAeval.c,v 1.10 2008/10/09 07:08:40 ivo Exp $";

PRIVATE char *costring(char *string);
PRIVATE char *tokenize(char *line);

int main(int argc, char *argv[]){
  struct RNAeval_args_info  args_info;
  unsigned int              input_type;
  char                      *line, *string, *structure, *input_string, *orig_sequence;
  char                      fname[80];
  char                      *ParamFile;
  int                       i, l, length1, length2;
  float                     energy;
  int                       istty;
  int                       circular=0;
  int                       noconv=0;
  int                       verbose = 0;

  string = orig_sequence = ParamFile = NULL;

  /*
  #############################################
  # check the command line parameters
  #############################################
  */
  if(RNAeval_cmdline_parser (argc, argv, &args_info) != 0) exit(1);
  /* temperature */
  if(args_info.temp_given)        temperature = args_info.temp_arg;
  /* do not take special tetra loop energies into account */
  if(args_info.noTetra_given)     tetra_loop=0;
  /* set dangle model */
  if(args_info.dangles_given)     dangles = args_info.dangles_arg;
  /* do not convert DNA nucleotide "T" to appropriate RNA "U" */
  if(args_info.noconv_given)      noconv = 1;
  /* set energy model */
  if(args_info.energyModel_given) energy_set = args_info.energyModel_arg;
  /* take another energy parameter set */
  if(args_info.paramFile_given)   ParamFile = strdup(args_info.paramFile_arg);
  /* assume RNA sequence to be circular */
  if(args_info.circ_given)        circular=1;
  /* logarithmic multiloop energies */
  if(args_info.logML_given)       logML = 1;
  /* be verbose */
  if(args_info.verbose_given)     verbose = 1;

  /* free allocated memory of command line data structure */
  RNAeval_cmdline_parser_free (&args_info);

  /*
  #############################################
  # begin initializing
  #############################################
  */

   if (ParamFile!=NULL) read_parameter_file(ParamFile);
   update_fold_params();

   istty = isatty(fileno(stdout))&&isatty(fileno(stdin));


  /*
  #############################################
  # main loop: continue until end of file
  #############################################
  */
  do {
    cut_point = -1;
    /*
    ########################################################
    # handle user input from 'stdin'
    ########################################################
    */
    if(istty){
      printf("Use '&' to connect 2 sequences that shall form a complex.\n");
      print_tty_input_seq_str("Input sequence & structure;   @ to quit");
    }

    /* extract filename from fasta header if available */
    fname[0] = '\0';
    while((input_type = get_input_line(&input_string, 0)) == VRNA_INPUT_FASTA_HEADER){
      printf(">%s\n", input_string);
      (void) sscanf(input_string, "%42s", fname);
      free(input_string);
    }

    /* break on any error, EOF or quit request */
    if(input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR)){ break;}
    /* else assume a proper sequence of letters of a certain alphabet (RNA, DNA, etc.) */
    else{
      string    = tokenize(input_string); /* frees input_string */
      length2   = (int) strlen(string);
      free(input_string);
    }

    /* get the structure */
    input_type = get_input_line(&input_string, 0);
    if(input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR)){ break;}
    else{
      structure = tokenize(input_string);
      length1   = (int) strlen(structure);
      free(input_string);
    }

    if(length1 != length2)
      nrerror("Sequence and Structure have unequal length.");

    /* convert DNA alphabet to RNA if not explicitely switched off */
    if(!noconv) str_DNA2RNA(string);
    /* store case-unmodified sequence */
    orig_sequence = strdup(string);
    /* convert sequence to uppercase letters only */
    str_uppercase(string);

    if(istty){
      if (cut_point == -1)
        printf("length = %d\n", length1);
      else
        printf("length1 = %d\nlength2 = %d\n", cut_point-1, length1-cut_point+1);
    }

    energy = (circular) ? energy_of_circ_structure(string, structure, verbose) : energy_of_structure(string, structure, verbose);

    if (cut_point == -1)
      printf("%s\n%s", orig_sequence, structure);
    else {
      char *pstring, *pstruct;
      pstring = costring(orig_sequence);
      pstruct = costring(structure);
      printf("%s\n%s", pstring,  pstruct);
      free(pstring);
      free(pstruct);
    }
    if (istty)
      printf("\n energy = %6.2f\n", energy);
    else
      printf(" (%6.2f)\n", energy);

    free(string);
    free(orig_sequence);
    free(structure);
    string = structure = orig_sequence = NULL;
    (void) fflush(stdout);
  } while (1);
  return 0;
}

PRIVATE char *tokenize(char *line)
{
  char *token, *copy, *ctmp;
  int cut = -1;

  copy = (char *) space(strlen(line)+1);
  ctmp = (char *) space(strlen(line)+1);
  (void) sscanf(line, "%s", copy);
  ctmp[0] = '\0';
  token = strtok(copy, "&");
  cut = strlen(token)+1;
  while (token) {
    strcat(ctmp, token);
    token = strtok(NULL, "&");
  }
  if (cut > strlen(ctmp)) cut = -1;
  if (cut > -1) {
    if (cut_point==-1) cut_point = cut;
    else if (cut_point != cut) {
      fprintf(stderr,"cut_point = %d cut = %d\n", cut_point, cut);
      nrerror("Sequence and Structure have different cut points.");
    }
  }
  free(copy);

  return ctmp;
}

PRIVATE char *costring(char *string)
{
  char *ctmp;
  int len;

  len = strlen(string);
  ctmp = (char *)space((len+2) * sizeof(char));
  /* first sequence */
  (void) strncpy(ctmp, string, cut_point-1);
  /* spacer */
  ctmp[cut_point-1] = '&';
  /* second sequence */
  (void) strcat(ctmp, string+cut_point-1);

  return ctmp;
}
