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
  char                      *string, *structure, *orig_sequence, *tmp;
  char                      *rec_sequence, *rec_id, **rec_rest;
  char                      fname[FILENAME_MAX_LENGTH];
  char                      *ParamFile;
  int                       i, length1, length2;
  float                     energy;
  int                       istty;
  int                       circular=0;
  int                       noconv=0;
  int                       verbose = 0;
  unsigned int              rec_type, read_opt;

  string  = orig_sequence = ParamFile = NULL;
  gquad   = 0;
  dangles = 2;

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
  if(args_info.dangles_given){
    if((args_info.dangles_arg < 0) || (args_info.dangles_arg > 3))
      warn_user("required dangle model not implemented, falling back to default dangles=2");
    else
      dangles = args_info.dangles_arg;
  }
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
  /* gquadruplex support */
  if(args_info.gquad_given)       gquad = 1;

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
    nrerror("G-Quadruplex support is currently not available for circular RNA structures");
  }

  /* set options we wanna pass to read_record */
  if(istty){
    read_opt |= VRNA_INPUT_NOSKIP_BLANK_LINES;
    print_tty_input_seq_str("Use '&' to connect 2 sequences that shall form a complex.\n"
                            "Input sequence (upper or lower case) followed by structure");
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
      if(!istty) printf("%s\n", rec_id);
      (void) sscanf(rec_id, ">%" XSTR(FILENAME_ID_LENGTH) "s", fname);
    }
    else fname[0] = '\0';

    cut_point = -1;

    string    = tokenize(rec_sequence);
    length2   = (int) strlen(string);
    tmp       = extract_record_rest_structure((const char **)rec_rest, 0, (rec_id) ? VRNA_OPTION_MULTILINE : 0);

    if(!tmp)
      nrerror("structure missing");

    structure = tokenize(tmp);
    length1   = (int) strlen(structure);
    if(length1 != length2)
      nrerror("structure and sequence differ in length!");

    free(tmp);

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

    if(gquad)
      energy = energy_of_gquad_structure(string, structure, verbose);
    else
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

    /* print user help for the next round if we get input from tty */
    if(istty){
      print_tty_input_seq_str("Use '&' to connect 2 sequences that shall form a complex.\n"
                              "Input sequence (upper or lower case) followed by structure");
    }
  }
  return EXIT_SUCCESS;
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
