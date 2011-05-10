/*
            Interactive access to inverse folding routines
                    c Ivo Hofacker, Peter Stadler
                          Vienna RNA Package
*/
/* Last changed Time-stamp: <2000-09-28 14:58:05 ivo> */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include "inverse.h"
#include "fold_vars.h"
#include "fold.h"
#include "part_func.h"
#include "utils.h"
#include "read_epars.h"
#include "RNAinverse_cmdl.h"

#ifdef dmalloc
#include  "/usr/local/include/dmalloc.h"
#define space(X) calloc(1,(X))
#endif

/*@unused@*/
static char rcsid[] = "$Id: RNAinverse.c,v 1.11 2000/09/28 12:59:21 ivo Rel $";

#define  REPEAT_DEFAULT  100

PRIVATE void usage(void);
extern int inv_verbose;

int main(int argc, char *argv[]){
  struct  RNAinverse_args_info args_info;
  int     input_type;
  char    *input_string, *start, *structure, *rstart, *str2, *line;
  char    *ParamFile=NULL, *c, *ns_bases;
  int     i,j, length, l, hd, sym;
  double  energy=0., kT;
  int     pf, mfe, istty;
  int     repeat, found;

  do_backtrack = 0; pf = 0; mfe = 1;
  repeat = 0;
  input_type = 0;
  input_string = ns_bases = NULL;
  init_rand();

  /*
  #############################################
  # check the command line parameters
  #############################################
  */
  if(RNAinverse_cmdline_parser (argc, argv, &args_info) != 0) exit(1);
  /* temperature */
  if(args_info.temp_given)        temperature = args_info.temp_arg;
  /* do not take special tetra loop energies into account */
  if(args_info.noTetra_given)     tetra_loop=0;
  /* set dangle model */
  if(args_info.dangles_given)     dangles = args_info.dangles_arg;
  /* do not allow wobble pairs (GU) */
  if(args_info.noGU_given)        noGU = 1;
  /* do not allow weak closing pairs (AU,GU) */
  if(args_info.noClosingGU_given) no_closingGU = 1;
  /* set energy model */
  if(args_info.energyModel_given) energy_set = args_info.energyModel_arg;
  /* take another energy parameter set */
  if(args_info.paramFile_given)   ParamFile = strdup(args_info.paramFile_arg);
  /* Allow other pairs in addition to the usual AU,GC,and GU pairs */
  if(args_info.nsp_given)         ns_bases = strdup(args_info.nsp_arg);
  /* alter the alphabet */
  if(args_info.alphabet_given){
    symbolset=args_info.alphabet_arg;
    /* symbolset should only have uppercase characters */
    for (l = 0; l < (int)strlen(symbolset); l++)
      symbolset[l] = toupper(symbolset[l]);
  }
  /* set function for optimization */
  if(args_info.function_given){
    if(strlen(args_info.function_arg) > 2){
      RNAinverse_cmdline_parser_print_help(); exit(EXIT_FAILURE);
    }
    else{
      if((*args_info.function_arg == 'm') || (*(args_info.function_arg+1) == 'm')) mfe = 1;
      if((*args_info.function_arg == 'p') || (*(args_info.function_arg+1) == 'p')) pf = 1;
    }
  }
  /* set repeat */
  if(args_info.repeat_given)      repeat = args_info.repeat_arg;
  /* set final cost */
  if(args_info.final_given)       final_cost = args_info.final_arg;
  /* do we wannabe verbose */
  if(args_info.verbose_given)     inv_verbose = 1;

  /* free allocated memory of command line data structure */
  RNAinverse_cmdline_parser_free (&args_info);

  kT = (temperature+273.15)*1.98717/1000.0;

  istty = (isatty(fileno(stdout))&&isatty(fileno(stdin)));

  if (ParamFile!=NULL)
    read_parameter_file(ParamFile);

  give_up = (repeat<0);

  do {
    /*
    ########################################################
    # handle user input from 'stdin'
    ########################################################
    */
    if(istty)
      print_tty_input_seq_str("Input structure & start string\n"
                              "(lower case letters for const positions) and 0 or empty line for random start string\n");

    input_type = get_multi_input_line(&input_string, 0);
    /* we are waiting for a structure (i.e. something like a constraint) so we skip all sequences, fasta-headers and misc lines */
    while(input_type & (VRNA_INPUT_SEQUENCE | VRNA_INPUT_MISC | VRNA_INPUT_FASTA_HEADER)){
      if(!istty && (input_type & VRNA_INPUT_FASTA_HEADER)) printf(">%s\n", input_string);
      free(input_string); input_string = NULL;
      input_type = get_multi_input_line(&input_string, 0);
    }
    if(input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR)) break;

    if(input_type & (VRNA_INPUT_CONSTRAINT)){
      structure = (char *)space(sizeof(char) * (strlen(input_string) + 1));
      (void)sscanf(input_string, "%s", structure); /* scanf gets rid of trailing junk */
      length = (int)strlen(structure);
      free(input_string); input_string = NULL;
      input_type = get_multi_input_line(&input_string, VRNA_INPUT_NOSKIP_BLANK_LINES | VRNA_INPUT_NOSKIP_COMMENTS);
    }
    if(input_type & VRNA_INPUT_QUIT) break;

    start = (char *)space(sizeof(char) * (length+1));
    /* now we assume to get a sequence (input_string may be empty as well) */
    if(input_type & VRNA_INPUT_SEQUENCE){
      (void)strncpy(start, input_string, length);
      start[length] = '\0';
      free(input_string); input_string = NULL;
    }
    /* fallback to empty start sequence */
    else start[0] = '\0';

    /*
    ########################################################
    # done with 'stdin' handling
    ########################################################
    */

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

    str2 = (char *) space((unsigned)length+1);
    if (istty) printf("length = %d\n", length);

    if (repeat!=0) found = (repeat>0)? repeat : (-repeat);
    else found = 1;

    /* initialize_fold(length); <- obsolete (hopefully commenting this out does not affect anything crucial ;) */

    rstart = (char *) space((unsigned)length+1);
    while(found>0) {
      char *string;
      string = (char *) space((unsigned)length+1);
      strcpy(string, start);
      for (i=0; i<length; i++) {
        /* lower case characters are kept fixed, any other character
           not in symbolset is replaced by a random character */
        if (islower(string[i])) continue;

        if (string[i]=='\0' || (strchr(symbolset,string[i])==NULL))
          string[i]=symbolset[int_urn(0,strlen(symbolset)-1)];
      }
      strcpy(rstart, string); /* remember start string */

      if (mfe) {
        energy = inverse_fold(string, structure);
        if( (repeat>=0) || (energy<=0.0) ) {
          found--;
          hd = hamming(rstart, string);
          printf("%s  %3d", string, hd);
          if (energy>0) { /* no solution found */
            printf("   d= %g\n", energy);
            if(istty) {
              energy = fold(string,str2);
              printf("%s\n", str2);
            }
          } else printf("\n");
        }
      }
      if (pf) {
        if (!(mfe && give_up && (energy>0))) {
          /* unless we gave up in the mfe part */
          double prob, min_en, sfact=1.07;

          /* get a reasonable pf_scale */
          min_en = fold(string,str2);
          pf_scale = exp(-(sfact*min_en)/kT/length);
          /* init_pf_fold(length); <- obsolete (hopefully commenting this out does not affect anything crucial ;) */

          energy = inverse_pf_fold(string, structure);
          prob = exp(-energy/kT);
          hd = hamming(rstart, string);
          printf("%s  %3d  (%g)\n", string, hd, prob);
          free_pf_arrays();
        }
        if (!mfe) found--;
      }
      (void) fflush(stdout);
      free(string);
    }
    free(rstart);
    free_arrays();

    free(structure);
    free(str2);
    free(start);
    (void) fflush(stdout);
  } while (1);
  return 0;
}


PRIVATE void usage(void)
{
  nrerror("usage: RNAinverse [-F[mp]] [-a ALPHABET] [-R [repeats]] [-f final]\n"
          "                  [-T temp] [-4] [-d[2]] [-noGU] [-P paramfile] [-e e_set] [-v]");
}
