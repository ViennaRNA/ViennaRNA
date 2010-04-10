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
#include "fold.h"
#include "part_func.h"
#include "fold_vars.h"
#include "utils.h"
#include "read_epars.h"
#include "Lfold.h"
#include "RNALfold_cmdl.h"

/*@unused@*/
static char rcsid[] = "$Id: RNALfold.c,v 1.2 2003/07/14 13:38:47 ivo Exp $";

int main(int argc, char *argv[]){
  struct  RNALfold_args_info args_info;
  char    *line, *c, *string=NULL, *structure=NULL, *ParamFile=NULL, *ns_bases=NULL;
  int     i, length, l, sym, r, istty, noconv=0, maxdist=150;
  double  energy, min_en;
   
  do_backtrack  = 1;

  /*
  #############################################
  # check the command line prameters
  #############################################
  */
  if(RNALfold_cmdline_parser (argc, argv, &args_info) != 0) exit(1);
  /* temperature */
  if(args_info.temp_given)        temperature = args_info.temp_arg;
  /* do not take special tetra loop energies into account */
  if(args_info.noTetra_given)     tetra_loop=0;
  /* set dangle model */
  if(args_info.dangles_given)     dangles = args_info.dangles_arg;
  /* do not allow weak pairs */
  if(args_info.noLP_given)        noLonelyPairs = 1;
  /* do not allow wobble pairs (GU) */
  if(args_info.noGU_given)        noGU = 1;
  /* do not allow weak closing pairs (AU,GU) */
  if(args_info.noClosingGU_given) no_closingGU = 1;
  /* do not convert DNA nucleotide "T" to appropriate RNA "U" */
  if(args_info.noconv_given)      noconv = 1;
  /* set energy model */
  if(args_info.energyModel_given) energy_set = args_info.energyModel_arg;
  /* take another energy parameter set */
  if(args_info.paramFile_given)   ParamFile = strdup(args_info.paramFile_arg);
  /* Allow other pairs in addition to the usual AU,GC,and GU pairs */
  if(args_info.nsp_given)         ns_bases = strdup(args_info.nsp_arg);
  /* set the maximum base pair span */
  if(args_info.span_given)        maxdist = args_info.span_arg;
  
  /* free allocated memory of command line data structure */
  RNALfold_cmdline_parser_free (&args_info);

  /*
  #############################################
  # begin initializing
  #############################################
  */
  if (ParamFile != NULL)
    read_parameter_file(ParamFile);
   
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

  istty = isatty(fileno(stdout))&&isatty(fileno(stdin));
        
  /*
  #############################################
  # main loop: continue until end of file
  #############################################
  */
  do{
    if(istty) print_tty_input_seq();

    /* skip comment lines and break on EOF */
    if(skip_comment_lines(&line) < 0) break;

    /* break on @ */
    if(strcmp(line, "@") == 0)        break;

    /* skip fasta header */
    if(*line == '>'){
      printf("%s\n", line); free(line);
      if(skip_comment_lines(&line) < 0) break;
    }

    string = (char *) space(strlen(line)+1);
    (void) sscanf(line,"%s",string);
    free(line);

    length    = (int) strlen(string);
    structure = (char *) space((unsigned) length+1);

    if(noconv)  str_RNA2RNA(string);
    else        str_DNA2RNA(string);

    if (istty)  printf("length = %d\n", length);

    /* initialize_fold(length); */
    update_fold_params();
    min_en = Lfold((const char *)string, structure, maxdist);
    printf("%s\n%s", string, structure);

    if (istty)
      printf("\n minimum free energy = %6.2f kcal/mol\n", min_en);
    else
      printf(" (%6.2f)\n", min_en);

    (void) fflush(stdout);
    free(string);
    free(structure); 
  } while (1);
  return 0;
}

