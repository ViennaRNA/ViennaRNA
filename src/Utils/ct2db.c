

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/data_structures.h"
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/MEA.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/structure_utils.h"
#include "ViennaRNA/params.h"
#include "ct2db_cmdl.h"

int main(int argc, char *argv[]){
  struct  ct2db_args_info args_info;

  char          *line       = NULL;

  int           pkfree      = 0;
  int           verbose     = 0;

  char          *seq        = NULL;
  char          *structure  = NULL;
  plist         *pairs      = NULL;
  int           num_pairs   = 0;
  unsigned int  seq_pos     = 0;
  unsigned int  length      = 0;

  /*
  #############################################
  # check the command line prameters
  #############################################
  */
  if(ct2db_cmdline_parser (argc, argv, &args_info) != 0) exit(1);
  /* remove Pseudknots? */
  if(args_info.removePK_given)
    pkfree = 1;
  /* be verbose ? */
  if(args_info.verbose_given)
    verbose = 1;

  /* free allocated memory of command line data structure */
  ct2db_cmdline_parser_free (&args_info);

  do{
    if ((line = vrna_read_line(stdin))==NULL) break;
    /* skip comment lines */
    while ((*line=='*')||(*line=='\0')||(*line=='>')||(*line=='#')) {
      free(line);
      if ((line = vrna_read_line(stdin))==NULL) break;
    }

    if ((line ==NULL) || (strcmp(line, "@") == 0)) break;

    if(length != 0){
      if(!seq){
        /*
        printf("sequence length: %u\n", length);
        */
        seq   = (char *)vrna_alloc(sizeof(char) * (length + 1));
        pairs = (plist *)vrna_alloc(sizeof(plist) * (2*length));
      }
    }

    char      *ptr, c;
    int       row       = 0;
    long int  curr_int  = -1;
    long int  last_int  = -1;
    int       i         = -1;
    int       j         = -1;

    ptr = strtok(line, " \t");

    while(ptr != NULL){
/*
      printf("%d. %s\n", row, ptr);
*/

      char* end;
      curr_int          = strtol(ptr, &end, 10);
      if(*end) curr_int = -1; /* not read an int */

      switch(row){
        
        case 0:   if( curr_int != -1)
                    i = (int) curr_int;
                  break;
        case 4:   if( curr_int != -1 )
                    if((i != -1) && (curr_int > i)){
                      j = curr_int;
                      pairs[num_pairs].i      = i;
                      pairs[num_pairs].j      = j;
                      pairs[num_pairs].p      = 1.;
                      pairs[num_pairs++].type = 0;
                    }
                  break;
        case 1:   if(sscanf(ptr, "%c", &c)){
                    switch(c){
                      case 'E': if(!strcmp("ENERGY", ptr)){
                                  if(last_int != -1)
                                  length = (unsigned int)last_int;
                                }
                                break;
                      case 'A':
                      case 'U':
                      case 'T':
                      case 'C':
                      case 'G': if(seq_pos < length){
                                  seq[seq_pos++] = c;
                                }
                                break;
                      default:  break;
                    }
                  }
                  break;
        default:  break;
      }
      ptr = strtok(NULL, " \t");
      if(curr_int != -1)
        last_int = curr_int;
      row++;
    }

    free(line);

  } while(1);

  if(seq){
    pairs = vrna_realloc(pairs, sizeof(plist) * (num_pairs+1));
    pairs[num_pairs].i    = 0;
    pairs[num_pairs].j    = 0;
    pairs[num_pairs].p    = 0.;
    pairs[num_pairs].type = 0;


/*
    plist *p = pairs;
    for(p = pairs; (*p).i; p++)
      printf("pair (%d,%d)\n", (*p).i, (*p).j);
*/

    if(pkfree){
      float mea, MEAgamma;
      MEAgamma = 2.0;
      vrna_exp_param_t *params = vrna_exp_params(NULL);

      structure = (char *)vrna_alloc(sizeof(char) * (length + 1));
      strcpy(structure, seq);

      mea = MEA_seq(pairs, seq, structure, MEAgamma, params);

      char *structure_tmp = vrna_db_from_plist(pairs, length);
      int d =  vrna_bp_distance((const char *)structure, (const char *)structure_tmp);
      if(verbose && (d > 0))
        fprintf(stderr, "removed %d pairs from pseudoknotted structure\n", d);

      free(structure_tmp);
      free(params);
    } else {
      structure = vrna_db_from_plist(pairs, length);
    }

    printf("%s\n%s\n", seq, structure);
    free(pairs);
    free(seq);
  }

  return(EXIT_SUCCESS);
}
