/*
    file_formats_aln.c

    Various functions dealing with file formats for RNA sequence alingments

    (c) 2016 Ronny Lorenz

    ViennaRNA package
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/aln_util.h"
#include "ViennaRNA/file_formats_aln.h"

#define MAX_NUM_NAMES    500

/*
#################################
# PRIVATE VARIABLES             #
#################################
*/

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/

typedef int (aln_parser_function)(FILE *fp, char ***names, char ***aln, char **id, char **structure);

PRIVATE aln_parser_function parse_aln_stockholm;

PRIVATE aln_parser_function parse_aln_clustal;

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PUBLIC int
vrna_file_alignment_read( const char *filename,
                          char ***names,
                          char ***aln,
                          char  **id,
                          char  **structure,
                          unsigned int options){

  FILE  *fp;
  char  *line = NULL;
  int   i, n, seq_num;

  seq_num   = 0;

  if(!(fp = fopen(filename, "r"))){
    vrna_message_warning("Alignment file could not be opened!");
    return seq_num;
  }

  if(names && aln){
    *names  = NULL;
    *aln    = NULL;
  } else {
    return seq_num;
  }

  if(id)
    *id = NULL;

  if(structure)
    *structure = NULL;

  /* if no alignment file format was specified, lets try to guess it */
  if(options == 0)
    options = VRNA_FILE_FORMAT_ALN_DEFAULT;

#define known_aln_formats_num   2

  unsigned int known_formats[known_aln_formats_num] = {
                                                        VRNA_FILE_FORMAT_ALN_STOCKHOLM,
                                                        VRNA_FILE_FORMAT_ALN_CLUSTAL
                                                      };

  aln_parser_function *parsers[known_aln_formats_num] = {
                                                          parse_aln_stockholm,
                                                          parse_aln_clustal
                                                        };

  int r = -1;
  long int fp_position = ftell(fp);

  for(i = 0; i < known_aln_formats_num; i++){
    if(options & known_formats[i]){
      /* go back to beginning of file */
      if(!fseek(fp, fp_position, SEEK_SET)){
        r = parsers[i](fp, names, aln, id, structure);
        if(r > 0)
          break;
      } else {
        fprintf(stderr, "ERROR: Something unexpected happened while parsing the alignment file");
      }
    }
  }

  fclose(fp);

  if(r == -1){
    vrna_message_warning("Alignment file parser is unknown (or not specified?)");
  } else {
    seq_num = r;
  }

  return seq_num;
}

PRIVATE int
parse_aln_stockholm(FILE *fp,
                    char ***names,
                    char ***aln,
                    char **id,
                    char **structure){

  return vrna_file_stockholm_read_record(fp, names, aln, id, structure, 0);
}

PRIVATE int
parse_aln_clustal(FILE *fp,
                  char ***names,
                  char ***aln,
                  char **id,
                  char **structure){

  int r;
  *id         = NULL;
  *structure  = NULL;

  /* read_clustal is limited to MAX_NUM_NAMES sequences in the alignment */
  *names  = (char **)vrna_alloc(sizeof(char *) * MAX_NUM_NAMES);
  *aln    = (char **)vrna_alloc(sizeof(char *) * MAX_NUM_NAMES);

  r = read_clustal(fp, *aln, *names);

  /* resize according to actual sequence number */
  *names  = (char **)vrna_realloc(*names, sizeof(char *) * (r + 1));
  *aln    = (char **)vrna_realloc(*aln, sizeof(char *) * (r + 1));

  return r;
}

PUBLIC int
vrna_file_stockholm_read_record(  FILE  *fp,
                                  char  ***names,
                                  char  ***aln,
                                  char  **id,
                                  char  **structure,
                                  int   verbosity){

  char  *line = NULL;
  int   i, n, seq_num, seq_length;

  seq_num     = 0;
  seq_length  = 0;

  if(!fp){
    if(verbosity >= 0)
      vrna_message_warning("can't read from filepointer while parsing Stockholm formatted sequence alignment!");
    return seq_num;
  }

  if(names && aln){
    *names  = NULL;
    *aln    = NULL;
  } else {
    return seq_num;
  }

  if(id)
    *id = NULL;

  if(structure)
    *structure = NULL;

  int inrecord = 0;
  while(line = get_line(fp)){
    if(strstr(line, "STOCKHOLM 1.0")){
      inrecord = 1;
      break;
    }
  }

  if(inrecord){
    while((line=get_line(fp))){

      if(strncmp(line, "//", 2) == 0){
        /* end of alignment */
        free(line);
        line = NULL;
        break;
      }

      n = (int)strlen(line);

      switch(*line){
        /* we skip lines that start with whitespace */
        case ' ': case '\0':
          goto stockholm_next_line;

        /* Stockholm markup, or comment */
        case '#':
          if(strstr(line, "STOCKHOLM 1.0")){
            if(verbosity >= 0)
              vrna_message_warning("Malformatted Stockholm record, missing // ?");
            /* drop everything we've read so far and start new, blank record */
            if(id != NULL){
              free(*id);
              *id = NULL;
            }
            if(structure != NULL){
              free(*structure);
              *structure = NULL;
            }
            for(i = 0; i < seq_num; i++){
              free((*names)[i]);
              free((*aln)[i]);
            }
            free(*names);
            *names = NULL;
            free(*aln);
            *aln = NULL;
            seq_num = 0;
          } else if(strncmp(line, "#=GF", 4) == 0){
            /* found feature markup */
            if((id != NULL) && (strncmp(line, "#=GF ID", 7) == 0)){
              *id = (char *)vrna_alloc(sizeof(char) * n);
              if(sscanf(line, "#=GF ID %s", *id) == 1){
                *id = (char *)vrna_realloc(*id, sizeof(char) * (strlen(*id) + 1));
              } else {
                free(*id);
                *id = NULL;
              }
            }
          } else if(strncmp(line, "#=GC", 4) == 0){
            /* found per-column annotation */
            if((structure != NULL) && (strncmp(line, "#=GC SS_cons", 12) == 0)){
              *structure = (char *)vrna_alloc(sizeof(char) * n);
              if(sscanf(line, "#=GC SS_cons %s", *structure) == 1){
                *structure = (char *)vrna_realloc(*structure, sizeof(char) * (strlen(*structure) + 1));
              } else {
                free(*structure);
                *structure = NULL;
              }
            }
          } else if(strncmp(line, "#=GS", 4) == 0){
            /* found generic per-sequence annotation */
          } else if(strncmp(line, "#=GR", 4) == 0){
            /* found generic per-Residue annotation */
          } else {
            /* may be comment? */
          }
          break;

        /* should be sequence */
        default:
          {
            int tmp_l;
            char *tmp_name  = (char *)vrna_alloc(sizeof(char) * (n + 1));
            char *tmp_seq   = (char *)vrna_alloc(sizeof(char) * (n + 1));
            if(sscanf(line, "%s %s", tmp_name, tmp_seq) == 2){
              seq_num++;
              tmp_l = (int)strlen(tmp_seq);

              if(seq_num == 1){
                seq_length = tmp_l;
              } else {  /* check sequence length against first */
                if(seq_length != tmp_l){
                  if(verbosity >= 0)
                    vrna_message_warning("Discarding Stockholm record! Sequence lengths do not match.");

                  /* drop everything we've read so far and abort parsing */
                  if(id != NULL){
                    free(*id);
                    *id = NULL;
                  }

                  if(structure != NULL){
                    free(*structure);
                    *structure = NULL;
                  }

                  for(i = 0; i < seq_num; i++){
                    free((*names)[i]);
                    free((*aln)[i]);
                  }

                  free(*names);
                  *names = NULL;

                  free(*aln);
                  *aln = NULL;

                  seq_num = 0;

                  free(tmp_name);
                  free(tmp_seq);
                  free(line);
                  line = NULL;

                  goto stockholm_exit;
                }
              }

              /* store sequence name */
              (*names) = (char **)vrna_realloc(*names, sizeof(char *) * seq_num);
              (*names)[seq_num - 1] = (char *)vrna_alloc(sizeof(char) * (strlen(tmp_name) + 1));
              strcpy((*names)[seq_num - 1], tmp_name);

              /* store aligned sequence */
              (*aln) = (char **)vrna_realloc(*aln, sizeof(char *) * seq_num);
              (*aln)[seq_num - 1] = (char *)vrna_alloc(sizeof(char) * (tmp_l + 1));
              strcpy((*aln)[seq_num - 1], tmp_seq);

            }
            free(tmp_name);
            free(tmp_seq);
          }
          break;
      }

stockholm_next_line:

      free(line);
    }
  } else {
    if(verbosity > 0)
      vrna_message_warning("Did not find any Stockholm formatted record\n");
  }

stockholm_exit:

  free(line);

  /*
    append additional entry in 'aln' and 'names' pointing to NULL (this may be
    used as an indication for the end of the sequence list)
  */
  if(seq_num > 0){
    (*aln)            = (char **)vrna_realloc(*aln, sizeof(char *) * (seq_num + 1));
    (*names)          = (char **)vrna_realloc(*names, sizeof(char *) * (seq_num + 1));
    (*aln)[seq_num]   = NULL;
    (*names)[seq_num] = NULL;
    if(verbosity >= 0)
      fprintf(stderr, "%d sequences; length of alignment %d.\n", seq_num, strlen((*aln)[0]));
  }

  return seq_num;
}
