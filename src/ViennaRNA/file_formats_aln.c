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
#include "ViennaRNA/file_formats.h"
#include "ViennaRNA/file_formats_aln.h"

#define MAX_NUM_NAMES    500

/*
#################################
# PRIVATE VARIABLES             #
#################################
*/

/*
#################################
# STATIC DECLARATIONS           #
#################################
*/

typedef int (aln_parser_function)(FILE *fp, char ***names, char ***aln, char **id, char **structure);

typedef struct {
  unsigned int        code;
  aln_parser_function *parser;
} parsable;

/* number of known alignment parsers */
#define NUM_PARSERS 3

PRIVATE aln_parser_function parse_aln_stockholm;

PRIVATE aln_parser_function parse_aln_clustal;

PRIVATE aln_parser_function parse_aln_fasta;

PRIVATE int
parse_fasta_alignment(FILE *fp,
                      char ***names,
                      char ***aln,
                      int verbosity);

PRIVATE int
parse_clustal_alignment(FILE *clust,
                        char ***names,
                        char ***aln,
                        int verbosity);

PRIVATE int
check_alignment(const char **names,
                const char **aln,
                int seq_num);

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
  parsable  known_parsers[NUM_PARSERS] =  {
    /* option, parser */
    { VRNA_FILE_FORMAT_ALN_STOCKHOLM, parse_aln_stockholm },
    { VRNA_FILE_FORMAT_ALN_CLUSTAL,   parse_aln_clustal },
    { VRNA_FILE_FORMAT_ALN_FASTA,     parse_aln_fasta }
  };

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

  int r = -1;
  long int fp_position = ftell(fp);

  for(i = 0; i < NUM_PARSERS; i++){
    if((options & known_parsers[i].code) && (known_parsers[i].parser)){
      /* go back to beginning of file */
      if(!fseek(fp, fp_position, SEEK_SET)){
        r = known_parsers[i].parser(fp, names, aln, id, structure);
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

    if((seq_num > 0) && (!(options & VRNA_FILE_FORMAT_ALN_NOCHECK))){
      if(!check_alignment((const char **)(*names), (const char **)(*aln), seq_num)){
        vrna_message_warning("Alignment did not pass sanity checks!");

        /* discard the data we've read! */
        if(id)
          *id = NULL;

        if(structure)
          *structure = NULL;

        for(i = 0; i < seq_num; i++){
          free((*aln)[i]);
          free((*names)[i]);
        }
        free(*aln);
        free(*names);

        *aln    = NULL;
        *names  = NULL;
        seq_num = 0;
      }
    }
  }

  return seq_num;
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
      free(line);
      break;
    }
    free(line);
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

  /* clustal format doesn't contain id's or structure information */
  if(id)
    *id = NULL;
  if(structure)
    *structure = NULL;

  return parse_clustal_alignment(fp, names, aln, 0);
}

PRIVATE int
parse_aln_fasta(FILE *fp,
                char ***names,
                char ***aln,
                char **id,
                char **structure){

  /* fasta alignments do not contain an id, or structure information */
  if(id)
    *id = NULL;
  if(structure)
    *structure = NULL;

  return parse_fasta_alignment(fp, names, aln, 0);
}

PRIVATE int
parse_fasta_alignment(FILE *fp,
                      char ***names,
                      char ***aln,
                      int verbosity){

  unsigned int  read_opt, rec_type;
  int           num_seq;
  char          *rec_id, *rec_sequence, **rec_rest;

  rec_id        = NULL;
  rec_sequence  = NULL;
  rec_rest      = NULL;
  num_seq       = 0;
  read_opt      = VRNA_INPUT_NO_REST; /* read sequence and header information only */

  /* read until EOF or user abort */
  while(
    !((rec_type = vrna_file_fasta_read_record(&rec_id, &rec_sequence, &rec_rest, fp, read_opt))
        & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT))){

    if(rec_id){ /* valid FASTA entry */
      num_seq++;

      char *id = (char *)vrna_alloc(sizeof(char) * strlen(rec_id));
      (void) sscanf(rec_id, ">%s", id);

      (*names) = (char **)vrna_realloc(*names, sizeof(char *) * num_seq);
      (*names)[num_seq - 1] = id;

      vrna_seq_toupper(rec_sequence);
      (*aln) = (char **)vrna_realloc(*aln, sizeof(char *) * num_seq);
      (*aln)[num_seq - 1] = strdup(rec_sequence);
    }

    free(rec_id);
    free(rec_sequence);
    free(rec_rest);
  }

  free(rec_id);
  free(rec_sequence);
  free(rec_rest);

  /*
    append additional entry in 'aln' and 'names' pointing to NULL (this may be
    used as an indication for the end of the sequence list)
  */
  if(num_seq > 0){
    (*aln)            = (char **)vrna_realloc(*aln, sizeof(char *) * (num_seq + 1));
    (*names)          = (char **)vrna_realloc(*names, sizeof(char *) * (num_seq + 1));
    (*aln)[num_seq]   = NULL;
    (*names)[num_seq] = NULL;
    if(verbosity >= 0)
      fprintf(stderr, "%d sequences; length of alignment %d.\n", num_seq, strlen((*aln)[0]));
  }

  return num_seq;
}


PRIVATE int
parse_clustal_alignment(FILE *clust,
                        char ***names,
                        char ***aln,
                        int verbosity){

  char *line, *name, *seq;
  int  n, r, nn=0, num_seq = 0, i;

  if((line=get_line(clust)) == NULL){
    if(verbosity >= 0)
      vrna_message_warning("Empty CLUSTALW file\n");
    return 0;
  }

  if(strncmp(line,"CLUSTAL", 7) != 0){
    if(verbosity > 0)
      vrna_message_warning("This doesn't look like a CLUSTAL/STOCKHOLM file, sorry\n");

    free(line);
    return 0;
  }

  free(line);
  line = get_line(clust);

  while (line!=NULL) {
    n = strlen(line);

    if((n < 4) || isspace((int)line[0])) {
      /* skip non-sequence line */
      free(line);
      line  = get_line(clust);
      nn    = 0;  /* reset sequence number */
      continue;
    }

    /* skip comments */
    if(line[0] == '#'){
      free(line);
      line = get_line(clust);
      continue;
    }

    seq   = (char *)vrna_alloc(sizeof(char) * (n + 1));
    name  = (char *)vrna_alloc(sizeof(char) * (n + 1));
    if(sscanf(line, "%s %s", name, seq) == 2){
      /* realloc to actual sizes */
      seq   = (char *)vrna_realloc(seq, sizeof(char) * (strlen(seq) + 1));
      name  = (char *)vrna_realloc(name, sizeof(char) * (strlen(name) + 1));
      for(i = 0; i < strlen(seq); i++){
        if(seq[i] == '.')
          seq[i] = '-'; /* replace '.' gaps with '-' */
      }

      /* convert sequence to uppercase letters */
      vrna_seq_toupper(seq);

      if(nn == num_seq){ /* first time */
        (*names)      = (char **)vrna_realloc(*names, sizeof(char *) * (nn + 1));
        (*names)[nn]  = strdup(name);
        (*aln)        = (char **)vrna_realloc(*aln, sizeof(char *) * (nn + 1));
        (*aln)[nn]    = strdup(seq);
      } else {
        if (strcmp(name, (*names)[nn]) != 0) {
          /* name doesn't match */
           fprintf(stderr,
                   "Sorry, your file is messed up (inconsitent seq-names)\n");
           free(line); free(seq);
           return 0;
        }
        (*aln)[nn] = (char *)vrna_realloc((*aln)[nn], strlen(seq) + strlen((*aln)[nn])+1);
        strcat((*aln)[nn], seq);
      }
      nn++;
      if(nn > num_seq)
        num_seq = nn;

      free(seq);
      free(name);
    }
    free(line);

    line = get_line(clust);
  }

  /*
    append additional entry in 'aln' and 'names' pointing to NULL (this may be
    used as an indication for the end of the sequence list)
  */
  if(num_seq > 0){
    (*aln)            = (char **)vrna_realloc(*aln, sizeof(char *) * (num_seq + 1));
    (*names)          = (char **)vrna_realloc(*names, sizeof(char *) * (num_seq + 1));
    (*aln)[num_seq]   = NULL;
    (*names)[num_seq] = NULL;
    if(verbosity >= 0)
      fprintf(stderr, "%d sequences; length of alignment %d.\n", num_seq, strlen((*aln)[0]));
  }

  return num_seq;
}

PRIVATE int
check_alignment(const char **names,
                const char **aln,
                int seq_num){

  int i, j, l, pass = 1;

  /* check for unique names */
  for(i = 0; i < seq_num; i++){
    for(j = i + 1; j < seq_num; j++){
      if(!strcmp(names[i], names[j])){
        vrna_message_warning("Sequence IDs in input alignment are not unique!");
        pass = 0;
      }
    }
  }

  /* check for equal lengths of sequences */
  l = (int)strlen(aln[0]);
  for(i = 1; i < seq_num; i++)
    if((int)strlen(aln[i]) != l){
      vrna_message_warning("Sequence lengths in input alignment do not match!");
      pass = 0;
    }

  return pass;
}
