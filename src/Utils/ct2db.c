#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/part_func.h"
#include "ViennaRNA/MEA.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/structures.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/io/utils.h"
#include "ct2db_cmdl.h"

typedef struct _parameters {
  int pkfree;
  int convertToRNA;
  int verbose;
} parameters;


typedef int
(structure_callback)(char       *sequence,
                     plist      *pairs,
                     int        pairs_length,
                     parameters *ct2db_params);


typedef struct _ct_file_sequence_header {
  int     seq_length;
  double  energy;
  char    *sequence_name;
} ct_file_sequence_header;


typedef struct _ct_file_index_line {
  int   index_i;
  char  base;
  int   index_minus_one;
  int   index_plus_one;
  int   pair_index_j;
  int   natural_number_index;
} ct_file_index_line;


int
read_index(char *ptr)
{
  char  *end;
  int   i = strtol(ptr, &end, 10);

  if (end == ptr)
    /* index could not be parsed. */
    i = -1;

  return i;
}


int
read_double(char *ptr)
{
  char    *end;
  double  d = strtod(ptr, &end);

  if (end == ptr) {
    /* double could not be parsed. */
  }

  return d;
}


/**
 * @param line - in
 * @param head - in & out
 * @param index_line - in & out
 */
int
read_line(char                    *line,
          ct_file_sequence_header *head,
          ct_file_index_line      *index_line)
{
  char  *ptr, c;
  int   col       = 0;
  int   first_int = -1;

  ptr = strtok(line, " \t");
  int   read_energy = 0;

  while (ptr != NULL) {
    switch (col) {
      case 0:
        first_int         = read_index(ptr);
        head->seq_length  = first_int;
        /* strcpy(head->sequence_name, ptr); //TODO; parse this if you need it. */
        /* or */
        index_line->index_i = first_int;
        /* we know it only if the second line (first index line) has been parsed. */
        break;
      case 1:
        if (strcmp("ENERGY", ptr) == 0 || strcmp("Energy", ptr) == 0 || strcmp("energy", ptr) == 0)
          read_energy = 1;

        if (sscanf(ptr, "%c", &c))
          index_line->base = c;

        break;
      case 2:
        index_line->index_minus_one = read_index(ptr);
        break;
      case 3:
        if (read_energy) {
          head->energy  = read_double(ptr);
          read_energy   = 0;
        }

        index_line->index_plus_one = read_index(ptr);
        break;
      case 4:
        index_line->pair_index_j = read_index(ptr);
        break;
      case 5:
        index_line->natural_number_index = read_index(ptr);
        break;
      default:
        break;
    }
    ptr = strtok(NULL, " \t");
    col++;
  }
  return 0;
}


int
read_ct_file(FILE               *file,
             structure_callback *report_structure,
             parameters         *ct2db_params)
{
  int                     num_pairs = 0;
  unsigned int            seq_pos   = 0;
  unsigned int            length    = 0;
  char                    *line     = NULL;
  vrna_ep_t               *pairs    = NULL;
  char                    *sequence = NULL;
  ct_file_sequence_header previous_head;

  do {
    if ((line = vrna_read_line(stdin)) == NULL)
      break;

    /* skip comment lines */
    while ((*line == '*') || (*line == '\0') || (*line == '>') || (*line == '#')) {
      free(line);
      if ((line = vrna_read_line(stdin)) == NULL)
        break;
    }

    if ((line == NULL) || (strcmp(line, "@") == 0))
      break;

    /*   We know the sequence length either from the first line (comments are not defined in ct) or from the line with the entry 'ENERGY'.
     * actually a ct file can contain more than one sequence with structure. After the last line with (index, char, 4*index) the new sequence could start.
     */
    ct_file_sequence_header head = {
      -1, 0, NULL
    };
    ct_file_index_line      index_line = {
      -1, '\0', -1, -1, -1, -1
    };
    read_line(line, &head, &index_line);

    if ((index_line.index_i != -1) &&
        (index_line.base != '\0') &&
        (index_line.index_minus_one != -1) &&
        (index_line.index_plus_one != -1) &&
        (index_line.pair_index_j != -1) &&
        (index_line.natural_number_index != -1)) {
      /* now we can be sure that a index line has been pared. */
      if (index_line.index_i == 1) {
        /* if it is the first line, the previous line was the header. */
        length = previous_head.seq_length;
        /* previous_head.energy; */

        /* report previous structure */
        if (sequence) {
          pairs                 = vrna_realloc(pairs, sizeof(plist) * (num_pairs + 1));
          pairs[num_pairs].i    = 0;
          pairs[num_pairs].j    = 0;
          pairs[num_pairs].p    = 0.;
          pairs[num_pairs].type = 0;

          report_structure(sequence, pairs, num_pairs, ct2db_params);
        }

        /* free previous pairs and sequence */
        free(pairs);
        free(sequence);
        seq_pos   = 0;
        num_pairs = 0;

        /* allocate new memory. */
        sequence          = (char *)vrna_alloc(sizeof(char) * (length + 1));
        sequence          = memset(sequence, 'N', sizeof(char) * length);
        sequence[length]  = '\0';
        /* printf("strlen: %d\n",strlen(sequence)); */
        pairs = (plist *)vrna_alloc(sizeof(plist) * (2 * length));
      }

      /* continue with sequence and add pairs. */
      /* replace all nucleotide chars to upper and replace T to U and non-RNA chars to N*/
      char c = index_line.base;
      /* printf("%c\n",c); */
      if (ct2db_params->convertToRNA == 1) {
        c = (char)toupper(c);
        if (c == 'T')
          c = 'U';

        if (c != 'A' && c != 'C' && c != 'G' && c != 'U')
          c = 'N';
      }

      /*  accept all characters (could be a modified base) */
      if (seq_pos < length)
        sequence[seq_pos++] = c;

      int i = index_line.index_i;
      int j = index_line.pair_index_j;
      if (i > length || j > length) {
        fprintf(stderr,
                "ERROR: The index pair (%d, %d) is not valid! The sequence Length is %d!\n",
                i,
                j,
                (int)length);
        exit(EXIT_FAILURE);
      }

      if (i > 0 && j > 0) {
        /* printf("%d %d\n",i,j); */
        if (j > i) {
          pairs[num_pairs].i      = i;
          pairs[num_pairs].j      = j;
          pairs[num_pairs].p      = 1.;
          pairs[num_pairs++].type = 0;
          /* printf("%d %d\n",i,j); */
        } else {
          int not_in_list = 1;
          for (int k = 0; k < num_pairs; k++) {
            if (pairs[k].i == j) {
              not_in_list = 0;
              break;
            }
          }
          if (not_in_list == 1) {
            /* printf("%d %d\n",i,j); */
            pairs[num_pairs].i      = j;
            pairs[num_pairs].j      = i;
            pairs[num_pairs].p      = 1.;
            pairs[num_pairs++].type = 0;
          }
        }
      }
    }

    previous_head = head;

    free(line);
  } while (1);

  if (sequence) {
    pairs                 = vrna_realloc(pairs, sizeof(plist) * (num_pairs + 1));
    pairs[num_pairs].i    = 0;
    pairs[num_pairs].j    = 0;
    pairs[num_pairs].p    = 0.;
    pairs[num_pairs].type = 0;

    report_structure(sequence, pairs, num_pairs, ct2db_params);
    free(pairs);
    free(sequence);
  }

  return 0;
}


int
convert_ct_structure_to_db(char       *sequence,
                           plist      *pairs,
                           int        pairs_length,
                           parameters *ct2db_params)
{
  char  *structure;
  int   length = strlen(sequence);

  if (ct2db_params->pkfree) {
    float             MEAgamma;
    MEAgamma = 2.0;
    vrna_exp_param_t  *params = vrna_exp_params(NULL);
    /* printf("length %d %d\n", pairs_length,strlen(sequence)); */
    structure         = (char *)vrna_alloc(sizeof(char) * (length + 1));
    structure         = memset(structure, '.', length * sizeof(char));
    structure[length] = '\0';
    /* strcpy(structure, seq); */

    (void) MEA_seq(pairs, sequence, structure, MEAgamma, params);

    char  *structure_tmp  = vrna_db_from_plist(pairs, length);
    int   d               = vrna_bp_distance((const char *)structure, (const char *)structure_tmp);
    if (ct2db_params->verbose && (d > 0))
      fprintf(stderr, "removed %d pairs from pseudoknotted structure\n", d);

    free(structure_tmp);
    free(params);
  } else {
    structure = vrna_db_from_plist(pairs, length);
  }

  printf("%s\n%s\n", sequence, structure);
  free(structure);

  return 0;
}


int
main(int  argc,
     char *argv[])
{
  struct ct2db_args_info  args_info;
  parameters              ct2db_params;

  ct2db_params.pkfree       = 0;
  ct2db_params.convertToRNA = 0;
  ct2db_params.verbose      = 0;

  /*
   #############################################
   # check the command line prameters
   #############################################
   */
  if (ct2db_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  /* remove Pseudknots? */
  if (args_info.removePK_given)
    ct2db_params.pkfree = 1;

  /* replace all characters by U or N */
  if (args_info.convertToRNA_given)
    ct2db_params.convertToRNA = 1;

  /* be verbose ? */
  if (args_info.verbose_given)
    ct2db_params.verbose = 1;

  /* free allocated memory of command line data structure */
  ct2db_cmdline_parser_free(&args_info);

  read_ct_file(stdin, convert_ct_structure_to_db, &ct2db_params);

  return EXIT_SUCCESS;
}
