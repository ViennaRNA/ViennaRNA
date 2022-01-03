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
#include "ViennaRNA/io/file_formats.h"
#include "ct2db_cmdl.h"

typedef struct _parameters {
  int no_pk;
  int no_mod;
  int verbose;
} parameters;


int
process_input(FILE        *fp,
              const char  *filename,
              parameters  *opt);


static void
remove_noncanonical_bases(char *sequence)
{
  if (sequence)
    for (char *ptr = sequence; *ptr; ptr++) {
      switch (*ptr) {
        case  'A': /* fall through */
        case  'a': /* fall through */
        case  'C': /* fall through */
        case  'c': /* fall through */
        case  'G': /* fall through */
        case  'g': /* fall through */
        case  'U': /* fall through */
        case  'u': /* fall through */
        case  'T': /* fall through */
        case  't':
          break;
        default:
          *ptr = 'N';
          break;
      }
    }
}

int
main(int  argc,
     char *argv[])
{
  struct ct2db_args_info  args_info;
  parameters              ct2db_params;

  ct2db_params.no_pk    = 0;
  ct2db_params.no_mod   = 0;
  ct2db_params.verbose  = 0;

  /*
   #############################################
   # check the command line prameters
   #############################################
   */
  if (ct2db_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  /* remove Pseudknots? */
  if (args_info.no_pk_given)
    ct2db_params.no_pk = 1;

  /* replace all characters by U or N */
  if (args_info.no_modified_given)
    ct2db_params.no_mod = 1;

  /* be verbose ? */
  if (args_info.verbose_given)
    ct2db_params.verbose = 1;

  char  **input_files = NULL;
  int   i, num_files = 0;

  /* collect all unnamed options */
  if (args_info.inputs_num > 0) {
    input_files = (char **)vrna_realloc(input_files, sizeof(char *) * args_info.inputs_num);
    for (i = 0; i < args_info.inputs_num; i++)
      input_files[num_files++] = strdup(args_info.inputs[i]);
  }

  /* free allocated memory of command line data structure */
  ct2db_cmdline_parser_free(&args_info);

  if (num_files > 0) {
    int i, skip;
    for (skip = i = 0; i < num_files; i++) {
      if (!skip) {
        FILE *input_stream = fopen((const char *)input_files[i], "r");

        if (!input_stream)
          vrna_message_error("Unable to open %d. input file \"%s\" for reading", i + 1,
                             input_files[i]);

        if (ct2db_params.verbose) {
          vrna_message_info(stderr,
                            "Processing %d. input file \"%s\"",
                            i + 1,
                            input_files[i]);
        }

        if (process_input(input_stream, (const char *)input_files[i], &ct2db_params) == 0)
          skip = 1;

        fclose(input_stream);
      }

      free(input_files[i]);
    }
  } else {
    (void)process_input(stdin, NULL, &ct2db_params);
  }

  free(input_files);

  return EXIT_SUCCESS;
}

int
process_input(FILE        *fp,
              const char  *filename,
              parameters  *opt)
{
  char          *id, *sequence, *structure, *ss, *remainder;
  unsigned int  parse_options;

  parse_options = 0;
  remainder     = NULL;

  if (opt->verbose)
    parse_options |= VRNA_INPUT_VERBOSE;

  while (vrna_file_connect_read_record(fp, &id, &sequence, &structure, &remainder, parse_options)) {
    if (opt->no_pk) {
      ss = vrna_db_pk_remove((const char *)structure, VRNA_BRACKETS_ANY);
      free(structure);
    } else {
      ss = structure;
    }

    if (opt->no_mod)
      remove_noncanonical_bases(sequence);

    printf(">%s\n%s\n%s\n", (id) ? id : "seq", sequence, ss);
    free(id);
    free(sequence);
    free(ss);
  }

  return 1;
}
