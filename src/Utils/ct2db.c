#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/partfunc/global.h"
#include "ViennaRNA/structures/mea.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/structures/pairtable.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/io/utils.h"
#include "ViennaRNA/io/file_formats.h"

#include "../bin/gengetopt_helpers.h"

#include "ct2db_cmdl.h"

typedef struct _parameters {
  int           no_pk;
  int           no_mod;
  int           verbose;
  unsigned int  id_from_filename;
  char          *filename_suffix;
  char          *id;
  char          *id_prefix;
  char          *id_suffix;
  unsigned int  id_num;
  unsigned int  id_num_start;
} parameters;


int
process_input(FILE        *fp,
              const char  *filename,
              parameters  *opt);


static void
remove_noncanonical_bases(char *sequence)
{
  if (sequence) {
    for (char *ptr = sequence; *ptr; ptr++) {
      switch (*ptr) {
        case  'A':  /* fall through */
        case  'a':  /* fall through */
        case  'C':  /* fall through */
        case  'c':  /* fall through */
        case  'G':  /* fall through */
        case  'g':  /* fall through */
        case  'U':  /* fall through */
        case  'u':  /* fall through */
        case  'T':  /* fall through */
        case  't':
          break;
        default:
          *ptr = 'N';
          break;
      }
    }
  }
}


static const char *
find_suffix(const char  *string,
            const char  *suffix)
{
  size_t  l = 0;
  size_t  n = 0;

  l = strlen(suffix);
  n = strlen(string);

  if ((l) &&
      (n) &&
      (n >= l)) {
    if (strcmp(string + (n - l), suffix) == 0)
      return string + n - l;
  }

  return NULL;
}


int
main(int  argc,
     char *argv[])
{
  struct ct2db_args_info  args_info;
  parameters              ct2db_params;

  ct2db_params.no_pk            = 0;
  ct2db_params.no_mod           = 0;
  ct2db_params.verbose          = 0;
  ct2db_params.id               = NULL;
  ct2db_params.id_from_filename = 0;
  ct2db_params.filename_suffix  = NULL;
  ct2db_params.id_num           = 0;
  ct2db_params.id_num_start     = 0;
  ct2db_params.id_prefix        = NULL;
  ct2db_params.id_suffix        = NULL;

  /*
   #############################################
   # check the command line prameters
   #############################################
   */
  if (ct2db_cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  /* prepare logging system and verbose mode */
  ggo_log_settings(args_info, ct2db_params.verbose);

  /* remove Pseudknots? */
  if (args_info.no_pk_given)
    ct2db_params.no_pk = 1;

  /* replace all characters by U or N */
  if (args_info.no_modified_given)
    ct2db_params.no_mod = 1;

  if (args_info.fasta_header_given)
    ct2db_params.id = strdup(args_info.fasta_header_arg);

  if (args_info.fh_from_filename_given) {
    switch (args_info.fh_from_filename_arg[0]) {
      case 'a':
        ct2db_params.id_from_filename = 2;
        break;
      case 'p':
        ct2db_params.id_from_filename = 3;
        break;
      default:
        ct2db_params.id_from_filename = 1;
        break;
    }
  }

  ct2db_params.filename_suffix = strdup(args_info.filename_suffix_arg);

  if ((args_info.filename_suffix_given) &&
      (strlen(args_info.filename_suffix_arg) == 0)) {
    free(ct2db_params.filename_suffix);
    ct2db_params.filename_suffix = NULL;
  }

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
          vrna_log_error("Unable to open %d. input file \"%s\" for reading", i + 1,
                         input_files[i]);

        vrna_log_info("Processing %d. input file \"%s\"",
                      i + 1,
                      input_files[i]);

        if (ct2db_params.id_from_filename) {
          char *fn = vrna_basename(input_files[i]);

          if (ct2db_params.filename_suffix) {
            const char *p = find_suffix(fn, ct2db_params.filename_suffix);

            if ((p) &&
                (p > fn))
              fn[p - fn] = '\0';
          }

          free(ct2db_params.id_prefix);
          free(ct2db_params.id_suffix);
          ct2db_params.id_prefix  = NULL;
          ct2db_params.id_suffix  = NULL;

          switch (ct2db_params.id_from_filename) {
            case 1:
              free(ct2db_params.id);
              ct2db_params.id = fn;
              break;
            case 2:
              ct2db_params.id_suffix = fn;
              break;
            case 3:
              ct2db_params.id_prefix = fn;
              break;
            default:
              break;
          }

          ct2db_params.id_num = ct2db_params.id_num_start;
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

  free(ct2db_params.filename_suffix);
  free(ct2db_params.id);
  free(ct2db_params.id_suffix);
  free(ct2db_params.id_prefix);

  free(input_files);

  if (vrna_log_fp() != stderr)
    fclose(vrna_log_fp());

  return EXIT_SUCCESS;
}


int
process_input(FILE        *fp,
              const char  *filename,
              parameters  *opt)
{
  char          *id, *sequence, *structure, *ss, *remainder;
  short         *pt_pk, *pt;
  size_t        bp_pk, bp;
  unsigned int  parse_options;

  parse_options = 0;
  remainder     = NULL;

  if (opt->verbose)
    parse_options |= VRNA_INPUT_VERBOSE;

  while (vrna_file_connect_read_record(fp, &id, &sequence, &structure, &remainder, parse_options)) {
    if (opt->no_pk) {
      if (structure) {
        pt_pk = vrna_ptable_from_string(structure, VRNA_BRACKETS_ANY);
        pt    = vrna_pt_pk_remove(pt_pk, 0);
        ss    = vrna_db_from_ptable(pt);

        bp_pk = 0;
        bp    = 0;

        for (size_t i = 1; i <= pt_pk[0]; i++) {
          if (pt_pk[i] > i)
            bp_pk++;

          if (pt[i] > i)
            bp++;
        }

        vrna_log_info("Removed %u of %u base pairs to make structure pseudo-knot free",
                      bp_pk - bp,
                      bp_pk);

        free(pt_pk);
        free(pt);
        free(structure);
      } else {
        ss = NULL;
      }
    } else {
      ss = structure;
    }

    if (opt->no_mod)
      remove_noncanonical_bases(sequence);

    if (opt->id) {
      free(id);
      id = strdup(opt->id);
    } else if (id == NULL) {
      id = strdup("seq");
    }

    if (opt->id_num > opt->id_num_start) {
      char *p = vrna_strdup_printf("%s_%u", id, opt->id_num);
      free(id);
      id = p;
    }

    if (opt->id_prefix) {
      char *p = vrna_strdup_printf("%s %s", opt->id_prefix, id);
      free(id);
      id = p;
    } else if (opt->id_suffix) {
      char *p = vrna_strdup_printf("%s %s", id, opt->id_suffix);
      free(id);
      id = p;
    }

    printf(">%s\n%s\n%s\n", id, sequence, ss);
    free(id);
    free(sequence);
    free(ss);
    opt->id_num++;
  }

  return 1;
}
