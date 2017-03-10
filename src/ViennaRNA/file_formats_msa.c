/*
 *  file_formats_msa.c
 *
 *  Various functions dealing with file formats for Multiple Sequence Alignments (MSA)
 *
 *  (c) 2016 Ronny Lorenz
 *
 *  ViennaRNA package
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
#include "ViennaRNA/file_utils.h"
#include "ViennaRNA/file_formats.h"
#include "ViennaRNA/file_formats_msa.h"

/*
 #################################
 # STATIC DECLARATIONS           #
 #################################
 */

typedef int (aln_parser_function)(FILE  *fp,
                                  char  ***names,
                                  char  ***aln,
                                  char  **id,
                                  char  **structure,
                                  int   verbosity);
typedef int (aln_writer_function)(FILE          *fp,
                                  const char    **names,
                                  const char    **aln,
                                  const char    *id,
                                  const char    *structure,
                                  const char    *source,
                                  unsigned int  options,
                                  int           verbosity);

typedef struct {
  unsigned int        code;
  aln_parser_function *parser;
  const char          *name;
} parsable;

typedef struct {
  unsigned int        code;
  aln_writer_function *writer;
  const char          *name;
} writable;

PRIVATE aln_parser_function parse_aln_stockholm;

PRIVATE aln_parser_function parse_aln_clustal;

PRIVATE aln_parser_function parse_aln_fasta;

PRIVATE aln_parser_function parse_aln_maf;

PRIVATE aln_writer_function write_aln_stockholm;

PRIVATE int
parse_fasta_alignment(FILE  *fp,
                      char  ***names,
                      char  ***aln,
                      int   verbosity);


PRIVATE int
parse_clustal_alignment(FILE  *clust,
                        char  ***names,
                        char  ***aln,
                        int   verbosity);


PRIVATE int
parse_stockholm_alignment(FILE  *fp,
                          char  ***aln,
                          char  ***names,
                          char  **id,
                          char  **structure,
                          int   verbosity);


PRIVATE int
parse_maf_alignment(FILE  *fp,
                    char  ***aln,
                    char  ***names,
                    int   verbosity);


PRIVATE int
write_stockholm_alignment(FILE          *fp,
                          const char    **names,
                          const char    **aln,
                          const char    *id,
                          const char    *structure,
                          const char    *source,
                          unsigned int  options,
                          int           verbosity);


PRIVATE int
check_alignment(const char  **names,
                const char  **aln,
                int         seq_num,
                int         verbosity);


PRIVATE void
free_msa_record(char  ***names,
                char  ***aln,
                char  **id,
                char  **structure);


PRIVATE void
add_sequence(const char *id,
             const char *seq,
             char       ***names,
             char       ***aln,
             int        seq_num);


PRIVATE void
endmarker_msa_record(char ***names,
                     char ***aln,
                     int  seq_num);


/*
 #################################
 # STATIC VARIABLES              #
 #################################
 */

/* number of known alignment parsers */
#define NUM_PARSERS 4
#define NUM_WRITERS 1

static parsable known_parsers[NUM_PARSERS] = {
  /* option, parser, name */
  { VRNA_FILE_FORMAT_MSA_STOCKHOLM, parse_aln_stockholm, "Stockholm 1.0 format" },
  { VRNA_FILE_FORMAT_MSA_CLUSTAL,   parse_aln_clustal,   "ClustalW format"      },
  { VRNA_FILE_FORMAT_MSA_FASTA,     parse_aln_fasta,     "FASTA format"         },
  { VRNA_FILE_FORMAT_MSA_MAF,       parse_aln_maf,       "MAF format"           }
};

static writable known_writers[NUM_WRITERS] = {
  /* option, parser, name */
  { VRNA_FILE_FORMAT_MSA_STOCKHOLM, write_aln_stockholm, "Stockholm 1.0 format" }
};

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC unsigned int
vrna_file_msa_detect_format(const char    *filename,
                            unsigned int  options)
{
  FILE          *fp;
  char          **names, **aln;
  unsigned int  format;
  int           i, r;
  long int      fp_position;

  names   = NULL;
  aln     = NULL;
  format  = VRNA_FILE_FORMAT_MSA_UNKNOWN;

  /* if no alignment file format(s) were specified we probe for all of them */
  if (options == 0)
    options = VRNA_FILE_FORMAT_MSA_DEFAULT;

  if (!(fp = fopen(filename, "r"))) {
    if (!(options & VRNA_FILE_FORMAT_MSA_SILENT))
      vrna_message_warning("Alignment file could not be opened!");

    return format;
  }

  r           = -1;
  fp_position = ftell(fp);

  for (i = 0; i < NUM_PARSERS; i++) {
    if ((options & known_parsers[i].code) && (known_parsers[i].parser)) {
      /* go back to beginning of file */
      if (!fseek(fp, fp_position, SEEK_SET)) {
        r = known_parsers[i].parser(fp, &names, &aln, NULL, NULL, -1);
        free_msa_record(&names, &aln, NULL, NULL);
        if (r > 0) {
          format = known_parsers[i].code;
          break;
        }
      } else {
        vrna_message_error("Something unexpected happened while parsing the alignment file");
      }
    }
  }

  fclose(fp);

  return format;
}


PUBLIC int
vrna_file_msa_read(const char   *filename,
                   char         ***names,
                   char         ***aln,
                   char         **id,
                   char         **structure,
                   unsigned int options)
{
  FILE      *fp;
  char      *line = NULL;
  int       i, n, seq_num, r, verb_level;
  long int  fp_position;

  verb_level  = 1; /* we default to be very verbose */
  seq_num     = 0;

  if (options & VRNA_FILE_FORMAT_MSA_QUIET)
    verb_level = 0;

  if (options & VRNA_FILE_FORMAT_MSA_SILENT)
    verb_level = -1;

  if (!(fp = fopen(filename, "r"))) {
    if (verb_level >= 0)
      vrna_message_warning("Alignment file could not be opened!");

    return seq_num;
  }

  if (names && aln) {
    *names  = NULL;
    *aln    = NULL;
  } else {
    return seq_num;
  }

  if (id)
    *id = NULL;

  if (structure)
    *structure = NULL;

  /* if no alignment file format was specified, lets try to guess it */
  if (options == 0)
    options = VRNA_FILE_FORMAT_MSA_DEFAULT;

  r           = -1;
  fp_position = ftell(fp);

  for (i = 0; i < NUM_PARSERS; i++) {
    if ((options & known_parsers[i].code) && (known_parsers[i].parser)) {
      /* go back to beginning of file */
      if (!fseek(fp, fp_position, SEEK_SET)) {
        r = known_parsers[i].parser(fp, names, aln, id, structure, verb_level);
        if (r > 0)
          break;
      } else {
        vrna_message_error("Something unexpected happened while parsing the alignment file");
      }
    }
  }

  fclose(fp);

  if (r == -1) {
    if (verb_level >= 0)
      vrna_message_warning("Alignment file parser is unknown (or not specified?)");
  } else {
    seq_num = r;

    if ((seq_num > 0) && (!(options & VRNA_FILE_FORMAT_MSA_NOCHECK))) {
      if (!check_alignment((const char **)(*names), (const char **)(*aln), seq_num, verb_level)) {
        if (verb_level >= 0)
          vrna_message_warning("Alignment did not pass sanity checks!");

        /* discard the data we've read! */
        free_msa_record(names, aln, id, structure);

        seq_num = 0;
      }
    }
  }

  return seq_num;
}


PUBLIC int
vrna_file_msa_read_record(FILE          *fp,
                          char          ***names,
                          char          ***aln,
                          char          **id,
                          char          **structure,
                          unsigned int  options)
{
  const char          *parser_name;
  int                 i, r, n, seq_num, verb_level;
  aln_parser_function *parser;

  verb_level  = 1; /* we default to be very verbose */
  seq_num     = 0;
  parser_name = NULL;
  parser      = NULL;

  if (options & VRNA_FILE_FORMAT_MSA_QUIET)
    verb_level = 0;

  if (options & VRNA_FILE_FORMAT_MSA_SILENT)
    verb_level = -1;

  if (!fp) {
    if (verb_level >= 0)
      vrna_message_warning("Can't read alignment from file pointer!");

    return seq_num;
  }

  if (names && aln) {
    *names  = NULL;
    *aln    = NULL;
  } else {
    return seq_num;
  }

  if (id)
    *id = NULL;

  if (structure)
    *structure = NULL;

  for (r = i = 0; i < NUM_PARSERS; i++) {
    if ((options & known_parsers[i].code) && (known_parsers[i].parser)) {
      if (!parser) {
        parser      = known_parsers[i].parser;
        parser_name = known_parsers[i].name;
      }

      r++;
    }
  }

  if (r == 0) {
    if (verb_level >= 0)
      vrna_message_warning("Did not find parser for specified MSA format!");
  } else {
    if ((r > 1) && (verb_level > 0))
      vrna_message_warning("More than one MSA format parser specified!\n"
                           "Using parser for %s", parser_name);

    seq_num = parser(fp, names, aln, id, structure, verb_level);

    if ((seq_num > 0) && (!(options & VRNA_FILE_FORMAT_MSA_NOCHECK))) {
      if (!check_alignment((const char **)(*names), (const char **)(*aln), seq_num, verb_level)) {
        if (verb_level >= 0)
          vrna_message_warning("Alignment did not pass sanity checks!");

        /* discard the data we've read! */
        free_msa_record(names, aln, id, structure);

        seq_num = -1;
      }
    }
  }

  return seq_num;
}


PUBLIC int
vrna_file_msa_write(const char    *filename,
                    const char    **names,
                    const char    **aln,
                    const char    *id,
                    const char    *structure,
                    const char    *source,
                    unsigned int  options)
{
  int ret, verb_level;

  verb_level  = 1;  /* we default to be very verbose */
  ret         = 0;  /* failure */

  if (options & VRNA_FILE_FORMAT_MSA_QUIET)
    verb_level = 0;

  if (options & VRNA_FILE_FORMAT_MSA_SILENT)
    verb_level = -1;

  if (filename && names && aln) {
    /* we require at least a filename, the sequence identifiers (names), and the alignment */
    FILE                *fp;
    int                 i, r, n, seq_num;
    const char          *writer_name;
    aln_writer_function *writer;

    seq_num     = 0;
    writer_name = NULL;
    writer      = NULL;

    /* find out the number of sequences in the alignment */
    for (seq_num = 0; aln[seq_num]; seq_num++) ;

    if (seq_num == 0) {
      if (verb_level >= 0)
        vrna_message_warning("Alignment did not pass sanity checks!");

      return ret;
    }

    /* check the alignment itself for consistency */
    if ((seq_num > 0) && (!(options & VRNA_FILE_FORMAT_MSA_NOCHECK))) {
      if (!check_alignment(names, aln, seq_num, verb_level)) {
        if (verb_level >= 0)
          vrna_message_warning("Alignment did not pass sanity checks!");

        return ret;
      }
    }

    /* detect output format */
    for (i = 0; i < NUM_WRITERS; i++) {
      if ((options & known_writers[i].code) && (known_writers[i].writer)) {
        if (!writer) {
          writer      = known_writers[i].writer;
          writer_name = known_writers[i].name;
        }

        r++;
      }
    }

    if (r == 0) {
      if (verb_level >= 0)
        vrna_message_warning("Did not find writer for specified MSA format!");
    } else {
      if ((r > 1) && (verb_level > 0))
        vrna_message_warning("More than one MSA format writer specified!\n"
                             "Using writer for %s",
                             writer_name);

      /* try opening the output stream */
      if (options & VRNA_FILE_FORMAT_MSA_APPEND)
        fp = fopen(filename, "a");
      else
        fp = fopen(filename, "w");

      if (!fp) {
        if (verb_level >= 0)
          vrna_message_warning("Alignment file could not be opened for writing!");

        return ret;
      }

      /* write output alignment to file */
      ret = writer(fp, names, aln, id, structure, source, options, verb_level);

      /* close output stream */
      fclose(fp);
    }
  } else if (verb_level >= 0) {
    vrna_message_warning("vrna_file_msa_write: insufficient input for writing anything!");
  }

  return ret;
}


PRIVATE int
parse_stockholm_alignment(FILE  *fp,
                          char  ***names,
                          char  ***aln,
                          char  **id,
                          char  **structure,
                          int   verbosity)
{
  char  *line = NULL;
  int   i, n, seq_num, seq_length, has_record;

  seq_num     = 0;
  seq_length  = 0;

  if (!fp) {
    if (verbosity >= 0)
      vrna_message_warning("Can't read from filepointer while parsing Stockholm formatted sequence alignment!");

    return -1;
  }

  if (names && aln) {
    *names  = NULL;
    *aln    = NULL;
  } else {
    return -1;
  }

  if (id)
    *id = NULL;

  if (structure)
    *structure = NULL;

  int inrecord = 0;
  while ((line = vrna_read_line(fp))) {
    if (strstr(line, "STOCKHOLM 1.0")) {
      inrecord    = 1;
      has_record  = 1;
      free(line);
      break;
    }

    free(line);
  }

  if (inrecord) {
    while ((line = vrna_read_line(fp))) {
      if (strncmp(line, "//", 2) == 0) {
        /* end of alignment */
        free(line);
        line = NULL;
        break;
      }

      n = (int)strlen(line);

      switch (*line) {
        /* we skip lines that start with whitespace */
        case ' ':
        case '\0':
          goto stockholm_next_line;

        /* Stockholm markup, or comment */
        case '#':
          if (strstr(line, "STOCKHOLM 1.0")) {
            if (verbosity >= 0)
              vrna_message_warning("Malformatted Stockholm record, missing // ?");

            /* drop everything we've read so far and start new, blank record */
            free_msa_record(names, aln, id, structure);

            seq_num = 0;
          } else if (strncmp(line, "#=GF", 4) == 0) {
            /* found feature markup */
            if ((id != NULL) && (strncmp(line, "#=GF ID", 7) == 0)) {
              *id = (char *)vrna_alloc(sizeof(char) * n);
              if (sscanf(line, "#=GF ID %s", *id) == 1) {
                *id = (char *)vrna_realloc(*id, sizeof(char) * (strlen(*id) + 1));
              } else {
                free(*id);
                *id = NULL;
              }
            }
          } else if (strncmp(line, "#=GC", 4) == 0) {
            /* found per-column annotation */
            if ((structure != NULL) && (strncmp(line, "#=GC SS_cons", 12) == 0)) {
              char *ss = (char *)vrna_alloc(sizeof(char) * n);
              if (sscanf(line, "#=GC SS_cons %s", ss) == 1) {
                *structure = (char *)vrna_realloc(*structure, sizeof(char) * (strlen(ss) + 1));
                strcpy(*structure, ss);
              }

              free(ss);
            }
          } else if (strncmp(line, "#=GS", 4) == 0) {
            /* found generic per-sequence annotation */
          } else if (strncmp(line, "#=GR", 4) == 0) {
            /* found generic per-Residue annotation */
          } else {
            /* may be comment? */
          }

          break;

        /* should be sequence */
        default:
        {
          int   tmp_l;
          char  *tmp_name = (char *)vrna_alloc(sizeof(char) * (n + 1));
          char  *tmp_seq  = (char *)vrna_alloc(sizeof(char) * (n + 1));
          if (sscanf(line, "%s %s", tmp_name, tmp_seq) == 2) {
            seq_num++;
            tmp_l = (int)strlen(tmp_seq);

            if (seq_num == 1) {
              seq_length = tmp_l;
            } else {
              /* check sequence length against first */
              if (seq_length != tmp_l) {
                if (verbosity >= 0)
                  vrna_message_warning("Discarding Stockholm record! Sequence lengths do not match.");

                /* drop everything we've read so far and abort parsing */
                free_msa_record(names, aln, id, structure);

                seq_num = 0;

                free(tmp_name);
                free(tmp_seq);
                free(line);
                line = NULL;

                goto stockholm_exit;
              }
            }

            add_sequence(tmp_name, tmp_seq,
                         names, aln,
                         seq_num);
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
    /*
     *  if (verbosity >= 0)
     *    vrna_message_warning("Did not find any Stockholm 1.0 formatted record!");
     */
    return -1;
  }

stockholm_exit:

  free(line);

  endmarker_msa_record(names, aln, seq_num);

  if ((seq_num > 0) && (verbosity > 0))
    vrna_message_info(stderr, "%d sequences; length of alignment %d.", seq_num, (int)strlen((*aln)[0]));

  return seq_num;
}


PRIVATE int
write_stockholm_alignment(FILE          *fp,
                          const char    **names,
                          const char    **aln,
                          const char    *id,
                          const char    *structure,
                          const char    *source,
                          unsigned int  options,
                          int           verbosity)
{
  int ret;

  ret = 1; /* success */

  if (fp) {
    int s, seq_num, l, longest_name;

    /* detect sequence number and longest sequence name in alignment */
    for (longest_name = seq_num = 0; names[seq_num]; seq_num++) {
      l = (int)strlen(names[seq_num]);
      if (l > longest_name)
        longest_name = l;
    }

    if (seq_num > 0) {
      /* print header */
      fprintf(fp, "# STOCKHOLM 1.0\n");

      if (id)   /* print the ID if available */
        fprintf(fp, "#=GF ID %s\n", id);

      if (structure) {
        /* print structure source if structure is present */
        if (!source)
          source = "ViennaRNA Package prediction";

        fprintf(fp, "#=GF SS %s\n", source);
        /*
         * in case we print a consensus structure, reset the longest_name
         * variable if longest name is shorter than '#=GC SS_cons' which
         * has 12 characters
         */
        if (longest_name < 12)
          longest_name = 12;
      }

      /* print actual alignment */
      for (s = 0; s < seq_num; s++)
        fprintf(fp, "%-*s  %s\n", longest_name, names[s], aln[s]);

      /* output reference sequence */
      char *cons = (options & VRNA_FILE_FORMAT_MSA_MIS) ? consens_mis(aln) : consensus(aln);
      fprintf(fp, "%-*s  %s\n", longest_name, "#=GC RF", cons);
      free(cons);

      /* print consensus structure */
      if (structure)
        fprintf(fp, "%-*s  %s\n", longest_name, "#=GC SS_cons", structure);

      /* print record-end marker */
      fprintf(fp, "//\n");
    }
  }

  return ret;
}


PRIVATE int
parse_fasta_alignment(FILE  *fp,
                      char  ***names,
                      char  ***aln,
                      int   verbosity)
{
  unsigned int  read_opt, rec_type;
  int           seq_num;
  char          *rec_id, *rec_sequence, **rec_rest;

  rec_id        = NULL;
  rec_sequence  = NULL;
  rec_rest      = NULL;
  seq_num       = 0;
  read_opt      = VRNA_INPUT_NO_REST; /* read sequence and header information only */

  /* read until EOF or user abort */
  while (
    !((rec_type = vrna_file_fasta_read_record(&rec_id, &rec_sequence, &rec_rest, fp, read_opt))
      & (VRNA_INPUT_ERROR | VRNA_INPUT_QUIT))) {
    if (rec_id) {
      /* valid FASTA entry */
      seq_num++;

      char *id = (char *)vrna_alloc(sizeof(char) * strlen(rec_id));
      (void)sscanf(rec_id, ">%s", id);
      vrna_seq_toupper(rec_sequence);

      add_sequence(id, rec_sequence,
                   names, aln,
                   seq_num);

      free(id);
    }

    free(rec_id);
    free(rec_sequence);
    free(rec_rest);
  }

  free(rec_id);
  free(rec_sequence);
  free(rec_rest);

  endmarker_msa_record(names, aln, seq_num);

  if (seq_num > 0) {
    if (verbosity > 0)
      vrna_message_info(stderr, "%d sequences; length of alignment %d.", seq_num, (int)strlen((*aln)[0]));
  } else {
    /*
     *  if (verbosity >= 0)
     *    vrna_message_warning("Did not find any FASTA formatted record!");
     */
    return -1;
  }

  return seq_num;
}


PRIVATE int
parse_clustal_alignment(FILE  *clust,
                        char  ***names,
                        char  ***aln,
                        int   verbosity)
{
  char  *line, *name, *seq;
  int   n, r, nn = 0, seq_num = 0, i;

  if ((line = vrna_read_line(clust)) == NULL)
    return -1;

  if (strncmp(line, "CLUSTAL", 7) != 0) {
    if (verbosity >= 0)
      vrna_message_warning("This doesn't look like a CLUSTALW file, sorry");

    free(line);
    return -1;
  }

  free(line);
  line = vrna_read_line(clust);

  while (line != NULL) {
    n = strlen(line);

    if ((n < 4) || isspace((int)line[0])) {
      /* skip non-sequence line */
      free(line);
      line  = vrna_read_line(clust);
      nn    = 0;  /* reset sequence number */
      continue;
    }

    /* skip comments */
    if (line[0] == '#') {
      free(line);
      line = vrna_read_line(clust);
      continue;
    }

    seq   = (char *)vrna_alloc(sizeof(char) * (n + 1));
    name  = (char *)vrna_alloc(sizeof(char) * (n + 1));
    if (sscanf(line, "%s %s", name, seq) == 2) {
      /* realloc to actual sizes */
      seq   = (char *)vrna_realloc(seq, sizeof(char) * (strlen(seq) + 1));
      name  = (char *)vrna_realloc(name, sizeof(char) * (strlen(name) + 1));
      for (i = 0; i < strlen(seq); i++)
        if (seq[i] == '.') /* replace '.' gaps with '-' */
          seq[i] = '-';

      /* convert sequence to uppercase letters */
      vrna_seq_toupper(seq);

      if (nn == seq_num) {
        /* first time */
        add_sequence(name, seq,
                     names, aln,
                     nn + 1);
      } else {
        if (strcmp(name, (*names)[nn]) != 0) {
          /* name doesn't match */
          if (verbosity >= 0)
            vrna_message_warning("Sorry, your file is messed up (inconsitent seq-names)");

          free(line);
          free(seq);
          return 0;
        }

        (*aln)[nn] = (char *)vrna_realloc((*aln)[nn], strlen(seq) + strlen((*aln)[nn]) + 1);
        strcat((*aln)[nn], seq);
      }

      nn++;
      if (nn > seq_num)
        seq_num = nn;

      free(seq);
      free(name);
    }

    free(line);

    line = vrna_read_line(clust);
  }

  endmarker_msa_record(names, aln, seq_num);

  if ((seq_num > 0) && (verbosity > 0))
    vrna_message_info(stderr, "%d sequences; length of alignment %d.", seq_num, (int)strlen((*aln)[0]));

  return seq_num;
}


PRIVATE int
parse_maf_alignment(FILE  *fp,
                    char  ***names,
                    char  ***aln,
                    int   verbosity)
{
  char  *line = NULL, *tmp_name, *tmp_sequence, strand;
  int   i, n, seq_num, seq_length, start, length, src_length;

  seq_num     = 0;
  seq_length  = 0;

  if (!fp) {
    if (verbosity >= 0)
      vrna_message_warning("Can't read from filepointer while parsing MAF formatted sequence alignment!");

    return -1;
  }

  if (names && aln) {
    *names  = NULL;
    *aln    = NULL;
  } else {
    return -1;
  }

  int inrecord = 0;
  while ((line = vrna_read_line(fp))) {
    if (*line == 'a') {
      if ((line[1] == '\0') || isspace(line[1])) {
        inrecord = 1;
        free(line);
        break;
      }
    }

    free(line);
  }

  if (inrecord) {
    while ((line = vrna_read_line(fp))) {
      n = (int)strlen(line);

      switch (*line) {
        case '#': /* comment */
          break;

        case 's': /* a sequence within the alignment block */
          tmp_name      = (char *)vrna_alloc(sizeof(char) * n);
          tmp_sequence  = (char *)vrna_alloc(sizeof(char) * n);
          if (sscanf(line, "s %s %d %d %c %d %s",
                     tmp_name,
                     &start,
                     &length,
                     &strand,
                     &src_length,
                     tmp_sequence) == 6) {
            seq_num++;
            tmp_name      = (char *)vrna_realloc(tmp_name, sizeof(char) * (strlen(tmp_name) + 1));
            tmp_sequence  = (char *)vrna_realloc(tmp_sequence, sizeof(char) * (strlen(tmp_sequence) + 1));

            vrna_seq_toupper(tmp_sequence);

            add_sequence(tmp_name, tmp_sequence,
                         names, aln,
                         seq_num);

            free(tmp_name);
            free(tmp_sequence);
            break;
          }

          free(tmp_name);
          free(tmp_sequence);
        /* all through */

        default: /* something else that ends the block */
          free(line);
          goto maf_exit;
      }

      free(line);
    }
  } else {
    /*
     *  if (verbosity >= 0)
     *    vrna_message_warning("Did not find any MAF formatted record!");
     */
    return -1;
  }

maf_exit:

  endmarker_msa_record(names, aln, seq_num);

  if ((seq_num > 0) && (verbosity > 0))
    vrna_message_info(stderr, "%d sequences; length of alignment %d.", seq_num, (int)strlen((*aln)[0]));

  return seq_num;
}


PRIVATE void
free_msa_record(char  ***names,
                char  ***aln,
                char  **id,
                char  **structure)
{
  int s, i;

  s = 0;
  if (aln && (*aln))
    for (; (*aln)[s]; s++) ;

  if (id != NULL) {
    free(*id);
    *id = NULL;
  }

  if (structure != NULL) {
    free(*structure);
    *structure = NULL;
  }

  for (i = 0; i < s; i++) {
    free((*names)[i]);
    free((*aln)[i]);
  }

  if (names && (*names)) {
    free(*names);
    *names = NULL;
  }

  if (aln && (*aln)) {
    free(*aln);
    *aln = NULL;
  }
}


PRIVATE int
parse_aln_stockholm(FILE  *fp,
                    char  ***names,
                    char  ***aln,
                    char  **id,
                    char  **structure,
                    int   verbosity)
{
  return parse_stockholm_alignment(fp, names, aln, id, structure, verbosity);
}


PRIVATE int
parse_aln_clustal(FILE  *fp,
                  char  ***names,
                  char  ***aln,
                  char  **id,
                  char  **structure,
                  int   verbosity)
{
  /* clustal format doesn't contain id's or structure information */
  if (id)
    *id = NULL;

  if (structure)
    *structure = NULL;

  return parse_clustal_alignment(fp, names, aln, verbosity);
}


PRIVATE int
parse_aln_fasta(FILE  *fp,
                char  ***names,
                char  ***aln,
                char  **id,
                char  **structure,
                int   verbosity)
{
  /* fasta alignments do not contain an id, or structure information */
  if (id)
    *id = NULL;

  if (structure)
    *structure = NULL;

  return parse_fasta_alignment(fp, names, aln, verbosity);
}


PRIVATE int
parse_aln_maf(FILE  *fp,
              char  ***names,
              char  ***aln,
              char  **id,
              char  **structure,
              int   verbosity)
{
  /* MAF alignments do not contain an id, or structure information */
  if (id)
    *id = NULL;

  if (structure)
    *structure = NULL;

  return parse_maf_alignment(fp, names, aln, verbosity);
}


PRIVATE int
write_aln_stockholm(FILE          *fp,
                    const char    **names,
                    const char    **aln,
                    const char    *id,
                    const char    *structure,
                    const char    *source,
                    unsigned int  options,
                    int           verbosity)
{
  return write_stockholm_alignment(fp, names, aln, id, structure, source, options, verbosity);
}


PRIVATE void
add_sequence(const char *id,
             const char *seq,
             char       ***names,
             char       ***aln,
             int        seq_num)
{
  (*names)              = (char **)vrna_realloc(*names, sizeof(char *) * (seq_num));
  (*names)[seq_num - 1] = strdup(id);
  (*aln)                = (char **)vrna_realloc(*aln, sizeof(char *) * (seq_num));
  (*aln)[seq_num - 1]   = strdup(seq);
}


PRIVATE void
append_sequence(char  *seq,
                char  **aln,
                int   seq_num)
{
}


PRIVATE void
endmarker_msa_record(char ***names,
                     char ***aln,
                     int  seq_num)
{
  /*
   * append additional entry in 'aln' and 'names' pointing to NULL (this may be
   * used as an indication for the end of the sequence list)
   */
  if (seq_num > 0) {
    (*aln)            = (char **)vrna_realloc(*aln, sizeof(char *) * (seq_num + 1));
    (*names)          = (char **)vrna_realloc(*names, sizeof(char *) * (seq_num + 1));
    (*aln)[seq_num]   = NULL;
    (*names)[seq_num] = NULL;
  }
}


PRIVATE int
check_alignment(const char  **names,
                const char  **aln,
                int         seq_num,
                int         verbosity)
{
  int i, j, l, pass = 1;

  /* check for unique names */
  for (i = 0; i < seq_num; i++) {
    for (j = i + 1; j < seq_num; j++) {
      if (!strcmp(names[i], names[j])) {
        if (verbosity >= 0)
          vrna_message_warning("Sequence IDs in input alignment are not unique!");

        pass = 0;
      }
    }
  }

  /* check for equal lengths of sequences */
  l = (int)strlen(aln[0]);
  for (i = 1; i < seq_num; i++)
    if ((int)strlen(aln[i]) != l) {
      if (verbosity >= 0)
        vrna_message_warning("Sequence lengths in input alignment do not match!");

      pass = 0;
    }

  return pass;
}
