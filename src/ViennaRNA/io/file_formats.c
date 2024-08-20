/*
 *  io/file_formats.c
 *
 *  Various functions dealing with file formats for RNA sequences, structures, and alignments
 *
 *  (c) 2014 Ronny Lorenz
 *
 *  Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "json/json.h"

#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/io/utils.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/io/file_formats.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */

typedef struct {
  unsigned int  effective_length;
  unsigned int  length;
  char          *id;
  char          *sequence;
  unsigned int  sequence_pos;
  short         *pt;
  unsigned int  strands;
  unsigned int  *actual_pos;
} ct_data;


PRIVATE char          *inbuf  = NULL;
PRIVATE char          *inbuf2 = NULL;
PRIVATE unsigned int  typebuf = 0;

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE unsigned int
read_multiple_input_lines(char          **string,
                          FILE          *file,
                          unsigned int  option);


PRIVATE void
elim_trailing_ws(char *string);


PRIVATE INLINE ct_data *
init_ct_data(unsigned int n);


PRIVATE INLINE void
resize_ct_data(ct_data      *dat,
               unsigned int n);


PRIVATE INLINE int
finalize_ct_data(ct_data *dat);


PRIVATE INLINE ct_data*
process_ct_header(unsigned int  n,
                  size_t        num_tok,
                  const char    **tok);

PRIVATE INLINE int
process_ct_nt_line(ct_data       *data,
                   unsigned int  i,
                   char          nucleotide,
                   unsigned int  predecessor,
                   unsigned int  j,
                   unsigned int  actual_i);

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */

/* eliminate whitespaces/non-printable characters at the end of a character string */
PRIVATE void
elim_trailing_ws(char *string)
{
  int i, l = strlen(string);

  for (i = l - 1; i >= 0; i--) {
    if (isspace(string[i]) || (!isprint(string[i])))
      continue;

    break;
  }
  string[(i >= 0) ? (i + 1) : 0] = '\0';
}


PUBLIC void
vrna_file_helixlist(const char  *seq,
                    const char  *db,
                    float       energy,
                    FILE        *file)
{
  int       s;
  short     *pt;
  vrna_hx_t *list;
  FILE      *out;

  if (strlen(seq) != strlen(db)) {
    vrna_log_warning("vrna_file_helixlist: "
                         "sequence and structure have unequal length (%d vs. %d)!",
                         strlen(seq),
                         strlen(db));
    return;
  }

  out   = (file) ? file : stdout;
  pt    = vrna_ptable(db);
  list  = vrna_hx_from_ptable(pt);

  fprintf(out, "%s\t%6.2f\n", seq, energy);
  for (s = 0; list[s].length > 0; s++)
    fprintf(out, "%d\t%d\t%d\n", list[s].start, list[s].end, list[s].length);

  free(pt);
  free(list);
}


PUBLIC void
vrna_file_connect(const char  *seq,
                  const char  *db,
                  float       energy,
                  const char  *identifier,
                  FILE        *file)
{
  int   i, power_d;
  FILE  *out = (file) ? file : stdout;

  if (strlen(seq) != strlen(db)) {
    vrna_log_warning("vrna_file_connect: "
                         "sequence and structure have unequal length (%d vs. %d)!",
                         strlen(seq),
                         strlen(db));
    return;
  }

  short *pt = vrna_ptable(db);

  for (power_d = 0; pow(10, power_d) <= (int)strlen(seq); power_d++);

  /*
   * Connect table file format looks like this:
   *
   * 300  ENERGY = 7.0  example
   *  1 G       0    2   22    1
   *  2 G       1    3   21    2
   *
   * where the headerline is followed by 6 columns with:
   * 1. Base number: index n
   * 2. Base (A, C, G, T, U, X)
   * 3. Index n-1  (0 if first nucleotide)
   * 4. Index n+1  (0 if last nucleotide)
   * 5. Number of the base to which n is paired. No pairing is indicated by 0 (zero).
   * 6. Natural numbering.
   */

  /* print header */
  fprintf(out, "%d  ENERGY = %6.2f", (int)strlen(seq), energy);
  if (identifier)
    fprintf(out, "  %s\n", identifier);

  /*
   * print structure information except for last line
   * TODO: modify the structure information for cofold
   */
  for (i = 0; i < strlen(seq) - 1; i++) {
    fprintf(out, "%*d %c %*d %*d %*d %*d\n",
            power_d, i + 1,               /* nucleotide index */
            (char)toupper(seq[i]),        /* nucleotide char */
            power_d, i,                   /* nucleotide predecessor index */
            power_d, i + 2,               /* nucleotide successor index */
            power_d, pt[i + 1],           /* pairing partner index */
            power_d, i + 1);              /* nucleotide natural numbering */
  }
  /* print last line */
  fprintf(out, "%*d %c %*d %*d %*d %*d\n",
          power_d, i + 1,
          (char)toupper(seq[i]),
          power_d, i,
          power_d, 0,
          power_d, pt[i + 1],
          power_d, i + 1);

  /* clean up */
  free(pt);
  fflush(out);
}


PUBLIC void
vrna_file_bpseq(const char  *seq,
                const char  *db,
                FILE        *file)
{
  int   i;
  FILE  *out = (file) ? file : stdout;

  if (strlen(seq) != strlen(db)) {
    vrna_log_warning("vrna_file_bpseq: "
                         "sequence and structure have unequal length (%d vs. %d)!",
                         strlen(seq),
                         strlen(db));
    return;
  }

  short *pt = vrna_ptable(db);

  for (i = 1; i <= pt[0]; i++)
    fprintf(out, "%d %c %d\n", i, (char)toupper(seq[i - 1]), pt[i]);

  /* clean up */
  free(pt);
  fflush(out);
}


PUBLIC void
vrna_file_json(const char *seq,
               const char *db,
               double     energy,
               const char *identifier,
               FILE       *file)
{
  FILE      *out = (file) ? file : stdout;

  JsonNode  *data = json_mkobject();
  JsonNode  *value;

  if (identifier) {
    value = json_mkstring(identifier);
    json_append_member(data, "id", value);
  }

  value = json_mkstring(seq);
  json_append_member(data, "sequence", value);

  value = json_mknumber(energy);
  json_append_member(data, "mfe", value);

  value = json_mkstring(db);
  json_append_member(data, "structure", value);


  fprintf(out, "%s\n", json_stringify(data, "\t"));

  fflush(out);
}


PRIVATE unsigned int
read_multiple_input_lines(char          **string,
                          FILE          *file,
                          unsigned int  option)
{
  char  *line;
  int   i, l;
  int   state       = 0;
  int   str_length  = 0;
  FILE  *in         = (file) ? file : stdin;

  line    = (inbuf2) ? inbuf2 : vrna_read_line(in);
  inbuf2  = NULL;
  do{
    /*
     * read lines until informative data appears or
     * report an error if anything goes wrong
     */
    if (!line)
      return VRNA_INPUT_ERROR;

    l = (int)strlen(line);

    /* eliminate whitespaces at the end of the line read */
    if (!(option & VRNA_INPUT_NO_TRUNCATION))
      elim_trailing_ws(line);

    l           = (int)strlen(line);
    str_length  = (*string) ? (int)strlen(*string) : 0;

    switch (*line) {
      case  '@':    /* user abort */
        if (state)
          inbuf2 = line;
        else
          free(line);

        return (state == 2) ? VRNA_INPUT_CONSTRAINT : (state ==
                                                       1) ? VRNA_INPUT_SEQUENCE : VRNA_INPUT_QUIT;

      case  '\0':   /* empty line */
        if (option & VRNA_INPUT_NOSKIP_BLANK_LINES) {
          if (state)
            inbuf2 = line;
          else
            free(line);

          return (state == 2) ? VRNA_INPUT_CONSTRAINT : (state ==
                                                         1) ? VRNA_INPUT_SEQUENCE :
                 VRNA_INPUT_BLANK_LINE;
        }

        break;

      case  '#':  /* fall through */
      case  '%':  /* fall through */
      case  ';':  /* fall through */
      case  '/':  /* fall through */
      case  '*':  /* fall through */
      case ' ':   /* comments */
        if (option & VRNA_INPUT_NOSKIP_COMMENTS) {
          if (state)
            inbuf2 = line;
          else
            *string = line;

          return (state == 2) ? VRNA_INPUT_CONSTRAINT : (state ==
                                                         1) ? VRNA_INPUT_SEQUENCE :
                 VRNA_INPUT_COMMENT;
        }

        break;

      case  '>':  /* fasta header */
        if (state)
          inbuf2 = line;
        else
          *string = line;

        return (state == 2) ? VRNA_INPUT_CONSTRAINT : (state ==
                                                       1) ? VRNA_INPUT_SEQUENCE :
               VRNA_INPUT_FASTA_HEADER;

      case 'x':   /* fall through */
      case 'e':   /* fall through */
      case 'l':   /* fall through */
      case '&':   /* seems to be a constraint or line starting with second sequence for dimer calculations */
        i = 1;
        /* lets see if this assumption holds for the complete line */
        while ((line[i] == 'x') || (line[i] == 'e') || (line[i] == 'l'))
          i++;
        /* lines solely consisting of 'x's, 'e's or 'l's will be considered as structure constraint */

        if (
          ((line[i] > 64) && (line[i] < 91))                /* A-Z */
          || ((line[i] > 96) && (line[i] < 123))            /* a-z */
          ) {
          if (option & VRNA_INPUT_FASTA_HEADER) {
            /* are we in structure mode? Then we remember this line for the next round */
            if (state == 2) {
              inbuf2 = line;
              return VRNA_INPUT_CONSTRAINT;
            } else {
              *string = (char *)vrna_realloc(*string, sizeof(char) * (str_length + l + 1));
              memcpy(*string + str_length,
                     line,
                     sizeof(char) * l);
              (*string)[str_length + l] = '\0';
              state                     = 1;
            }

            break;
          }
          /* otherwise return line read */
          else {
            *string = line;
            return VRNA_INPUT_SEQUENCE;
          }
        }

      /*
       * mmmh? it really seems to be a constraint
       * fallthrough
       */
      case  '<':  /* fall through */
      case  '.':  /* fall through */
      case  '|':  /* fall through */
      case  '(':  /* fall through */
      case ')':   /* fall through */
      case '[':   /* fall through */
      case ']':   /* fall through */
      case '{':   /* fall through */
      case '}':   /* fall through */
      case ',':   /* fall through */
      case '+':   /* fall through */
        /*
         * seems to be a structure or a constraint
         * either we concatenate this line to one that we read previously
         */
        if (option & VRNA_INPUT_FASTA_HEADER) {
          if (state == 1) {
            inbuf2 = line;
            return VRNA_INPUT_SEQUENCE;
          } else {
            *string = (char *)vrna_realloc(*string, sizeof(char) * (str_length + l + 1));
            memcpy(*string + str_length,
                   line,
                   sizeof(char) * l);
            (*string)[str_length + l] = '\0';
            state                     = 2;
          }
        }
        /* or we return it as it is */
        else {
          *string = line;
          return VRNA_INPUT_CONSTRAINT;
        }

        break;
      default:
        if (option & VRNA_INPUT_FASTA_HEADER) {
          /* are we already in sequence mode? */
          if (state == 2) {
            inbuf2 = line;
            return VRNA_INPUT_CONSTRAINT;
          } else {
            *string = (char *)vrna_realloc(*string, sizeof(char) * (str_length + l + 1));
            memcpy(*string + str_length,
                   line,
                   sizeof(char) * l);
            (*string)[str_length + l] = '\0';
            state                     = 1;
          }
        }
        /* otherwise return line read */
        else {
          *string = line;
          return VRNA_INPUT_SEQUENCE;
        }
    }
    free(line);
    line = vrna_read_line(in);
  } while (line);

  return (state == 2) ? VRNA_INPUT_CONSTRAINT : (state ==
                                                 1) ? VRNA_INPUT_SEQUENCE : VRNA_INPUT_ERROR;
}


PUBLIC unsigned int
vrna_file_fasta_read_record(char          **header,
                            char          **sequence,
                            char          ***rest,
                            FILE          *file,
                            unsigned int  options)
{
  unsigned int  input_type, return_type, tmp_type;
  int           rest_count;
  char          *input_string;

  rest_count    = 0;
  return_type   = tmp_type = 0;
  input_string  = *header = *sequence = NULL;
  *rest         = (char **)vrna_alloc(sizeof(char *));

  /* remove unnecessary option flags from options variable... */
  options &= ~VRNA_INPUT_FASTA_HEADER;

  /* read first input or last buffered input */
  if (typebuf) {
    input_type    = typebuf;
    input_string  = inbuf;
    typebuf       = 0;
    inbuf         = NULL;
  } else {
    input_type = read_multiple_input_lines(&input_string, file, options);
  }

  if (input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR))
    return input_type;

  /* skip everything until we read either a fasta header or a sequence */
  while (input_type & (VRNA_INPUT_MISC | VRNA_INPUT_CONSTRAINT | VRNA_INPUT_BLANK_LINE)) {
    free(input_string);
    input_string  = NULL;
    input_type    = read_multiple_input_lines(&input_string, file, options);
    if (input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR))
      return input_type;
  }

  if (input_type & VRNA_INPUT_FASTA_HEADER) {
    return_type   |= VRNA_INPUT_FASTA_HEADER; /* remember that we've read a fasta header */
    *header       = input_string;
    input_string  = NULL;
    /* get next data-block with fasta support if not explicitely forbidden by VRNA_INPUT_NO_SPAN */
    input_type = read_multiple_input_lines(
      &input_string,
      file,
      ((options & VRNA_INPUT_NO_SPAN) ? 0 : VRNA_INPUT_FASTA_HEADER) | options
      );
    if (input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR))
      return return_type | input_type;
  }

  if (input_type & VRNA_INPUT_SEQUENCE) {
    return_type   |= VRNA_INPUT_SEQUENCE; /* remember that we've read a sequence */
    *sequence     = input_string;
    input_string  = NULL;
  } else {
    vrna_log_warning("vrna_file_fasta_read_record: "
                         "sequence input missing!");
    return VRNA_INPUT_ERROR;
  }

  /* read the rest until we find user abort, EOF, new sequence or new fasta header */
  if (!(options & VRNA_INPUT_NO_REST)) {
    options   |= VRNA_INPUT_NOSKIP_COMMENTS; /* allow commetns to appear in rest output */
    tmp_type  = VRNA_INPUT_QUIT | VRNA_INPUT_ERROR | VRNA_INPUT_SEQUENCE | VRNA_INPUT_FASTA_HEADER;
    if (options & VRNA_INPUT_NOSKIP_BLANK_LINES)
      tmp_type |= VRNA_INPUT_BLANK_LINE;

    while (!((input_type = read_multiple_input_lines(&input_string, file, options)) & tmp_type)) {
      *rest                   = vrna_realloc(*rest, sizeof(char **) * (++rest_count + 1));
      (*rest)[rest_count - 1] = input_string;
      input_string            = NULL;
    }
    /*
     * if(input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR)) return input_type;
     */

    /*  finished reading everything...
     *   we now put the last line into the buffer if necessary
     *   since it should belong to the next record
     */
    inbuf   = input_string;
    typebuf = input_type;
  }

  (*rest)[rest_count] = NULL;
  return return_type;
}


PUBLIC char *
vrna_extract_record_rest_structure(const char   **lines,
                                   unsigned int length,
                                   unsigned int options)
{
  char  *structure = NULL;
  int   r, i, l, cl, stop;
  char  *c;

  if (lines) {
    for (r = i = stop = 0; lines[i]; i++) {
      l = (int)strlen(lines[i]);
      c = (char *)vrna_alloc(sizeof(char) * (l + 1));
      (void)sscanf(lines[i], "%s", c);
      cl = (int)strlen(c);

      /* line commented out ? */
      if ((*c == '#') || (*c == '%') || (*c == ';') || (*c == '/') || (*c == '*' || (*c == '\0'))) {
        /* skip leading comments only, i.e. do not allow comments inside the constraint */
        if (!r)
          continue;
        else
          break;
      }

      /* append the structure part to the output */
      r         += cl + 1;
      structure = (char *)vrna_realloc(structure, r * sizeof(char));
      strcat(structure, c);
      free(c);
      /* stop if the assumed structure length has been reached */
      if ((length > 0) && (r - 1 == length))
        break;

      /* stop if not allowed to read from multiple lines */
      if (!(options & VRNA_OPTION_MULTILINE))
        break;
    }
  }

  return structure;
}


PUBLIC void
vrna_extract_record_rest_constraint(char          **cstruc,
                                    const char    **lines,
                                    unsigned int  option)
{
  *cstruc =
    vrna_extract_record_rest_structure(lines,
                                       0,
                                       option |
                                       (option &
                                        VRNA_CONSTRAINT_MULTILINE) ? VRNA_OPTION_MULTILINE : 0);
}


PUBLIC int
vrna_file_SHAPE_read(const char *file_name,
                     int        length,
                     double     default_value,
                     char       *sequence,
                     double     *values)
{
  FILE  *fp;
  char  *line;
  int   i;
  int   count = 0;

  if (!file_name)
    return 0;

  if (!(fp = fopen(file_name, "r"))) {
    vrna_log_warning("SHAPE data file could not be opened");
    return 0;
  }

  for (i = 0; i < length; ++i) {
    sequence[i]   = 'N';
    values[i + 1] = default_value;
  }

  sequence[length] = '\0';

  while ((line = vrna_read_line(fp))) {
    int           position;
    unsigned char nucleotide    = 'N';
    double        reactivity    = default_value;
    char          *second_entry = 0;
    char          *third_entry  = 0;
    char          *c;

    if (sscanf(line, "%d", &position) != 1) {
      free(line);
      continue;
    }

    if (position <= 0 || position > length) {
      vrna_log_warning("Provided SHAPE data outside of sequence scope");
      fclose(fp);
      free(line);
      return 0;
    }

    for (c = line + 1; *c; ++c) {
      if (isspace(*(c - 1)) && !isspace(*c)) {
        if (!second_entry) {
          second_entry = c;
        } else {
          third_entry = c;
          break;
        }
      }
    }

    if (second_entry) {
      if (third_entry) {
        sscanf(second_entry, "%c", &nucleotide);
        sscanf(third_entry, "%lf", &reactivity);
      } else if (sscanf(second_entry, "%lf", &reactivity) != 1) {
        sscanf(second_entry, "%c", &nucleotide);
      }
    }

    sequence[position - 1]  = nucleotide;
    values[position]        = reactivity;
    ++count;

    free(line);
  }

  fclose(fp);

  if (!count) {
    vrna_log_warning("SHAPE data file is empty");
    return 0;
  }

  return 1;
}


PUBLIC int
vrna_file_RNAstrand_db_read_record(FILE         *fp,
                                   char         **name_p,
                                   char         **sequence_p,
                                   char         **structure_p,
                                   char         **source_p,
                                   char         **fname_p,
                                   char         **id_p,
                                   unsigned int options)
{
  char          *ptr;
  unsigned int  state = 0;
  int           ret   = 0;
  size_t        line_length;
  size_t        seq_len = 0;
  size_t        struct_len = 0;

  char *line;

  *name_p = *sequence_p = *structure_p = *source_p = *fname_p = *id_p = NULL;

  while ((line = vrna_read_line(fp))) {
    /* skip lines starting with whitespace */
    if ((*line == '\0') || isspace(*line)) {
      /* whitespace should separate blocks from each other, but we allow for whitespace line before record */
      if (state > 0)
        state++;

      if (state > 3)
        break;

      continue;
    }

    if (state > 3)
      break; /* we should have read the entire record */

    line_length = strlen(line);

    if (*line == '#') {
      state = 1;
      /* still in header */
      if (!strncmp(line, "# File", 6)) {
        char *name = (char *)vrna_alloc(sizeof(char) * (line_length - 5));
        if (sscanf(line, "# File %s", name) == 1) {
          *name_p = name;
        } else {
          free(name);
          goto RNAstrand_parser_end;
        }
      } else if (!strncmp(line, "# External source:", 18)) {
        char *source = (char *)vrna_alloc(sizeof(char) * (line_length - 18));
        size_t pos        = 19;
        size_t source_len = 0;
        while (line[pos] != '\0') {
          /* only read until first comma */
          if (line[pos] == ',')
            break;
          source_len++;
          pos++;
        }
        if (source_len > 0) {
          source = (char *)vrna_realloc(source, sizeof(char) * (source_len + 1));
          strncpy(source, line + 19, sizeof(char) * source_len);
          source[source_len] = '\0';
          *source_p = source;

          /* try detecting the 'file name' */
          if ((ptr = strstr(line + 19, "file name:"))) {
            pos = 11;
            source_len = 0;
            while (ptr[pos] != '\0') {
              /* only read until next comma */
              if (ptr[pos] == ',')
                break;
              source_len++;
              pos++;
            }
            if (source_len > 0) {
              *fname_p = (char *)vrna_alloc(sizeof(char) * (source_len + 1));
              strncpy(*fname_p, ptr + 11, sizeof(char) * source_len);
              (*fname_p)[source_len] = '\0';
            }
          }

          /* try detecting 'ID' */
          if ((ptr = strstr(line + 19, "ID:"))) {
            pos = 4;
            source_len = 0;
            while (ptr[pos] != '\0') {
              /* only read until next comma */
              if (ptr[pos] == ',')
                break;
              source_len++;
              pos++;
            }
            if (source_len > 0) {
              *id_p = (char *)vrna_alloc(sizeof(char) * (source_len + 1));
              strncpy(*id_p, ptr + 4, sizeof(char) * source_len);
              (*id_p)[source_len] = '\0';
            }
          }
        } else {
          free(source);
          goto RNAstrand_parser_end;
        }
      }
    } else {
      /*
          from here on, we need some heuristic to actually decide whether we
          read the sequence, or a dot-parenthesis structure.
          In fact, the sequence may contain some special characters that are used
          to annotate modified nucleotides. On the other hand, the dot-parenthesis
          structure may consist of characters other than the usual parenthesis/bracket
          characters '()', '[]', '{}', and '.'. Pseudo-knots sometimes are annotated
          as matching upper-case/lowercase letters from the alphabet [A-Z].
      */

      if (state == 2) {
        /* lets first count the number of alphabetic characters and brackets/parenthesis */
        size_t alpha      = 0;
        size_t dot_parent = 0;
        for (size_t i = 0; i < line_length; i++) {
          if ((isalpha(line[i])) ||
              (line[i] == '~')) {
            alpha++;
          } else if ((line[i] == '.') ||
                     (line[i] == '(') ||
                     (line[i] == ')') ||
                     (line[i] == '[') ||
                     (line[i] == ']') ||
                     (line[i] == '{') ||
                     (line[i] == '}') ||
                     (line[i] == '<') ||
                     (line[i] == '>')) {
            dot_parent++;
          }
        }

        /*
            here, we simply assume that if the entire line
            looks like dot-parenthesis, or at least the majority
            of characters are more dot-parenthesis-like, we are
            actually looking at a structure.
        */
        if ((dot_parent == line_length) ||
            ((alpha != line_length) &&
              (dot_parent > alpha))) {
            state = 3;
        }

        if (state == 2) {
          /* still in sequence scan mode? */
          size_t tmp_len = seq_len + line_length + 1;
          *sequence_p = (char *)vrna_realloc(*sequence_p, sizeof(char) * tmp_len);
          memcpy(*sequence_p + seq_len, line, sizeof(char) * line_length);
          (*sequence_p)[tmp_len - 1] = '\0';
          seq_len += line_length;
        }
      }

      if (state == 3) {
        /* we are in state == 3, so this line must be structure */
        size_t tmp_len = struct_len + line_length + 1;
        *structure_p = (char *)vrna_realloc(*structure_p, sizeof(char) * tmp_len);
        memcpy(*structure_p + struct_len, line, sizeof(char) * line_length);
        (*structure_p)[tmp_len - 1] = '\0';
        struct_len += line_length;
      }
    }
  }

RNAstrand_parser_end:

  if (*name_p)
    ret++;
  if (*source_p)
    ret++;
  if (*sequence_p)
    ret++;
  if (*structure_p)
    ret++;
  if (*fname_p)
    ret++;
  if (*id_p)
    ret++;

  if ((*sequence_p) && (*structure_p))
    return ret;
  else {
    return 0;
  }
}



PUBLIC int
vrna_file_connect_read_record(FILE          *fp,
                              char          **id,
                              char          **sequence,
                              char          **structure,
                              char          **remainder,
                              unsigned int  options)
{
  char          c, *line, *end, tok2, **tok, **ptr;
  size_t        num_tok;
  int           is_nt_line, is_header;
  long          tok1, tok3, tok4, tok5, tok6;
  ct_data       *ct_entry = NULL;

  if (!fp) {
    if (options & VRNA_INPUT_VERBOSE)
      vrna_log_warning("vrna_file_connect_read_record@file_formats.c: "
                           "Can't read from file pointer while parsing connectivity table formatted sequence input!");
    return -1;
  }

  if (id)
    *id = NULL;

  if (sequence)
    *sequence = NULL;

  if (structure)
    *structure = NULL;

  if ((remainder) &&
      (*remainder)) {
    line        = *remainder;
    *remainder  = NULL;
  } else {
    line = vrna_read_line(fp);
  }

  if (!line)
    return 0;

  do {
    /* trim leading and trailing white spaces */
    vrna_strtrim(line, NULL, 0, VRNA_TRIM_DEFAULT);

    /* reduce all consecutive whitespaces to a single space */
    vrna_strtrim(line, NULL, 1, VRNA_TRIM_IN_BETWEEN | VRNA_TRIM_SUBST_BY_FIRST);

    /*
        For now, we do not extract any information from the comments
        in the file, so lets skip empty and comment lines all together
    */
    c = *line;

    if ((c == '\0') ||
        (c == '*') ||
        (c == '>') ||
        (c == '#') ||
        (c == ';')) {
      free(line);
      continue;
    }

    /* if this is a non-comment and non-empty line, tokenize it */
    tok = vrna_strsplit(line, " ");

    /* count number of tokens */
    for (num_tok = 0; tok[num_tok]; num_tok++);


    /*
        First, use some heuristic to judge what kind of
        line we are currently faced with. We will distinguish
        two types of data: headers and nucleotide information
    */
    is_nt_line = is_header = 0;

    /*
        1st, if we have more than 5 whitespace separated blocks, this
        seems like a nucleotide information line. Check, if this guess
        is valid
    */
    if (num_tok > 5) {
      tok2 = tok[1][0];

      /* check 1st token, i.e. nt position within segment */
      tok1 = strtol(tok[0], &end, 10);
      if (tok[0] != end) {
        /* check 3rd token, i.e. 5' connecting base */
        tok3 = strtol(tok[2], &end, 10);
        if (tok[2] != end) {
          /* check 4th token, i.e. 3' connecting base */
          tok4 = strtol(tok[3], &end, 10);
          if (tok[3] != end) {
            /* check 5th token, i.e. pairing partner */
            tok5 = strtol(tok[4], &end, 10);
            if (tok[4] != end) {
              /* check 6th token, i.e. historical numbering */
              tok6 = strtol(tok[5], &end, 10);
              if (tok[6] != end)
                is_nt_line = 1;
            }
          }
        }
      }
    }

    /* 2nd, we might have read a header line */
    if ((!is_nt_line) &&
        (num_tok > 0)) {
      tok1 = strtol(tok[0], &end, 10);
      if (tok[0] != end) {
        /* first token is a number, should be sequence length */
        is_header = 1;
      }
    }

    /* now, do something with our guess and extract the data */
    if (is_header) {
      /* return result and save header for next round if we already have a complete data set */
      if (ct_entry) {
        if ((finalize_ct_data(ct_entry) != 0) &&
            (options & VRNA_INPUT_VERBOSE))
          vrna_log_warning("vrna_file_connect_read_record@file_formats.c: "
                               "Malformed input file! Sequence length stated: %u, actual length: %u\n",
                               ct_entry->length,
                               ct_entry->effective_length);

        *id         = ct_entry->id;
        *sequence   = ct_entry->sequence;
        *structure  = vrna_db_from_ptable(ct_entry->pt);
        *remainder  = line;

        free(ct_entry->pt);
        free(ct_entry->actual_pos);
        free(ct_entry);

        for (ptr = tok; *ptr; ptr++)
          free(*ptr);
        free(tok);

        return 1;
      }

      ct_entry = process_ct_header((unsigned int)tok1, num_tok, (const char **)tok);
    } else if ((is_nt_line) &&
               (ct_entry)) {
      /* this is a nucleotide information line, so let's extract the relevant data */
      if (!process_ct_nt_line(ct_entry, tok1, tok2, tok3, tok5, tok6))
        printf("Something went wrong with storing nucleotide information\n");
    } else if (options & VRNA_INPUT_VERBOSE) {
      vrna_log_warning("vrna_file_connect_read_record@file_formats.c: "
                           "Unusal line in input:\n%s\n", line);
    }

    free(line);
    for (ptr = tok; *ptr; ptr++)
      free(*ptr);
    free(tok);
  } while ((line = vrna_read_line(fp)));

  /* end of file */
  if (ct_entry) {
    if ((finalize_ct_data(ct_entry) != 0) &&
        (options & VRNA_INPUT_VERBOSE))
      vrna_log_warning("vrna_file_connect_read_record@file_formats.c: "
                           "Malformed input file! Sequence length stated: %u, actual length: %u\n",
                           ct_entry->length,
                           ct_entry->effective_length);

    *id         = ct_entry->id;
    *sequence   = ct_entry->sequence;
    *structure  = vrna_db_from_ptable(ct_entry->pt);
    *remainder  = NULL;

    free(ct_entry->pt);
    free(ct_entry->actual_pos);
    free(ct_entry);

    return 1;
  }

  return 0;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */

PRIVATE INLINE ct_data *
init_ct_data(unsigned int n)
{
  ct_data *dat      = (ct_data *)vrna_alloc(sizeof(ct_data));

  dat->length           = n;
  dat->effective_length = n;
  dat->strands          = 1;
  dat->sequence         = (char *)vrna_alloc(sizeof(char) * (2 * n + 1));
  dat->sequence_pos     = 0;
  dat->pt               = (short *)vrna_alloc(sizeof(short) * (n + 1));
  dat->actual_pos       = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n + 1));
  dat->pt[0]            = (short )n;
  dat->id               = NULL;

  return dat;
}


PRIVATE INLINE void
resize_ct_data(ct_data      *dat,
               unsigned int n)
{
  dat->effective_length = n;
  dat->sequence         = (char *)vrna_realloc(dat->sequence, sizeof(char) * (2 * n + 1));
  dat->pt               = (short *)vrna_realloc(dat->pt, sizeof(short) * (n + 1));
  dat->actual_pos       = (unsigned int *)vrna_realloc(dat->actual_pos, sizeof(unsigned int) * (n + 1));
}


PRIVATE INLINE int
finalize_ct_data(ct_data *dat)
{
  dat->sequence[dat->sequence_pos] = '\0';
  if (strlen(dat->sequence) < dat->effective_length) {
    memset(dat->sequence, 'N', sizeof(char) * (dat->effective_length - strlen(dat->sequence)));
    dat->sequence[dat->effective_length] = '\0';
  }
  dat->pt[0] = (short)dat->effective_length;

  if (dat->length != dat->effective_length)
    return 1;

  return 0;
}


PRIVATE INLINE ct_data*
process_ct_header(unsigned int  n,
                  size_t        num_tok,
                  const char    **tok)
{
  char    *ptr;
  size_t  id_pos;
  float   energy;
  ct_data *dat = init_ct_data(n);

  if (num_tok > 1) {
    id_pos = 1;

    /* find out where the sequence id actually starts */
    /* sometimes it is preceeded by an ENERGY = x.y entry
     * denoting the free energy of the structure encoded in
     * the connectivity data.
     */
    ptr = strdup(tok[1]);
    vrna_seq_toupper(ptr);

    if (sscanf(ptr, "ENERGY = %f", &energy) == 1) {
      id_pos++;
    } else if (num_tok > 2) {
      free(ptr);
      ptr = vrna_strdup_printf("%s %s",
                               tok[1],
                               tok[2]);
      vrna_seq_toupper(ptr);

      if (sscanf(ptr, "ENERGY = %f", &energy) == 1) {
        id_pos += 2;
      } else if (num_tok > 3) {
        free(ptr);
        ptr = vrna_strdup_printf("%s %s %s",
                                 tok[1],
                                 tok[2],
                                 tok[3]);
        vrna_seq_toupper(ptr);

        if (sscanf(ptr, "ENERGY = %f", &energy) == 1) {
          id_pos += 3;
        }

        free(ptr);
      }
    }

    /*
        we've also seen .ct files that only contain the
        ENERGY keyword, without any further values, so
        let's also skip this keyword if necessary
        usage of strstr() is save here, since we converted
        ptr to uppercase before
    */
    if ((id_pos == 1) &&
        (strstr(ptr, "ENERGY") == ptr)) {
      id_pos++;
    }

    /*
        finally, if there is anything left that might be a
        sequence identifier, keep it
    */
    if (id_pos < num_tok)
      dat->id = vrna_strjoin((const char **)(tok + id_pos), " ");
  }

  return dat;
}


PRIVATE INLINE int
process_ct_nt_line(ct_data       *data,
                   unsigned int  i,
                   char          nucleotide,
                   unsigned int  predecessor,
                   unsigned int  j,
                   unsigned int  actual_i)
{
  /* check for malformed input line */
  unsigned int  max = i;
  max = MAX2(max, j);

  if (max > data->effective_length)
    resize_ct_data(data, max);

  if (i <= data->effective_length) {
    if ((i > 1) &&
        (predecessor == 0)) {
      data->strands++;
      data->sequence[data->sequence_pos++] = '&';
    }

    data->pt[i]                           = j;
    data->sequence[data->sequence_pos++]  = nucleotide;
    data->actual_pos[i]                   = actual_i;

    return 1;
  }

  return 0;
}


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
 * ###########################################
 * # deprecated functions below              #
 *###########################################
 */

PUBLIC unsigned int
get_multi_input_line(char         **string,
                     unsigned int option)
{
  return read_multiple_input_lines(string, NULL, option);
}


PUBLIC unsigned int
read_record(char          **header,
            char          **sequence,
            char          ***rest,
            unsigned int  options)
{
  return vrna_file_fasta_read_record(header, sequence, rest, NULL, options);
}


PUBLIC char *
extract_record_rest_structure(const char    **lines,
                              unsigned int  length,
                              unsigned int  options)
{
  return vrna_extract_record_rest_structure(lines, length, options);
}


#endif
