/*
    file_formats.c

    Various functions dealing with file formats for RNA sequences, structures, and alignments

    (c) 2014 Ronny Lorenz

    Vienna RNA package
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/file_formats.h"

/*
#################################
# PRIVATE VARIABLES             #
#################################
*/

PRIVATE char          *inbuf  = NULL;
PRIVATE char          *inbuf2 = NULL;
PRIVATE unsigned int  typebuf = 0;

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/

PRIVATE void
find_helices(short *pt, int i, int j, FILE *file);

PRIVATE unsigned int
read_multiple_input_lines(char **string, FILE *file, unsigned int option);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PRIVATE void
find_helices(short *pt, int i, int j, FILE *file){

  FILE *out = (file) ? file : stdout;
  int h_start, h_length, h_end;

  h_start = h_length = h_end = 0;

  for(;i < j; i++){
    if(i > pt[i]) continue;
    h_start = i;
    h_end   = pt[i];
    h_length = 1;
    while(pt[i+1] == (pt[i]-1)){
      h_length++;
      i++;
    }
    if(i < h_end){
      find_helices(pt, i+1, h_end, file);
    }
    if(h_length > 1){
      fprintf(out, "%d %d %d\n", h_start, h_end, h_length);
    }
    i = pt[h_start] - 1;
  }
}

PUBLIC void
vrna_structure_print_helix_list(const char *db,
                                FILE *file){

  short *pt = vrna_pt_get(db);

  find_helices(pt, 1, pt[0], file);
  free(pt);
}

PUBLIC void
vrna_structure_print_ct(const char *seq,
                        const char *db,
                        float energy,
                        const char *identifier,
                        FILE *file){

  int i, power_d;
  FILE *out = (file) ? file : stdout;

  if(strlen(seq) != strlen(db))
    nrerror("vrna_structure_print_ct: sequence and structure have unequal length!");

  short *pt = vrna_pt_get(db);

  for(power_d=0;pow(10,power_d) <= (int)strlen(seq);power_d++);

  /*
    Connect table file format looks like this:

    300  ENERGY = 7.0  example
      1 G       0    2   22    1
      2 G       1    3   21    2

    where the headerline is followed by 6 columns with:
    1. Base number: index n
    2. Base (A, C, G, T, U, X)
    3. Index n-1  (0 if first nucleotide)
    4. Index n+1  (0 if last nucleotide)
    5. Number of the base to which n is paired. No pairing is indicated by 0 (zero).
    6. Natural numbering.
  */

  /* print header */
  fprintf(out, "%d  ENERGY = %6.2f", (int)strlen(seq), energy);
  if(identifier)
    fprintf(out, "  %s\n", identifier);

  /* print structure information except for last line */
  /* TODO: modify the structure information for cofold */
  for(i = 0; i < strlen(seq) - 1; i++){
    fprintf(out, "%*d %c %*d %*d %*d %*d\n",
                  power_d, i+1,           /* nucleotide index */
                  (char)toupper(seq[i]),  /* nucleotide char */
                  power_d, i,             /* nucleotide predecessor index */
                  power_d, i+2,           /* nucleotide successor index */
                  power_d, pt[i+1],       /* pairing partner index */
                  power_d, i+1);          /* nucleotide natural numbering */
  }
  /* print last line */
  fprintf(out, "%*d %c %*d %*d %*d %*d\n",
                power_d, i+1,
                (char)toupper(seq[i]),
                power_d, i,
                power_d, 0,
                power_d, pt[i+1],
                power_d, i+1);

  /* clean up */
  free(pt);
  fflush(out);
}

PUBLIC void
vrna_structure_print_bpseq( const char *seq,
                            const char *db,
                            FILE *file){

  int i;
  FILE *out = (file) ? file : stdout;

  if(strlen(seq) != strlen(db))
    nrerror("vrna_structure_print_bpseq: sequence and structure have unequal length!");

  short *pt = vrna_pt_get(db);

  for(i = 1; i <= pt[0]; i++){
    fprintf(out, "%d %c %d\n", i, (char)toupper(seq[i-1]), pt[i]);
  }

  /* clean up */
  free(pt);
  fflush(out);
}

PRIVATE  unsigned int
read_multiple_input_lines(char **string,
                          FILE *file,
                          unsigned int option){

  char  *line;
  int   i, l;
  int   state = 0;
  int   str_length = 0;
  FILE  *in = (file) ? file : stdin;

  line = (inbuf2) ? inbuf2 : get_line(in);
  inbuf2 = NULL;
  do{

    /*
    * read lines until informative data appears or
    * report an error if anything goes wrong
    */
    if(!line) return VRNA_INPUT_ERROR;

    l = (int)strlen(line);

    /* eliminate whitespaces at the end of the line read */
    if(!(option & VRNA_INPUT_NO_TRUNCATION)){
      for(i = l-1; i >= 0; i--){
        if      (line[i] == ' ')  continue;
        else if (line[i] == '\t') continue;
        else                      break;
      }
      line[(i >= 0) ? (i+1) : 0] = '\0';
    }

    l           = (int)strlen(line);
    str_length  = (*string) ? (int) strlen(*string) : 0;

    switch(*line){
      case  '@':    /* user abort */
                    if(state) inbuf2 = line;
                    else      free(line);
                    return (state==2) ? VRNA_INPUT_CONSTRAINT : (state==1) ? VRNA_INPUT_SEQUENCE : VRNA_INPUT_QUIT;

      case  '\0':   /* empty line */
                    if(option & VRNA_INPUT_NOSKIP_BLANK_LINES){
                      if(state) inbuf2 = line;
                      else      free(line);
                      return (state==2) ? VRNA_INPUT_CONSTRAINT : (state==1) ? VRNA_INPUT_SEQUENCE : VRNA_INPUT_BLANK_LINE;
                    }
                    break;

      case  '#': case  '%': case  ';': case  '/': case  '*': case ' ':
                    /* comments */
                    if(option & VRNA_INPUT_NOSKIP_COMMENTS){
                      if(state) inbuf2   = line;
                      else      *string = line;
                      return (state == 2) ? VRNA_INPUT_CONSTRAINT : (state==1) ? VRNA_INPUT_SEQUENCE : VRNA_INPUT_COMMENT;
                    }
                    break;

      case  '>':    /* fasta header */
                    if(state) inbuf2   = line;
                    else      *string = line;
                    return (state==2) ? VRNA_INPUT_CONSTRAINT : (state==1) ? VRNA_INPUT_SEQUENCE : VRNA_INPUT_FASTA_HEADER;

      case  'x': case 'e': case 'l': case '&':   /* seems to be a constraint or line starting with second sequence for dimer calculations */
                    i = 1;
                    /* lets see if this assumption holds for the complete line */
                    while((line[i] == 'x') || (line[i] == 'e') || (line[i] == 'l')) i++;
                    /* lines solely consisting of 'x's, 'e's or 'l's will be considered as structure constraint */
                    
                    if(
                            ((line[i]>64) && (line[i]<91))  /* A-Z */
                        ||  ((line[i]>96) && (line[i]<123)) /* a-z */
                      ){
                      if(option & VRNA_INPUT_FASTA_HEADER){
                        /* are we in structure mode? Then we remember this line for the next round */
                        if(state == 2){ inbuf2 = line; return VRNA_INPUT_CONSTRAINT;}
                        else{
                          *string = (char *)xrealloc(*string, sizeof(char) * (str_length + l + 1));
                          strcpy(*string + str_length, line);
                          state = 1;
                        }
                        break;
                      }
                      /* otherwise return line read */
                      else{ *string = line; return VRNA_INPUT_SEQUENCE;}
                    }
                    /* mmmh? it really seems to be a constraint */
                    /* fallthrough */
      case  '<': case  '.': case  '|': case  '(': case ')': case '[': case ']': case '{': case '}': case ',': case '+':
                    /* seems to be a structure or a constraint */
                    /* either we concatenate this line to one that we read previously */
                    if(option & VRNA_INPUT_FASTA_HEADER){
                      if(state == 1){
                        inbuf2 = line;
                        return VRNA_INPUT_SEQUENCE;
                      }
                      else{
                        *string = (char *)xrealloc(*string, sizeof(char) * (str_length + l + 1));
                        strcpy(*string + str_length, line);
                        state = 2;
                      }
                    }
                    /* or we return it as it is */
                    else{
                      *string = line;
                      return VRNA_INPUT_CONSTRAINT;
                    }
                    break;
      default:      if(option & VRNA_INPUT_FASTA_HEADER){
                      /* are we already in sequence mode? */
                      if(state == 2){
                        inbuf2 = line;
                        return VRNA_INPUT_CONSTRAINT;
                      }
                      else{
                        *string = (char *)xrealloc(*string, sizeof(char) * (str_length + l + 1));
                        strcpy(*string + str_length, line);
                        state = 1;
                      }
                    }
                    /* otherwise return line read */
                    else{
                      *string = line;
                      return VRNA_INPUT_SEQUENCE;
                    }
    }
    free(line);
    line = get_line(in);
  }while(line);

  return (state==2) ? VRNA_INPUT_CONSTRAINT : (state==1) ? VRNA_INPUT_SEQUENCE : VRNA_INPUT_ERROR;
}

PUBLIC  unsigned int
vrna_read_fasta_record( char **header,
                        char **sequence,
                        char ***rest,
                        FILE *file,
                        unsigned int options){

  unsigned int  input_type, return_type, tmp_type;
  int           rest_count;
  char          *input_string;

  rest_count    = 0;
  return_type   = tmp_type = 0;
  input_string  = *header = *sequence = NULL;
  *rest         = (char **)space(sizeof(char *));

  /* remove unnecessary option flags from options variable... */
  options &= ~VRNA_INPUT_FASTA_HEADER;

  /* read first input or last buffered input */
  if(typebuf){
    input_type    = typebuf;
    input_string  = inbuf;
    typebuf       = 0;
    inbuf         = NULL;
  }
  else input_type  = read_multiple_input_lines(&input_string, file, options);

  if(input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR)) return input_type;

  /* skip everything until we read either a fasta header or a sequence */
  while(input_type & (VRNA_INPUT_MISC | VRNA_INPUT_CONSTRAINT | VRNA_INPUT_BLANK_LINE)){
    free(input_string); input_string = NULL;
    input_type    = read_multiple_input_lines(&input_string, file, options);
    if(input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR)) return input_type;
  }

  if(input_type & VRNA_INPUT_FASTA_HEADER){
    return_type  |= VRNA_INPUT_FASTA_HEADER; /* remember that we've read a fasta header */
    *header       = input_string;
    input_string  = NULL;
    /* get next data-block with fasta support if not explicitely forbidden by VRNA_INPUT_NO_SPAN */
    input_type  = read_multiple_input_lines(
                    &input_string,
                    file,
                    ((options & VRNA_INPUT_NO_SPAN) ? 0 : VRNA_INPUT_FASTA_HEADER) | options
                  );
    if(input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR)) return (return_type | input_type);
  }

  if(input_type & VRNA_INPUT_SEQUENCE){
    return_type  |= VRNA_INPUT_SEQUENCE; /* remember that we've read a sequence */
    *sequence     = input_string;
    input_string  = NULL;
  } else nrerror("sequence input missing");

  /* read the rest until we find user abort, EOF, new sequence or new fasta header */
  if(!(options & VRNA_INPUT_NO_REST)){
    options |= VRNA_INPUT_NOSKIP_COMMENTS; /* allow commetns to appear in rest output */
    tmp_type = VRNA_INPUT_QUIT | VRNA_INPUT_ERROR | VRNA_INPUT_SEQUENCE | VRNA_INPUT_FASTA_HEADER;
    if(options & VRNA_INPUT_NOSKIP_BLANK_LINES) tmp_type |= VRNA_INPUT_BLANK_LINE;
    while(!((input_type = read_multiple_input_lines(&input_string, file, options)) & tmp_type)){
      *rest = xrealloc(*rest, sizeof(char **)*(++rest_count + 1));
      (*rest)[rest_count-1] = input_string;
      input_string = NULL;
    }
    /*
    if(input_type & (VRNA_INPUT_QUIT | VRNA_INPUT_ERROR)) return input_type;
    */

    /*  finished reading everything...
    *   we now put the last line into the buffer if necessary
    *   since it should belong to the next record
    */
    inbuf = input_string;
    typebuf = input_type;
  }
  (*rest)[rest_count] = NULL;
  return (return_type);
}

PUBLIC char *
vrna_extract_record_rest_structure( const char **lines,
                                    unsigned int length,
                                    unsigned int options){

  char *structure = NULL;
  int r, i, l, cl, stop;
  char *c;

  if(lines){
    for(r = i = stop = 0; lines[i]; i++){
      l   = (int)strlen(lines[i]);
      c   = (char *) space(sizeof(char) * (l+1));
      (void) sscanf(lines[i], "%s", c);
      cl  = (int)strlen(c);

      /* line commented out ? */
      if((*c == '#') || (*c == '%') || (*c == ';') || (*c == '/') || (*c == '*' || (*c == '\0'))){
        /* skip leading comments only, i.e. do not allow comments inside the constraint */
        if(!r)  continue;
        else    break;
      }

      /* append the structure part to the output */
      r += cl+1;
      structure = (char *)xrealloc(structure, r*sizeof(char));
      strcat(structure, c);
      free(c);
      /* stop if the assumed structure length has been reached */
      if((length > 0) && (r-1 == length)) break;
      /* stop if not allowed to read from multiple lines */
      if(!(options & VRNA_OPTION_MULTILINE)) break;
    }
  }
  return structure;
}

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

PUBLIC unsigned int
get_multi_input_line( char **string,
                      unsigned int option){

  return read_multiple_input_lines(string, NULL, option);
}

PUBLIC unsigned int
read_record(char **header,
            char **sequence,
            char  ***rest,
            unsigned int options){

  return vrna_read_fasta_record(header, sequence, rest, NULL, options);
}

PUBLIC char *
extract_record_rest_structure(const char **lines,
                              unsigned int length,
                              unsigned int options){

  return vrna_extract_record_rest_structure(lines, length, options);
}
