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
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/

PRIVATE void
find_helices(short *pt, int i, int j, FILE *file);

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
    nrerror("vrna_ct_from_dbsequence and ");

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
    nrerror("vrna_ct_from_dbsequence and ");

  short *pt = vrna_pt_get(db);

  for(i = 1; i <= pt[0]; i++){
    fprintf(out, "%d %c %d\n", i, (char)toupper(seq[i-1]), pt[i]);
  }

  /* clean up */
  free(pt);
  fflush(out);
}

