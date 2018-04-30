/*
 *  sequence.c
 *
 *  Code for handling nucleotide sequences
 *
 *  Part of the ViennaRNA Package
 *  (c) 2016 Ronny Lorenz
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ViennaRNA/data_structures.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/string_utils.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/sequence.h"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE void set_sequence(vrna_seq_t    *obj,
                          const char    *string,
                          unsigned int  options);


PRIVATE void free_sequence_data(vrna_seq_t *obj);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_seq_t *
vrna_sequence(const char    *string,
              unsigned int  options)
{
  vrna_seq_t *data = NULL;

  if (string) {
    data = (vrna_seq_t *)vrna_alloc(sizeof(vrna_seq_t));
    set_sequence(data, string, options);
  }

  return data;
}


PUBLIC int
vrna_sequence_add(vrna_fold_compound_t  *vc,
                  const char            *string,
                  unsigned int          options)
{
  int ret = 0;

  if ((vc) && (vc->type == VRNA_FC_TYPE_SINGLE) && (string)) {
    vc->nucleotides =
      (vrna_seq_t *)vrna_realloc(vc->nucleotides, sizeof(vrna_seq_t) * (vc->strands + 1));
    set_sequence(&(vc->nucleotides[vc->strands]), string, options);
    vc->strands++;

    ret = 1;
  }

  return ret;
}


PUBLIC int
vrna_sequence_remove(vrna_fold_compound_t *vc,
                     unsigned int         i)
{
  unsigned int  size;
  int           ret = 0;

  if (vc) {
    if (i < vc->strands) {
      free_sequence_data(&(vc->nucleotides[i]));
      /* roll all nucleotide sequences behind the deleted one to the front */
      size = vc->strands - i - 1;
      if (size > 0)
        memmove(vc->nucleotides + i, vc->nucleotides + i + 1, sizeof(vrna_seq_t) * size);

      vc->strands--;
      vc->nucleotides =
        (vrna_seq_t *)vrna_realloc(vc->nucleotides, sizeof(vrna_seq_t) * vc->strands);

      ret = 1;
    }
  }

  return ret;
}


PUBLIC void
vrna_sequence_remove_all(vrna_fold_compound_t *vc)
{
  unsigned int i;

  if (vc) {
    for (i = 0; i < vc->strands; i++)
      free_sequence_data(&(vc->nucleotides[i]));

    free(vc->nucleotides);

    vc->nucleotides = NULL;
    vc->strands     = 0;
  }
}


PUBLIC void
vrna_sequence_prepare(vrna_fold_compound_t *fc)
{
  unsigned int cnt, i;

  if (fc) {
    free(fc->strand_number);
    free(fc->strand_order);
    free(fc->strand_start);
    free(fc->strand_end);

    fc->strand_order  = NULL;
    fc->strand_start  = NULL;
    fc->strand_end    = NULL;

    fc->strand_number = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (fc->length + 2));

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        /* 1. store initial strand order */
        fc->strand_order = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (fc->strands + 1));
        for (cnt = 0; cnt < fc->strands; cnt++)
          fc->strand_order[cnt] = cnt;

        /* 2. mark start and end positions of sequences */
        fc->strand_start  = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (fc->strands + 1));
        fc->strand_end    = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (fc->strands + 1));

        fc->strand_start[0] = 1;
        fc->strand_end[0]   = fc->strand_start[0] + fc->nucleotides[0].length - 1;

        for (cnt = 1; cnt < fc->strands; cnt++) {
          fc->strand_start[cnt] = fc->strand_end[cnt - 1] + 1;
          fc->strand_end[cnt]   = fc->strand_start[cnt] + fc->nucleotides[cnt].length - 1;
          for (i = fc->strand_start[cnt]; i <= fc->strand_end[cnt]; i++)
            fc->strand_number[i] = cnt;
        }
        /* this sets pos. n + 1 as well */
        fc->strand_number[fc->length + 1] = fc->strands - 1;

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        /*
         *  for now, comparative structure prediction does not allow for RNA-RNA interactions,
         *  so we pretend working on a single strand
         */
        fc->strands     = 1;
        fc->nucleotides = (vrna_seq_t *)vrna_realloc(fc->nucleotides,
                                                     sizeof(vrna_seq_t) * (fc->strands + 1));
        fc->nucleotides[0].string = NULL;
        fc->nucleotides[0].type   = VRNA_SEQ_RNA;
        fc->nucleotides[0].length = fc->length;

        /* 1. store initial strand order */
        fc->strand_order = (unsigned int *)vrna_alloc(sizeof(unsigned int) * 2);

        /* 2. mark start and end positions of sequences */
        fc->strand_start    = (unsigned int *)vrna_alloc(sizeof(unsigned int) * 2);
        fc->strand_end      = (unsigned int *)vrna_alloc(sizeof(unsigned int) * 2);
        fc->strand_start[0] = 1;
        fc->strand_end[0]   = fc->strand_start[0] + fc->length - 1;

        break;
    }
  }
}


PRIVATE void
set_sequence(vrna_seq_t   *obj,
             const char   *string,
             unsigned int options)
{
  obj->string = strdup(string);
  vrna_seq_toupper(obj->string);
  obj->length = strlen(obj->string);

  switch (options) {
    default:
      obj->type = VRNA_SEQ_RNA;
  }

  obj->encoding = vrna_seq_encode(obj->string, NULL);
}


PRIVATE void
free_sequence_data(vrna_seq_t *obj)
{
  free(obj->string);
  free(obj->encoding);
  obj->type   = VRNA_SEQ_UNKNOWN;
  obj->length = 0;
}
