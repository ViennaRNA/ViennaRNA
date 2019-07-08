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

#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/sequence.h"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE void set_sequence(vrna_seq_t    *obj,
                          const char    *string,
                          const char    *name,
                          vrna_md_t     *md,
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
    set_sequence(data, string, NULL, NULL, options);
  }

  return data;
}


PUBLIC int
vrna_sequence_add(vrna_fold_compound_t  *vc,
                  const char            *string,
                  unsigned int          options)
{
  unsigned int add_length;
  int ret = 0;

  if ((vc) && (vc->type == VRNA_FC_TYPE_SINGLE) && (string)) {
    add_length = strlen(string);
    
    /* add the sequence to the nucleotides container */
    vc->nucleotides = (vrna_seq_t *)vrna_realloc(vc->nucleotides,
                                                 sizeof(vrna_seq_t) *
                                                 (vc->strands + 1));
    set_sequence(&(vc->nucleotides[vc->strands]),
                 string,
                 NULL,
                 &(vc->params->model_details),
                 options);

    /* increase strands counter */
    vc->strands++;

    /* add new sequence to initial order of all strands */
    vc->sequence = (char *)vrna_realloc(vc->sequence,
                                        sizeof(char) *
                                        (vc->length + add_length + 1));
    memcpy(vc->sequence + vc->length,
           vc->nucleotides[vc->strands - 1].string,
           add_length * sizeof(char));
    vc->sequence[vc->length + add_length] = '\0';

    /* add encoding for new strand */
    vc->sequence_encoding = (short *)vrna_realloc(vc->sequence_encoding,
                                                  sizeof(short) *
                                                  (vc->length + add_length + 2));

    memcpy(vc->sequence_encoding + vc->length + 1,
           vc->nucleotides[vc->strands - 1].encoding + 1,
           add_length * sizeof(short));

    /* restore circular encoding */
    vc->sequence_encoding[vc->length + add_length + 1]  = vc->sequence_encoding[1];
    vc->sequence_encoding[0]                            = vc->sequence_encoding[vc->length + add_length];

    /* add encoding2 (simple encoding) for new strand */
    vc->sequence_encoding2 = (short *)vrna_realloc(vc->sequence_encoding2,
                                                   sizeof(short) *
                                                   (vc->length + add_length + 2));
    short *enc = vrna_seq_encode_simple(vc->nucleotides[vc->strands - 1].string,
                                        &(vc->params->model_details));
    memcpy(vc->sequence_encoding2 + vc->length + 1,
           enc + 1,
           add_length * sizeof(short));
    free(enc);
    vc->sequence_encoding2[vc->length + add_length + 1] = vc->sequence_encoding2[1];
    vc->sequence_encoding2[0] = (short)(vc->length + add_length);

    /* finally, increase length property of the fold compound */
    vc->length = vc->length + add_length;

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
vrna_sequence_remove_all(vrna_fold_compound_t *fc)
{
  unsigned int i;

  if (fc) {
    for (i = 0; i < fc->strands; i++)
      free_sequence_data(&(fc->nucleotides[i]));

    free(fc->nucleotides);
    free(fc->strand_number);
    free(fc->strand_order);
    free(fc->strand_start);
    free(fc->strand_end);

    fc->strands       = 0;
    fc->nucleotides   = NULL;
    fc->strand_number = NULL;
    fc->strand_order  = NULL;
    fc->strand_start  = NULL;
    fc->strand_end    = NULL;
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
             const char   *name,
             vrna_md_t    *md,
             unsigned int options)
{
  obj->name   = name ? strdup(name) : NULL;
  obj->string = strdup(string);
  vrna_seq_toupper(obj->string);
  obj->length = strlen(obj->string);

  switch (options) {
    default:
      obj->type = VRNA_SEQ_RNA;
  }

  obj->encoding   = vrna_seq_encode(obj->string, md);
  obj->encoding5  = (short *)vrna_alloc(sizeof(short) * (obj->length + 1));
  obj->encoding3  = (short *)vrna_alloc(sizeof(short) * (obj->length + 1));

  if (md->circ) {
    for (size_t i = obj->length; i > 0; i--) {
      if (obj->encoding[i] == 0) /* no nucleotide, i.e. gap */
        continue;

      obj->encoding5[1] = obj->encoding[i];
      break;
    }
    for (size_t i = 1; i <= obj->length; i++) {
      if (obj->encoding[i] == 0) /* no nucleotide, i.e. gap */
        continue;

      obj->encoding3[obj->length] = obj->encoding[i];
      break;
    }
  } else {
    obj->encoding5[1] = obj->encoding3[obj->length] = 0;
  }

  for (size_t i = 1, p = 0; i < obj->length; i++) {
    if (obj->encoding[i] == 0) {
      obj->encoding5[i + 1] = obj->encoding5[i];
    } else {
      obj->encoding5[i + 1]  = obj->encoding[i];
    }
  }

  for (size_t i = obj->length; i > 1; i--) {
    if (obj->encoding[i] == 0)
      obj->encoding3[i - 1] = obj->encoding3[i];
    else
      obj->encoding3[i - 1] = obj->encoding[i];
  }

}


PRIVATE void
free_sequence_data(vrna_seq_t *obj)
{
  free(obj->string);
  free(obj->name);
  free(obj->encoding);
  free(obj->encoding5);
  free(obj->encoding3);
  obj->string = NULL;
  obj->name = NULL;
  obj->encoding = NULL;
  obj->encoding5 = NULL;
  obj->encoding3 = NULL;
  obj->type   = VRNA_SEQ_UNKNOWN;
  obj->length = 0;
}
