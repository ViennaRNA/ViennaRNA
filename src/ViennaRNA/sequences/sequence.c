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
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/sequences/sequence.h"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE void
set_sequence(vrna_seq_t   *obj,
             const char   *string,
             const char   *name,
             vrna_md_t    *md,
             unsigned int options);


PRIVATE void
free_sequence_data(vrna_seq_t *obj);


PRIVATE void
update_strand_positions(vrna_fold_compound_t *fc);


PRIVATE void
update_sequence(vrna_fold_compound_t *fc);


PRIVATE void
update_encodings(vrna_fold_compound_t *fc);


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
vrna_sequences_add(vrna_fold_compound_t *fc,
                   const char           **sequences,
                   const unsigned int   *order,
                   unsigned int         options)
{
  if ((fc) &&
      (sequences)) {
    size_t        s, s_new;
    unsigned int  strands_old, n_new;
    vrna_md_t     *md;

    md            = &(fc->params->model_details);
    strands_old   = fc->strands;

    /* count the number of sequences to add */
    for (s_new = 0; sequences[s_new]; s_new++);

    /* append sequences to storage container */
    fc->nucleotides = (vrna_seq_t *)vrna_realloc(fc->nucleotides,
                                                 sizeof(vrna_seq_t) *
                                                 (strands_old + s_new));
    for (s = 0, n_new = 0; s < s_new; s++) {
      set_sequence(&(fc->nucleotides[strands_old + s]),
                   sequences[s],
                   NULL,
                   md,
                   options);
      n_new += fc->nucleotides[strands_old + s].length;
    }

    /* adjust strands counter */
    fc->strands += s_new;

    /* adjust total length of concatenated sequences */
    fc->length += n_new;

    /* adjust sequence order */
    fc->strand_order = (unsigned int *)vrna_realloc(fc->strand_order,
                                                    sizeof(unsigned int) *
                                                    (fc->strands + 1));
    if (order) {
      memcpy(fc->strand_order + strands_old + 1,
             order,
             sizeof(unsigned int) * s_new);
    } else {
      /* default order is same as input order */
      for (s = 0; s < s_new; s++)
        fc->strand_order[strands_old + s + 1] = s;
    }

    for (s = 0; s < s_new; s++)
      fc->strand_order[strands_old + s + 1] += strands_old;

    /* adjust strand positions */
    fc->strand_start = (unsigned int *)vrna_realloc(fc->strand_start,
                                                    sizeof(unsigned int) *
                                                    (fc->strands + 1));
    fc->strand_end   = (unsigned int *)vrna_realloc(fc->strand_end,
                                                    sizeof(unsigned int) *
                                                    (fc->strands + 1));
    fc->strand_number = (unsigned int *)vrna_realloc(fc->strand_number,
                                                     sizeof(unsigned int) *
                                                     (fc->length + 2));

    update_strand_positions(fc);

    /* adjust sequence array */
    fc->sequence  = (char *)vrna_realloc(fc->sequence,
                                         sizeof(char) *
                                         (fc->length + 1));

    update_sequence(fc);
    fc->sequence[fc->length] = '\0';

    /* adjust encodings arrays */
    fc->sequence_encoding = (short *)vrna_realloc(fc->sequence_encoding,
                                                  sizeof(short) *
                                                  (fc->length + 2));
    fc->sequence_encoding2 = (short *)vrna_realloc(fc->sequence_encoding2,
                                                  sizeof(short) *
                                                  (fc->length + 2));
    fc->encoding5 = (short *)vrna_realloc(fc->encoding5,
                                          sizeof(short) *
                                          (fc->length + 2));
    fc->encoding3 = (short *)vrna_realloc(fc->encoding3,
                                          sizeof(short) *
                                          (fc->length + 2));

    update_encodings(fc);

  }

  return 0;
}


PUBLIC int
vrna_sequence_add(vrna_fold_compound_t  *vc,
                  const char            *string,
                  unsigned int          options)
{
  unsigned int  add_length;
  int           ret = 0;

  if ((vc) &&
      (vc->type == VRNA_FC_TYPE_SINGLE) &&
      (string)) {
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
    vc->sequence_encoding[0]                            =
      vc->sequence_encoding[vc->length + add_length];

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
    vc->sequence_encoding2[0]                           = (short)(vc->length + add_length);

    /* finally, increase length property of the fold compound */
    vc->length = vc->length + add_length;

    ret = 1;
  }

  return ret;
}


PUBLIC int
vrna_msa_add(vrna_fold_compound_t     *fc,
             const char               **alignment,
             const char               **names,
             const unsigned char      *orientation,
             const unsigned long long *start,
             const unsigned long long *genome_size,
             unsigned int             options)
{
  int ret;

  ret = 0;

  if ((fc) &&
      (fc->type == VRNA_FC_TYPE_COMPARATIVE) &&
      (alignment)) {
    size_t      s, ss, cnt, num_names, num_orientations, num_starts, num_genome_sizes;
    vrna_msa_t  *msa;

    num_names = num_orientations = num_starts = num_genome_sizes = 0;

    /* add the sequence to the nucleotides container */
    fc->alignment = (vrna_msa_t *)vrna_realloc(fc->alignment,
                                               sizeof(vrna_msa_t) *
                                               (fc->strands + 1));

    /* count number of sequences in alignment */
    for (s = 0; alignment[s]; s++);

    msa = &(fc->alignment[fc->strands]);

    msa->n_seq        = s;
    msa->sequences    = (vrna_seq_t *)vrna_alloc(sizeof(vrna_seq_t) * s);
    msa->orientation  = NULL;
    msa->start        = NULL;
    msa->genome_size  = NULL;
    msa->a2s          = NULL; /* alignment column to nt number mapping */
    msa->gapfree_seq  = NULL; /* gap-free sequence */
    msa->gapfree_size = NULL; /* gap-free sequence length */

    if (names) {
      for (s = 0; s < msa->n_seq; s++) {
        if (!names[s])
          break;

        num_names++;
      }

      if (num_names != msa->n_seq) {
        vrna_log_warning("vrna_msa_add(): "
                             "Too few names provided for sequences in MSA input! "
                             "Expected %u but received %u ",
                             msa->n_seq,
                             num_names);
      }
    }

    for (s = 0; alignment[s]; s++) {
      set_sequence(&(msa->sequences[s]),
                   alignment[s],
                   (s < num_names) ? names[s] : NULL,
                   &(fc->params->model_details),
                   options);
    }

    if (orientation) {
      /* check whether the all orientations are provided */
      for (s = 0; s < msa->n_seq; s++) {
        if (!orientation[s])
          break;

        num_orientations++;
      }

      if (s != msa->n_seq) {
        vrna_log_warning("vrna_msa_add(): "
                             "Too few orientations provided for sequences in MSA input! "
                             "Expected %u but received %u ",
                             msa->n_seq,
                             num_orientations);
      }

      msa->orientation = (unsigned char *)vrna_alloc(sizeof(unsigned char) * msa->n_seq);

      memcpy(msa->orientation, orientation, sizeof(unsigned char) * num_orientations);
    }

    if (start) {
      /* check whether the all orientations are provided */
      for (s = 0; s < msa->n_seq; s++) {
        if (!start[s])
          break;

        num_starts++;
      }

      if (s != msa->n_seq) {
        vrna_log_warning("vrna_msa_add(): "
                             "Too few start positions provided for sequences in MSA input! "
                             "Expected %u but received %u ",
                             msa->n_seq,
                             num_starts);
      }

      msa->start = (unsigned long long *)vrna_alloc(sizeof(unsigned long long) * msa->n_seq);

      memcpy(msa->start, start, sizeof(unsigned long long) * num_starts);
    }

    if (genome_size) {
      /* check whether the all orientations are provided */
      for (s = 0; s < msa->n_seq; s++) {
        if (!genome_size[s])
          break;

        num_genome_sizes++;
      }

      if (s != msa->n_seq) {
        vrna_log_warning("vrna_msa_add(): "
                             "Too few genome sizes provided for sequences in MSA input! "
                             "Expected %u but received %u ",
                             msa->n_seq,
                             num_genome_sizes);
      }

      msa->genome_size = (unsigned long long *)vrna_alloc(sizeof(unsigned long long) * msa->n_seq);

      memcpy(msa->genome_size, genome_size, sizeof(unsigned long long) * num_genome_sizes);
    }

    /* now for the gap-free sequence properties */
    msa->gapfree_seq  = (char **)vrna_alloc(sizeof(char *) * msa->n_seq);
    msa->gapfree_size = (unsigned int *)vrna_alloc(sizeof(unsigned int) * msa->n_seq);
    msa->a2s          = (unsigned int **)vrna_alloc(sizeof(unsigned int *) * msa->n_seq);

    for (s = 0; s < msa->n_seq; s++) {
      msa->gapfree_seq[s]   = vrna_seq_ungapped(msa->sequences[s].string);
      msa->gapfree_size[s]  = strlen(msa->gapfree_seq[s]);
      msa->a2s[s]           =
        (unsigned int *)vrna_alloc(sizeof(unsigned int) * (msa->sequences[s].length + 1));

      for (ss = 1, cnt = 0; ss <= msa->sequences[s].length; ss++) {
        if (msa->sequences[s].encoding[ss])
          cnt++;

        msa->a2s[s][ss] = cnt;
      }
    }

    /* increase strands counter */
    fc->strands++;
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
  unsigned int i, s;

  if (fc) {
    if (fc->type == VRNA_FC_TYPE_SINGLE) {
      for (i = 0; i < fc->strands; i++)
        free_sequence_data(&(fc->nucleotides[i]));

      free(fc->nucleotides);
      fc->nucleotides = NULL;
    } else {
      for (i = 0; i < fc->strands; i++) {
        for (s = 0; s < fc->alignment[i].n_seq; s++) {
          free_sequence_data(&(fc->alignment[i].sequences[s]));
          free(fc->alignment[i].gapfree_seq[s]);
          free(fc->alignment[i].a2s[s]);
        }

        free(fc->alignment[i].sequences);
        free(fc->alignment[i].gapfree_seq);
        free(fc->alignment[i].a2s);
        free(fc->alignment[i].gapfree_size);
        free(fc->alignment[i].genome_size);
        free(fc->alignment[i].start);
        free(fc->alignment[i].orientation);
      }
      free(fc->alignment);
      fc->alignment = NULL;
      /* free memory occupied by temporary hack in vrna_sequence_prepare() */
      free_sequence_data(&(fc->nucleotides[0]));
    }

    free(fc->strand_number);
    free(fc->strand_order);
    free(fc->strand_order_uniq);
    free(fc->strand_start);
    free(fc->strand_end);

    fc->strands           = 0;
    fc->strand_number     = NULL;
    fc->strand_order      = NULL;
    fc->strand_order_uniq = NULL;
    fc->strand_start      = NULL;
    fc->strand_end        = NULL;
  }
}


PUBLIC void
vrna_sequence_prepare(vrna_fold_compound_t *fc)
{
  unsigned int cnt, i;

  if (fc) {
    free(fc->strand_number);
    free(fc->strand_order);
    free(fc->strand_order_uniq);
    free(fc->strand_start);
    free(fc->strand_end);

    fc->strand_order      = NULL;
    fc->strand_order_uniq = NULL;
    fc->strand_start      = NULL;
    fc->strand_end        = NULL;

    fc->strand_number = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (fc->length + 2));

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        /* 1. store initial strand order */
        fc->strand_order_uniq =
          (unsigned int *)vrna_alloc(sizeof(unsigned int) * (fc->strands + 1));
        fc->strand_order =
          (unsigned int *)vrna_alloc(sizeof(unsigned int) * (fc->strands + 1));
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
        /* this sets pos. 0 and n + 1 as well for convenience reasons */
        fc->strand_number[0]              = fc->strand_number[1];
        fc->strand_number[fc->length + 1] = fc->strand_number[fc->length];

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        /*
         *  for now, comparative structure prediction does not allow for RNA-RNA interactions,
         *  so we pretend working on a single strand
         */
        fc->nucleotides = (vrna_seq_t *)vrna_realloc(fc->nucleotides,
                                                     sizeof(vrna_seq_t) * (fc->strands + 1));
        fc->nucleotides[0].string = NULL;
        fc->nucleotides[0].type   = VRNA_SEQ_RNA;
        fc->nucleotides[0].length = fc->length;

        /* 1. store initial strand order */
        fc->strand_order_uniq = (unsigned int *)vrna_alloc(sizeof(unsigned int) * 2);
        fc->strand_order      = (unsigned int *)vrna_alloc(sizeof(unsigned int) * 2);

        /* 2. mark start and end positions of sequences */
        fc->strand_start    = (unsigned int *)vrna_alloc(sizeof(unsigned int) * 2);
        fc->strand_end      = (unsigned int *)vrna_alloc(sizeof(unsigned int) * 2);
        fc->strand_start[0] = 1;
        fc->strand_end[0]   = fc->strand_start[0] + fc->length - 1;

        break;
    }
  }
}


PUBLIC int
vrna_sequence_order_update(vrna_fold_compound_t *fc,
                           const unsigned int   *order)
{
  if ((fc) &&
      (order)) {
    /* assign new order to strand_order arrays */
    memcpy(fc->strand_order_uniq, order, sizeof(unsigned int) * fc->strands);
    memcpy(fc->strand_order, order, sizeof(unsigned int) * fc->strands);

    update_strand_positions(fc);
    update_sequence(fc);
    update_encodings(fc);

    return 1;
  }

  return 0;
}


PRIVATE void
update_strand_positions(vrna_fold_compound_t *fc)
{
  unsigned int *order = fc->strand_order;

  /* now, update strand_start/end positions and strand_number association */
  fc->strand_start[order[0]]  = 1;
  fc->strand_end[order[0]]    = fc->strand_start[order[0]] +
                                fc->nucleotides[order[0]].length -
                                1;

  for (size_t j = fc->strand_start[order[0]]; j <= fc->strand_end[order[0]]; j++)
    fc->strand_number[j] = order[0];

  for (size_t i = 1; i < fc->strands; i++) {
    fc->strand_start[order[i]]  = fc->strand_end[order[i - 1]] +
                                  1;
    fc->strand_end[order[i]]    = fc->strand_start[order[i]] +
                                  fc->nucleotides[order[i]].length -
                                  1;

    for (size_t j = fc->strand_start[order[i]]; j <= fc->strand_end[order[i]]; j++)
      fc->strand_number[j] = order[i];
  }

  /* also set pos. 0 and n + 1 for convenience reasons */
  fc->strand_number[0]              = fc->strand_number[1];
  fc->strand_number[fc->length + 1] = fc->strand_number[fc->length];
}

PRIVATE void
update_sequence(vrna_fold_compound_t *fc)
{
  unsigned int *order = fc->strand_order;

  /* update the global concatenated sequence string */
  for (size_t i = 0; i < fc->strands; i++)
    memcpy(fc->sequence + fc->strand_start[order[i]] - 1,
           fc->nucleotides[order[i]].string,
           sizeof(char) * fc->nucleotides[order[i]].length);
}


PRIVATE void
update_encodings(vrna_fold_compound_t *fc)
{
  unsigned int *order = fc->strand_order;

  /* Update global sequence encoding(s) for current strand order */
  for (size_t i = 0; i < fc->strands; i++)
    memcpy(fc->sequence_encoding + fc->strand_start[order[i]],
           fc->nucleotides[order[i]].encoding + 1,
           sizeof(short) * fc->nucleotides[order[i]].length);

  fc->sequence_encoding[0]              = fc->sequence_encoding[fc->length];
  fc->sequence_encoding[fc->length + 1] = fc->sequence_encoding[1];

  for (size_t i = 0; i < fc->strands; i++) {
    short *enc = vrna_seq_encode_simple(fc->nucleotides[order[i]].string,
                                        &(fc->params->model_details));
    memcpy(fc->sequence_encoding2 + fc->strand_start[order[i]],
           enc + 1,
           sizeof(short) * fc->nucleotides[order[i]].length);
    free(enc);
  }

  fc->sequence_encoding2[0]               = (short)fc->length;
  fc->sequence_encoding2[fc->length + 1]  = fc->sequence_encoding2[1];

  /* Update 5' and 3' neighbor encodings for current strand order */

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

  for (size_t i = 1; i < obj->length; i++) {
    if (obj->encoding[i] == 0)
      obj->encoding5[i + 1] = obj->encoding5[i];
    else
      obj->encoding5[i + 1] = obj->encoding[i];
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
  obj->string     = NULL;
  obj->name       = NULL;
  obj->encoding   = NULL;
  obj->encoding5  = NULL;
  obj->encoding3  = NULL;
  obj->type       = VRNA_SEQ_UNKNOWN;
  obj->length     = 0;
}
