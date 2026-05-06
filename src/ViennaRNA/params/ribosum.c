#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/io/utils.h"

#include "ViennaRNA/params/ribosum.h"
#include "ViennaRNA/params/ribosum_scores.inc"

static float (*ribosum_chooser[21][21])[7][7] = {
  { 0 }, /* NA */
  { 0 }, /* NA */
  { 0 }, /* NA */
  { 0 }, /* NA */
  { 0 }, /* NA */
  { 0 }, /* NA */
  { 0 }, /* NA */
  { 0 }, /* NA */
  { 0 }, /* NA */
  { 0 }, /* NA */
  { 0 }, /* NA */
  { 0 }, /* NA */
  { 0, 0, 0, 0, 0, &dm_12_5, &dm_12_6, &dm_12_7, &dm_12_8, &dm_12_9, &dm_12_10, &dm_12_11 },
  { 0, 0, 0, 0, 0, &dm_13_5, &dm_13_6, &dm_13_7, &dm_13_8, &dm_13_9, &dm_13_10, &dm_13_11, &dm_13_12 },
  { 0, 0, 0, 0, 0, &dm_14_5, &dm_14_6, &dm_14_7, &dm_14_8, &dm_14_9, &dm_14_10, &dm_14_11, &dm_14_12, &dm_14_13 },
  { 0, 0, 0, 0, 0, &dm_15_5, &dm_15_6, &dm_15_7, &dm_15_8, &dm_15_9, &dm_15_10, &dm_15_11, &dm_15_12, &dm_15_13, &dm_15_14 },
  { 0, 0, 0, 0, 0, &dm_16_5, &dm_16_6, &dm_16_7, &dm_16_8, &dm_16_9, &dm_16_10, &dm_16_11, &dm_16_12, &dm_16_13, &dm_16_14, &dm_16_15 },
  { 0, 0, 0, 0, 0, &dm_17_5, &dm_17_6, &dm_17_7, &dm_17_8, &dm_17_9, &dm_17_10, &dm_17_11, &dm_17_12, &dm_17_13, &dm_17_14, &dm_17_15, &dm_17_16 },
  { 0, 0, 0, 0, 0, &dm_18_5, &dm_18_6, &dm_18_7, &dm_18_8, &dm_18_9, &dm_18_10, &dm_18_11, &dm_18_12, &dm_18_13, &dm_18_14, &dm_18_15, &dm_18_16, &dm_18_17 },
  { 0, 0, 0, 0, 0, &dm_19_5, &dm_19_6, &dm_19_7, &dm_19_8, &dm_19_9, &dm_19_10, &dm_19_11, &dm_19_12, &dm_19_13, &dm_19_14, &dm_19_15, &dm_19_16, &dm_19_17, &dm_19_18 },
  { 0, 0, 0, 0, 0, &dm_20_5, &dm_20_6, &dm_20_7, &dm_20_8, &dm_20_9, &dm_20_10, &dm_20_11, &dm_20_12, &dm_20_13, &dm_20_14, &dm_20_15, &dm_20_16, &dm_20_17, &dm_20_18, &dm_20_19 },
};


PUBLIC double **
vrna_ribosum_mx(double min_identity_pct,
                double max_identity_pct)
{
  int     min, max;
  float   (*scores)[7][7], dm_identity[7][7] = {0};
  double  **rscore;
  
  if (min_identity_pct > 100.)
    min_identity_pct = 100.;

  if (min_identity_pct < 0.)
    min_identity_pct = 0.;

  if (max_identity_pct > 100.)
    max_identity_pct = 100.;

  if (max_identity_pct < 0.)
    max_identity_pct = 0.;

  min = (int)(min_identity_pct + 0.5) / 5;
  max = (int)(max_identity_pct + 0.5) / 5;

  if (max < 12)
    max = 12;

  if (min < 5)
    min = 5;

  if (min == 20) {
    scores = &dm_identity;
  } else {
    if (min >= max)
      min = max - 1;

    scores = ribosum_chooser[max][min];
  }

  rscore = (double **)vrna_alloc(7 * sizeof(double *));
  for (size_t i = 0; i < 7; i++)
    rscore[i] = (double *)vrna_alloc(7 * sizeof(double));

  for (size_t i = 0; i < 7; i++)
    for (size_t j = 0; j < 7; j++)
      rscore[i][j] = (double)(*scores)[i][j];

  return rscore;
}

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

float **
get_ribosum(const char  **Alseq,
            int         n_seq,
            int         length)
{
  int   i, j, k;
  float ident   = 0;
  float minimum = 1;
  float maximum = 0.;
  float **ribo;
  double **rmx;

  ribo = (float **)vrna_alloc(7 * sizeof(float *));
  for (i = 0; i < 7; i++)
    ribo[i] = (float *)vrna_alloc(7 * sizeof(float));

  /* identify min/max pairwise distances, normalized per nucleotide (aln column) */
  for (j = 0; j < n_seq - 1; j++)
    for (k = j + 1; k < n_seq; k++) {
      ident = length - vrna_hamming_distance(Alseq[k], Alseq[j]);
      if ((ident / (length)) < minimum)
        minimum = ident / (float)(length);

      if ((ident / (length)) > maximum)
        maximum = ident / (float)(length);
    }

  /*  express min/max in terms of percent identity
   *  +0.5 to ALWAYS round up
   */
  minimum *= 100;
  maximum *= 100;
  minimum += 0.5;
  maximum += 0.5;


  if ((n_seq == 1) ||
      (minimum > 100.45)) {
    for (i = 0; i < 7; i++)
      for (j = 0; j < 7; j++)
        ribo[i][j] = 0.;

    return ribo;
  }

  rmx = vrna_ribosum_mx((double)((int)minimum),
                        (double)((int)maximum));

  for (i = 0; i < 7; i++) {
    for (j = 0; j < 7; j++)
      ribo[i][j] = (float)rmx[i][j];

    free(rmx[i]);
  }

  free(rmx);

  return ribo;
}


PUBLIC float **
readribosum(char *name)
{
  float **dm;
  char  *line;
  FILE  *fp;
  int   i = 0;
  int   who = 0;
  float a, b, c, d, e, f;
  int   translator[7] = {
    0, 5, 1, 2, 3, 6, 4
  };

  fp  = fopen(name, "r");
  dm  = (float **)vrna_alloc(7 * sizeof(float *));
  for (i = 0; i < 7; i++)
    dm[i] = (float *)vrna_alloc(7 * sizeof(float));
  while (1) {
    /*bisma hoit fertisch san*/
    line = vrna_read_line(fp);
    if (*line == '#')
      continue;

    i = 0;
    i = sscanf(line, "%f %f %f %f %f %f", &a, &b, &c, &d, &e, &f);
    if (i == 0)
      break;

    dm[translator[++who]][translator[1]]  = a;
    dm[translator[who]][translator[2]]    = b;
    dm[translator[who]][translator[3]]    = c;
    dm[translator[who]][translator[4]]    = d;
    dm[translator[who]][translator[5]]    = e;
    dm[translator[who]][translator[6]]    = f;
    free(line);
    if (who == 6)
      break;
  }
  fclose(fp);
  return dm;
}

#endif
