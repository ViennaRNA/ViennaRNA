/* Crosslinking data handling */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/params/default.h"
#include "ViennaRNA/params/constants.h" /* defines MINPSCORE */
#include "ViennaRNA/datastructures/array.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/sequences/alignments.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/io/utils.h"
#include "ViennaRNA/io/file_formats.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/probing/basic.h"


typedef struct {
  unsigned int  a;      /* start of 5' interval */
  unsigned int  b;      /* end of 5' interval */
  unsigned int  c;      /* start of 3' interval */
  unsigned int  d;      /* end of 3' interval */
  double        score;  /* interval score */
} interact_region_s;


struct vrna_crosslink_data_s {
  vrna_array(interact_region_s)   regions;  /* list of interacting regions */
  vrna_array(interact_region_s **) overlap5;
  vrna_array(interact_region_s **) overlap3;
};



PUBLIC int
vrna_sc_xl_prepare(vrna_fold_compound_t  *fc)
{
  if ((fc) &&
      (fc->sc) &&
      (fc->sc->data)) {
    /* here we assume that sc->data actually points to a vrna_crosslink_data_s */
     
  }
  
  return 0;
}




PRIVATE int
xlink_hp_en(vrna_fold_compound_t  *fc,
            int                   i,
            int                   j,
            int                   k,
            int                   l,
            void                  *data)
{
  int                           energy = 0;
  struct vrna_crosslink_data_s  *dat = (struct vrna_crosslink_data_s *)data;

  /* go through each interval that overlaps with position i to check whether j connects both segments */
  if (dat->overlap5[i])
    for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5[i]); cnt++) {
      interact_region_s *r = dat->overlap5[i][cnt];
      if ((r->b < j) &&
          (r->c <= j) &&
          (r->d >= j))
      energy += r->score;
    }

  return energy;

}


PRIVATE int
xlink_il_en(vrna_fold_compound_t  *fc,
            int                   i,
            int                   j,
            int                   k,
            int                   l,
            void                  *data)
{
  int                           energy = 0;
  struct vrna_crosslink_data_s  *dat = (struct vrna_crosslink_data_s *)data;

  /*  go through each interval that overlaps with position i to check whether
   *  k does not overlap with the 5' interval and j connects both segments */
  if (dat->overlap5[i])
    for (size_t cnt = 0; cnt < vrna_array_size(dat->overlap5[i]); cnt++) {
      interact_region_s *r = dat->overlap5[i][cnt];
      if ((r->b < k) &&
          (r->c <= j) &&
          (r->d >= j))
      energy += r->score;
    }

  return energy;

}
