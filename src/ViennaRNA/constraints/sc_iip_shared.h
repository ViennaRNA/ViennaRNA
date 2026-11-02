#ifndef VIENNA_RNA_PACKAGE_CONSTRAINTS_SOFT_IIP_SHARED_H
#define VIENNA_RNA_PACKAGE_CONSTRAINTS_SOFT_IIP_SHARED_H

#include "ViennaRNA/datastructures/array.h"

typedef struct {
  unsigned int  a;          /* start of 5' interval */
  unsigned int  b;          /* end of 5' interval */
  unsigned int  c;          /* start of 3' interval */
  unsigned int  d;          /* end of 3' interval */
  int           score;      /* interval score */
  FLT_OR_DBL    exp_score;  /* interval score as Boltzmann factor */
  double        energy;     /* energy for the interaction between the intervals */
} interact_region_s;


struct vrna_crosslink_data_s {
  vrna_array(interact_region_s)   regions;  /* list of interacting regions */
  vrna_array(interact_region_s **) overlap5;  /* for each nucleotide position i we store the corresponding interval pair (I, I') where I overlaps with i */
  vrna_array(interact_region_s **) overlap3;  /* for each nucleotide position i we store the corresponding interval pair (I, I') where I' overlaps with i */
  int         *A;
  int         *B;
  FLT_OR_DBL  *qA;
  FLT_OR_DBL  *qB;
};

#endif
