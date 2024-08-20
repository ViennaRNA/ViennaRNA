#ifndef VIENNA_RNA_PACKAGE_CONSTRAINTS_SOFT_INTERN_H
#define VIENNA_RNA_PACKAGE_CONSTRAINTS_SOFT_INTERN_H

#include "ViennaRNA/datastructures/array.h"

#define MOD_PARAMS_STACK_dG     (1 << 0)
#define MOD_PARAMS_STACK_dH     (1 << 1)
#define MOD_PARAMS_MISMATCH_dG  (1 << 2)
#define MOD_PARAMS_MISMATCH_dH  (1 << 3)
#define MOD_PARAMS_TERMINAL_dG  (1 << 4)
#define MOD_PARAMS_TERMINAL_dH  (1 << 5)
#define MOD_PARAMS_DANGLES_dG   (1 << 6)
#define MOD_PARAMS_DANGLES_dH   (1 << 7)

#define MAX_ALPHABET  (6)
#define MAX_PAIRS     (NBPAIRS + 1 + 25)


/* a container to store the data read from a json parameter file */
struct vrna_sc_mod_param_s {
  unsigned int  available;

  char          *name;
  char          one_letter_code;
  char          unmodified;
  char          fallback;
  char          pairing_partners[7];
  unsigned int  pairing_partners_encoding[7];
  unsigned int  unmodified_encoding;
  unsigned int  fallback_encoding;

  size_t        num_ptypes;
  size_t        ptypes[MAX_ALPHABET][MAX_ALPHABET];

  int           stack_dG[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET];
  int           stack_dH[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET];

  int           dangle5_dG[MAX_PAIRS][MAX_ALPHABET];
  int           dangle5_dH[MAX_PAIRS][MAX_ALPHABET];
  int           dangle3_dG[MAX_PAIRS][MAX_ALPHABET];
  int           dangle3_dH[MAX_PAIRS][MAX_ALPHABET];

  int           mismatch_dG[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET];
  int           mismatch_dH[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET];

  int           terminal_dG[MAX_PAIRS];
  int           terminal_dH[MAX_PAIRS];
};

/* the actual data structure passed around while evaluating */
typedef struct {
  short   *enc;

  size_t  strands;
  vrna_array(unsigned int  *)  modification_sites;

  size_t  ptypes[MAX_ALPHABET][MAX_ALPHABET];

  int     stack_diff[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET];

  int     dangle5_diff[MAX_PAIRS][MAX_ALPHABET];
  int     dangle3_diff[MAX_PAIRS][MAX_ALPHABET];

  int     mismatch_diff[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET];

  int     terminal_diff[MAX_PAIRS];
} energy_corrections;


#endif
