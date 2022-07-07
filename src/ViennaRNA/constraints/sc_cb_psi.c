#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/fold_compound.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/params/default.h"

#ifndef INLINE
# ifdef __GNUC__
#   define INLINE inline
# else
#   define INLINE
# endif
#endif

#define DEBUG

/*  Absolute stacking energy parameters for Psi-A pairs */
/*  taken from Hudson et al. 2013 */
PRIVATE int stacking_psi_A_37[6][6] = {
  /*  N   A     C   G     U     P */
  {   0,  0,    0,  0,    0,    0 }, /*  N */
  {   0,  0,    0,  0,    -210, 0 }, /*  A */
  {   0,  0,    0,  -249, 0,    0 }, /*  C */
  {   0,  0,    -220, 0,  0,    0 }, /*  G */
  {   0,  -274, 0,  0,    0,    0 }, /*  U */
  {   0,  0,    0,  0,    0,    0 }, /*  P */
};

PRIVATE int stacking_psi_A_dH[6][6] = {
  /*  N   A     C   G     U     P */
  {   0,  0,    0,  0,    0,    0 }, /*  N */
  {   0,  0,    0,  0,    -1247, 0 }, /*  A */
  {   0,  0,    0,  -1729, 0,    0 }, /*  C */
  {   0,  0,    -1119, 0,  0,    0 }, /*  G */
  {   0,  -2694, 0,  0,    0,    0 }, /*  U */
  {   0,  0,    0,  0,    0,    0 }, /*  P */
};

PRIVATE int stacking_A_psi_37[6][6] = {
  /*  N   A     C     G     U     P */
  {   0,  0,    0,    0,    0,    0 }, /*  N */
  {   0,  0,    0,    0,    -162, 0 }, /*  A */
  {   0,  0,    0,    -329, 0,    0 }, /*  C */
  {   0,  0,    -277, 0,    0,    0 }, /*  G */
  {   0,  -280, 0,    0,    0,    0 }, /*  U */
  {   0,  0,    0,    0,    0,    0 }, /*  P */
};

PRIVATE int stacking_A_psi_dH[6][6] = {
  /*  N   A     C     G     U     P */
  {   0,  0,    0,    0,    0,    0 }, /*  N */
  {   0,  0,    0,    0,    -2081, 0 }, /*  A */
  {   0,  0,    0,    -2407, 0,    0 }, /*  C */
  {   0,  0,    -1623, 0,    0,    0 }, /*  G */
  {   0,  -2208, 0,    0,    0,    0 }, /*  U */
  {   0,  0,    0,    0,    0,    0 }, /*  P */
};


PRIVATE int term_psi_A_37 = 31;
PRIVATE int term_psi_A_dH = -204;


PRIVATE int stack_psi_A_diff[6][6] = {
  0
};
PRIVATE int stack_A_psi_diff[6][6] = {
  0
};

PRIVATE int delta_terminal_psi_A  = 0; /*  Terminal Psi-A pairs at helix ends more stable than their AU counterpart */

PRIVATE INLINE void
init_psi_A_params(vrna_param_t *P)
{
  vrna_md_t *md     = &(P->model_details);
  double    tempf   = (md->temperature + K0) / (37. + K0);
  unsigned int  pair_AU = md->pair[1][4];
  unsigned int  pair_UA = md->pair[4][1];
  char nt[6] = {'\0', 'A', 'C','G','U'};
  int e;

  /* according to Hudson et al. 2013, terminal psi-A is de-stabilizing by 0.31 kcal/mol */
  if (term_psi_A_37 != 0)
    delta_terminal_psi_A = ((term_psi_A_dH) - ((term_psi_A_dH) - (term_psi_A_37)) * tempf) - P->TerminalAU;
#ifdef DEBUG
  printf("delta TerminalPsi-A, TerminalUA: %d\n", delta_terminal_psi_A);
#endif
  /* compute differences to 'regular' stacking energies */
  for (size_t si = 1; si < 5; si++)
    for (size_t sj = 1; sj < 5; sj++) {
      unsigned int tt = md->pair[sj][si];
      if (tt) {
        if (stacking_psi_A_37[si][sj] != 0) {
          e = ((stacking_psi_A_dH[si][sj]) - ((stacking_psi_A_dH[si][sj]) - (stacking_psi_A_37[si][sj])) * tempf);
          stack_psi_A_diff[si][sj] = e - P->stack[pair_UA][tt];
#ifdef DEBUG
          printf("d(psi_A, %c%c) = %d = %d - %d\n", nt[si], nt[sj], stack_psi_A_diff[si][sj], e, P->stack[pair_UA][tt]);
#endif
        }

        if (stacking_A_psi_37[si][sj] != 0) {
          e = ((stacking_A_psi_dH[si][sj]) - ((stacking_A_psi_dH[si][sj]) - (stacking_A_psi_37[si][sj])) * tempf);
          stack_A_psi_diff[si][sj] = e - P->stack[pair_AU][tt];
#ifdef DEBUG
          printf("d(A_psi, %c%c) = %d = %d - %d\n", nt[si], nt[sj], stack_A_psi_diff[si][sj], e, P->stack[pair_AU][tt]);
#endif
        }
      }
    }
}


PRIVATE INLINE int
mismatch_psi_A(int  i,
               int  j,
               int  i1,
               int  j1)
{
  /* Return correction for terminal psi-A pairs */
  if ((i == 5 && j == 1) ||
      (i == 1 && j == 5))
    return delta_terminal_psi_A;

  return 0;
}


PRIVATE int
sc_psi_A_IL(int                  i,
            int                  j,
            int                  k,
            int                  l,
            short                *enc)
{
  /* stacks */
  int e     = 0;
  int enc_i = enc[i];
  int enc_j = enc[j];
  int enc_k = enc[k];
  int enc_l = enc[l];

  if ((i + 1 == k) &&
      (l + 1 == j)) {
    /* correct for known stacks that involve at least one pseudouridine */

    /* 1. P-A pair enclosing other pair (k,l) */
    if ((enc_i == 5) && (enc_j == 1))
      return stack_psi_A_diff[enc_k][enc_l];
    /*  2. A-P pair enclosed by other pair (i,j) */
    else if ((enc_k == 1) && (enc_l == 5))
      return stack_psi_A_diff[enc_j][enc_i];
    /*  3. A-P pair enclosing other pair (k, l) */
    else if ((enc_i == 1) && (enc_j == 5))
      return stack_A_psi_diff[enc_k][enc_l];
    /*  4. P-A pair enclosed by other pair (i,j) */
    else if ((enc_k == 5) && (enc_l == 1))
      return stack_A_psi_diff[enc_j][enc_i];
  } else {
    if ((enc_i == 1 && enc_j == 5) ||
        (enc_j == 1 && enc_i == 5))
      e += delta_terminal_psi_A;

    if ((enc_k == 1 && enc_l == 5) ||
        (enc_l == 1 && enc_k == 5))
      e += delta_terminal_psi_A;
  }

  return e;
}

PRIVATE int
sc_psi_A_hp(vrna_fold_compound_t  *fc,
            int                   i,
            int                   j,
            int                   k,
            int                   l,
            void                  *data)
{
  short *enc = (short *)data;

  return mismatch_psi_A(enc[i], enc[j], enc[i + 1], enc[j - 1]);
}


PRIVATE int
sc_psi_A_int(vrna_fold_compound_t *fc,
             int                  i,
             int                  j,
             int                  k,
             int                  l,
             void                 *data)
{
  short *enc = (short *)data;

  return sc_psi_A_IL(i, j, k, l, enc);
}


PRIVATE int
sc_psi_A_ml(vrna_fold_compound_t  *fc,
         int                   i,
         int                   j,
         int                   k,
         int                   l,
         void                  *data)
{
  short *enc = (short *)data;

  return mismatch_psi_A(enc[i], enc[j], enc[i + 1], enc[j - 1]);
}


PRIVATE int
sc_psi_A_stem(vrna_fold_compound_t  *fc,
            int                   i,
            int                   j,
            int                   k,
            int                   l,
            void                  *data)
{
  short *enc = (short *)data;

  return mismatch_psi_A(enc[l], enc[k], enc[l + 1], enc[k - 1]);
}


PRIVATE int
sc_psi_A_ext_stem_ext(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   k,
                    int                   l,
                    void                  *data)
{
  short *enc = (short *)data;

  return mismatch_psi_A(enc[k], enc[i], enc[k + 1], enc[i - 1]);
}


PRIVATE int
sc_psi_A_ext_ext_stem(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   k,
                    int                   l,
                    void                  *data)
{
  short *enc = (short *)data;

  return mismatch_psi_A(enc[j], enc[l], enc[j + 1], enc[l - 1]);
}


PRIVATE int
sc_psi_A_ext_stem_outside(vrna_fold_compound_t  *fc,
                        int                   i,
                        int                   j,
                        int                   k,
                        int                   l,
                        void                  *data)
{
  short *enc = (short *)data;

  return mismatch_psi_A(enc[l], enc[k], enc[l + 1], enc[k - 1]);
}


PRIVATE int
sc_psi_A_ml_ml_stem(vrna_fold_compound_t  *fc,
                  int                   i,
                  int                   j,
                  int                   k,
                  int                   l,
                  void                  *data)
{
  short *enc = (short *)data;

  return mismatch_psi_A(enc[j], enc[l], enc[j + 1], enc[l - 1]);
}


PUBLIC void
vrna_sc_psi(vrna_fold_compound_t  *fc,
            const unsigned int    *modification_sites)
{
  if ((fc) &&
      (modification_sites)) {
    short *enc = (short *)vrna_alloc(sizeof(short) * (fc->length + 2));
    memcpy(enc, fc->sequence_encoding, sizeof(short) * (fc->length + 1));

    for (unsigned int i = 0; modification_sites[i]; i++)
      if (modification_sites[i] <= fc->length)
        enc[modification_sites[i]] = 5;

    init_psi_A_params(fc->params);
    vrna_sc_multi_cb_add(fc, &sc_psi_A_hp, NULL, (void *)enc, &free, VRNA_DECOMP_PAIR_HP);
    vrna_sc_multi_cb_add(fc, &sc_psi_A_int, NULL, (void *)enc, NULL, VRNA_DECOMP_PAIR_IL);
    vrna_sc_multi_cb_add(fc, &sc_psi_A_ml, NULL, (void *)enc, NULL, VRNA_DECOMP_PAIR_ML);
    vrna_sc_multi_cb_add(fc, &sc_psi_A_stem, NULL, (void *)enc, NULL, VRNA_DECOMP_EXT_STEM);
    vrna_sc_multi_cb_add(fc, &sc_psi_A_ext_stem_ext, NULL, (void *)enc, NULL, VRNA_DECOMP_EXT_STEM_EXT);
    vrna_sc_multi_cb_add(fc, &sc_psi_A_ext_ext_stem, NULL, (void *)enc, NULL, VRNA_DECOMP_EXT_EXT_STEM);
    vrna_sc_multi_cb_add(fc,
                         &sc_psi_A_ext_stem_outside, NULL,
                         (void *)enc,
                         NULL,
                         VRNA_DECOMP_EXT_STEM_OUTSIDE);
    vrna_sc_multi_cb_add(fc, &sc_psi_A_stem, NULL, (void *)enc, NULL, VRNA_DECOMP_ML_STEM);
    vrna_sc_multi_cb_add(fc, &sc_psi_A_ml_ml_stem, NULL, (void *)enc, NULL, VRNA_DECOMP_ML_ML_STEM);
  }
}
