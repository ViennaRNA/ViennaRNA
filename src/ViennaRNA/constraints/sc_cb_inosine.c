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

/*  Absolute stacking energy parameters for I-U pairs */
/*  taken from Wright et al. 2007 */
PRIVATE int stacking_I_U_37[6][6] = {
  /*  N   A     C     G     U     I */
  { 0,  0,    0,    0,    0,    0 }, /*  N */
  { 0,  0,    0,    0,    43,   0 }, /*  A */
  { 0,  0,    0,    -103, 0,    0 }, /*  C */
  { 0,  0,    -122, 0,    0,    0 }, /*  G */
  { 0,  -50,  0,    0,    0,    358 }, /*  U */
  { 0,  0,    0,    0,    266,  0 } /*  I */
};

PRIVATE int stacking_I_U_dH[6][6] = {
  /*  N   A       C       G       U     I */
  { 0,  0,      0,      0,      0,    0 }, /*  N */
  { 0,  0,      0,      0,      -822, 0 }, /*  A */
  { 0,  0,      0,      -1156,  0,    0 }, /*  C */
  { 0,  0,      -1338,  0,      0,    0 }, /*  G */
  { 0,  -1583,  0,      0,      0,    1700 }, /*  U */
  { 0,  0,      0,      0,      953,  0 } /*  I */
};

PRIVATE int stacking_U_I_37[6][6] = {
  /*  N   A     C     G     U     I */
  { 0,  0,    0,    0,    0,    0 }, /*  N */
  { 0,  0,    0,    0,    37,   0 }, /*  A */
  { 0,  0,    0,    -134, 0,    0 }, /*  C */
  { 0,  0,    -77,  0,    0,    0 }, /*  G */
  { 0,  -41,  0,    0,    0,    266 }, /*  U */
  { 0,  0,    0,    0,    223,  0 } /*  I */
};

PRIVATE int stacking_U_I_dH[6][6] = {
  /*  N   A       C       G     U       I */
  { 0,  0,      0,      0,    0,      0 }, /*  N */
  { 0,  0,      0,      0,    -1008,  0 }, /*  A */
  { 0,  0,      0,      -981, 0,      0 }, /*  C */
  { 0,  0,      -1199,  0,    0,      0 }, /*  G */
  { 0,  -1168,  0,      0,    0,      953 }, /*  U */
  { 0,  0,      0,      0,    841,    0 } /*  I */
};

PRIVATE int term_IU_37 = -133;
PRIVATE int term_IU_dH = -8;


PRIVATE int stack_IU_diff[6][6] = {
  0
};
PRIVATE int stack_UI_diff[6][6] = {
  0
};

PRIVATE int delta_terminalIU  = 0; /*  Terminal IU pairs more stable than their AU counterpart */


PRIVATE INLINE void
init_IU_params(vrna_param_t *P)
{
  vrna_md_t *md     = &(P->model_details);
  double    tempf   = (md->temperature + K0) / 37.;
  unsigned int  pair_AU = md->pair[1][4];
  unsigned int  pair_UA = md->pair[4][1];

  delta_terminalIU = ((term_IU_dH) - ((term_IU_dH) - (term_IU_37)) * tempf) - P->TerminalAU; /*  according to Wright et al. 2007, terminal IU is stabilizing by -1.33 kcal/mol */

  /* compute differences to 'regular' stacking energies */
  for (size_t si = 1; si < 5; si++)
    for (size_t sj = 1; sj < 5; sj++) {
      unsigned int tt = md->pair[sj][si];
      if (tt) {
        stack_IU_diff[si][sj] = ((stacking_I_U_dH[si][sj]) - ((stacking_I_U_dH[si][sj]) - (stacking_I_U_37[si][sj])) * tempf) - P->stack[pair_AU][tt];
        stack_UI_diff[si][sj] = ((stacking_U_I_dH[si][sj]) - ((stacking_U_I_dH[si][sj]) - (stacking_U_I_37[si][sj])) * tempf) - P->stack[pair_UA][tt];
      }
    }
}

PRIVATE INLINE int
mismatch_IU(int  i,
            int  j,
            int  i1,
            int  j1)
{
  /* Return correction for terminal IU pairs */
  if ((i == 5 && j == 4) ||
      (i == 4 && j == 5))
    return delta_terminalIU;

  return 0;
}


PRIVATE INLINE int
sc_IU_IL(int         i,
         int         j,
         int         k,
         int         l,
         const short *enc)
{
  /* stacks */
  int e     = 0;
  int enc_i = enc[i];
  int enc_j = enc[j];
  int enc_k = enc[k];
  int enc_l = enc[l];

  if ((i + 1 == k) && (l == j - 1)) {
    if ((enc_i == 5) && (enc_j == 4))
      e += stack_IU_diff[enc_k][enc_l];
    else if ((enc_l == 5) && (enc_k == 4))
      e += stack_IU_diff[enc_j][enc_i];
    else if ((enc_i == 4) && (enc_j == 5))
      e += stack_UI_diff[enc_k][enc_l];
    else if ((enc_l == 4) && (enc_k == 5))
      e += stack_UI_diff[enc_j][enc_i];
  } else {
    if ((enc_i == 4 && enc_j == 5) ||
        (enc_j == 4 && enc_i == 5))
      e += delta_terminalIU;

    if ((enc_k == 4 && enc_l == 5) ||
        (enc_l == 4 && enc_k == 5))
      e += delta_terminalIU;
  }

  return e;
}


PRIVATE int
sc_IU_hp(vrna_fold_compound_t  *fc,
          int                   i,
          int                   j,
          int                   k,
          int                   l,
          void                  *data)
{
  short *enc = (short *)data;

  return mismatch_IU(enc[i], enc[j], enc[i + 1], enc[j - 1]);
}


PRIVATE int
sc_IU_int(vrna_fold_compound_t *fc,
           int                  i,
           int                  j,
           int                  k,
           int                  l,
           void                 *data)
{
  short *enc = (short *)data;

  return sc_IU_IL(i, j, k, l, enc);
}


PRIVATE int
sc_IU_ml(vrna_fold_compound_t  *fc,
         int                   i,
         int                   j,
         int                   k,
         int                   l,
         void                  *data)
{
  short *enc = (short *)data;

  return mismatch_IU(enc[i], enc[j], enc[i + 1], enc[j - 1]);
}


PRIVATE int
sc_IU_stem(vrna_fold_compound_t  *fc,
            int                   i,
            int                   j,
            int                   k,
            int                   l,
            void                  *data)
{
  short *enc = (short *)data;

  return mismatch_IU(enc[l], enc[k], enc[l + 1], enc[k - 1]);
}


PRIVATE int
sc_IU_ext_stem_ext(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   k,
                    int                   l,
                    void                  *data)
{
  short *enc = (short *)data;

  return mismatch_IU(enc[k], enc[i], enc[k + 1], enc[i - 1]);
}


PRIVATE int
sc_IU_ext_ext_stem(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   k,
                    int                   l,
                    void                  *data)
{
  short *enc = (short *)data;

  return mismatch_IU(enc[j], enc[l], enc[j + 1], enc[l - 1]);
}


PRIVATE int
sc_IU_ext_stem_outside(vrna_fold_compound_t  *fc,
                        int                   i,
                        int                   j,
                        int                   k,
                        int                   l,
                        void                  *data)
{
  short *enc = (short *)data;

  return mismatch_IU(enc[l], enc[k], enc[l + 1], enc[k - 1]);
}


PRIVATE int
sc_IU_ml_ml_stem(vrna_fold_compound_t  *fc,
                  int                   i,
                  int                   j,
                  int                   k,
                  int                   l,
                  void                  *data)
{
  short *enc = (short *)data;

  return mismatch_IU(enc[j], enc[l], enc[j + 1], enc[l - 1]);
}



PUBLIC void
vrna_sc_inosine(vrna_fold_compound_t  *fc,
                const unsigned int    *modification_sites)
{
  if ((fc) &&
      (modification_sites)) {
    short *enc = (short *)vrna_alloc(sizeof(short) * (fc->length + 2));
    memcpy(enc, fc->sequence_encoding, sizeof(short) * (fc->length + 1));

    for (unsigned int i = 0; modification_sites[i]; i++)
      if (modification_sites[i] <= fc->length)
        enc[modification_sites[i]] = 5;

    init_IU_params(fc->params);
    vrna_sc_multi_cb_add(fc, &sc_IU_hp, NULL, (void *)enc, &free, VRNA_DECOMP_PAIR_HP);
    vrna_sc_multi_cb_add(fc, &sc_IU_int, NULL, (void *)enc, NULL, VRNA_DECOMP_PAIR_IL);
    vrna_sc_multi_cb_add(fc, &sc_IU_ml, NULL, (void *)enc, NULL, VRNA_DECOMP_PAIR_ML);
    vrna_sc_multi_cb_add(fc, &sc_IU_stem, NULL, (void *)enc, NULL, VRNA_DECOMP_EXT_STEM);
    vrna_sc_multi_cb_add(fc, &sc_IU_ext_stem_ext, NULL, (void *)enc, NULL, VRNA_DECOMP_EXT_STEM_EXT);
    vrna_sc_multi_cb_add(fc, &sc_IU_ext_ext_stem, NULL, (void *)enc, NULL, VRNA_DECOMP_EXT_EXT_STEM);
    vrna_sc_multi_cb_add(fc,
                         &sc_IU_ext_stem_outside, NULL,
                         (void *)enc,
                         NULL,
                         VRNA_DECOMP_EXT_STEM_OUTSIDE);
    vrna_sc_multi_cb_add(fc, &sc_IU_stem, NULL, (void *)enc, NULL, VRNA_DECOMP_ML_STEM);
    vrna_sc_multi_cb_add(fc, &sc_IU_ml_ml_stem, NULL, (void *)enc, NULL, VRNA_DECOMP_ML_ML_STEM);
  }
}
