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

/* m6A data taken from Kierzek et al. 2022 */

/*  absolute stacking energies */
PRIVATE int stack_6U[6][6] = {
  /* N   A     C      G     U     6 */
  { 0,   0,    0,     0,    0,    0     }, /*  N */
  { 0,   0,    0,     0,    -92,  0     }, /*  A */
  { 0,   0,    0,     -179, 0,    0     }, /*  C */
  { 0,   0,    -156,  0,    -3,   0     }, /*  G */
  { 0,   -110, 0,     -69,  0,    -46   }, /*  U */
  { 0,   0,    0,     0,    -21,  0     } /*  6 */
};

PRIVATE int stack_U6[6][6] = {
  /*  N  A    C     G     U     6 */
  { 0,   0,   0,    0,    0,    0     }, /*  N */
  { 0,   0,   0,    0,    -73,  0     }, /*  A */
  { 0,   0,   0,    -172, 0,    0     }, /*  C */
  { 0,   0,   -124, 0,    -32,  0     }, /*  G */
  { 0,   -83, 0,    -32,  0,    -21   }, /*  U */
  { 0,   0,   0,    0,    145,  0     } /*  6 */
};

PRIVATE int stack_6U_diff[6][6] = {
  0
};
PRIVATE int stack_U6_diff[6][6] = {
  0
};

/*  dangle differences to unmethylated construct */
PRIVATE int dangle_5[6][6] = {
  /*  N   A    C    G    U    6 */
  {   0,  0,   0,   0,   0,   0   },  /*  N */
  {   0,  0,   0,   0,   -56, 0   },  /*  A */
  {   0,  0,   0,   0,   0,   0   },  /*  C */
  {   0,  0,   0,   0,   0,   0   },  /*  G */
  {   0,  0,   0,   0,   0,   0   },  /*  U */
  {   0,  0,   0,   0,   0,   0   }   /*  6 */
};

PRIVATE int dangle_3[6][6] = {
  /*  N   A    C    G     U     6 */
  {   0,  0,   0,   0,    0,    0   },  /*  N */
  {   0,  0,   0,   0,    -22,  0   },  /*  A */
  {   0,  0,   0,   -47,  0,    0   },  /*  C */
  {   0,  0,   -49, 0,    0,    0   },  /*  G */
  {   0,  0,   0,   0,    0,    0   },  /*  U */
  {   0,  0,   0,   0,    0,    0   }   /*  6 */
};

/*
 *  mismatch differences to unmethylated construct
 *  mismatches with unpaired m6A receive an average
 *  of -0.3 kcal/mol unless experimentally measured.
 *  mismatches with closing m6A-U pair receive an additional
 *  -0.3 kcal/mol
 */
PRIVATE int mm_CG[6][6] = {
  /*  N   A     C     G     U     6 */
  { 0, 0,   0,   0,   0,   0     },     /*  N */
  { 0, 0,   0,   0,   0,   -30   },     /*  A */
  { 0, 0,   0,   0,   0,   -30   },     /*  C */
  { 0, 0,   0,   0,   0,   -30   },     /*  G */
  { 0, 0,   0,   0,   0,   -30   },     /*  U */
  { 0, -30, -38, -30, -30, -29   }      /*  6 */
};

PRIVATE int mm_AU[6][6] = {
  /*  N   A     C     G     U     6 */
  { 0, 0,   0,   0,   0,   0     },     /*  N */
  { 0, 0,   0,   0,   0,   -30   },     /*  A */
  { 0, 0,   0,   0,   0,   -30   },     /*  C */
  { 0, 0,   0,   0,   0,   2     },     /*  G */
  { 0, 0,   0,   0,   0,   -30   },     /*  U */
  { 0, -30, -30, -30, -30, -12   }      /*  6 */
};

PRIVATE int mm_6U[6][6] = {
  /*  N   A     C     G     U     6 */
  { 0, 0,   0,   0,   0,   0     },     /*  N */
  { 0, -30, -30, -30, -30, -60   },     /*  A */
  { 0, -30, -30, -30, -30, -60   },     /*  C */
  { 0, -30, -30, -30, -30, -60   },     /*  G */
  { 0, -30, -30, -30, -30, -60   },     /*  U */
  { 0, -60, -60, -60, -60, -54   }      /*  6 */
};

PRIVATE int mm_UG[6][6] = {
  /*  N   A     C     G     U     6 */
  { 0, 0,   0,   0,   0,   0     },     /*  N */
  { 0, 0,   0,   0,   0,   -30   },     /*  A */
  { 0, 0,   0,   0,   0,   -30   },     /*  C */
  { 0, 0,   0,   0,   0,   -74   },     /*  G */
  { 0, 0,   0,   0,   0,   -30   },     /*  U */
  { 0, -30, -30, -30, -30, -30   }      /*  6 */
};

PRIVATE int mm_UA[6][6] = {
  /*  N   A     C     G     U     6 */
  { 0, 0,   0,   0,   0,   0     },     /*  N */
  { 0, 0,   0,   0,   0,   -30   },     /*  A */
  { 0, 0,   0,   0,   0,   -30   },     /*  C */
  { 0, 0,   0,   0,   0,   -30   },     /*  G */
  { 0, 0,   0,   0,   0,   -30   },     /*  U */
  { 0, -30, -30, -18, -30, -30   }      /*  6 */
};

PRIVATE int mm_m6A_avg        = -30;
PRIVATE int int_G6_G6         = 33;
PRIVATE int hp_CG_GU6AUA      = 23;
PRIVATE int delta_terminalAU  = -40; /*  Terminal 6U pairs are about 0.4 kcal/mol more stable than their AU counterpart */

PRIVATE INLINE void
init_m6A_stacks(vrna_param_t *P)
{
  vrna_md_t *md     = &(P->model_details);
  double    tempf   = (md->temperature + K0) / 37.;
  unsigned int  pair_AU = md->pair[1][4];
  unsigned int  pair_UA = md->pair[4][1];

  delta_terminalAU = -P->TerminalAU; /*  according to Kierzek et al. paper, we simply remove terminal AU as applied in ViennaRNA parameter set */

  /* compute differences to 'regular' stacking energies */
  for (size_t si = 1; si < 5; si++)
    for (size_t sj = 1; sj < 5; sj++) {
      unsigned int tt = md->pair[sj][si];
      if (tt) {
        stack_6U_diff[si][sj] = stack_6U[si][sj] - P->stack[pair_AU][tt];
        stack_U6_diff[si][sj] = stack_U6[si][sj] - P->stack[pair_UA][tt];
      }
    }
}


PRIVATE INLINE int
default_mm(int  k,
           int  l)
{
  if ((k == 5) || (l == 5))
    return mm_m6A_avg;

  return 0;
}


PRIVATE INLINE int
mismatch_m6A(int  i,
             int  j,
             int  i1,
             int  j1)
{
  /*
   * Return mismatch energy correction for pair (i,j) with i3 dangling from 3' to i
   * and j5 dangling from 5' to j, so i3 = i + 1, j5 + 1 = j
   */
  int e = default_mm(i1, j1);

  if ((i == 1) && (j == 4)) {
    e = mm_AU[i1][j1];
  } else if ((i == 2) && (j == 3)) {
    e = mm_CG[i1][j1];
  } else if (i == 4) {
    if (j == 1)
      e = mm_UA[i1][j1];
    else if (j == 3)
      e = mm_UG[i1][j1];
    else if (j == 5)
      e += mm_m6A_avg;
  } else if ((i == 5) && (j == 4)) {
    e = mm_6U[i1][j1];
  }

  if ((i == 5 && j == 4) ||
      (i == 4 && j == 5))
    e += delta_terminalAU;

  return e;
}


PRIVATE INLINE int
sc_IL(int         i,
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
      e += stack_6U_diff[enc_k][enc_l];
    else if ((enc_l == 5) && (enc_k == 4))
      e += stack_6U_diff[enc_j][enc_i];
    else if ((enc_i == 4) && (enc_j == 5))
      e += stack_U6_diff[enc_k][enc_l];
    else if ((enc_l == 4) && (enc_k == 5))
      e += stack_U6_diff[enc_j][enc_i];
  }
  /* special GG6U/GG6U 2x2 internal loops */
  else if ((i + 3 == k) &&
           (l + 3 == j) &&
           (enc[i] == 3) &&
           (enc[i + 1] == 3) &&
           (enc[i + 2] == 5) &&
           (enc[k] == 4) &&
           (enc[l] == 3) &&
           (enc[l + 1] == 3) &&
           (enc[l + 2] == 5) &&
           (enc[j] == 4)) {
    e += int_G6_G6;

    if ((enc_k == 4 && enc_l == 5) ||
        (enc_l == 4 && enc_k == 5))
      e += delta_terminalAU;

    if ((enc_i == 4 && enc_j == 5) ||
        (enc_j == 4 && enc_i == 5))
      e += delta_terminalAU;
  }
  /* other loops that receive mismatch energies, i.e. > 2x2 loops */
  else if (k - i > 2 && j - l > 2) {
    e += mismatch_m6A(enc_i, enc_j, enc[i + 1], enc[j - 1]);
    e += mismatch_m6A(enc_l, enc_k, enc[l + 1], enc[k - 1]);
  }
  /* everything else, i.e. bulges, 1xn && 2x2 loops */
  else {
    if ((enc_i == 4 && enc_j == 5) ||
        (enc_j == 4 && enc_i == 5))
      e += delta_terminalAU;

    if ((enc_k == 4 && enc_l == 5) ||
        (enc_l == 4 && enc_k == 5))
      e += delta_terminalAU;
  }

  return e;
}


PRIVATE int
sc_m6a_hp(vrna_fold_compound_t  *fc,
          int                   i,
          int                   j,
          int                   k,
          int                   l,
          void                  *data)
{
  short *enc = (short *)data;

  return mismatch_m6A(enc[i], enc[j], enc[i + 1], enc[j - 1]);
}


PRIVATE int
sc_m6a_int(vrna_fold_compound_t *fc,
           int                  i,
           int                  j,
           int                  k,
           int                  l,
           void                 *data)
{
  short *enc = (short *)data;

  return sc_IL(i, j, k, l, enc);
}


PRIVATE int
sc_m6a_ml(vrna_fold_compound_t  *fc,
          int                   i,
          int                   j,
          int                   k,
          int                   l,
          void                  *data)
{
  short *enc = (short *)data;

  return mismatch_m6A(enc[i], enc[j], enc[i + 1], enc[j - 1]);
}


PRIVATE int
sc_m6a_stem(vrna_fold_compound_t  *fc,
            int                   i,
            int                   j,
            int                   k,
            int                   l,
            void                  *data)
{
  short *enc = (short *)data;

  return mismatch_m6A(enc[l], enc[k], enc[l + 1], enc[k - 1]);
}


PRIVATE int
sc_m6a_ext_stem_ext(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   k,
                    int                   l,
                    void                  *data)
{
  short *enc = (short *)data;

  return mismatch_m6A(enc[k], enc[i], enc[k + 1], enc[i - 1]);
}


PRIVATE int
sc_m6a_ext_ext_stem(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   k,
                    int                   l,
                    void                  *data)
{
  short *enc = (short *)data;

  return mismatch_m6A(enc[j], enc[l], enc[j + 1], enc[l - 1]);
}


PRIVATE int
sc_m6a_ext_stem_outside(vrna_fold_compound_t  *fc,
                        int                   i,
                        int                   j,
                        int                   k,
                        int                   l,
                        void                  *data)
{
  short *enc = (short *)data;

  return mismatch_m6A(enc[l], enc[k], enc[l + 1], enc[k - 1]);
}


PRIVATE int
sc_m6a_ml_ml_stem(vrna_fold_compound_t  *fc,
                  int                   i,
                  int                   j,
                  int                   k,
                  int                   l,
                  void                  *data)
{
  short *enc = (short *)data;

  return mismatch_m6A(enc[j], enc[l], enc[j + 1], enc[l - 1]);
}


#if 0

PRIVATE FLT_OR_DBL
sc_m6a_mm_pf(int            i,
             int            j,
             int            k,
             int            l,
             unsigned char  d,
             void           *data)
{
  int e = sc_m6a_mm(i, j, k, l, d, data) * 10.;

  return (FLT_OR_DBL)exp(-(FLT_OR_DBL)e / ((37. + K0) * GASCONST));
}


#endif


PUBLIC void
vrna_sc_m6A(vrna_fold_compound_t  *fc,
            const unsigned int    *modification_sites)
{
  if ((fc) &&
      (modification_sites)) {
    short *enc = (short *)vrna_alloc(sizeof(short) * (fc->length + 2));
    memcpy(enc, fc->sequence_encoding, sizeof(short) * (fc->length + 1));

    for (unsigned int i = 0; modification_sites[i]; i++)
      if (modification_sites[i] <= fc->length)
        enc[modification_sites[i]] = 5;

    init_m6A_stacks(fc->params);
    vrna_sc_multi_cb_add(fc, &sc_m6a_hp, NULL, (void *)enc, &free, VRNA_DECOMP_PAIR_HP);
    vrna_sc_multi_cb_add(fc, &sc_m6a_int, NULL, (void *)enc, NULL, VRNA_DECOMP_PAIR_IL);
    vrna_sc_multi_cb_add(fc, &sc_m6a_ml, NULL, (void *)enc, NULL, VRNA_DECOMP_PAIR_ML);
    vrna_sc_multi_cb_add(fc, &sc_m6a_stem, NULL, (void *)enc, NULL, VRNA_DECOMP_EXT_STEM);
    vrna_sc_multi_cb_add(fc, &sc_m6a_ext_stem_ext, NULL, (void *)enc, NULL, VRNA_DECOMP_EXT_STEM_EXT);
    vrna_sc_multi_cb_add(fc, &sc_m6a_ext_ext_stem, NULL, (void *)enc, NULL, VRNA_DECOMP_EXT_EXT_STEM);
    vrna_sc_multi_cb_add(fc,
                         &sc_m6a_ext_stem_outside, NULL,
                         (void *)enc,
                         NULL,
                         VRNA_DECOMP_EXT_STEM_OUTSIDE);
    vrna_sc_multi_cb_add(fc, &sc_m6a_stem, NULL, (void *)enc, NULL, VRNA_DECOMP_ML_STEM);
    vrna_sc_multi_cb_add(fc, &sc_m6a_ml_ml_stem, NULL, (void *)enc, NULL, VRNA_DECOMP_ML_ML_STEM);
  }
}
