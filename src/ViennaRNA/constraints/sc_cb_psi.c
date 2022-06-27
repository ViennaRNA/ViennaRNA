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

/*  energy parameters for Psi-A pairs */
/*  taken from Hudson et al. 2013 */
PRIVATE int stacking_psi_A[6][6] = {
  /*  N   A     C   G     U     P */
  { 0,  0,    0,  0,    0,    0 }, /*  N */
  { 0,  0,    0,  0,    -77,  0 }, /*  A */
  { 0,  0,    0,  -14,  0,    0 }, /*  C */
  { 0,  0,    -9, 0,    0,    0 }, /*  G */
  { 0,  -181, 0,  0,    0,    0 }, /*  U */
  { 0,  0,    0,  0,    0,    0 }, /*  P */
};

PRIVATE int stacking_A_psi[6][6] = {
  /*  N   A     C     G     U     P */
  { 0,  0,    0,    0,    0,    0 }, /*  N */
  { 0,  0,    0,    0,    -69,  0 }, /*  A */
  { 0,  0,    0,    -105, 0,    0 }, /*  C */
  { 0,  0,    -69,  0,    0,    0 }, /*  G */
  { 0,  -170, 0,    0,    0,    0 }, /*  U */
  { 0,  0,    0,    0,    0,    0 }, /*  P */
};


/* soft constraints stacking energy correction callback for P-A stacks */
PRIVATE int
psi_A_stack_correction(vrna_fold_compound_t *fc,
                       int                  i,
                       int                  j,
                       int                  k,
                       int                  l,
                       void                 *data)
{
  short *enc = (short *)data;

  if ((i + 1 == k) &&
      (l + 1 == j)) {
    /* correct for known stacks that involve at least one pseudouridine */

    /* 1. P-A pair enclosing other pair (k,l) */
    if ((enc[i] == 5) && (enc[j] == 1))
      return stacking_psi_A[enc[k]][enc[l]];
    /*  2. A-P pair enclosed by other pair (i,j) */
    else if ((enc[k] == 1) && (enc[l] == 5))
      return stacking_psi_A[enc[j]][enc[i]];
    /*  3. A-P pair enclosing other pair (k, l) */
    else if ((enc[i] == 1) && (enc[j] == 5))
      return stacking_A_psi[enc[k]][enc[l]];
    /*  4. P-A pair enclosed by other pair (i,j) */
    else if ((enc[k] == 5) && (enc[l] == 1))
      return stacking_A_psi[enc[j]][enc[i]];

    return 0; /* return 0 by default */
  }
}

PUBLIC void
vrna_sc_psi(vrna_fold_compound_t  *fc,
            unsigned int          *modification_sites)
{
  if ((fc) &&
      (modification_sites)) {
    short *enc = (short *)vrna_alloc(sizeof(short) * (fc->length + 2));
    memcpy(enc, fc->sequence_encoding, sizeof(short) * (fc->length + 1));

    for (unsigned int i = 0; modification_sites[i]; i++)
      if (modification_sites[i] <= fc->length)
        enc[modification_sites[i]] = 5;

    vrna_sc_multi_cb_add(fc, &psi_A_stack_correction, (void *)enc, &free, VRNA_DECOMP_PAIR_IL);
  }
}
