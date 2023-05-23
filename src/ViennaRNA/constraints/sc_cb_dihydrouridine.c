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


/* dihydrouridines de-stabilize stacking by at least 1.5kcal/mol, Dalluge et al. 1996 */
PRIVATE int
D_stack_correction(vrna_fold_compound_t *fc,
                   int                  i,
                   int                  j,
                   int                  k,
                   int                  l,
                   void                 *data)
{
  short *enc = (short *)data;

  if ((i + 1 == k) &&
      (l + 1 == j) &&
      ((enc[i] == 5) ||
       (enc[j] == 5) ||
       (enc[k] == 5) ||
       (enc[l] == 5)))
    return 150;

  return 0; /* return 0 by default */
}


PUBLIC int
vrna_sc_mod_dihydrouridine(vrna_fold_compound_t *fc,
                           const unsigned int   *modification_sites,
                           unsigned int         options)
{
  unsigned int ret = 0;

  if ((fc) &&
      (modification_sites)) {
    short *enc = (short *)vrna_alloc(sizeof(short) * (fc->length + 2));
    memcpy(enc, fc->sequence_encoding, sizeof(short) * (fc->length + 1));

    for (unsigned int i = 0; modification_sites[i]; i++)
      if (modification_sites[i] <= fc->length) {
        enc[modification_sites[i]] = 5;
        ret++;
      }

    if (vrna_sc_multi_cb_add(fc,
                             &D_stack_correction,
                             NULL,
                             (void *)enc,
                             &free,
                             VRNA_DECOMP_PAIR_IL))
      return ret;
  }

  return 0;
}
