#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/backtrack/exterior.h"

#include "ViennaRNA/intern/grammar_dat.h"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC unsigned int
vrna_bt_f(vrna_fold_compound_t  *fc,
          unsigned int          i,
          unsigned int          j,
          vrna_bps_t            bp_stack,
          vrna_bts_t            bt_stack)
{
  unsigned int  ret = 0;
  int           e;

  e = INF;

  if ((fc) &&
      (bp_stack) &&
      (bt_stack) &&
      (fc->matrices)) {
    if ((i == 1) &&
        (fc->matrices->type == VRNA_MX_DEFAULT) &&
        (fc->matrices->f5)) {
      e   = fc->matrices->f5[j];
      ret = vrna_bt_exterior_f5(fc, j, bp_stack, bt_stack);
    } else if ((fc->matrices->type == VRNA_MX_WINDOW) &&
               (fc->matrices->f3_local)) {
      e   = fc->matrices->f3_local[i];
      ret = vrna_bt_exterior_f3(fc, i, j, bp_stack, bt_stack);
    }

    if ((!ret) &&
        (fc->aux_grammar)) {
      for (size_t c = 0; c < vrna_array_size(fc->aux_grammar->f); c++) {
        if ((fc->aux_grammar->f[c].cb_bt) &&
            ((ret =
                fc->aux_grammar->f[c].cb_bt(fc, i, j, e, bp_stack, bt_stack,
                                            fc->aux_grammar->f[c].data)) != 0))
          break;
      }
    }
  }

  return ret;
}
