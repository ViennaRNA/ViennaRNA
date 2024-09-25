#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/eval/hairpin.h"
#include "ViennaRNA/backtrack/hairpin.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "ViennaRNA/constraints/hairpin_hc.inc"

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

/**
 *  @brief Backtrack a hairpin loop closed by @f$ (i,j) @f$
 *
 *  @note This function is polymorphic! The provided #vrna_fold_compound_t may be of type
 *  #VRNA_FC_TYPE_SINGLE or #VRNA_FC_TYPE_COMPARATIVE
 *
 */
PUBLIC int
vrna_bt_hairpin(vrna_fold_compound_t  *fc,
                unsigned int          i,
                unsigned int          j,
                int                   en,
                vrna_bps_t            bp_stack,
                vrna_bts_t            bt_stack)
{
  unsigned int  u;
  int           e;
  vrna_sc_t     *sc;

  sc = NULL;

  if ((fc) &&
      (bp_stack) &&
      (bt_stack)) {
    u = j - i - 1;

    if (fc->hc->up_hp[i + 1] < u)
      return 0;

    e = vrna_eval_hairpin(fc, i, j, VRNA_EVAL_LOOP_DEFAULT);

    if (e == en) {
      switch (fc->type) {
        case  VRNA_FC_TYPE_SINGLE:
          sc = fc->sc;
          break;

        case  VRNA_FC_TYPE_COMPARATIVE:
          if (fc->scs)
            sc = fc->scs[0];

          break;

        default:
          break;
      }

      if (sc) {
        if (sc->bt) {
          vrna_basepair_t *ptr, *aux_bps;
          aux_bps = sc->bt(i, j, i, j, VRNA_DECOMP_PAIR_HP, sc->data);
          for (ptr = aux_bps; ptr && ptr->i != 0; ptr++) {
            vrna_bps_push(bp_stack,
                          (vrna_bp_t){
                            .i = ptr->i,
                            .j = ptr->j
                          });
          }
          free(aux_bps);
        }
      }

      return 1;
    }
  }

  return 0;
}


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

PUBLIC int
vrna_BT_hp_loop(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j,
                int                   en,
                vrna_bp_stack_t       *bp_stack,
                unsigned int          *stack_count)
{
  int r = 0;

  if ((fc) &&
      (bp_stack) &&
      (stack_count)) {
    vrna_bps_t  bps = vrna_bps_init(0);
    vrna_bts_t  bts = vrna_bts_init(0);

    r = vrna_bt_hairpin(fc, i, j, en, bps, bts);

    while (vrna_bps_size(bps) > 0) {
      vrna_bp_t bp = vrna_bps_pop(bps);
      bp_stack[++(*stack_count)].i  = bp.i;
      bp_stack[*stack_count].j      = bp.j;
    }

    vrna_bps_free(bps);
    vrna_bts_free(bts);
  }

  return r;
}


#endif
