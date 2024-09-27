#include <stdio.h>
#include <stdlib.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/grammar/partfunc.h"
#include "ViennaRNA/partfunc/multifold.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "ViennaRNA/constraints/exterior_hc.inc"

PRIVATE FLT_OR_DBL
mf_rule_pair(vrna_fold_compound_t *fc,
             unsigned int         i,
             unsigned int         j,
             void                 *data);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_pf_multifold_prepare(vrna_fold_compound_t *fc)
{
  if (fc)
    return vrna_gr_add_aux_exp_c(fc,
                                 &mf_rule_pair,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL);

  return 0;
}


/*
 #################################
 # STATIC helper functions below #
 #################################
 */
PRIVATE FLT_OR_DBL
mf_rule_pair(vrna_fold_compound_t *fc,
             unsigned int         i,
             unsigned int         j,
             void                 *data VRNA_UNUSED)
{
  short                 *S1, *S2, s5, s3;
  unsigned int          *sn, *ends, type, nick;
  int                   *my_iindx;
  FLT_OR_DBL            contribution, *q, *scale, qbase, tmp, tmp2;
  vrna_exp_param_t      *pf_params;
  vrna_md_t             *md;
  vrna_hc_eval_f        evaluate;
  struct hc_ext_def_dat hc_dat_local;
  vrna_sc_t             *sc;

  contribution  = 0;
  S1            = fc->sequence_encoding;
  S2            = fc->sequence_encoding2;
  pf_params     = fc->exp_params;
  md            = &(pf_params->model_details);
  sn            = fc->strand_number;
  ends          = fc->strand_end;
  q             = fc->exp_matrices->q;
  scale         = fc->exp_matrices->scale;
  my_iindx      = fc->iindx;
  sc            = fc->sc;
  evaluate      = prepare_hc_ext_def(fc, &hc_dat_local);

  if ((sn[i] != sn[j]) &&
      (evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local))) {
    /* most obious strand nick is at end of sn[i] and start of sn[j] */
    type  = vrna_get_ptype_md(S2[j], S2[i], md);
    s5    = (sn[j] == sn[j - 1]) ? S1[j - 1] : -1;
    s3    = (sn[i] == sn[i + 1]) ? S1[i + 1] : -1;
    qbase = vrna_exp_E_exterior_stem(type, s5, s3, pf_params) *
            scale[2];

    if (sc) {
      if (sc->exp_f)
        qbase *= sc->exp_f(j, i, j, i, VRNA_DECOMP_EXT_STEM, sc->data);
    }

    tmp = 0.;

    /*
     *  if (evaluate(i + 1,
     *               j - 1,
     *               ends[sn[i]],
     *               ends[sn[i]] + 1,
     *               VRNA_DECOMP_EXT_EXT_EXT,
     *               &hc_dat_local))
     */
    if (sn[i] != sn[i + 1]) {
      if ((sn[j - 1] != sn[j]) &&
          (i + 1 == j))
        tmp = 1.;
      else if (sn[j - 1] == sn[j])
        tmp = q[my_iindx[i + 1] - j + 1];
    } else if (sn[j - 1] != sn[j]) {
      tmp = q[my_iindx[i + 1] - j + 1];
    } else {
      tmp = q[my_iindx[i + 1] - ends[sn[i]]] *
            q[my_iindx[ends[sn[i]] + 1] - j + 1];

      /* check whether we find more strand nicks between i and j */
      nick = ends[sn[i]] + 1;
      while (sn[nick] != sn[j]) {
        /*
         *      if (evaluate(i + 1,
         *                   j - 1,
         *                   ends[sn[nick]],
         *                   ends[sn[nick]] + 1,
         *                   VRNA_DECOMP_EXT_EXT_EXT,
         *                   &hc_dat_local))
         */
        tmp2 = 1.;
        if (i + 1 <= ends[sn[nick]])
          tmp2 *= q[my_iindx[i + 1] - ends[sn[nick]]];

        if (ends[sn[nick]] + 1 <= j - 1)
          tmp2 *= q[my_iindx[ends[sn[nick]] + 1] - j + 1];

        tmp += tmp2;


        nick = ends[sn[nick]] + 1;
      }
    }

    contribution = qbase *
                   tmp;
  }

  return contribution;
}
