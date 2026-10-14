/*
 *                Helix Probabilities
 *
 *                Ronny Lorenz
 *                Vienna RNA package
 */

/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE void
compute_hx_prob(vrna_fold_compound_t  *fc,
                unsigned int          l,
                helper_arrays_t       dp_helpers,
                void                  *data);


PRIVATE void
init_hxp(vrna_fold_compound_t *fc,
         helper_arrays_t      d,
         void                 *data);

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */

/*
 #################################
 # STATIC helper functions below #
 #################################
 */

/*
 *  Compute helix probabilities for all helices ending at l
 *  This function assumes probs[kl] contains the outside probabilities for 
 *  basepair (k,l)
 *  Helix probabilities are returned by callback(fc, pr, l, k, s, VRNA_EQM_PROBS_HX, data)
 *  for all helices of size s ending at l and starting at k, if any.
 */
PRIVATE void
compute_hx_prob(vrna_fold_compound_t  *fc,
                unsigned int          l,
                helper_arrays_t       dp_helpers,
                void                  *data)
{
  vrna_eqm_probs_f       cb       = dp_helpers->cb;
  vrna_eqm_probs_data_t  cb_data  = dp_helpers->cb_data;

  unsigned int      k, s, start, end;
  int               *my_iindx, kl;
  FLT_OR_DBL        *qb, *probs, p, *hxp, qbt1;
  vrna_mx_pf_t      *matrices;

  my_iindx  = fc->iindx;
  matrices  = fc->exp_matrices;
  qb        = matrices->qb;
  probs     = matrices->probs;
  hxp       = dp_helpers->hxp;

  for (k = 1; k < l; k++) {
    /* reset probabilities */
    hxp = memset(hxp, 0, sizeof(FLT_OR_DBL) * (dp_helpers->options.hx_length + 1));

    kl = my_iindx[k] - l;
    s  = 1;

    /* get initial probs, i.e. base pair probabilities */
    p = probs[kl];

    hxp[s++] = p *
               qb[kl];

    if (p > FLT_MIN) {
      /* loop over all helix sizes */
      for (s = 2; s <= dp_helpers->options.hx_length; s++) {
        start = k + s - 2;
        end   = l + 2 - s;

        if (start + 2 >= end)
          break;

        qbt1 = qb[my_iindx[start + 1] - (end - 1)];

        if (qbt1 < FLT_MIN)
          break;

        p *= vrna_exp_eval_stack(fc, start, end, VRNA_EVAL_LOOP_DEFAULT);

        hxp[s] = qbt1 *
                 p;
      }
    }

    cb(fc, hxp, l, k, s - 1, VRNA_EQM_PROBS_HX, cb_data);
  }
}


PRIVATE void
init_hxp(vrna_fold_compound_t *fc,
         helper_arrays_t      d,
         void                 *data)
{
  unsigned int    n = fc->length;

  if (d->options.hx_length == 0)
    d->options.hx_length = (n + 1) / 2;

  /* reserve memory for probability output */
  d->hxp = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (d->options.hx_length + 1));
}
