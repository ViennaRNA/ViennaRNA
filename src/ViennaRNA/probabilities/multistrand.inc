/*
 *                Multistrand base pair probabilities
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
typedef struct {
  unsigned int  strands;
  FLT_OR_DBL    *Y5;
  FLT_OR_DBL    **Y5p;
  FLT_OR_DBL    **Y3;
  FLT_OR_DBL    **Y3p;
} multistrand_t;


PRIVATE void
init_multistrand(vrna_fold_compound_t *fc,
                 helper_arrays_t      d,
                 void                 *data);


PRIVATE void
update_multistrand(vrna_fold_compound_t *fc,
                   unsigned int         l,
                   helper_arrays_t      d,
                   void                 *data);


PRIVATE void
multistrand_update_Y5(vrna_fold_compound_t  *fc,
                      unsigned int          l,
                      helper_arrays_t       dp_helpers,
                      void                  *data);


PRIVATE void
multistrand_update_Y3(vrna_fold_compound_t  *fc,
                      unsigned int          l,
                      helper_arrays_t       dp_helpers,
                      void                  *data);


PRIVATE void
multistrand_contrib(vrna_fold_compound_t  *fc,
                    unsigned int          l,
                    helper_arrays_t       dp_helpers,
                    void                  *data);


PRIVATE multistrand_t *
get_multistrand_helpers(unsigned int strands);


PRIVATE void
free_multistrand_helpers(void *data);

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

/*#####################################
 *# Multistrand implementation        #
 *#####################################
 */
PRIVATE multistrand_t *
get_multistrand_helpers(unsigned int strands)
{
  multistrand_t *data;

  data          = (multistrand_t *)vrna_alloc(sizeof(multistrand_t));
  data->strands = strands;

  return data;
}


PRIVATE void
free_multistrand_helpers(void *data)
{
  multistrand_t *d = (multistrand_t *)data;

  free(d->Y5);

  if (d->Y5p)
    for (unsigned int s = 0; s < d->strands; s++)
      free(d->Y5p[s]);

  free(d->Y5p);

  if (d->Y3)
    for (unsigned int s = 0; s < d->strands; s++)
      free(d->Y3[s]);

  free(d->Y3);

  if (d->Y3p)
    for (unsigned int s = 0; s < d->strands; s++)
      free(d->Y3p[s]);

  free(d->Y3p);

  free(d);
}


PRIVATE void
init_multistrand(vrna_fold_compound_t *fc,
                 helper_arrays_t      d,
                 void                 *data)
{
  unsigned int  n       = fc->length;
  unsigned int  strands = fc->strands;
  multistrand_t *dd     = (multistrand_t *)data;

  dd->strands = strands;

  /* reserve memory for probability output */
  dd->Y5  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * strands);
  dd->Y5p = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * strands);
  dd->Y3  = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * strands);
  dd->Y3p = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * strands);

  for (unsigned int s = 0; s < strands; s++) {
    dd->Y5p[s] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
    dd->Y3[s]  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
    dd->Y3p[s] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
  }
}


PRIVATE void
update_multistrand(vrna_fold_compound_t *fc,
                   unsigned int         l,
                   helper_arrays_t      d,
                   void                 *data)
{
  if (l > 1) {
    multistrand_update_Y5(fc, l - 1, d, data);
    multistrand_update_Y3(fc, l - 1, d, data);
  }
}


PRIVATE void
multistrand_update_Y5(vrna_fold_compound_t  *fc,
                      unsigned int          l,
                      helper_arrays_t       dp_helpers,
                      void                  *data)
{
  short             *S, *S1;
  unsigned int      i, s, *se, *sn, type, n, end;
  int               *my_iindx;
  FLT_OR_DBL        *probs, *q, *scale, qtmp;
  FLT_OR_DBL        *Y5;
  FLT_OR_DBL        **Y5p;
  vrna_md_t         *md;
  vrna_exp_param_t  *pf_params;
  struct sc_ext_exp_dat *sc_wrapper;
  sc_ext_exp_cb         sc_red_stem;
  sc_ext_exp_split      sc_split;
  multistrand_t         *dd;
  
  dd          = (multistrand_t *)data;
  n           = fc->length;
  sn          = fc->strand_number;
  se          = fc->strand_end;
  my_iindx    = fc->iindx;
  q           = fc->exp_matrices->q;
  probs       = fc->exp_matrices->probs;
  scale       = fc->exp_matrices->scale;
  pf_params   = fc->exp_params;
  md          = &(pf_params->model_details);
  S           = fc->sequence_encoding2;
  S1          = fc->sequence_encoding;
  sc_wrapper  = &(dp_helpers->sc_wrapper_ext);
  sc_red_stem = sc_wrapper->red_stem;
  sc_split    = sc_wrapper->split;
  Y5          = dd->Y5;
  Y5p         = dd->Y5p;

  /* compute Y5 for all strands */
  for (s = 0; s < fc->strands; s++) {
    unsigned int j;

    Y5[s] = 0;

    if ((se[s] < l) &&
        (sn[l] == sn[l + 1])) {
      /* pre-compute newly available Y5p[s][j] with j == l + 1 */
      end = se[s];
      j   = l + 1;

      Y5p[s][j] = 0.;

      if (probs[my_iindx[end] - j] > 0) {
        type  = vrna_get_ptype_md(S[j], S[end], md);
        qtmp  = probs[my_iindx[end] - j] *
                vrna_exp_E_exterior_stem(type,
                                    S1[j - 1],
                                    -1,
                                    pf_params) *
                scale[2];

        if (sc_red_stem)
          qtmp *= sc_red_stem(j, end, j, end, sc_wrapper);

        Y5p[s][j] += qtmp;
      }

      for (i = 1; i < end; i++) {
        if ((probs[my_iindx[i] - j] > 0) &&
            (sn[i] == sn[i + 1])) {
          type  = vrna_get_ptype_md(S[j], S[i], md);
          qtmp  = probs[my_iindx[i] - j] *
                  vrna_exp_E_exterior_stem(type,
                  S1[j - 1],
                  S1[i + 1],
                  pf_params) *
                  q[my_iindx[i + 1] - end] *
                  scale[2];

          if (sc_red_stem)
            qtmp *= sc_red_stem(j, i, j, i, sc_wrapper);
          if (sc_split)
            qtmp *= sc_split(i, end, i + 1, sc_wrapper);

          Y5p[s][j] += qtmp;
        }
      }

      if ((probs[my_iindx[i] - j] > 0) &&
          (sn[i] == sn[i + 1])) {
        type  = vrna_get_ptype_md(S[j], S[i], md);
        qtmp  = probs[my_iindx[i] - j] *
                vrna_exp_E_exterior_stem(type,
                                    S1[j - 1],
                                    S1[i + 1],
                                    pf_params) *
                scale[2];

        if (sc_red_stem)
          qtmp *= sc_red_stem(j, i, j, i, sc_wrapper);

        Y5p[s][j] += qtmp;
      }

      /* recompute Y5[s] */
      Y5[s] += Y5p[s][l + 1];
      for (j = l + 2; j <= n; j++) {
        qtmp = q[my_iindx[l + 1] - (j - 1)] *
               Y5p[s][j];

        if (sc_split)
          qtmp *= sc_split(l + 1, j, j, sc_wrapper);

        Y5[s] += qtmp;
      }
    }
  }
}


PRIVATE void
multistrand_update_Y3(vrna_fold_compound_t  *fc,
                      unsigned int          l,
                      helper_arrays_t       dp_helpers,
                      void                  *data)
{
  short             *S, *S1;
  unsigned int      i, j, k, s, n, start, type, *ss, *sn;
  int               *my_iindx;
  FLT_OR_DBL        *q, *probs, *scale, qtmp;
  FLT_OR_DBL        **Y3;
  FLT_OR_DBL        **Y3p;
  vrna_md_t         *md;
  vrna_exp_param_t  *pf_params;
  struct sc_ext_exp_dat   *sc_wrapper;
  sc_ext_exp_cb           sc_red_stem;
  sc_ext_exp_split        sc_split;
  multistrand_t         *dd;

  dd          = (multistrand_t *)data;
  n           = fc->length;
  sn          = fc->strand_number;
  ss          = fc->strand_start;
  S           = fc->sequence_encoding2;
  S1          = fc->sequence_encoding;
  my_iindx    = fc->iindx;
  q           = fc->exp_matrices->q;
  probs       = fc->exp_matrices->probs;
  scale       = fc->exp_matrices->scale;
  pf_params   = fc->exp_params;
  md          = &(pf_params->model_details);
  sc_wrapper  = &(dp_helpers->sc_wrapper_ext);
  sc_red_stem = sc_wrapper->red_stem;
  sc_split    = sc_wrapper->split;
  Y3          = dd->Y3;
  Y3p         = dd->Y3p;

  for (s = 0; s < fc->strands; s++) {
    start = ss[s];
    if (start == l + 1) {
      for (i = 1; i < start; i++) {
        Y3p[s][i] = 0;

        if (sn[i] == sn[i + 1]) {
          if (probs[my_iindx[i] - start] > 0) {
            type = vrna_get_ptype_md(S[start], S[i], md);
            qtmp = probs[my_iindx[i] - start] *
                         vrna_exp_E_exterior_stem(type,
                                             -1,
                                             S1[i + 1],
                                             pf_params) *
                         scale[2];

            if (sc_red_stem)
              qtmp *= sc_red_stem(start, i, start, i, sc_wrapper);

            Y3p[s][i] += qtmp;
          }

          for (j = start + 1; j <= n; j++) {
            if ((probs[my_iindx[i] - j] > 0) &&
                (sn[j - 1] == sn[j])) {
              type = vrna_get_ptype_md(S[j], S[i], md);
              qtmp = probs[my_iindx[i] - j] *
                           vrna_exp_E_exterior_stem(type,
                                               S1[j - 1],
                                               S1[i + 1],
                                               pf_params) *
                           q[my_iindx[start] - (j - 1)] *
                           scale[2];

              if (sc_red_stem)
                qtmp *= sc_red_stem(j, i, j, i, sc_wrapper);
              if (sc_split)
                qtmp *= sc_split(start, j, j, sc_wrapper);

              Y3p[s][i] += qtmp;
            }
          }
        }
      }

      for (k = 1; k < start; k++) {
        Y3[s][k] = 0.;

        if (sn[k - 1] == sn[k]) {
          for (i = 1; i < k - 1; i++) {
            if (sn[i] == sn[i + 1]) {
              qtmp = q[my_iindx[i + 1] - (k - 1)] *
                     Y3p[s][i];

              if (sc_split)
                qtmp *= sc_split(i, k - 1, i + 1, sc_wrapper);

              Y3[s][k] += qtmp;
            }
          }

          Y3[s][k] += Y3p[s][k - 1];
        }
      }
    }
  }
}


PRIVATE void
multistrand_contrib(vrna_fold_compound_t  *fc,
                    unsigned int          l,
                    helper_arrays_t       dp_helpers,
                    void                  *data)
{
  unsigned char     report_bpp;
  short             *S, *S1, s5, s3;
  unsigned int      k, s, *sn, *se, *ss, end, start, type;
  int               *my_iindx, kl;
  FLT_OR_DBL        *q, *qb, *probs, tmp, qtmp, *probs_ext;
  FLT_OR_DBL        *Y5;
  FLT_OR_DBL        **Y3;
  vrna_md_t         *md;
  vrna_exp_param_t  *pf_params;
  struct sc_ext_exp_dat   *sc_wrapper;
  sc_ext_exp_cb           sc_red_stem;
  multistrand_t         *dd;

  if ((l <= 1) ||
      (l == fc->length))
    return;

  dd          = (multistrand_t *)data;
  sn          = fc->strand_number;
  ss          = fc->strand_start;
  se          = fc->strand_end;
  S           = fc->sequence_encoding2;
  S1          = fc->sequence_encoding;
  pf_params   = fc->exp_params;
  md          = &(pf_params->model_details);
  my_iindx    = fc->iindx;
  q           = fc->exp_matrices->q;
  qb          = fc->exp_matrices->qb;
  probs       = fc->exp_matrices->probs;
  sc_wrapper  = &(dp_helpers->sc_wrapper_ext);
  sc_red_stem = sc_wrapper->red_stem;
  Y5          = dd->Y5;
  Y3          = dd->Y3;
  probs_ext   = dp_helpers->probs_ext;
  report_bpp  = ((dp_helpers->options.mode & VRNA_EQM_PROBS_BP) && (probs_ext)) ? 1 : 0;

  for (k = l - 1; k > 1; k--) {
    kl = my_iindx[k] - l;
    if (qb[kl] > 0) {
      tmp = 0.;
      for (s = 0; s < fc->strands; s++) {
        end   = se[s];
        start = ss[s];
        if (end == k - 1)
          tmp += Y5[s];
        else if ((end < k - 1) &&
                 (sn[k - 1] == sn[k]))
          tmp += Y5[s] *
                 q[my_iindx[end + 1] - (k - 1)];
        else if (start == l + 1)
          tmp += Y3[s][k];
        else if ((start > l + 1) &&
                 (sn[l] == sn[l + 1]))
          tmp += Y3[s][k] *
                 q[my_iindx[l + 1] - (start - 1)];
      }

      type      = vrna_get_ptype_md(S[k], S[l], md);
      s5        = (sn[k - 1] == sn[k]) ? S1[k - 1] : -1;
      s3        = (sn[l] == sn[l + 1]) ? S1[l + 1] : -1;
      qtmp      = vrna_exp_E_exterior_stem(type,
                                      s5,
                                      s3,
                                      pf_params);

      if (sc_red_stem)
        qtmp *= sc_red_stem(k, l, k, l, sc_wrapper);

      probs[kl] += tmp * qtmp;

      if (report_bpp)
        probs_ext[k] += tmp *
                        qtmp;
    }
  }
}
