/*
 *                Accessibility Computations
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
init_accessibility(vrna_fold_compound_t *fc,
                   helper_arrays_t      eqm_prob_helpers,
                   void                 *data);


PRIVATE void
report_n_clear_accessibility(vrna_fold_compound_t *fc,
                             unsigned int         l,
                             helper_arrays_t      eqm_prob_helpers,
                             void                 *data);


PRIVATE void
compute_accessibility(vrna_fold_compound_t  *fc,
                      unsigned int          l,
                      helper_arrays_t       dp_helpers,
                      void                  *data);


PRIVATE void
compute_up_ext(vrna_fold_compound_t *fc,
               unsigned int         l,
               helper_arrays_t      dp_helpers,
               void                 *data);


PRIVATE void
compute_up_hp(vrna_fold_compound_t *fc,
              unsigned int         l,
              helper_arrays_t      dp_helpers,
              void                 *data);


PRIVATE void
compute_up_int(vrna_fold_compound_t *fc,
               unsigned int         l,
               helper_arrays_t      dp_helpers,
               void                 *data);


PRIVATE void
compute_up_mb(vrna_fold_compound_t *fc,
              unsigned int         l,
              helper_arrays_t      dp_helpers,
              void                 *data);


PRIVATE FLT_OR_DBL
contrib_ext_up(vrna_fold_compound_t *fc,
               unsigned int         i,
               unsigned int         j,
               helper_arrays_t      dp_helpers);


PRIVATE FLT_OR_DBL
contrib_ext_up_comparative(vrna_fold_compound_t *fc,
                           unsigned int         i,
                           unsigned int         j,
                           helper_arrays_t      dp_helpers);


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
PRIVATE void
compute_accessibility(vrna_fold_compound_t  *fc,
                      unsigned int          l,
                      helper_arrays_t       dp_helpers,
                      void                  *data)
{
  /* actual computation for all subsequences of up to size u ending at position l */
  compute_up_ext(fc, l, dp_helpers, data);
  compute_up_hp(fc, l, dp_helpers, data);
  compute_up_int(fc, l, dp_helpers, data);
  compute_up_mb(fc, l, dp_helpers, data);
}


PRIVATE void
init_accessibility(vrna_fold_compound_t *fc,
                   helper_arrays_t      d,
                   void                 *data)
{
  unsigned int    n = fc->length;

  if (d->options.up_length == 0)
    d->options.up_length = n;

  d->qbh = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));

  /* reserve memory for accessibility output */
  d->access_ext  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (d->options.up_length + 1));
  d->access_hp   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (d->options.up_length + 1));
  d->access_int  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (d->options.up_length + 1));
  d->access_mb   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (d->options.up_length + 1));

  /* recompute qm2 if necessary */
  if (fc->exp_matrices->qm2_real == NULL) {
    short int             *S, *S1;
    unsigned int          turn, type, i, j, u;
    int                   *my_iindx, ii, ij;
    FLT_OR_DBL            *qm, *qb, temp, *expMLbase;
    vrna_exp_param_t      *params;
    vrna_md_t             *md;
    vrna_hc_t             *hc;
    struct sc_mb_exp_dat  *sc_wrapper;

    S           = fc->sequence_encoding2;
    S1          = fc->sequence_encoding;
    hc          = fc->hc;
    qm          = fc->exp_matrices->qm;
    qb          = fc->exp_matrices->qb;
    expMLbase   = fc->exp_matrices->expMLbase;
    my_iindx    = fc->iindx;
    params      = fc->exp_params;
    md          = &(params->model_details);
    turn        = md->min_loop_size;
    my_iindx    = fc->iindx;
    sc_wrapper  = &(d->sc_wrapper_mb);

    d->qm2_aux = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * ((n + 1) * (n + 2) / 2 + 1));

    for (i = 1; i < n; i++) {
      ii = my_iindx[i];
      for (j = i + 1; j <= n; j++) {
        ij = ii - j;

        /* j is unpaired */
        if (hc->eval_mb(i, j, i, j - 1, VRNA_DECOMP_ML_ML, hc)) {
          temp = d->qm2_aux[ij + 1] *
                 expMLbase[1];

          if (sc_wrapper->red_ml)
            temp *= sc_wrapper->red_ml(i, j, i, j - 1, sc_wrapper);

          d->qm2_aux[ij] = temp;
        }

        /* j pairs with some u */
        for (u = i + 1; u + turn < j; u++) {
          if (hc->eval_mb(i, j, u - 1, u, VRNA_DECOMP_ML_ML_ML, hc)) {
            type = vrna_get_ptype_md(S[u], S[j], md);

            temp = qm[my_iindx[i] - (u - 1)] *
                   qb[my_iindx[u] - j] *
                   vrna_exp_E_multibranch_stem(type,
                                               S1[u - 1],
                                               S1[j + 1],
                                               params);

            if (sc_wrapper->decomp_ml)
              temp *= sc_wrapper->decomp_ml(i, j, u - 1, u, sc_wrapper);

            d->qm2_aux[ij] += temp;
          }
        }
      }
    }
  }
}


PRIVATE void
report_n_clear_accessibility(vrna_fold_compound_t *fc,
                             unsigned int         l,
                             helper_arrays_t      d,
                             void                 *data)
{
  vrna_eqm_probs_f      cb      = d->cb;
  vrna_eqm_probs_data_t cb_data = d->cb_data;
  unsigned int          end     = (l < d->options.up_length) ? l : d->options.up_length;

  /* loop-type specific output */
  if ((d->cb) &&
      (d->options.mode & VRNA_EQM_PROBS_UP)) {
    if (d->options.mode & VRNA_EQM_PROBS_SPLIT) {
      cb(fc, d->access_ext, l, 1, end, VRNA_EQM_PROBS_UP | VRNA_LOOP_TYPE_EXT, cb_data);
      cb(fc, d->access_hp, l, 1, end, VRNA_EQM_PROBS_UP | VRNA_LOOP_TYPE_HP, cb_data);
      cb(fc, d->access_int, l, 1, end, VRNA_EQM_PROBS_UP | VRNA_LOOP_TYPE_INT, cb_data);
      cb(fc, d->access_mb, l, 1, end, VRNA_EQM_PROBS_UP | VRNA_LOOP_TYPE_MB, cb_data);
    }

    /* combined probabilities */
    for (unsigned int u = 1; u <= end; u++)
      d->access_ext[u] += d->access_hp[u] +
                                   d->access_int[u] +
                                   d->access_mb[u];

    cb(fc, d->access_ext, l, 1, end, VRNA_EQM_PROBS_UP | VRNA_LOOP_TYPE_ANY, cb_data);
  }

  /* reset data */
  memset(d->access_ext, 0, sizeof(FLT_OR_DBL) * (d->options.up_length + 1));
  memset(d->access_hp, 0, sizeof(FLT_OR_DBL) * (d->options.up_length + 1));
  memset(d->access_int, 0, sizeof(FLT_OR_DBL) * (d->options.up_length + 1));
  memset(d->access_mb, 0, sizeof(FLT_OR_DBL) * (d->options.up_length + 1));
}


PRIVATE void
compute_up_ext(vrna_fold_compound_t *fc,
               unsigned int         l,
               helper_arrays_t      dp_helpers,
               void                 *data)
{
  unsigned int          i, u, max_u, n;
  FLT_OR_DBL            *q1k, *qln, temp, *acc_ext;
  vrna_mx_pf_t          *matrices;
  vrna_hc_t             *hc;
  FLT_OR_DBL            (*contrib_up)(vrna_fold_compound_t *,
                                      unsigned int,
                                      unsigned int,
                                      helper_arrays_t);

  n         = fc->length;
  matrices  = fc->exp_matrices;
  q1k       = matrices->q1k;
  qln       = matrices->qln;
  hc        = fc->hc;
  max_u     = dp_helpers->options.up_length;
  acc_ext   = dp_helpers->access_ext;

  contrib_up = (fc->type == VRNA_FC_TYPE_SINGLE) ? &contrib_ext_up : &contrib_ext_up_comparative;

  /* compute unpaired probabilities for stretches of size u ending at position l */
  if (l < n) {
    for (u = 1; u <= max_u; u++) {
      i = l + 1 - u;
      if (i > 1) {
        if (hc->eval_ext(1, l, 1, i - 1, VRNA_DECOMP_EXT_EXT, hc)) {
          temp = contrib_up(fc, i, l, dp_helpers) *
                 q1k[i - 1] *
                 qln[l + 1] /
                 q1k[n];

          acc_ext[u] += temp;
        }
      } else {
        if (hc->eval_ext(1, l, 1, l, VRNA_DECOMP_EXT_UP, hc)) {
          temp = contrib_up(fc, 1, l, dp_helpers) *
                 qln[l + 1] /
                 q1k[n];

          acc_ext[u] += temp;
        }
        break;
      }
    }
  } else {
    for (u = 1; u <= max_u; u++) {
      i = l + 1 - u;
      if (i > 1) {
        if (hc->eval_ext(1, l, 1, i - 1, VRNA_DECOMP_EXT_EXT, hc)) {
          temp = contrib_up(fc, i, l, dp_helpers) *
                 q1k[i - 1] /
                 q1k[n];

          acc_ext[u] += temp;
        }
      } else {
        if (hc->eval_ext(1, l, 1, l, VRNA_DECOMP_EXT_UP, hc)) {
          temp = contrib_up(fc, 1, l, dp_helpers) /
                 q1k[n];

          acc_ext[u] += temp;
        }
        break;
      }
    }
  }
}


PRIVATE FLT_OR_DBL
contrib_ext_up(vrna_fold_compound_t *fc,
                 unsigned int         i,
                 unsigned int         j,
                 helper_arrays_t      dp_helpers)
{
  unsigned int      n;
  FLT_OR_DBL        contribution, *scale;
  vrna_sc_t         *sc;

  n             = fc->length;
  scale         = fc->exp_matrices->scale;
  sc            = fc->sc;
  contribution  = scale[j - i + 1];

  if ((sc) &&
      (sc->exp_f)) {
    if (i > 1) {
      if (j < n) {
        contribution *= sc->exp_f(1, n, i - 1, j + 1, VRNA_DECOMP_EXT_EXT_EXT, sc->data);
      } else {
        contribution *= sc->exp_f(1, n, 1, i - 1, VRNA_DECOMP_EXT_EXT, sc->data);
      }
    } else if (j < n) {
      contribution *= sc->exp_f(1, n, j + 1, n, VRNA_DECOMP_EXT_EXT, sc->data);
    } else {
      contribution *= sc->exp_f(1, n, 1, n, VRNA_DECOMP_EXT_UP, sc->data);
    }
  }

  return contribution;
}


PRIVATE FLT_OR_DBL
contrib_ext_up_comparative(vrna_fold_compound_t *fc,
                             unsigned int         i,
                             unsigned int         j,
                             helper_arrays_t      dp_helpers)
{
  unsigned int      n, s, n_seq;
  FLT_OR_DBL        contribution, *scale;
  vrna_sc_t         **scs;

  n         = fc->length;
  n_seq     = fc->n_seq;
  scale     = fc->exp_matrices->scale;
  scs       = fc->scs;

  contribution = scale[j - i + 1];

  if (scs) {
    if (i > 1) {
      if (j < n) {
        for (s = 0; s < n_seq; s++)
          if (scs[s]->exp_f)
            contribution *= scs[s]->exp_f(1, n, i - 1, j + 1, VRNA_DECOMP_EXT_EXT_EXT, scs[s]->data);
      } else {
        for (s = 0; s < n_seq; s++)
          if (scs[s]->exp_f)
            contribution *= scs[s]->exp_f(1, n, 1, i - 1, VRNA_DECOMP_EXT_EXT, scs[s]->data);
      }
    } else if (j < n) {
      for (s = 0; s < n_seq; s++)
        if (scs[s]->exp_f)
          contribution *= scs[s]->exp_f(1, n, j + 1, n, VRNA_DECOMP_EXT_EXT, scs[s]->data);
    } else {
      for (s = 0; s < n_seq; s++)
        if (scs[s]->exp_f)
          contribution *= scs[s]->exp_f(1, n, 1, n, VRNA_DECOMP_EXT_UP, scs[s]->data);
    }
  }

  return contribution;
}


PRIVATE void
compute_up_hp(vrna_fold_compound_t *fc,
              unsigned int         l,
              helper_arrays_t      dp_helpers,
              void                 *data)
{
  char                  **Ss;
  unsigned int          **a2s;
  short                 *S, *S2, **SS, **S5, **S3;
  int                   u, type, n_seq, s, max_u, ij;
  FLT_OR_DBL            q, qbt1, *probs, *scale, *puh, *qbh;
  vrna_exp_param_t      *P;
  vrna_md_t             *md;
  vrna_ud_t             *domains_up;
  unsigned int          i, j, n;
  int                   *my_iindx;
  vrna_hc_t             *hc;
  struct sc_hp_exp_dat  *sc_wrapper_hp;

  n           = fc->length;
  P           = fc->exp_params;
  md          = &(P->model_details);
  scale       = fc->exp_matrices->scale;
  domains_up  = fc->domains_up;
  probs       = fc->exp_matrices->probs;
  my_iindx  = fc->iindx;
  max_u     = dp_helpers->options.up_length;
  puh       = dp_helpers->access_hp;
  qbh       = dp_helpers->qbh;
  hc        = fc->hc;

  sc_wrapper_hp = &(dp_helpers->sc_wrapper_hp);

  if (l == n)
    return;

  j = l + 1;
  for (i = 1; i < l; i++) {
    ij = my_iindx[i] - j;
  
    if (probs[ij] < FLT_MIN)
      continue;

    if (hc->eval_hp(i, j, hc)) {
      q = 0.;
      switch (fc->type) {
        case VRNA_FC_TYPE_SINGLE:
          S     = fc->sequence_encoding;
          S2    = fc->sequence_encoding2;
          u     = j - i - 1;
          type  = vrna_get_ptype_md(S2[i], S2[j], md);

          q = exp_E_Hairpin(u, type, S[i + 1], S[j - 1], fc->sequence + i - 1, P);
          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          SS    = fc->S;
          S5    = fc->S5;   /* S5[s][i] holds next base 5' of i in sequence s */
          S3    = fc->S3;   /* Sl[s][i] holds next base 3' of i in sequence s */
          Ss    = fc->Ss;
          a2s   = fc->a2s;
          n_seq = fc->n_seq;
          qbt1  = 1.;

          for (s = 0; s < n_seq; s++) {
            u = a2s[s][j - 1] - a2s[s][i];
            if (a2s[s][i] < 1)
              continue;

            type  = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
            qbt1  *= exp_E_Hairpin(u, type, S3[s][i], S5[s][j], Ss[s] + a2s[s][i] - 1, P);
          }

          q = qbt1;
          break;

        default:
          break;
      }

      /* add soft constraints */
      if (sc_wrapper_hp->pair)
        q *= sc_wrapper_hp->pair(i, j, sc_wrapper_hp);

      if (domains_up && domains_up->exp_energy_cb) {
        /* we always consider both, bound and unbound state */
        q += q * domains_up->exp_energy_cb(fc,
                                           i + 1, j - 1,
                                           VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
                                           domains_up->data);
      }

      q *= scale[j - i + 1];
      q *= probs[ij];
      qbh[i] += q;
    }
  }

  /* now, for all unpaired stretches of size u ending at l */
  for (u = 1; u <= MIN2(max_u, l - 1); u++)
    for (i = 1; i < l + 1 - u; i++)
      puh[u] += qbh[i];
}


PRIVATE void
compute_up_int(vrna_fold_compound_t *fc,
               unsigned int         l,
               helper_arrays_t      dp_helpers,
               void                 *data)
{
  short                 *S, *S1;
  unsigned int          n, u, max_u, k, i, j, p, u1, u2, turn, type, type2, with_ud;
  int                   ij, *my_iindx;
  FLT_OR_DBL            *pui, *probs, *scale, *qb, tmp;
  vrna_exp_param_t      *params;
  vrna_md_t             *md;
  vrna_hc_t             *hc;
  vrna_ud_t             *domains_up;
  struct sc_int_exp_dat *sc_wrapper_int;


  n           = fc->length;
  S           = fc->sequence_encoding2;
  S1          = fc->sequence_encoding;
  my_iindx    = fc->iindx;
  qb          = fc->exp_matrices->qb;
  probs       = fc->exp_matrices->probs;
  scale       = fc->exp_matrices->scale;
  max_u     = dp_helpers->options.up_length;
  pui       = dp_helpers->access_int;
  params    = fc->exp_params;
  md        = &(params->model_details);
  turn      = md->min_loop_size;
  domains_up  = fc->domains_up;
  with_ud     = (domains_up && domains_up->exp_energy_cb) ? 1 : 0;
  hc        = fc->hc;
  sc_wrapper_int  = &(dp_helpers->sc_wrapper_int);

  if (l == n)
    return;

  FLT_OR_DBL *q3 = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
  FLT_OR_DBL *q5 = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));

  for (k = l - 1; k + MAXLOOP >= l; k--) {
    if (k == 0)
      break;

    for (j = l + 1; k + MAXLOOP + 1 >= j; j++) {
      if (j > n)
        break;

      u2 = j - k - 1;
      for (p = 2; p + turn < k; p++) {
        type2 = vrna_get_ptype_md(S[k], S[p], md);

        for (i = p - 1; i + MAXLOOP + 1 >= p + u2; i--) {
          if (i == 0)
            break;

          ij = my_iindx[i] - j;

          if (probs[ij] == 0.)
            continue;


          if (hc->eval_int(i, j, p, k, hc)) {
            u1   = p - i - 1;
            type = vrna_get_ptype_md(S[i], S[j], md);
            tmp  = probs[ij] *
                   vrna_exp_E_internal(u1,
                                  u2,
                                  type,
                                  type2,
                                  S1[i + 1],
                                  S1[j - 1],
                                  S1[p - 1],
                                  S1[k + 1],
                                  params) *
                    scale[u1 + u2 + 2];

            if (sc_wrapper_int->pair)
              tmp *= sc_wrapper_int->pair(i, j, p, k, sc_wrapper_int);

            if (with_ud) {
              /* here some unstructured domains stuff */
            }

            q3[k] += tmp *
                     qb[my_iindx[p] - k];
          } 
        }
      }
    }
  }

  for (i = l - 1; i + MAXLOOP >= l; i--) {
    if (i == 0)
      break;

    for (j = l + turn + 3; j <= n; j++) {
      ij = my_iindx[i] - j;

      if (probs[ij] < FLT_MIN)
        continue;

      type = vrna_get_ptype_md(S[i], S[j], md);

      for (k = j - 1; k + MAXLOOP >= j; k--) {
        if (l + turn + 2 > k)
          break;

        u2 = j - k - 1;

        for (p = l + 1; p + turn < k; p++) {
          u1 = p - i - 1;

          if (u1 + u2 > MAXLOOP)
            break;

          if (hc->eval_int(i, j, p, k, hc)) {
            type2 = vrna_get_ptype_md(S[k], S[p], md);
            tmp  = probs[ij] *
                   vrna_exp_E_internal(u1,
                                  u2,
                                  type,
                                  type2,
                                  S1[i + 1],
                                  S1[j - 1],
                                  S1[p - 1],
                                  S1[k + 1],
                                  params) *
                    scale[u1 + u2 + 2];

            if (sc_wrapper_int->pair)
              tmp *= sc_wrapper_int->pair(i, j, p, k, sc_wrapper_int);

            if (with_ud) {
              /* here some unstructured domains stuff */
            }

            q5[i] += tmp *
                     qb[my_iindx[p] - k];
          } 
        }
      }
    }  
  }

  for (u = 1; u <= max_u; u++) {
    unsigned int mink = (l > MAXLOOP) ? l - MAXLOOP - 1: 1;
    unsigned int maxk = (l > u) ? l - u : 0;

    for (k = mink; k <= maxk; k++)
      pui[u] += q3[k];

    unsigned int mini = (l > MAXLOOP) ? l - MAXLOOP - 1 : 1;
    unsigned int maxi = (l > u) ? l - u : 0;
    for (i = mini; i <= maxi; i++)
      pui[u] += q5[i];
      
  }
   
  free(q3);
  free(q5);
}


PRIVATE void
compute_up_mb(vrna_fold_compound_t *fc,
              unsigned int         l,
              helper_arrays_t      dp_helpers,
              void                 *data)
{
  short                 *S, *S1;
  unsigned int          n, u, max_u, i, j, turn, start, type;
  int                   *my_iindx, ij;
  FLT_OR_DBL            *pum, *qm, *qm2, *probs, *X, *Y, *Z, temp, up, qqm, *expMLbase, expMLclosing, *scale;
  vrna_exp_param_t      *params;
  vrna_md_t             *md;
  vrna_hc_t             *hc;
  struct sc_mb_exp_dat  *sc_wrapper;

  n         = fc->length;
  S         = fc->sequence_encoding2;
  S1        = fc->sequence_encoding;
  max_u     = dp_helpers->options.up_length;
  pum       = dp_helpers->access_mb;
  X         = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
  Y         = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
  Z         = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
  probs     = fc->exp_matrices->probs;
  qm        = fc->exp_matrices->qm;
  qm2           = (fc->exp_matrices->qm2_real) ? fc->exp_matrices->qm2_real : dp_helpers->qm2_aux;
  params    = fc->exp_params;
  md        = &(params->model_details);
  turn          = md->min_loop_size;
  expMLbase = fc->exp_matrices->expMLbase;
  scale         = fc->exp_matrices->scale;
  expMLclosing  = params->expMLclosing;
  my_iindx  = fc->iindx;
  hc        = fc->hc;
  sc_wrapper  = &(dp_helpers->sc_wrapper_mb);

  if (l == n)
    return;

  /* fill X, Y, Z */
  for (i = 1; i < l; i++) {
    X[i] = Y[i] = Z[i] = 0.;

    for (j = l + 1; j <= n; j++) {
      if (hc->eval_mb(i + 1, j - 1, i + 1, l, VRNA_DECOMP_ML_ML, hc)) {
        ij    = my_iindx[i] - j;
        type  = vrna_get_ptype_md(S[j], S[i], md);
        temp  = probs[ij] *
                expMLbase[j - l - 1] *
                vrna_exp_E_multibranch_stem(type, S1[j - 1], S1[i + 1], params);

        if (sc_wrapper->red_ml)
          temp *= sc_wrapper->red_ml(i + 1, j - 1, i + 1, l, sc_wrapper);

        X[i] += temp;
      }
    }

    for (j = l + turn + 1; j <= n; j++) {
      if (hc->eval_mb(i + 1, j - 1, l, l + 1, VRNA_DECOMP_ML_ML_ML, hc)) {
        ij    = my_iindx[i] - j;
        type  = vrna_get_ptype_md(S[j], S[i], md);
        temp  = probs[ij] *
                qm[my_iindx[l + 1] - j + 1] *
                vrna_exp_E_multibranch_stem(type, S1[j - 1], S1[i + 1], params);

        if (sc_wrapper->decomp_ml)
          temp *= sc_wrapper->decomp_ml(i + 1, j - 1, l, l + 1, sc_wrapper);

        Y[i] += temp;
      }
    }

    for (j = l + 2 * turn + 1; j <= n; j++) {
      if (hc->eval_mb(i + 1, j - 1, l + 1, j - 1, VRNA_DECOMP_ML_ML, hc)) {
        ij    = my_iindx[i] - j;
        type  = vrna_get_ptype_md(S[j], S[i], md);
        temp  = probs[ij] *
                qm2[my_iindx[l + 1] - j + 1] *
                vrna_exp_E_multibranch_stem(type, S1[j - 1], S1[i + 1], params);

        if (sc_wrapper->red_ml)
          temp *= sc_wrapper->red_ml(i + 1, j - 1, l + 1, j - 1, sc_wrapper);

        Z[i] += temp;
      }
    }
  }

  for (u = 1; u <= max_u; u++) {
    /* process unpaired stretches of size u ending at l */
    if (l >= u) {
      start = l + 1 - u;
      temp = 0.;
      up = expMLbase[u];

      for (i = 1; i + turn < start; i++) {
        if (hc->eval_mb(i + 1, l, i + 1, start - 1, VRNA_DECOMP_ML_ML, hc)) {
          /* unpaired stretch is somewhere in the middle */
          qqm = qm[my_iindx[i + 1] - start + 1] *
                Y[i] *
                up;

          if (sc_wrapper->red_ml)
            qqm *= sc_wrapper->red_ml(i + 1, l, i + 1, start - 1, sc_wrapper);

          temp += qqm;

          /* unpaired stretch is right-most */
          qqm = qm2[my_iindx[i + 1] - start + 1] *
                X[i] *
                up;

          if (sc_wrapper->red_ml)
            qqm *= sc_wrapper->red_ml(i + 1, l, i + 1, start - 1, sc_wrapper);

          temp += qqm;
        }
      }

      for (i = 1; i < start; i++) {
        /* unpaired stretch is left-most */
        qqm = expMLbase[l - i] *
              Z[i];

        if (sc_wrapper->red_ml)
          qqm *= sc_wrapper->red_ml(i + 1, l, start, 1, sc_wrapper);

        temp += qqm;
      }

      pum[u] += temp *
                scale[2] *
                expMLclosing;
    }
  }

  free(X);
  free(Y);
  free(Z);
}
