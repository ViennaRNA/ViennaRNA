
PRIVATE int
outside_bp_circ(vrna_fold_compound_t  *fc,
                unsigned int          l,
                FLT_OR_DBL            *Qmax,
                int                   *ov,
                helper_arrays_t       dp_helpers)
{
  char                      *ptype;
  unsigned char             *hard_constraints, eval, with_gquad;
  short                     *S, *S1, **SS, **S5, **S3, si, sj;
  unsigned int              s, n_seq, type, rt, *tt, **a2s, turn, u1, u2, u3, us1, us2, us3,
                            imin, imax, jmin, jmax, kmin, p, q, n, i, j, k;

  int                       kl, *rtype, *my_iindx, *jindx;
  FLT_OR_DBL                tmp, tmp2, tmp3, expMLclosing, *qb, *qm, *probs, *scale, *expMLbase, qo,
                            *qm2_real, q_g, qbt1;
  double                    *expintern;
  vrna_hc_t                 *hc;
  vrna_exp_param_t          *pf_params;
  vrna_mx_pf_t              *matrices;
  vrna_md_t                 *md;

  vrna_smx_csr(FLT_OR_DBL)  *q_gq = fc->exp_matrices->q_gq;

  FLT_OR_DBL                (*numerator_f)(vrna_fold_compound_t *fc,
                                           int                  i,
                                           int                  j);
  struct sc_mb_exp_dat  *sc_mb_wrapper;
  struct sc_hp_exp_dat  *sc_hp_wrapper;
  struct sc_ext_exp_dat *sc_ext_wrapper;
  struct sc_int_exp_dat *sc_int_wrapper;


  n                 = fc->length;
  n_seq             = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  SS                = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5                = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3                = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  a2s               = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->a2s;
  pf_params         = fc->exp_params;
  md                = &(pf_params->model_details);
  S                 = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding2 : NULL;
  S1                = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  my_iindx          = fc->iindx;
  jindx             = fc->jindx;
  ptype             = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->ptype : NULL;
  turn              = md->min_loop_size;
  type              = 0;
  rt                = 0;
  tt                = NULL;
  hc                = fc->hc;
  matrices          = fc->exp_matrices;
  qb                = matrices->qb;
  qm                = matrices->qm;
  qm2_real          = matrices->qm2_real;
  probs             = matrices->probs;
  scale             = matrices->scale;
  expMLbase         = matrices->expMLbase;
  qo                = matrices->qo;
  hard_constraints  = hc->mx;

  expMLclosing  = pf_params->expMLclosing;
  expintern     = &(pf_params->expinternal[0]);
  rtype         = &(pf_params->model_details.rtype[0]);
  with_gquad    = (unsigned char)md->gquad;

  sc_mb_wrapper = &(dp_helpers->sc_wrapper_mb);
  sc_hp_wrapper = &(dp_helpers->sc_wrapper_hp);
  sc_ext_wrapper = &(dp_helpers->sc_wrapper_ext);
  sc_int_wrapper = &(dp_helpers->sc_wrapper_int);

  switch (fc->type) {
    case  VRNA_FC_TYPE_SINGLE:
      numerator_f = numerator_single;
      tt          = NULL;
      break;
    case  VRNA_FC_TYPE_COMPARATIVE:
      numerator_f = numerator_comparative;
      tt          = (unsigned int *)vrna_alloc(sizeof(unsigned int) * n_seq);
      break;
    default:
      numerator_f = NULL;
      break;
  }

  /* reset all probabilities */
  probs -= l;

  for (k = 1; k <= l; k++)
    probs[my_iindx[k]] = 0;

  probs += l;


  /* 1. exterior pair (k,l) */
  for (k = 1; k + turn < l; k++) {
    kl = my_iindx[k] - l;

    if (qb[kl] > 0.) {
      probs[kl] = numerator_f(fc, k, l) / qo;

      if (fc->type == VRNA_FC_TYPE_SINGLE) {
        type  = vrna_get_ptype_md(S[k], S[l], md);
        rt    = vrna_get_ptype_md(S[l], S[k], md);
      } else {
        for (s = 0; s < n_seq; s++)
          tt[s] = vrna_get_ptype_md(SS[s][l], SS[s][k], md);
      }

      /* 1.1. Exterior Hairpin Contribution */
      tmp2 = vrna_exp_eval_hairpin(fc, l, k, VRNA_EVAL_LOOP_DEFAULT);

      /* 1.2. Exterior Interior Loop Contribution */
      if (hard_constraints[k * n + l] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
        /*
         * 1.2.1. (k,l)  delimts the "left" part of the internal loop
         * (l,k) is "outer pair", (p,q) before (k,l)
         */
        if ((turn + 2 < k) &&     /* does another base pair (p,q) fit before k? */
            (l + MAXLOOP + 1 > n)) {
          unsigned int ln3    = n - l;
          unsigned int max_p  = 1 + MAXLOOP - ln3;
          if (max_p + turn + 2 >= k)
            max_p = k - turn - 2;

          for (p = 1; p <= max_p; p++) {
            unsigned int ln1 = p - 1;
          
            if (hc->up_int[l + 1] < (ln1 + ln3))
              break;

            unsigned int min_q = p + turn + 1;
            if ((ln1 + ln3 + (k - min_q - 1)) > MAXLOOP)
              min_q = k - 1 - (MAXLOOP - ln1 - ln3);

            for (q = min_q; q < k; q++) {
              unsigned int type_2, ln1a, ln2a, ln2;
              ln2 = k - q - 1;
            
              if (hc->up_int[q + 1] < ln2)
                continue;

              if (qb[my_iindx[p] - q] == 0.)
                continue;

              eval = (hard_constraints[p * n + q] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) ? 1 : 0;
              if (hc->f)
                eval = hc->f(p, q, k, l, VRNA_DECOMP_PAIR_IL, hc->data);

              if (eval) {
                tmp = qb[my_iindx[p] - q];

                if (fc->type == VRNA_FC_TYPE_SINGLE) {
                  type_2 = vrna_get_ptype(jindx[q] + p, ptype);

                  tmp *= vrna_exp_E_internal(ln1 + ln3,
                                       ln2,
                                       rt,
                                       rtype[type_2],
                                       S1[l + 1],
                                       S1[k - 1],
                                       S1[p - 1],
                                       S1[q + 1],
                                       pf_params);
                } else {
                  for (s = 0; s < n_seq; s++) {
                    ln2a    = a2s[s][k - 1];
                    ln2a    -= a2s[s][q];
                    ln1a    = a2s[s][n] - a2s[s][l];
                    ln1a    += a2s[s][p - 1];
                    type_2  = vrna_get_ptype_md(SS[s][q], SS[s][p], md);
                    tmp     *= vrna_exp_E_internal(ln1a, ln2a, tt[s], type_2,
                                             S3[s][l],
                                             S5[s][k],
                                             S5[s][p],
                                             S3[s][q],
                                             pf_params);
                  }
                }

                if (sc_int_wrapper->pair_ext)
                  tmp *= sc_int_wrapper->pair_ext(p, q, k, l, sc_int_wrapper);

                tmp2 += tmp *
                        scale[ln1 + ln2 + ln3];
              }
            }
          }
        } /* end (p,q) before (k,l) internal loop */

        /*
         * 1.2.2. i,j  delimtis the "right" part of the internal loop
         * (p,q) after (k,l)
         */
        if ((l + turn + 2 <= n) &&  /* does another base pair (p,q) fit after l? */
            ((k - 1) <= MAXLOOP)) {
          unsigned int ln2    = k - 1;
          unsigned int max_p  = l + 1 + (MAXLOOP - ln2);

          if (max_p + turn + 1 > n)
            max_p = n - turn - 1;

          for (p = l + 1; p <= max_p; p++) {
            unsigned int ln1 = p - l - 1;

            if (hc->up_int[l + 1] < ln1)
              break;

            unsigned int min_q = p + turn + 1;

            if ((ln1 + ln2 + (n - min_q)) > MAXLOOP)
              min_q = n + ln1 + ln2 - MAXLOOP;

            for (q = min_q; q <= n; q++) {
              unsigned int ln3 = n - q;
              unsigned int ln2a, ln1a, type_2;

              if (hc->up_int[q + 1] < (ln2 + ln3))
                continue;

              if (qb[my_iindx[p] - q] == 0.)
                continue;

              eval = (hard_constraints[p * n + q] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) ? 1 : 0;
              if (hc->f)
                eval = hc->f(k, l, p, q, VRNA_DECOMP_PAIR_IL, hc->data) ? eval : 0;

              if (eval) {
                tmp = qb[my_iindx[p] - q];

                if (fc->type == VRNA_FC_TYPE_SINGLE) {
                  type_2  = vrna_get_ptype(jindx[q] + p, ptype);
                  tmp     *= vrna_exp_E_internal(ln2 + ln3,
                                           ln1,
                                           rtype[type_2],
                                           rt,
                                           S1[q + 1],
                                           S1[p - 1],
                                           S1[k - 1],
                                           S1[l + 1],
                                           pf_params);
                } else {
                  for (s = 0; s < n_seq; s++) {
                    ln1a    = a2s[s][p] - a2s[s][l + 1];
                    ln2a    = a2s[s][k - 1] + a2s[s][n] - a2s[s][q];
                    type_2  = vrna_get_ptype_md(SS[s][q], SS[s][p], md);
                    tmp     *= vrna_exp_E_internal(ln2a, ln1a, type_2, tt[s],
                                             S3[s][q],
                                             S5[s][p],
                                             S5[s][k],
                                             S3[s][l],
                                             pf_params);
                  }
                }

                if (sc_int_wrapper->pair_ext)
                  tmp *= sc_int_wrapper->pair_ext(k, l, p, q, sc_int_wrapper);

                tmp2 += tmp *
                        scale[ln1 + ln2 + ln3];
              }
            }
          }
        } /* end (p,q) after (k,l) internal loop */

        /* 1.3 Exterior multiloop decomposition */
        if (hc->eval_mb(k, l, k - 1, l + 1, VRNA_DECOMP_PAIR_ML_EXT, hc)) {
          FLT_OR_DBL sc_contrib = 1.;

          if (sc_mb_wrapper->pair_ext)
            sc_contrib = sc_mb_wrapper->pair_ext(k, l, sc_mb_wrapper);

          /* 1.3.1 Middle part                    */
          if ((k > 2) &&
              (l < n - 1)) {
            tmp = qm[my_iindx[1] - k + 1] *
                  qm[my_iindx[l + 1] - n];

            if (fc->type == VRNA_FC_TYPE_SINGLE) {
              tmp *= vrna_exp_E_multibranch_stem(type,
                                  S1[k - 1],
                                  S1[l + 1],
                                  pf_params) *
                     expMLclosing;
            } else {
              for (s = 0; s < n_seq; s++)
                tmp *= vrna_exp_E_multibranch_stem(rtype[tt[s]],
                                    S5[s][k],
                                    S3[s][l],
                                    pf_params);

              tmp *= pow(expMLclosing, n_seq);
            }

            tmp2 += tmp *
                    sc_contrib;
          }

          /* 1.3.2 right-most part  */
          if (hc->up_ml[l + 1] >= (n - l)) {
            if (k > 1) {
              if ((hc->eval_mb(k, n, k, l, VRNA_DECOMP_ML_ML, hc)) &&
                  (hc->eval_mb(1, l, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc))) {
                tmp = qm2_real[my_iindx[1] - k + 1] *
                      expMLbase[n - l];

                if (fc->type == VRNA_FC_TYPE_SINGLE) {
                  tmp *= vrna_exp_E_multibranch_stem(type,
                                      S1[k - 1],
                                      S1[l + 1],
                                      pf_params) *
                         expMLclosing;
                } else {
                  for (s = 0; s < n_seq; s++)
                    tmp *= vrna_exp_E_multibranch_stem(rtype[tt[s]],
                                        S5[s][k],
                                        S3[s][l],
                                        pf_params);

                  tmp *= pow(expMLclosing, n_seq);
                }

                if (sc_mb_wrapper->red_ml)
                  tmp *= sc_mb_wrapper->red_ml(k, n, k, l, sc_mb_wrapper);

                if (sc_mb_wrapper->decomp_ml)
                  tmp *= sc_mb_wrapper->decomp_ml(1, l, k - 1, k, sc_mb_wrapper);

                tmp2 += tmp *
                        sc_contrib;

              }
            }
          }

          /* 1.3.3 left-most part */
          if (hc->up_ml[1] >= (k - 1)) {
            if (l + 1 < n) {
              if ((hc->eval_mb(1, l, k, l, VRNA_DECOMP_ML_ML, hc)) &&
                  (hc->eval_mb(k, n, l, l + 1, VRNA_DECOMP_ML_ML_ML, hc))) {
                tmp = qm2_real[my_iindx[l + 1] - n] *
                      expMLbase[k - 1];

                if (fc->type == VRNA_FC_TYPE_SINGLE) {
                  tmp *= vrna_exp_E_multibranch_stem(type,
                                      S1[k - 1],
                                      S1[l + 1],
                                      pf_params) *
                         expMLclosing;
                } else {
                  for (s = 0; s < n_seq; s++)
                    tmp *= vrna_exp_E_multibranch_stem(rtype[tt[s]],
                                        S5[s][k],
                                        S3[s][l],
                                        pf_params);

                  tmp *= pow(expMLclosing, n_seq);
                }

                if (sc_mb_wrapper->red_ml)
                  tmp *= sc_mb_wrapper->red_ml(1, l, k, l, sc_mb_wrapper);

                if (sc_mb_wrapper->decomp_ml)
                  tmp *= sc_mb_wrapper->decomp_ml(k, n, l, l + 1, sc_mb_wrapper);

                tmp2 += tmp *
                        sc_contrib;
              }
            }
          }
        } /* end multibranch loop part */

        /* all exterior loop decompositions for pair (k,l) done  */

        if (with_gquad) {
          unsigned int u1, u2, u3,  us1, us2, us3;
          /* consider all the other possibilities where at least one gquad is enclosed in the same loop as (k,l) */
          /* 1. internal-loop like structure */
          if (hard_constraints[k * n + l] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
            /* 1.1 gquad is 5' of base pair (k,l) */

            FLT_OR_DBL qint = 1.;
            short sk, sl;

            switch (fc->type) {
              case VRNA_FC_TYPE_SINGLE:
                sl = S1[l + 1];
                sk = S1[k - 1];

                if (md->dangles == 2)
                  qint *= pf_params->expmismatchI[rt][sl][sk];

                if (rt > 2)
                  qint *= pf_params->expTermAU;

                break;

              case VRNA_FC_TYPE_COMPARATIVE:
                for (s = 0; s < n_seq; s++) {
                  if (md->dangles == 2)
                    qint *= pf_params->expmismatchI[tt[s]][S3[s][l]][S5[s][k]];

                  if (tt[s] > 2)
                    qint *= pf_params->expTermAU;
                }
                break;
            }

            /* int loop with [p,q] gquad followed by (k,l) base pair */
            if (n - l <= MAXLOOP) {
              unsigned int u1     = n - l;
              unsigned int max_p  = 1 + MAXLOOP - u1;
              for (p = 1; p <= max_p; p++) {
                u2 = p - 1;
                unsigned int max_q = k - 1;
                if (p + VRNA_GQUAD_MAX_BOX_SIZE - 1 < max_q)
                  max_q = p + VRNA_GQUAD_MAX_BOX_SIZE - 1;

                unsigned int min_q = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
                if (min_q >= k)
                  min_q = k;
                else if ((k - min_q) + u1 + u2 > MAXLOOP + 1)
                  min_q = k - (MAXLOOP - u1 - u2) - 1;
                   
                for (q = min_q; q <= max_q; q++) {
                  u3 = k - q - 1;
                  if (((u1 + u2 == 0) && (u3 < 3)) ||
                      ((u1 + u2 < 3) && (u3 == 0)))
                    continue;

#ifndef VRNA_DISABLE_C11_FEATURES
                  q_g = vrna_smx_csr_get(q_gq, p, q, 0.);
#else
                  q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, p, q, 0.);
#endif
                  if (q_g != 0.) {
                    tmp = scale[u1 + u2 + u3];

                    switch (fc->type) {
                      case VRNA_FC_TYPE_SINGLE:
                        tmp *= (FLT_OR_DBL)expintern[u1 + u2 + u3];
                        break;

                      case VRNA_FC_TYPE_COMPARATIVE:
                        for (s = 0; s < n_seq; s++) {
                          us1   = (p > 1) ? a2s[s][p - 1] - a2s[s][1] : 0;
                          us2   = a2s[s][k - 1] - a2s[s][q];
                          us3   = a2s[s][n] - a2s[s][l];
                          tmp  *= (FLT_OR_DBL)expintern[us1 + us2 + us3];
                        }
                        break;
                    }

                    if (sc_int_wrapper->pair_ext)
                      tmp *= sc_int_wrapper->pair_ext(p, q, k, l, sc_int_wrapper);

                    /* store outside part for base pair */
                    tmp2 += qint *
                            q_g *
                            tmp;
                    /* store outside part for gquad */
                    probs[my_iindx[p] - q] += qint * tmp * qb[my_iindx[k] - l] / qo;
                  }
                }
              }
            }

            /* 1.2 gquad is 3' of base pair (k,l) */
            u1 = k - 1;
            unsigned int max_p = l + MAXLOOP - u1 + 1;
            if (l + VRNA_GQUAD_MIN_BOX_SIZE <= n)
              max_p = MIN2(max_p, n - VRNA_GQUAD_MIN_BOX_SIZE + 1);
            else
              max_p = l;

            for (p = l + 1; p <= max_p; p++) {
              u2 = p - l - 1;
              unsigned int min_q = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
              if (min_q > n)
                break;
              if (min_q + MAXLOOP - u1 - u2 < n)
                min_q = n + u1 + u2 - MAXLOOP;

              unsigned int max_q = p + VRNA_GQUAD_MAX_BOX_SIZE - 1;
              if (max_q > n)
                max_q = n;

              for (q = min_q; q <= max_q; q++) {
                u3  = n - q;
                if (((u1 + u3 == 0) && (u2 < 3)) ||
                    ((u1 + u3 < 3) && (u2 == 0)))
                  continue;
              
#ifndef VRNA_DISABLE_C11_FEATURES
                q_g = vrna_smx_csr_get(q_gq, p, q, 0.);
#else
                q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, p, q, 0.);
#endif
                if (q_g != 0.) {
                  tmp = scale[u1 + u2 + u3];

                  switch (fc->type) {
                    case VRNA_FC_TYPE_SINGLE:
                      tmp *= (FLT_OR_DBL)expintern[u1 + u2 + u3];
                      break;

                    case VRNA_FC_TYPE_COMPARATIVE:
                      for (s = 0; s < n_seq; s++) {
                        us1   = a2s[s][k - 1] - a2s[s][1];
                        us2   = a2s[s][p - 1] - a2s[s][l];
                        us3   = a2s[s][n] - a2s[s][q];
                        tmp  *= (FLT_OR_DBL)expintern[us1 + us2 + us3];
                      }
                      break;
                  }

                  if (sc_int_wrapper->pair_ext)
                    tmp *= sc_int_wrapper->pair_ext(k, l, p, q, sc_int_wrapper);

                  /* store outside part for base pair */
                  tmp2 +=  qint *
                          q_g *
                          tmp;
                  /* store outside part for gquad */
                  probs[my_iindx[p] - q] += tmp * qint * qb[my_iindx[k] - l] / qo;
                }
              }
            }
          }

          if (hard_constraints[k * n + l] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
            /* 3. all base pairs (k,l) where G-Quad spans across n,1 boundary */
            /* 3.2 multibranch-loop like configurations */
            FLT_OR_DBL qmb = 1.;
          
            switch (fc->type) {
              case VRNA_FC_TYPE_SINGLE:
                qmb *=  vrna_exp_E_multibranch_stem(type, S1[k - 1], S1[l + 1], pf_params);
                break;
              case VRNA_FC_TYPE_COMPARATIVE:
                for (s = 0; s < n_seq; s++)
                  qmb *= vrna_exp_E_multibranch_stem(rtype[tt[s]],
                                      S5[s][k],
                                      S3[s][l],
                                      pf_params);
                break;
            }

            qmb *= pow(expMLclosing, (double)n_seq);

            FLT_OR_DBL expGQstem = pow(vrna_exp_E_multibranch_stem(0, -1, -1, pf_params), (double)n_seq);

            unsigned int min_p = l + 1;
            if (min_p + VRNA_GQUAD_MAX_BOX_SIZE - 1 <= n)
              min_p = n - VRNA_GQUAD_MAX_BOX_SIZE + 1;

            for (p = min_p; p <= n; p++) {

              unsigned int min_q = 1;
              unsigned int max_q = VRNA_GQUAD_MAX_BOX_SIZE - 1 - (n - p);
              if (max_q >= k)
                max_q = k - 1;
              for (q = min_q; q <= max_q; q++) {
#ifndef VRNA_DISABLE_C11_FEATURES
                if ((q_g = vrna_smx_csr_get(q_gq, p, q, 0.)) != 0.) {
#else
                if ((q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, p, q, 0.)) != 0.) {
#endif
                  /* 3.2.1 (k,l) last branch before gquad, ie. [p,q] + ML + (k,l) */
                  if (q + turn < k) {
                    u1 = p - l - 1;
                    if (hc->up_ml[l + 1] >= u1) {
                      tmp = qmb *
                            qm[my_iindx[q + 1] - k + 1] *
                            expGQstem *
                            expMLbase[u1];

                      if (sc_mb_wrapper->pair_ext)
                        tmp *= sc_mb_wrapper->pair_ext(k, l, sc_mb_wrapper);

                      /* store outside part for base pair (k,l) */
                      tmp2 += q_g *
                              tmp;
                    }
                  }

                  /* 3.2.2 (k,l) first branch after gquad, ie. [p,q] + (k,l) + ML */
                  if (l + turn < p) {
                    u2 = k - q - 1;
                    if (hc->up_ml[q + 1] >= u2) {
                      tmp = qmb *
                            qm[my_iindx[l + 1] - p + 1] *
                            expGQstem *
                            expMLbase[u2];
                    
                      if (sc_mb_wrapper->pair_ext)
                        tmp *= sc_mb_wrapper->pair_ext(k, l, sc_mb_wrapper);

                      /* store outside part for base pair (k,l) */
                      tmp2 += q_g *
                              tmp;
                    }
                  }            

                  /* 3.2.3 (k,l) somewhere in between, ie. [p,q] + ML + (k,l) + ML */
                  if ((q + turn < k) &&
                      (l + turn < p)) {
                    tmp = qmb *
                          qm[my_iindx[q + 1] - k + 1] *
                          qm[my_iindx[l + 1] - p + 1] *
                          expGQstem;

                    if (sc_mb_wrapper->pair_ext)
                      tmp *= sc_mb_wrapper->pair_ext(k, l, sc_mb_wrapper);

                      /* store outside part for base pair (k,l) */
                    tmp2 += q_g *
                            tmp;
                  }
                }
              }
            }
          }
        }

        probs[kl] *= tmp2;
      }
    } /* end if qb[kl] > 0 */
  } /* end for k */

  /* Done with outside probabilities for base pairs (k,l) */

  /* Next follows outside probabilities for G-Quadruplexes delimited by [k : l] */
  if (with_gquad) {
    if (VRNA_GQUAD_MIN_BOX_SIZE <= l) {
      unsigned int min_k = 1;
      if (min_k + VRNA_GQUAD_MAX_BOX_SIZE - 1 < l)
        min_k = l + 1 - VRNA_GQUAD_MAX_BOX_SIZE;

      for (k = min_k; k + VRNA_GQUAD_MIN_BOX_SIZE - 1 <= l; k++) {
        kl = my_iindx[k] - l;
      
#ifndef VRNA_DISABLE_C11_FEATURES
        if ((q_g = vrna_smx_csr_get(q_gq, k, l, 0.)) != 0.) {
#else
        if ((q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.)) != 0.) {
#endif
          tmp2 = 0.;

          /* 1. all possible G-quads spanning from k to l */
          
          /* 1.1 hairpin-loop like configurations */
          if ((k - 1) + (n - l) >= 3) {
            eval = (hc->up_ext[1] >= (k - 1)) ? 1 : 0;
            if (l < n)
              eval = (hc->up_ext[l + 1] >= (n - l)) ? eval : 0;
            if (hc->f) {
              if (k > 1)
                eval = (hc->f(1, k - 1, 1, k - 1, VRNA_DECOMP_EXT_UP, hc->data)) ? eval : 0;
              if (l < n)
                eval = (hc->f(l + 1, n, l + 1, n, VRNA_DECOMP_EXT_UP, hc->data)) ? eval : 0;
            }
            
            if (eval) {
              tmp = scale[(k - 1) + (n - l)]; 

              if (md->circ_penalty)
                tmp *= pow(vrna_exp_E_hairpin((k - 1) + (n - l), 0, -1, -1, NULL, pf_params), (double)n_seq);

              if (sc_ext_wrapper->red_up)
                tmp *= sc_ext_wrapper->red_up(1, k - 1, sc_ext_wrapper) *
                       sc_ext_wrapper->red_up(l + 1, n, sc_ext_wrapper);
 
              tmp2 += tmp;
            }
          }

          /* 1.2 internal-loop like configurations */
          /* 1.2.3 [k,l] [p,q] everything else already handled above */
          u1 = k - 1;
          if ((l + VRNA_GQUAD_MIN_BOX_SIZE <= n) &&
              (u1 <= MAXLOOP) &&
              (hc->up_int[1] >= u1)) {
            unsigned int max_p = l + 1 + MAXLOOP - u1;
            if (max_p + VRNA_GQUAD_MIN_BOX_SIZE - 1 > n)
              max_p = n - VRNA_GQUAD_MIN_BOX_SIZE + 1;

            for (p = l + 1; p <= max_p; p++) {
              u2 = p - l - 1;
              if (hc->up_int[l + 1] < u2)
                break;

              unsigned int min_q = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
              if (min_q + MAXLOOP < n + u1 + u2)
                min_q = n + u1 + u2 - MAXLOOP;

              for (u3 = 0, q = n; q >= min_q; q--, u3++) {
                if (hc->up_int[q + 1] < u3)
                  break;

                if ((((u1 + u3) == 0) && (u2 < 3)) ||
                    (((u1 + u3) < 3) && (u2 == 0)))
                  continue;

#ifndef VRNA_DISABLE_C11_FEATURES
                FLT_OR_DBL tmp3 = vrna_smx_csr_get(q_gq, p, q, 0.);
#else
                FLT_OR_DBL tmp3 = vrna_smx_csr_FLT_OR_DBL_get(q_gq, p, q, 0.);
#endif
                
                if (tmp3 != 0.) {
                    tmp = scale[u1 + u2 + u3];

                  switch (fc->type) {
                    case VRNA_FC_TYPE_SINGLE:
                      tmp *= (FLT_OR_DBL)expintern[u1 + u2 + u3];
                      break;

                    case VRNA_FC_TYPE_COMPARATIVE:
                      for (s = 0; s < n_seq; s++) {
                        us1   = (k > 1) ? a2s[s][k - 1] - a2s[s][1] : 0;
                        us2   = a2s[s][p - 1] - a2s[s][l];
                        us3   = a2s[s][n] - a2s[s][q];
                        tmp  *= (FLT_OR_DBL)expintern[us1 + us2 + us3];
                      }
                      break;
                  }

                  /* soft constraints for at least unpaired bases in the internal-loop like conformation?! */

                  tmp2 += tmp * tmp3;
                  probs[my_iindx[p] - q] += tmp * q_g / qo; /* add contribution for 2nd gquad here instead of extra block below */
                }
              }
            }
          }

          /* 1.3 multibranch-loop like configurations with [k,l] gquad */
          FLT_OR_DBL expGQstem = pow(vrna_exp_E_multibranch_stem(0, -1, -1, pf_params) *
                                 expMLclosing, (double)n_seq);

          /* 1.3.1 no element spanning n,1 junction */
          /* 1.3.1.1 Middle part                    */
          if ((k > 2) &&
              (l < n - 1)) {
            tmp = qm[my_iindx[1] - k + 1] *
                  qm[my_iindx[l + 1] - n];


            tmp2 += tmp * expGQstem;
          }

          /* 1.3.1,2 right-most part  */
          if (hc->up_ml[l + 1] >= (n - l)) {
            if (k > 1) {
              if ((hc->eval_mb(k, n, k, l, VRNA_DECOMP_ML_ML, hc)) &&
                  (hc->eval_mb(1, l, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc))) {
                tmp = qm2_real[my_iindx[1] - k + 1] *
                      expMLbase[n - l];

                if (sc_mb_wrapper->red_ml)
                  tmp *= sc_mb_wrapper->red_ml(k, n, k, l, sc_mb_wrapper);

                if (sc_mb_wrapper->decomp_ml)
                  tmp *= sc_mb_wrapper->decomp_ml(1, l, k - 1, k, sc_mb_wrapper);

                tmp2 += tmp *
                        expGQstem;
              }
            }
          }

          /* 1.3.1.3 left-most part */
          if (hc->up_ml[1] >= (k - 1)) {
            if (l + 1 < n) {
              if ((hc->eval_mb(1, l, k, l, VRNA_DECOMP_ML_ML, hc)) &&
                  (hc->eval_mb(k, n, l, l + 1, VRNA_DECOMP_ML_ML_ML, hc))) {
                tmp = qm2_real[my_iindx[l + 1] - n] *
                      expMLbase[k - 1];

                if (sc_mb_wrapper->red_ml)
                  tmp *= sc_mb_wrapper->red_ml(1, l, k, l, sc_mb_wrapper);

                if (sc_mb_wrapper->decomp_ml)
                  tmp *= sc_mb_wrapper->decomp_ml(k, n, l, l + 1, sc_mb_wrapper);

                tmp2 += tmp *
                        expGQstem;
              }
            }
          }

          /* 1.3.2 one element spanning n,1 junction */
          FLT_OR_DBL elem = pow(vrna_exp_E_multibranch_stem(0, -1, -1, pf_params), (double)n_seq);

          unsigned int min_p = l + 1;
          if (min_p + VRNA_GQUAD_MAX_BOX_SIZE - 1 <= n)
            min_p = n - VRNA_GQUAD_MAX_BOX_SIZE + 1;

          for (p = min_p; p <= n; p++) {

            unsigned int min_q = 1;
            unsigned int max_q = VRNA_GQUAD_MAX_BOX_SIZE - 1 - (n - p);
            if (max_q >= k)
              max_q = k - 1;
            for (q = min_q; q <= max_q; q++) {
#ifndef VRNA_DISABLE_C11_FEATURES
              if ((q_g = vrna_smx_csr_get(q_gq, p, q, 0.)) != 0.) {
#else
              if ((q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, p, q, 0.)) != 0.) {
#endif
                /* 3.2.1 (k,l) last branch before gquad, ie. [p,q] + ML + (k,l) */
                if (q + turn < k) {
                  u1 = p - l - 1;
                  if (hc->up_ml[l + 1] >= u1) {
                    tmp = expGQstem *
                          qm[my_iindx[q + 1] - k + 1] *
                          elem *
                          expMLbase[u1];

                    if (sc_mb_wrapper->pair_ext)
                      tmp *= sc_mb_wrapper->pair_ext(k, l, sc_mb_wrapper);

                    /* store outside part for base pair (k,l) */
                    tmp2 += q_g *
                            tmp;
                  }
                }

                /* 3.2.2 (k,l) first branch after gquad, ie. [p,q] + (k,l) + ML */
                if (l + turn < p) {
                  u2 = k - q - 1;
                  if (hc->up_ml[q + 1] >= u2) {
                    tmp = expGQstem *
                          qm[my_iindx[l + 1] - p + 1] *
                          elem *
                          expMLbase[u2];
                  
                    if (sc_mb_wrapper->pair_ext)
                      tmp *= sc_mb_wrapper->pair_ext(k, l, sc_mb_wrapper);

                    /* store outside part for base pair (k,l) */
                    tmp2 += q_g *
                            tmp;
                  }
                }            

                /* 3.2.3 (k,l) somewhere in between, ie. [p,q] + ML + (k,l) + ML */
                if ((q + turn < k) &&
                    (l + turn < p)) {
                  tmp = expGQstem *
                        qm[my_iindx[q + 1] - k + 1] *
                        qm[my_iindx[l + 1] - p + 1] *
                        elem;

                  if (sc_mb_wrapper->pair_ext)
                    tmp *= sc_mb_wrapper->pair_ext(k, l, sc_mb_wrapper);

                    /* store outside part for base pair (k,l) */
                  tmp2 += q_g *
                          tmp;
                }
              }
            }
          }

          /* store total contribution */
          probs[kl] += tmp2 / qo;
        }
      }
    }

    /* finally, outside probabilities for all Gquads spanning the n,1 junction */
    if (VRNA_GQUAD_MAX_BOX_SIZE > l) {
      unsigned int min_k, max_k;

      if (n <= VRNA_GQUAD_MAX_BOX_SIZE + 3) {
        min_k = l + 1;
      } else {
        min_k = n - (VRNA_GQUAD_MAX_BOX_SIZE - l);
        if (min_k <= l)
          min_k = l + 1;
      }

      max_k = n;
      if (VRNA_GQUAD_MIN_BOX_SIZE > l)
        max_k = n + l + 1 - VRNA_GQUAD_MIN_BOX_SIZE;

      for (k = min_k; k <= max_k; k++) {

#ifndef VRNA_DISABLE_C11_FEATURES
        if ((q_g = vrna_smx_csr_get(q_gq, k, l, 0.)) != 0.) {
#else
        if ((q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.)) != 0.) {
#endif
          tmp2 = 0;

          /* 1. hairpin-loop like case, i.e. single gquad, rest unpaired */
          u1    = k - l - 1;
          eval  = (hc->up_hp[l + 1] >= u1) ? 1 : 0;

          if ((eval) &&
              (u1 >= 3)) {
            qbt1 = scale[u1];

            if (md->circ_penalty)
              qbt1 *= pow(vrna_exp_E_hairpin(u1, 0, -1, -1, NULL, pf_params), (double)n_seq);

            if (sc_hp_wrapper->pair_ext)
              qbt1 *= sc_hp_wrapper->pair_ext(k, l, sc_hp_wrapper);

            tmp2 += qbt1;
          }

          /* 2. internal-loop like case, i.e. single gquad plus base pair or second gquad */
          imin = l + 1;
          imax = l + MAXLOOP + 1;
          if (imax >= k)
            imax = k - 1;

          for (i = imin; i <= imax; i++) {
            u1 = i - l - 1;
            jmax = k - 1;
            jmin = i + turn + 1;
            if ((MAXLOOP <= n) && 
                (jmin + MAXLOOP < k + u1))
              jmin = k + u1 - MAXLOOP - 1;

            for (j = jmax; j >= jmin; j--) {
              u2 = k - j - 1;
              if (((u1 == 0) && (u2 < 3)) ||
                  ((u1 < 3) && (u2 == 0)))
                continue;

              eval = (hc->up_int[l + 1] >= u1) ? 1 : 0;
              eval = (hc->up_int[j + 1] >= u2) ? eval : 0;
              /* 2.1 second element is a base pair */
              if ((eval) &&
                  (hc->mx[n * i + j] & (VRNA_CONSTRAINT_CONTEXT_INT_LOOP | VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC)) && 
                  (qb[my_iindx[i] - j] > 0.)) {

                qbt1 = scale[u1 + u2];
                FLT_OR_DBL qint = 1.;

                switch (fc->type) {
                  case VRNA_FC_TYPE_SINGLE:
                    sj = S1[j + 1];
                    si = S1[i - 1];

                    type    = vrna_get_ptype_md(S[j], S[i], md);
                    if (md->dangles == 2)
                      qint *= pf_params->expmismatchI[type][sj][si];

                    if (type > 2)
                      qint *= pf_params->expTermAU;

                    qbt1 *= (FLT_OR_DBL)expintern[u1 + u2];
                    break;

                  case VRNA_FC_TYPE_COMPARATIVE:
                    for (s = 0; s < n_seq; s++) {
                      type = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
                      if (md->dangles == 2)
                        qbt1 *= pf_params->expmismatchI[type][S3[s][j]][S5[s][i]];

                      if (type > 2)
                        qbt1 *= pf_params->expTermAU;

                      us1   = a2s[s][i - 1] - a2s[s][l];
                      us2   = a2s[s][k - 1] - a2s[s][j];
                      qbt1  *= (FLT_OR_DBL)expintern[us1 + us2];
                    }
                    break;
                }

                if (sc_int_wrapper->pair_ext)
                  qbt1 *= sc_int_wrapper->pair_ext(i, j, k, l, sc_int_wrapper);

                /* store quad probs */
                tmp2 += qbt1 * qint *
                        qb[my_iindx[i] - j];

                /* store base pair probs */
                probs[my_iindx[i] - j] += qbt1 * qint *
                                          q_g / qo;
              }

              /* 2.2 second element is a gquad */
              if ((eval) &&
                  (j - i + 1 >= VRNA_GQUAD_MIN_BOX_SIZE) &&
                  (j - i + 1 <= VRNA_GQUAD_MAX_BOX_SIZE) &&
#ifndef VRNA_DISABLE_C11_FEATURES
                  ((tmp3 = vrna_smx_csr_get(q_gq, i, j, 0.)) != 0.)) {
#else
                  ((tmp3 = vrna_smx_csr_FLT_OR_DBL_get(q_gq, i, j, 0.)) != 0.)) {
#endif
                qbt1 = scale[u1 + u2];

                switch (fc->type) {
                  case VRNA_FC_TYPE_SINGLE:
                    qbt1 *= (FLT_OR_DBL)expintern[u1 + u2];
                    break;

                  case VRNA_FC_TYPE_COMPARATIVE:
                    for (s = 0; s < n_seq; s++) {
                      us1   = a2s[s][i - 1] - a2s[s][l];
                      us2   = a2s[s][k - 1] - a2s[s][j];
                      qbt1  *= (FLT_OR_DBL)expintern[us1 + us2];
                    }
                    break;
                }

                if (sc_int_wrapper->pair_ext)
                  qbt1 *= sc_int_wrapper->pair_ext(i, j, k, l, sc_int_wrapper);

                /* store quad probs */
                tmp2 += tmp3 *
                        qbt1;

                /* store other gquad probs */
                probs[my_iindx[i] - j] += qbt1 *
                                          q_g / qo;
              }
            }
          }

          /* 3. multibranch-loop like case, i.e. gquad plus at least two more stems */
          qbt1 = qm2_real[my_iindx[l + 1] - k + 1] *
                 pow(vrna_exp_E_multibranch_stem(0, -1, -1, pf_params) * expMLclosing, (double)n_seq);

          tmp2 += qbt1;

          /* 4. store gquad outside probs for all gquads spanning the n,1 junction */
#ifndef VRNA_DISABLE_C11_FEATURES
          vrna_smx_csr_insert(matrices->p_gq, k, l, tmp2 * q_g / qo);

#else
          vrna_smx_csr_FLT_OR_DBL_insert(matrices->p_gq, k, l, tmp2 * q_g / qo);
#endif
        } /* end if gquad [k:l] */
      } /* end for(k = ... */
    } /* end gquad spanning n1-junction */
  } /* end with_gquad == 1 */

  free(tt);
}
