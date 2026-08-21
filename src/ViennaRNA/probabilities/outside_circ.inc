PRIVATE void
outside_bp_circ(vrna_fold_compound_t  *fc,
                unsigned int          j,
                helper_arrays_t       dp_helpers,
                void                  *data)
{
  char                  *ptype;
  unsigned char         *hard_constraints, eval, with_gquad, report_bpp;
  short                 *S, *S1, **SS, **S5, **S3, si, sj;
  unsigned int          k, l, n, i, s, n_seq, type, type_2, rt, *tt, **a2s, turn,
                        u1, u2, u3, us1, us2, us3,
                        min_i, max_i, min_j, max_j, min_k, max_k, min_l, max_l;
  int                   ij, *rtype, *my_iindx, *jindx;
  FLT_OR_DBL            tmp, tmp2, tmp3, qbt1, qint, qmb,
                        expMLclosing, *expMLbase, expGQstem, norm,
                        *qb, *qm, *probs, *scale, qo, *qm2_real, q_g;
  FLT_OR_DBL            *probs_hp, *probs_int, *probs_mb;
  double                *expintern;
  vrna_hc_t             *hc;
  vrna_exp_param_t      *pf_params;
  vrna_mx_pf_t          *matrices;
  vrna_md_t             *md;

  vrna_smx_csr(FLT_OR_DBL) * q_gq = fc->exp_matrices->q_gq;

  FLT_OR_DBL            (*numerator_f)(vrna_fold_compound_t *fc,
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

  expMLclosing  = pow(pf_params->expMLclosing, (double)n_seq);
  expintern     = &(pf_params->expinternal[0]);
  rtype         = &(pf_params->model_details.rtype[0]);
  with_gquad    = (unsigned char)md->gquad;
  expGQstem     =
    with_gquad ? pow(vrna_exp_E_multibranch_stem(0, -1, -1, pf_params), (double)n_seq) : 1.;

  sc_mb_wrapper   = &(dp_helpers->sc_wrapper_mb);
  sc_hp_wrapper   = &(dp_helpers->sc_wrapper_hp);
  sc_ext_wrapper  = &(dp_helpers->sc_wrapper_ext);
  sc_int_wrapper  = &(dp_helpers->sc_wrapper_int);

  report_bpp  = ((dp_helpers->options.mode & VRNA_EQM_PROBS_BPP) && (dp_helpers->probs_hp)) ? 1 : 0;
  probs_hp    = dp_helpers->probs_hp;
  probs_int   = dp_helpers->probs_int;
  probs_mb    = dp_helpers->probs_mb;

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

  /* 1. exterior pair (i,j) */
  for (i = 1; i + turn < j; i++) {
    ij = my_iindx[i] - j;

    if (qb[ij] > 0.) {
      norm      = numerator_f(fc, i, j) / qo;
      probs[ij] = norm;


      if (fc->type == VRNA_FC_TYPE_SINGLE) {
        type  = vrna_get_ptype_md(S[i], S[j], md);
        rt    = vrna_get_ptype_md(S[j], S[i], md);
      } else {
        for (s = 0; s < n_seq; s++)
          tt[s] = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
      }

      /* 1.1. Exterior Hairpin Contribution */
      tmp2 = vrna_exp_eval_hairpin(fc, j, i, VRNA_EVAL_LOOP_DEFAULT);

      if (report_bpp)
        probs_hp[i] += tmp2 *
                       norm;

      /* 1.2. Exterior Interior Loop Contribution */
      if (hard_constraints[i * n + j] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
        /*
         * 1.2.1. (i,j)  delimts the "left" part of the internal loop
         * (j,i) is "outer pair", (k,l) before (i,j)
         */
        if ((turn + 2 < i) &&     /* does another base pair (k,l) fit before i? */
            (j + MAXLOOP + 1 > n)) {
          u3    = n - j;
          max_k = 1 + MAXLOOP - u3;
          if (max_k + turn + 2 >= i)
            max_k = i - turn - 2;

          for (k = 1; k <= max_k; k++) {
            u1 = k - 1;

            if (hc->up_int[j + 1] < (u1 + u3))
              break;

            min_l = k + turn + 1;
            if ((u1 + u3 + (i - min_l - 1)) > MAXLOOP)
              min_l = i - 1 - (MAXLOOP - u1 - u3);

            for (l = min_l; l < i; l++) {
              u2 = i - l - 1;

              if (hc->up_int[l + 1] < u2)
                continue;

              if (qb[my_iindx[k] - l] == 0.)
                continue;

              eval = (hard_constraints[k * n + l] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) ? 1 : 0;
              if (hc->f)
                eval = hc->f(k, l, i, j, VRNA_DECOMP_PAIR_IL, hc->data);

              if (eval) {
                tmp = qb[my_iindx[k] - l];

                if (fc->type == VRNA_FC_TYPE_SINGLE) {
                  type_2 = vrna_get_ptype(jindx[l] + k, ptype);

                  tmp *= vrna_exp_E_internal(u1 + u3,
                                             u2,
                                             rt,
                                             rtype[type_2],
                                             S1[j + 1],
                                             S1[i - 1],
                                             S1[k - 1],
                                             S1[l + 1],
                                             pf_params);
                } else {
                  for (s = 0; s < n_seq; s++) {
                    us2     = a2s[s][i - 1];
                    us2     -= a2s[s][l];
                    us1     = a2s[s][n] - a2s[s][j];
                    us1     += a2s[s][k - 1];
                    type_2  = vrna_get_ptype_md(SS[s][l], SS[s][k], md);
                    tmp     *= vrna_exp_E_internal(us1, us2, tt[s], type_2,
                                                   S3[s][j],
                                                   S5[s][i],
                                                   S5[s][k],
                                                   S3[s][l],
                                                   pf_params);
                  }
                }

                if (sc_int_wrapper->pair_ext)
                  tmp *= sc_int_wrapper->pair_ext(k, l, i, j, sc_int_wrapper);

                tmp2 += tmp *
                        scale[u1 + u2 + u3];

                if (report_bpp)
                  probs_int[i] += tmp *
                                  scale[u1 + u2 + u3] *
                                  norm;
              }
            }
          }
        } /* end (k,l) before (i,j) internal loop */

        /*
         * 1.2.2. i,j  delimtis the "right" part of the internal loop
         * (k,l) after (i,j)
         */
        if ((j + turn + 2 <= n) &&  /* does another base pair (k,l) fit after j? */
            ((i - 1) <= MAXLOOP)) {
          u2    = i - 1;
          max_k = j + 1 + (MAXLOOP - u2);

          if (max_k + turn + 1 > n)
            max_k = n - turn - 1;

          for (k = j + 1; k <= max_k; k++) {
            u1 = k - j - 1;

            if (hc->up_int[j + 1] < u1)
              break;

            min_l = k + turn + 1;

            if ((u1 + u2 + (n - min_l)) > MAXLOOP)
              min_l = n + u1 + u2 - MAXLOOP;

            for (l = min_l; l <= n; l++) {
              u3 = n - l;

              if (hc->up_int[l + 1] < (u2 + u3))
                continue;

              if (qb[my_iindx[k] - l] == 0.)
                continue;

              eval = (hard_constraints[k * n + l] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) ? 1 : 0;
              if (hc->f)
                eval = hc->f(i, j, k, l, VRNA_DECOMP_PAIR_IL, hc->data) ? eval : 0;

              if (eval) {
                tmp = qb[my_iindx[k] - l];

                if (fc->type == VRNA_FC_TYPE_SINGLE) {
                  type_2  = vrna_get_ptype(jindx[l] + k, ptype);
                  tmp     *= vrna_exp_E_internal(u2 + u3,
                                                 u1,
                                                 rtype[type_2],
                                                 rt,
                                                 S1[l + 1],
                                                 S1[k - 1],
                                                 S1[i - 1],
                                                 S1[j + 1],
                                                 pf_params);
                } else {
                  for (s = 0; s < n_seq; s++) {
                    us1     = a2s[s][k] - a2s[s][j + 1];
                    us2     = a2s[s][i - 1] + a2s[s][n] - a2s[s][l];
                    type_2  = vrna_get_ptype_md(SS[s][l], SS[s][k], md);
                    tmp     *= vrna_exp_E_internal(us2, us1, type_2, tt[s],
                                                   S3[s][l],
                                                   S5[s][k],
                                                   S5[s][i],
                                                   S3[s][j],
                                                   pf_params);
                  }
                }

                if (sc_int_wrapper->pair_ext)
                  tmp *= sc_int_wrapper->pair_ext(i, j, k, l, sc_int_wrapper);

                tmp2 += tmp *
                        scale[u1 + u2 + u3];

                if (report_bpp)
                  probs_int[i] += tmp *
                                  scale[u1 + u2 + u3] *
                                  norm;
              }
            }
          }
        } /* end (k,l) after (i,j) internal loop */

        /* 1.3 Exterior multiloop decomposition */
        if (hc->eval_mb(i, j, i - 1, j + 1, VRNA_DECOMP_PAIR_ML_EXT, hc)) {
          FLT_OR_DBL sc_contrib = 1.;

          if (sc_mb_wrapper->pair_ext)
            sc_contrib = sc_mb_wrapper->pair_ext(i, j, sc_mb_wrapper);

          /* 1.3.1 Middle part                    */
          if ((i > 2) &&
              (j < n - 1)) {
            tmp = qm[my_iindx[1] - i + 1] *
                  qm[my_iindx[j + 1] - n];

            if (fc->type == VRNA_FC_TYPE_SINGLE) {
              tmp *= vrna_exp_E_multibranch_stem(type,
                                                 S1[i - 1],
                                                 S1[j + 1],
                                                 pf_params);
            } else {
              for (s = 0; s < n_seq; s++)
                tmp *= vrna_exp_E_multibranch_stem(rtype[tt[s]],
                                                   S5[s][i],
                                                   S3[s][j],
                                                   pf_params);
            }

            tmp2 += tmp *
                    sc_contrib *
                    expMLclosing;

            if (report_bpp)
              probs_mb[i] += tmp *
                             sc_contrib *
                             expMLclosing *
                             norm;
          }

          /* 1.3.2 right-most part  */
          if (hc->up_ml[j + 1] >= (n - j)) {
            if (i > 1) {
              if ((hc->eval_mb(i, n, i, j, VRNA_DECOMP_ML_ML, hc)) &&
                  (hc->eval_mb(1, j, i - 1, i, VRNA_DECOMP_ML_ML_ML, hc))) {
                tmp = qm2_real[my_iindx[1] - i + 1] *
                      expMLbase[n - j];

                if (fc->type == VRNA_FC_TYPE_SINGLE) {
                  tmp *= vrna_exp_E_multibranch_stem(type,
                                                     S1[i - 1],
                                                     S1[j + 1],
                                                     pf_params);
                } else {
                  for (s = 0; s < n_seq; s++)
                    tmp *= vrna_exp_E_multibranch_stem(rtype[tt[s]],
                                                       S5[s][i],
                                                       S3[s][j],
                                                       pf_params);
                }

                if (sc_mb_wrapper->red_ml)
                  tmp *= sc_mb_wrapper->red_ml(i, n, i, j, sc_mb_wrapper);

                if (sc_mb_wrapper->decomp_ml)
                  tmp *= sc_mb_wrapper->decomp_ml(1, j, i - 1, i, sc_mb_wrapper);

                tmp2 += tmp *
                        sc_contrib *
                        expMLclosing;

                if (report_bpp)
                  probs_mb[i] += tmp *
                                 sc_contrib *
                                 expMLclosing *
                                 norm;
              }
            }
          }

          /* 1.3.3 left-most part */
          if (hc->up_ml[1] >= (i - 1)) {
            if (j + 1 < n) {
              if ((hc->eval_mb(1, j, i, j, VRNA_DECOMP_ML_ML, hc)) &&
                  (hc->eval_mb(i, n, j, j + 1, VRNA_DECOMP_ML_ML_ML, hc))) {
                tmp = qm2_real[my_iindx[j + 1] - n] *
                      expMLbase[i - 1];

                if (fc->type == VRNA_FC_TYPE_SINGLE) {
                  tmp *= vrna_exp_E_multibranch_stem(type,
                                                     S1[i - 1],
                                                     S1[j + 1],
                                                     pf_params);
                } else {
                  for (s = 0; s < n_seq; s++)
                    tmp *= vrna_exp_E_multibranch_stem(rtype[tt[s]],
                                                       S5[s][i],
                                                       S3[s][j],
                                                       pf_params);
                }

                if (sc_mb_wrapper->red_ml)
                  tmp *= sc_mb_wrapper->red_ml(1, j, i, j, sc_mb_wrapper);

                if (sc_mb_wrapper->decomp_ml)
                  tmp *= sc_mb_wrapper->decomp_ml(i, n, j, j + 1, sc_mb_wrapper);

                tmp2 += tmp *
                        sc_contrib *
                        expMLclosing;

                if (report_bpp)
                  probs_mb[i] += tmp *
                                 sc_contrib *
                                 expMLclosing *
                                 norm;
              }
            }
          }
        } /* end multibranch loop part */

        /* all exterior loop decompositions for pair (i,j) done  */

        if (with_gquad) {
          /*
           * consider all the other possibilities where at least one gquad is enclosed in the same loop as (i,j)
           * 1. internal-loop like structure
           */
          if (hard_constraints[i * n + j] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
            /* 1.1 gquad is 5' of base pair (i,j) */

            qint = 1.;

            switch (fc->type) {
              case VRNA_FC_TYPE_SINGLE:
                sj  = S1[j + 1];
                si  = S1[i - 1];

                if (md->dangles == 2)
                  qint *= pf_params->expmismatchI[rt][sj][si];

                if (rt > 2)
                  qint *= pf_params->expTermAU;

                break;

              case VRNA_FC_TYPE_COMPARATIVE:
                for (s = 0; s < n_seq; s++) {
                  if (md->dangles == 2)
                    qint *= pf_params->expmismatchI[tt[s]][S3[s][j]][S5[s][i]];

                  if (tt[s] > 2)
                    qint *= pf_params->expTermAU;
                }
                break;
            }

            /* int loop with [k,l] gquad followed by (i,j) base pair */
            if (n - j <= MAXLOOP) {
              u1    = n - j;
              max_k = 1 + MAXLOOP - u1;
              for (k = 1; k <= max_k; k++) {
                u2    = k - 1;
                max_l = i - 1;
                if (k + VRNA_GQUAD_MAX_BOX_SIZE - 1 < max_l)
                  max_l = k + VRNA_GQUAD_MAX_BOX_SIZE - 1;

                min_l = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
                if (min_l >= i)
                  min_l = i;
                else if ((i - min_l) + u1 + u2 > MAXLOOP + 1)
                  min_l = i - (MAXLOOP - u1 - u2) - 1;

                for (l = min_l; l <= max_l; l++) {
                  u3 = i - l - 1;
                  if (((u1 + u2 == 0) && (u3 < 3)) ||
                      ((u1 + u2 < 3) && (u3 == 0)))
                    continue;

#ifndef VRNA_DISABLE_C11_FEATURES
                  q_g = vrna_smx_csr_get(q_gq, k, l, 0.);
#else
                  q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.);
#endif
                  if (q_g != 0.) {
                    tmp = scale[u1 + u2 + u3];

                    switch (fc->type) {
                      case VRNA_FC_TYPE_SINGLE:
                        tmp *= (FLT_OR_DBL)expintern[u1 + u2 + u3];
                        break;

                      case VRNA_FC_TYPE_COMPARATIVE:
                        for (s = 0; s < n_seq; s++) {
                          us1 = (k > 1) ? a2s[s][k - 1] - a2s[s][1] : 0;
                          us2 = a2s[s][i - 1] - a2s[s][l];
                          us3 = a2s[s][n] - a2s[s][j];
                          tmp *= (FLT_OR_DBL)expintern[us1 + us2 + us3];
                        }
                        break;
                    }

                    if (sc_int_wrapper->pair_ext)
                      tmp *= sc_int_wrapper->pair_ext(k, l, i, j, sc_int_wrapper);

                    /* store outside part for base pair */
                    tmp2 += qint *
                            q_g *
                            tmp;

                    if (report_bpp)
                      probs_int[i] += qint *
                                      q_g *
                                      tmp *
                                      norm;
                  }
                }
              }
            }

            /* 1.2 gquad is 3' of base pair (i,j) */
            u1    = i - 1;
            max_k = j + MAXLOOP - u1 + 1;
            if (j + VRNA_GQUAD_MIN_BOX_SIZE <= n)
              max_k = MIN2(max_k, n - VRNA_GQUAD_MIN_BOX_SIZE + 1);
            else
              max_k = j;

            for (k = j + 1; k <= max_k; k++) {
              u2    = k - j - 1;
              min_l = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
              if (min_l > n)
                break;

              if (min_l + MAXLOOP - u1 - u2 < n)
                min_l = n + u1 + u2 - MAXLOOP;

              max_l = k + VRNA_GQUAD_MAX_BOX_SIZE - 1;
              if (max_l > n)
                max_l = n;

              for (l = min_l; l <= max_l; l++) {
                u3 = n - l;
                if (((u1 + u3 == 0) && (u2 < 3)) ||
                    ((u1 + u3 < 3) && (u2 == 0)))
                  continue;

#ifndef VRNA_DISABLE_C11_FEATURES
                q_g = vrna_smx_csr_get(q_gq, k, l, 0.);
#else
                q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.);
#endif
                if (q_g != 0.) {
                  tmp = scale[u1 + u2 + u3];

                  switch (fc->type) {
                    case VRNA_FC_TYPE_SINGLE:
                      tmp *= (FLT_OR_DBL)expintern[u1 + u2 + u3];
                      break;

                    case VRNA_FC_TYPE_COMPARATIVE:
                      for (s = 0; s < n_seq; s++) {
                        us1 = a2s[s][i - 1] - a2s[s][1];
                        us2 = a2s[s][k - 1] - a2s[s][j];
                        us3 = a2s[s][n] - a2s[s][l];
                        tmp *= (FLT_OR_DBL)expintern[us1 + us2 + us3];
                      }
                      break;
                  }

                  if (sc_int_wrapper->pair_ext)
                    tmp *= sc_int_wrapper->pair_ext(i, j, k, l, sc_int_wrapper);

                  /* store outside part for base pair */
                  tmp2 += qint *
                          q_g *
                          tmp;

                  if (report_bpp)
                    probs_int[i] += qint *
                                    q_g *
                                    tmp *
                                    norm;
                }
              }
            }

            /* gquad [k,l] spanning n,1-junction follows after base pair (i,j) */
            if ((i > 1) &&
                (j < n)) {
              min_l = 1;
              max_l = VRNA_GQUAD_MAX_BOX_SIZE - 1;
              if (max_l >= i)
                max_l = i - 1;

              u1 = i - max_l - 1;
              for (l = max_l; l >= min_l; l--, u1++) {
                if (hc->up_int[l + 1] < u1)
                  break;

                if (u1 > MAXLOOP)
                  break;

                min_k = j + 1;
                if (n - min_k + 1 + l > VRNA_GQUAD_MAX_BOX_SIZE)
                  min_k = n + l + 1 - VRNA_GQUAD_MAX_BOX_SIZE;

                if (min_k - j - 1 + u1 > MAXLOOP)
                  continue;

                max_k = n;
                if (l + n - max_k + 1 < VRNA_GQUAD_MIN_BOX_SIZE)
                  max_k = n + l + 1 - VRNA_GQUAD_MIN_BOX_SIZE;

                for (k = min_k; k <= max_k; k++) {
                  u2 = k - j - 1;

                  if (hc->up_int[j + 1] < u2)
                    break;

                  if (u1 + u2 > MAXLOOP)
                    break;

                  if (((u1 == 0) && (u2 < 3)) ||
                      ((u1 < 3) && (u2 == 0)))
                    continue;

#ifndef VRNA_DISABLE_C11_FEATURES
                  if ((q_g = vrna_smx_csr_get(q_gq, k, l, 0.)) != 0.) {
#else
                  if ((q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.)) != 0.) {
#endif
                    tmp = scale[u1 + u2];

                    switch (fc->type) {
                      case VRNA_FC_TYPE_SINGLE:
                        tmp *= (FLT_OR_DBL)expintern[u1 + u2];
                        break;

                      case VRNA_FC_TYPE_COMPARATIVE:
                        for (s = 0; s < n_seq; s++) {
                          us1 = a2s[s][i - 1] - a2s[s][l];
                          us2 = a2s[s][k - 1] - a2s[s][j];
                          tmp *= (FLT_OR_DBL)expintern[us1 + us2];
                        }
                        break;
                    }

                    if (sc_int_wrapper->pair_ext)
                      tmp *= sc_int_wrapper->pair_ext(i, j, k, l, sc_int_wrapper);

                    tmp2 += qint *
                            q_g *
                            tmp;

                    if (report_bpp)
                      probs_int[i] += qint *
                                      q_g *
                                      tmp *
                                      norm;
                  }
                }
              }
            }
          }

          if (hard_constraints[i * n + j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
            /*
             * 3. all base pairs (i,j) where G-Quad spans across n,1 boundary
             * 3.2 multibranch-loop like configurations
             */
            qmb = 1.;

            switch (fc->type) {
              case VRNA_FC_TYPE_SINGLE:
                qmb *= vrna_exp_E_multibranch_stem(type, S1[i - 1], S1[j + 1], pf_params);
                break;
              case VRNA_FC_TYPE_COMPARATIVE:
                for (s = 0; s < n_seq; s++)
                  qmb *= vrna_exp_E_multibranch_stem(rtype[tt[s]],
                                                     S5[s][i],
                                                     S3[s][j],
                                                     pf_params);
                break;
            }

            qmb *= expMLclosing;

            min_k = j + 1;
            if (min_k + VRNA_GQUAD_MAX_BOX_SIZE - 1 <= n)
              min_k = n - VRNA_GQUAD_MAX_BOX_SIZE + 1;

            for (k = min_k; k <= n; k++) {
              min_l = 1;
              max_l = VRNA_GQUAD_MAX_BOX_SIZE - 1 - (n - k);
              if (max_l >= i)
                max_l = i - 1;

              for (l = min_l; l <= max_l; l++) {
#ifndef VRNA_DISABLE_C11_FEATURES
                if ((q_g = vrna_smx_csr_get(q_gq, k, l, 0.)) != 0.) {
#else
                if ((q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.)) != 0.) {
#endif
                  /* 3.2.1 (i,j) last branch before gquad, ie. [k,l] + ML + (i,j) */
                  if (l + turn < i) {
                    u1 = k - j - 1;
                    if (hc->up_ml[j + 1] >= u1) {
                      tmp = qmb *
                            qm[my_iindx[l + 1] - i + 1] *
                            expGQstem *
                            expMLbase[u1];

                      if (sc_mb_wrapper->pair_ext)
                        tmp *= sc_mb_wrapper->pair_ext(i, j, sc_mb_wrapper);

                      /* store outside part for base pair (i,j) */
                      tmp2 += q_g *
                              tmp;

                      if (report_bpp)
                        probs_int[i] += q_g *
                                        tmp *
                                        norm;
                    }
                  }

                  /* 3.2.2 (i,j) first branch after gquad, ie. [k,l] + (i,j) + ML */
                  if (j + turn < k) {
                    u2 = i - l - 1;
                    if (hc->up_ml[l + 1] >= u2) {
                      tmp = qmb *
                            qm[my_iindx[j + 1] - k + 1] *
                            expGQstem *
                            expMLbase[u2];

                      if (sc_mb_wrapper->pair_ext)
                        tmp *= sc_mb_wrapper->pair_ext(i, j, sc_mb_wrapper);

                      /* store outside part for base pair (i,j) */
                      tmp2 += q_g *
                              tmp;

                      if (report_bpp)
                        probs_int[i] += q_g *
                                        tmp *
                                        norm;
                    }
                  }

                  /* 3.2.3 (i,j) somewhere in between, ie. [k,l] + ML + (i,j) + ML */
                  if ((l + turn < i) &&
                      (j + turn < k)) {
                    tmp = qmb *
                          qm[my_iindx[l + 1] - i + 1] *
                          qm[my_iindx[j + 1] - k + 1] *
                          expGQstem;

                    if (sc_mb_wrapper->pair_ext)
                      tmp *= sc_mb_wrapper->pair_ext(i, j, sc_mb_wrapper);

                    /* store outside part for base pair (i,j) */
                    tmp2 += q_g *
                            tmp;

                    if (report_bpp)
                      probs_int[i] += q_g *
                                      tmp *
                                      norm;
                  }
                }
              }
            }
          }
        }

        probs[ij] *= tmp2;
      }
    } /* end if qb[ij] > 0 */
  }   /* end for i */

  /* Done with outside probabilities for base pairs (i,j) */

  /* Next follows outside probabilities for G-Quadruplexes delimited by [i : j] */
  if (with_gquad) {
    if (VRNA_GQUAD_MIN_BOX_SIZE <= j) {
      min_k = 1;
      if (min_k + VRNA_GQUAD_MAX_BOX_SIZE - 1 < j)
        min_k = j + 1 - VRNA_GQUAD_MAX_BOX_SIZE;

      for (i = min_k; i + VRNA_GQUAD_MIN_BOX_SIZE - 1 <= j; i++) {
        ij = my_iindx[i] - j;

#ifndef VRNA_DISABLE_C11_FEATURES
        if ((q_g = vrna_smx_csr_get(q_gq, i, j, 0.)) != 0.) {
#else
        if ((q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, i, j, 0.)) != 0.) {
#endif
          tmp2 = 0.;

          /* 1. all possible G-quads spanning from i to j */

          /* 1.1 hairpin-loop like configurations */
          if ((i - 1) + (n - j) >= 3) {
            eval = (hc->up_ext[1] >= (i - 1)) ? 1 : 0;
            if (j < n)
              eval = (hc->up_ext[j + 1] >= (n - j)) ? eval : 0;

            if (hc->f) {
              if (i > 1)
                eval = (hc->f(1, i - 1, 1, i - 1, VRNA_DECOMP_EXT_UP, hc->data)) ? eval : 0;

              if (j < n)
                eval = (hc->f(j + 1, n, j + 1, n, VRNA_DECOMP_EXT_UP, hc->data)) ? eval : 0;
            }

            if (eval) {
              tmp = scale[(i - 1) + (n - j)];

              if (md->circ_penalty)
                tmp *=
                  pow(vrna_exp_E_hairpin((i - 1) + (n - j), 0, -1, -1, NULL, pf_params),
                      (double)n_seq);

              if (sc_ext_wrapper->red_up)
                tmp *= sc_ext_wrapper->red_up(1, i - 1, sc_ext_wrapper) *
                       sc_ext_wrapper->red_up(j + 1, n, sc_ext_wrapper);

              tmp2 += tmp;
            }
          }

          /* 1.2 internal-loop like configurations */
          u1 = i - 1;

          /* 1.2.1 gquad [i,j] precedes a base pair (k,l) */
          if ((u1 <= MAXLOOP) &&
              (j + turn + 1 < n) &&
              (hc->up_int[1] >= u1)) {
            max_k = j + 1 + MAXLOOP - u1;

            if (max_k + turn >= n)
              max_k = n - turn - 1;

            for (k = j + 1; k <= max_k; k++) {
              u2 = k - j - 1;
              if (hc->up_int[j + 1] < u2)
                break;

              min_l = k + turn + 1;
              if ((u1 + u2 + (n - min_l)) > MAXLOOP)
                min_l = n - (MAXLOOP - u1 - u2);

              for (l = min_l; l <= n; l++) {
                u3 = n - l;

                if (((u1 + u3 == 0) && (u2 < 3)) ||
                    ((u1 + u3 < 3) && (u2 == 0)))
                  continue;

                if (qb[my_iindx[k] - l] == 0.)
                  continue;

                if (hard_constraints[k * n + l] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
                  qint = 1.;

                  switch (fc->type) {
                    case VRNA_FC_TYPE_SINGLE:
                      sj  = S1[l + 1];
                      si  = S1[k - 1];
                      rt  = vrna_get_ptype_md(S[l], S[k], md);

                      if (md->dangles == 2)
                        qint *= pf_params->expmismatchI[rt][sj][si];

                      if (rt > 2)
                        qint *= pf_params->expTermAU;

                      break;

                    case VRNA_FC_TYPE_COMPARATIVE:
                      for (s = 0; s < n_seq; s++) {
                        rt = vrna_get_ptype_md(SS[s][l], SS[s][k], md);

                        if (md->dangles == 2)
                          qint *= pf_params->expmismatchI[rt][S3[s][l]][S5[s][k]];

                        if (rt > 2)
                          qint *= pf_params->expTermAU;
                      }
                      break;
                  }

                  tmp = scale[u1 + u2 + u3];

                  switch (fc->type) {
                    case VRNA_FC_TYPE_SINGLE:
                      tmp *= (FLT_OR_DBL)expintern[u1 + u2 + u3];
                      break;

                    case VRNA_FC_TYPE_COMPARATIVE:
                      for (s = 0; s < n_seq; s++) {
                        us1 = (i > 1) ? a2s[s][i - 1] - a2s[s][1] : 0;
                        us2 = a2s[s][k - 1] - a2s[s][j];
                        us3 = a2s[s][n] - a2s[s][l];
                        tmp *= (FLT_OR_DBL)expintern[us1 + us2 + us3];
                      }
                      break;
                  }

                  if (sc_int_wrapper->pair_ext)
                    tmp *= sc_int_wrapper->pair_ext(i, j, k, l, sc_int_wrapper);

                  probs[my_iindx[i] - j] += qint *
                                            tmp *
                                            qb[my_iindx[k] - l] /
                                            qo;
                }
              }
            }
          }

          /* 1.2.2 base pair (k,l) precedes the gquad [i,j] */
          u3 = n - j;

          if ((turn + 2 < i) &&
              (u3 <= MAXLOOP) &&
              (hc->up_int[j + 1] >= u3)) {
            min_l = turn + 2;
            if (j + u3 - min_l - 1 > MAXLOOP)
              min_l = j + 2 - MAXLOOP;

            for (l = i - 1; l >= min_l; l--) {
              u2 = i - l;
              if (hc->up_int[l + 1] < u2)
                break;

              max_k = l - turn - 1;
              if (max_k - 1 + u2 + u3 > MAXLOOP)
                max_k = l + u2 + u3 - turn - 1 - MAXLOOP;

              for (k = 1; k <= max_k; k++) {
                u1 = k - 1;

                if (((u1 + u3 == 0) && (u2 < 3)) ||
                    ((u1 + u3 < 3) && (u2 == 0)))
                  continue;

                if (qb[my_iindx[k] - l] == 0.)
                  continue;

                if (hard_constraints[k * n + l] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
                  qint = 1.;

                  switch (fc->type) {
                    case VRNA_FC_TYPE_SINGLE:
                      sj  = S1[l + 1];
                      si  = S1[k - 1];
                      rt  = vrna_get_ptype_md(S[l], S[k], md);

                      if (md->dangles == 2)
                        qint *= pf_params->expmismatchI[rt][sj][si];

                      if (rt > 2)
                        qint *= pf_params->expTermAU;

                      break;

                    case VRNA_FC_TYPE_COMPARATIVE:
                      for (s = 0; s < n_seq; s++) {
                        rt = vrna_get_ptype_md(SS[s][l], SS[s][k], md);

                        if (md->dangles == 2)
                          qint *= pf_params->expmismatchI[rt][S3[s][l]][S5[s][k]];

                        if (rt > 2)
                          qint *= pf_params->expTermAU;
                      }
                      break;
                  }

                  tmp = scale[u1 + u2 + u3];

                  switch (fc->type) {
                    case VRNA_FC_TYPE_SINGLE:
                      tmp *= (FLT_OR_DBL)expintern[u1 + u2 + u3];
                      break;

                    case VRNA_FC_TYPE_COMPARATIVE:
                      for (s = 0; s < n_seq; s++) {
                        us1 = (i > 1) ? a2s[s][i - 1] - a2s[s][1] : 0;
                        us2 = a2s[s][k - 1] - a2s[s][j];
                        us3 = a2s[s][n] - a2s[s][l];
                        tmp *= (FLT_OR_DBL)expintern[us1 + us2 + us3];
                      }
                      break;
                  }

                  if (sc_int_wrapper->pair_ext)
                    tmp *= sc_int_wrapper->pair_ext(i, j, k, l, sc_int_wrapper);

                  probs[my_iindx[i] - j] += qint *
                                            tmp *
                                            qb[my_iindx[k] - l] /
                                            qo;
                }
              }
            }
          }

          /* 1.2.3 gquad [i,j] preceded by gquad [k,l] */
          u3 = n - j;

          if ((u3 < MAXLOOP) &&
              (VRNA_GQUAD_MIN_BOX_SIZE < i) &&
              (hc->up_int[j + 1] >= u3)) {
            min_l = VRNA_GQUAD_MIN_BOX_SIZE;
            if (u3 + i - min_l - 1 > MAXLOOP)
              min_l = i + u3 - 1 - MAXLOOP;

            for (l = i - 1; l >= min_l; l--) {
              u2 = i - l - 1;
              if (hc->up_int[l + 1] < u2)
                break;

              min_k = 1;
              if (min_k + VRNA_GQUAD_MAX_BOX_SIZE - 1 < l)
                min_k = l - VRNA_GQUAD_MAX_BOX_SIZE + 1;

              for (k = min_k; k + VRNA_GQUAD_MIN_BOX_SIZE - 1 <= l; k++) {
                u1 = k - 1;

                if (hc->up_int[1] < u1)
                  break;

                if (((u1 + u3 == 0) && (u2 < 3)) ||
                    ((u1 + u3 < 3) && (u2 == 0)))
                  continue;

#ifndef VRNA_DISABLE_C11_FEATURES
                tmp3 = vrna_smx_csr_get(q_gq, k, l, 0.);
#else
                tmp3 = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.);
#endif

                if (tmp3 != 0.) {
                  tmp = scale[u1 + u2 + u3];

                  switch (fc->type) {
                    case VRNA_FC_TYPE_SINGLE:
                      tmp *= (FLT_OR_DBL)expintern[u1 + u2 + u3];
                      break;

                    case VRNA_FC_TYPE_COMPARATIVE:
                      for (s = 0; s < n_seq; s++) {
                        us1 = (i > 1) ? a2s[s][i - 1] - a2s[s][1] : 0;
                        us2 = a2s[s][k - 1] - a2s[s][j];
                        us3 = a2s[s][n] - a2s[s][l];
                        tmp *= (FLT_OR_DBL)expintern[us1 + us2 + us3];
                      }
                      break;
                  }

                  /* soft constraints for at least unpaired bases in the internal-loop like conformation?! */

                  probs[my_iindx[i] - j] += tmp *
                                            tmp3;
                }
              }
            }
          }

          u1 = i - 1;

          /* 1.2.4 gquad [i,j] followed by gquad [k,l] */
          if ((j + VRNA_GQUAD_MIN_BOX_SIZE <= n) &&
              (u1 <= MAXLOOP) &&
              (hc->up_int[1] >= u1)) {
            max_k = j + 1 + MAXLOOP - u1;
            if (max_k + VRNA_GQUAD_MIN_BOX_SIZE - 1 > n)
              max_k = n - VRNA_GQUAD_MIN_BOX_SIZE + 1;

            for (k = j + 1; k <= max_k; k++) {
              u2 = k - j - 1;
              if (hc->up_int[j + 1] < u2)
                break;

              min_l = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
              if (min_l + MAXLOOP < n + u1 + u2)
                min_l = n + u1 + u2 - MAXLOOP;

              for (u3 = 0, l = n; l >= min_l; l--, u3++) {
                if (hc->up_int[l + 1] < u3)
                  break;

                if ((((u1 + u3) == 0) && (u2 < 3)) ||
                    (((u1 + u3) < 3) && (u2 == 0)))
                  continue;

#ifndef VRNA_DISABLE_C11_FEATURES
                tmp3 = vrna_smx_csr_get(q_gq, k, l, 0.);
#else
                tmp3 = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.);
#endif

                if (tmp3 != 0.) {
                  tmp = scale[u1 + u2 + u3];

                  switch (fc->type) {
                    case VRNA_FC_TYPE_SINGLE:
                      tmp *= (FLT_OR_DBL)expintern[u1 + u2 + u3];
                      break;

                    case VRNA_FC_TYPE_COMPARATIVE:
                      for (s = 0; s < n_seq; s++) {
                        us1 = (i > 1) ? a2s[s][i - 1] - a2s[s][1] : 0;
                        us2 = a2s[s][k - 1] - a2s[s][j];
                        us3 = a2s[s][n] - a2s[s][l];
                        tmp *= (FLT_OR_DBL)expintern[us1 + us2 + us3];
                      }
                      break;
                  }

                  /* soft constraints for at least unpaired bases in the internal-loop like conformation?! */

                  tmp2 += tmp * tmp3;
                }
              }
            }
          }

          /* 1.2.5  guad [i,j] is followed by gquad [k,l] spanning n,1-junction */
          if ((i > 1) &&
              (j < n)) {
            min_l = 1;
            max_l = VRNA_GQUAD_MAX_BOX_SIZE - 1;
            if (max_l >= i)
              max_l = i - 1;

            u1 = i - max_l - 1;
            for (l = max_l; l >= min_l; l--, u1++) {
              if (hc->up_int[l + 1] < u1)
                break;

              if (u1 > MAXLOOP)
                break;

              min_k = j + 1;
              if (n - min_k + 1 + l > VRNA_GQUAD_MAX_BOX_SIZE)
                min_k = n + l + 1 - VRNA_GQUAD_MAX_BOX_SIZE;

              if (min_k - j - 1 + u1 > MAXLOOP)
                continue;

              max_k = n;
              if (l + n - max_k + 1 < VRNA_GQUAD_MIN_BOX_SIZE)
                max_k = n + l + 1 - VRNA_GQUAD_MIN_BOX_SIZE;

              for (k = min_k; k <= max_k; k++) {
                u2 = k - j - 1;

                if (hc->up_int[j + 1] < u2)
                  break;

                if (u1 + u2 > MAXLOOP)
                  break;

                if (((u1 == 0) && (u2 < 3)) ||
                    ((u1 < 3) && (u2 == 0)))
                  continue;

#ifndef VRNA_DISABLE_C11_FEATURES
                tmp3 = vrna_smx_csr_get(q_gq, k, l, 0.);
#else
                tmp3 = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.);
#endif

                if (tmp3 > 0) {
                  tmp = scale[u1 + u2];

                  switch (fc->type) {
                    case VRNA_FC_TYPE_SINGLE:
                      tmp *= (FLT_OR_DBL)expintern[u1 + u2];
                      break;

                    case VRNA_FC_TYPE_COMPARATIVE:
                      for (s = 0; s < n_seq; s++) {
                        us1 = a2s[s][i - 1] - a2s[s][l];
                        us2 = a2s[s][k - 1] - a2s[s][j];
                        tmp *= (FLT_OR_DBL)expintern[us1 + us2];
                      }
                      break;
                  }

                  if (sc_int_wrapper->pair_ext)
                    tmp *= sc_int_wrapper->pair_ext(i, j, k, l, sc_int_wrapper);

                  /* soft constraints for at least unpaired bases in the internal-loop like conformation?! */

                  tmp2 += tmp *
                          tmp3;
                }
              }
            }
          }

          /* 1.3 multibranch-loop like configurations with [i,j] gquad */

          /*
           * 1.3.1 no element spanning n,1 junction
           * 1.3.1.1 Middle part
           */
          if ((i > 2) &&
              (j < n - 1)) {
            tmp = qm[my_iindx[1] - i + 1] *
                  qm[my_iindx[j + 1] - n];


            tmp2 += tmp * expGQstem;
          }

          /* 1.3.1,2 right-most part  */
          if (hc->up_ml[j + 1] >= (n - j)) {
            if (i > 1) {
              if ((hc->eval_mb(i, n, i, j, VRNA_DECOMP_ML_ML, hc)) &&
                  (hc->eval_mb(1, j, i - 1, i, VRNA_DECOMP_ML_ML_ML, hc))) {
                tmp = qm2_real[my_iindx[1] - i + 1] *
                      expMLbase[n - j];

                if (sc_mb_wrapper->red_ml)
                  tmp *= sc_mb_wrapper->red_ml(i, n, i, j, sc_mb_wrapper);

                if (sc_mb_wrapper->decomp_ml)
                  tmp *= sc_mb_wrapper->decomp_ml(1, j, i - 1, i, sc_mb_wrapper);

                tmp2 += tmp *
                        expGQstem;
              }
            }
          }

          /* 1.3.1.3 left-most part */
          if (hc->up_ml[1] >= (i - 1)) {
            if (j + 1 < n) {
              if ((hc->eval_mb(1, j, i, j, VRNA_DECOMP_ML_ML, hc)) &&
                  (hc->eval_mb(i, n, j, j + 1, VRNA_DECOMP_ML_ML_ML, hc))) {
                tmp = qm2_real[my_iindx[j + 1] - n] *
                      expMLbase[i - 1];

                if (sc_mb_wrapper->red_ml)
                  tmp *= sc_mb_wrapper->red_ml(1, j, i, j, sc_mb_wrapper);

                if (sc_mb_wrapper->decomp_ml)
                  tmp *= sc_mb_wrapper->decomp_ml(i, n, j, j + 1, sc_mb_wrapper);

                tmp2 += tmp *
                        expGQstem;
              }
            }
          }

          /* 1.3.2 one element spanning n,1 junction */
          min_k = j + 1;

          if (min_k + VRNA_GQUAD_MAX_BOX_SIZE - 1 <= n)
            min_k = n - VRNA_GQUAD_MAX_BOX_SIZE + 1;

          for (k = min_k; k <= n; k++) {
            min_l = 1;
            max_l = VRNA_GQUAD_MAX_BOX_SIZE - 1 - (n - k);

            if (max_l >= i)
              max_l = i - 1;

            for (l = min_l; l <= max_l; l++) {
#ifndef VRNA_DISABLE_C11_FEATURES
              if ((q_g = vrna_smx_csr_get(q_gq, k, l, 0.)) != 0.) {
#else
              if ((q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.)) != 0.) {
#endif
                /* 3.2.1 (i,j) last branch before gquad, ie. [k,l] + ML + (i,j) */
                if (l + turn < i) {
                  u1 = k - j - 1;
                  if (hc->up_ml[j + 1] >= u1) {
                    tmp = expGQstem *
                          qm[my_iindx[l + 1] - i + 1] *
                          expGQstem *
                          expMLbase[u1];

                    if (sc_mb_wrapper->pair_ext)
                      tmp *= sc_mb_wrapper->pair_ext(i, j, sc_mb_wrapper);

                    /* store outside part for base pair (i,j) */
                    tmp2 += q_g *
                            tmp;
                  }
                }

                /* 3.2.2 (i,j) first branch after gquad, ie. [k,l] + (i,j) + ML */
                if (j + turn < k) {
                  u2 = i - l - 1;
                  if (hc->up_ml[l + 1] >= u2) {
                    tmp = expGQstem *
                          qm[my_iindx[j + 1] - k + 1] *
                          expGQstem *
                          expMLbase[u2];

                    if (sc_mb_wrapper->pair_ext)
                      tmp *= sc_mb_wrapper->pair_ext(i, j, sc_mb_wrapper);

                    /* store outside part for base pair (i,j) */
                    tmp2 += q_g *
                            tmp;
                  }
                }

                /* 3.2.3 (i,j) somewhere in between, ie. [k,l] + ML + (i,j) + ML */
                if ((l + turn < i) &&
                    (j + turn < k)) {
                  tmp = expGQstem *
                        qm[my_iindx[l + 1] - i + 1] *
                        qm[my_iindx[j + 1] - k + 1] *
                        expGQstem;

                  if (sc_mb_wrapper->pair_ext)
                    tmp *= sc_mb_wrapper->pair_ext(i, j, sc_mb_wrapper);

                  /* store outside part for base pair (i,j) */
                  tmp2 += q_g *
                          tmp;
                }
              }
            }
          }

          /* store total contribution */
          probs[ij] += tmp2 / qo;
        }
      }
    }

    /* finally, outside probabilities for all Gquads spanning the n,1 junction */
    if (VRNA_GQUAD_MAX_BOX_SIZE > j) {
      if (n <= VRNA_GQUAD_MAX_BOX_SIZE + 3) {
        min_k = j + 1;
      } else {
        min_k = n - (VRNA_GQUAD_MAX_BOX_SIZE - j);
        if (min_k <= j)
          min_k = j + 1;
      }

      max_k = n;
      if (VRNA_GQUAD_MIN_BOX_SIZE > j)
        max_k = n + j + 1 - VRNA_GQUAD_MIN_BOX_SIZE;

      for (i = min_k; i <= max_k; i++) {
#ifndef VRNA_DISABLE_C11_FEATURES
        if ((q_g = vrna_smx_csr_get(q_gq, i, j, 0.)) != 0.) {
#else
        if ((q_g = vrna_smx_csr_FLT_OR_DBL_get(q_gq, i, j, 0.)) != 0.) {
#endif
          tmp2 = 0;

          /* 1. hairpin-loop like case, k.e. single gquad, rest unpaired */
          u1    = i - j - 1;
          eval  = (hc->up_hp[j + 1] >= u1) ? 1 : 0;

          if ((eval) &&
              (u1 >= 3)) {
            qbt1 = scale[u1];

            if (md->circ_penalty)
              qbt1 *= pow(vrna_exp_E_hairpin(u1, 0, -1, -1, NULL, pf_params), (double)n_seq);

            if (sc_hp_wrapper->pair_ext)
              qbt1 *= sc_hp_wrapper->pair_ext(i, j, sc_hp_wrapper);

            tmp2 += qbt1;
          }

          /* 2. internal-loop like case, k.e. single gquad plus base pair or second gquad */
          min_i = j + 1;
          max_i = j + MAXLOOP + 1;
          if (max_i >= i)
            max_i = i - 1;

          for (k = min_i; k <= max_i; k++) {
            u1    = k - j - 1;
            max_j = i - 1;
            min_j = k + turn + 1;
            if ((MAXLOOP <= n) &&
                (min_j + MAXLOOP < i + u1))
              min_j = i + u1 - MAXLOOP - 1;

            for (l = max_j; l >= min_j; l--) {
              u2 = i - l - 1;
              if (((u1 == 0) && (u2 < 3)) ||
                  ((u1 < 3) && (u2 == 0)))
                continue;

              eval  = (hc->up_int[j + 1] >= u1) ? 1 : 0;
              eval  = (hc->up_int[l + 1] >= u2) ? eval : 0;
              /* 2.1 second element is a base pair */
              if ((eval) &&
                  (hc->mx[n * k + l] &
                   (VRNA_CONSTRAINT_CONTEXT_INT_LOOP | VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC)) &&
                  (qb[my_iindx[k] - l] > 0.)) {
                qbt1  = scale[u1 + u2];
                qint  = 1.;

                switch (fc->type) {
                  case VRNA_FC_TYPE_SINGLE:
                    sj  = S1[l + 1];
                    si  = S1[k - 1];

                    type = vrna_get_ptype_md(S[l], S[k], md);
                    if (md->dangles == 2)
                      qint *= pf_params->expmismatchI[type][sj][si];

                    if (type > 2)
                      qint *= pf_params->expTermAU;

                    qbt1 *= (FLT_OR_DBL)expintern[u1 + u2];
                    break;

                  case VRNA_FC_TYPE_COMPARATIVE:
                    for (s = 0; s < n_seq; s++) {
                      type = vrna_get_ptype_md(SS[s][l], SS[s][k], md);
                      if (md->dangles == 2)
                        qbt1 *= pf_params->expmismatchI[type][S3[s][l]][S5[s][k]];

                      if (type > 2)
                        qbt1 *= pf_params->expTermAU;

                      us1   = a2s[s][k - 1] - a2s[s][j];
                      us2   = a2s[s][i - 1] - a2s[s][l];
                      qbt1  *= (FLT_OR_DBL)expintern[us1 + us2];
                    }
                    break;
                }

                if (sc_int_wrapper->pair_ext)
                  qbt1 *= sc_int_wrapper->pair_ext(k, l, i, j, sc_int_wrapper);

                /* store quad probs */
                tmp2 += qbt1 *
                        qint *
                        qb[my_iindx[k] - l];
              }

              /* 2.2 second element is a gquad */
              if ((eval) &&
                  (l - k + 1 >= VRNA_GQUAD_MIN_BOX_SIZE) &&
                  (l - k + 1 <= VRNA_GQUAD_MAX_BOX_SIZE) &&
#ifndef VRNA_DISABLE_C11_FEATURES
                  ((tmp3 = vrna_smx_csr_get(q_gq, k, l, 0.)) != 0.)) {
#else
                  ((tmp3 = vrna_smx_csr_FLT_OR_DBL_get(q_gq, k, l, 0.)) != 0.)) {
#endif
                qbt1 = scale[u1 + u2];

                switch (fc->type) {
                  case VRNA_FC_TYPE_SINGLE:
                    qbt1 *= (FLT_OR_DBL)expintern[u1 + u2];
                    break;

                  case VRNA_FC_TYPE_COMPARATIVE:
                    for (s = 0; s < n_seq; s++) {
                      us1   = a2s[s][k - 1] - a2s[s][j];
                      us2   = a2s[s][i - 1] - a2s[s][l];
                      qbt1  *= (FLT_OR_DBL)expintern[us1 + us2];
                    }
                    break;
                }

                if (sc_int_wrapper->pair_ext)
                  qbt1 *= sc_int_wrapper->pair_ext(k, l, i, j, sc_int_wrapper);

                /* store quad probs */
                tmp2 += tmp3 *
                        qbt1;
              }
            }
          }

          /* 3. multibranch-loop like case, k.e. gquad plus at least two more stems */
          qbt1 = qm2_real[my_iindx[j + 1] - i + 1] *
                 expMLclosing *
                 expGQstem;

          tmp2 += qbt1;

          /* 4. store gquad outside probs for all gquads spanning the n,1 junction */
#ifndef VRNA_DISABLE_C11_FEATURES
          vrna_smx_csr_insert(matrices->p_gq, i, j, tmp2 * q_g / qo);

#else
          vrna_smx_csr_FLT_OR_DBL_insert(matrices->p_gq, i, j, tmp2 * q_g / qo);
#endif
        } /* end if gquad [i:j] */
      }   /* end for(i = ... */
    }     /* end gquad spanning n1-junction */
  }       /* end with_gquad == 1 */

  free(tt);
}
