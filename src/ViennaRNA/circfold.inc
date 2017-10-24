/* -*-C-*- */
/* this file contains code for folding circular RNAs */
/* it's #include'd into fold.c */

PRIVATE void
fill_arrays_circ( vrna_fold_compound_t *vc,
                  sect bt_stack[],
                  int *bt){

  /* variant of fold() for circular RNAs */
  /* auxiliarry arrays:
     fM2 = multiloop region with exactly two stems, extending to 3' end
     for stupid dangles=1 case we also need:
     fM_d3 = multiloop region with >= 2 stems, starting at pos 2
             (a pair (k,n) will form 3' dangle with pos 1)
     fM_d5 = multiloop region with >= 2 stems, extending to pos n-1
             (a pair (1,k) will form a 5' dangle with pos n)
  */
  int               Hi, Hj, Ii, Ij, Ip, Iq, ip, iq, Mi;
  int               *fM_d3, *fM_d5, Md3i, Md5i, FcMd3, FcMd5;
  int               FcH, FcI, FcM, Fc, *fM2;
  int               i,j, ij, u, length, new_c, fm, type;
  int               *my_c, *my_fML, *indx, dangle_model, turn;
  vrna_param_t      *P;
  char              *ptype;
  unsigned char     *hard_constraints, eval;
  short             *S1;
  vrna_hc_t         *hc;
  vrna_sc_t         *sc;

  P                 = vc->params;
  ptype             = vc->ptype;
  indx              = vc->jindx;
  S1                = vc->sequence_encoding;

  dangle_model      = P->model_details.dangles;
  turn              = P->model_details.min_loop_size;

  hc                = vc->hc;
  sc                = vc->sc;
  hard_constraints  = hc->matrix;

  my_c              = vc->matrices->c;
  my_fML            = vc->matrices->fML;
  fM2               = vc->matrices->fM2;

  length            = vc->length;

  FcH = FcI = FcM = FcMd3 = FcMd5 = INF;

  /* unfolded state */
  eval = (hc->up_ext[1] >= length) ? 1 : 0;
  if (hc->f)
    eval = (hc->f(1, length, 1, length, VRNA_DECOMP_EXT_UP, hc->data)) ? eval : 0;

  if (eval) {
    Fc = 0;
    if(sc){
      if(sc->energy_up)
        Fc += sc->energy_up[1][length];

      if (sc->f)
        Fc += sc->f(1, length, 1, length, VRNA_DECOMP_EXT_UP, sc->data);
    }
  } else {
    Fc = INF;
  }

  for(i = 1; i < length; i++)
    for(j = i + turn + 1; j <= length; j++){
      u = length-j + i-1;
      if (u<turn) continue;

      ij = indx[j] + i;

      if (!hard_constraints[ij])
        continue;

      /* exterior hairpin case */
      new_c =   vrna_E_hp_loop(vc, j, i)
              + my_c[ij];

      if (new_c<FcH) {
        FcH = new_c;
        Hi  = i;
        Hj  = j;
      }

      /* exterior interior loop case */
      ip = iq = 0;
      new_c =   vrna_E_ext_int_loop(vc, i, j, &ip, &iq)
              + my_c[ij];

      if(ip != 0){
        if(new_c < FcI){
          FcI = new_c;
          Ii  = i;
          Ij  = j;
          Ip  = ip;
          Iq  = iq;
        }
      }
    } /* end of i,j loop */
  Fc = MIN2(Fc, FcH);
  Fc = MIN2(Fc, FcI);
  /* compute the fM2 array (multi loops with exactly 2 helices) */
  /* to get a unique ML decomposition, just use fM1 instead of fML
     below. However, that will not work with dangle_model==1  */
  if ((sc) && (sc->f)) {
    if (hc->f) {
      for (i=1; i<length-turn; i++) {
        fM2[i] = INF;
        for (u=i+turn; u<length-turn; u++) {
          if (hc->f(i, length, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
            fm = my_fML[indx[u]+i] +
                 my_fML[indx[length]+u+1] +
                 sc->f(i, length, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
            fM2[i] = MIN2(fM2[i], fm);
          }
        }
      }
    } else {
      for (i=1; i<length-turn; i++) {
        fM2[i] = INF;
        for (u=i+turn; u<length-turn; u++) {
          fm = my_fML[indx[u]+i] +
               my_fML[indx[length]+u+1] +
               sc->f(i, length, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
          fM2[i] = MIN2(fM2[i], fm);
        }
      }
    }
  } else {
    if (hc->f) {
      for (i=1; i<length-turn; i++) {
        fM2[i] = INF;
        for (u=i+turn; u<length-turn; u++) {
          if (hc->f(i, length, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
            fm = my_fML[indx[u]+i] +
                 my_fML[indx[length]+u+1];
            fM2[i] = MIN2(fM2[i], fm);
          }
        }
      }
    } else {
      for (i=1; i<length-turn; i++) {
        fM2[i] = INF;
        for (u=i+turn; u<length-turn; u++) {
          fm = my_fML[indx[u]+i] +
               my_fML[indx[length]+u+1];
          fM2[i] = MIN2(fM2[i], fm);
        }
      }
    }
  }

  if ((sc) && (sc->f)) {
    if (hc->f) {
      for (i=turn+1; i<length-2*turn; i++) {
        if (hc->f(1, length, i, i + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
          fm = my_fML[indx[i]+1] +
               fM2[i+1] +
               P->MLclosing +
               sc->f(1, length, i, i + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
          if (fm<FcM) {
            FcM=fm; Mi=i;
          }
        }
      }
    } else {
      for (i=turn+1; i<length-2*turn; i++) {
        fm = my_fML[indx[i]+1] +
             fM2[i+1] +
             P->MLclosing +
             sc->f(1, length, i, i + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
        if (fm<FcM) {
          FcM=fm; Mi=i;
        }
      }
    }
  } else {
    if (hc->f) {
      for (i=turn+1; i<length-2*turn; i++) {
        if (hc->f(1, length, i, i + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
          fm = my_fML[indx[i]+1]+fM2[i+1]+P->MLclosing;
          if (fm<FcM) {
            FcM=fm; Mi=i;
          }
        }
      }
    } else {
      for (i=turn+1; i<length-2*turn; i++) {
        fm = my_fML[indx[i]+1]+fM2[i+1]+P->MLclosing;
        if (fm<FcM) {
          FcM=fm; Mi=i;
        }
      }
    }
  }
  Fc = MIN2(Fc, FcM);

  if ((dangle_model == 1) || (dangle_model == 3)){
    fM_d3 =  (int *) vrna_alloc(sizeof(int)*(length+2));
    fM_d5 =  (int *) vrna_alloc(sizeof(int)*(length+2));

    for (i=turn+1; i<length-turn; i++) {
      fM_d3[i] = INF;
      fM_d5[i] = INF;
    }


    if (hc->f) {
      if (hc->up_ml[1]) {
        for (i=turn+1; i<length-turn; i++) {
          if (hc->f(1, i, 2, i, VRNA_DECOMP_ML_ML, hc->data)) {
            for (u=2+turn; u<i-turn; u++) {
              if (hc->f(2, i, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
                fm = my_fML[indx[u]+2] +
                     my_fML[indx[i]+u+1];

                if (sc) {
                  if (sc->energy_up)
                    fm += sc->energy_up[1][1];

                  if (sc->f)
                    fm += sc->f(1, i, 2, i, VRNA_DECOMP_ML_ML, sc->data) +
                          sc->f(2, i, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
                }

                fM_d3[i] = MIN2(fM_d3[i], fm);
              }
            }
          }
        }

        for (i=2*turn+1; i<length-turn; i++) {
          type = ptype[indx[length]+i+1];
          eval = (hard_constraints[indx[length] + i + 1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) ? 1 : 0;
          eval = hc->f(i + 1, length, 2, i, VRNA_DECOMP_PAIR_ML_EXT, hc->data) ? eval : 0;
          if (eval) {
            fm = fM_d3[i] +
                 my_c[indx[length]+i+1]+
                 E_MLstem(type, -1, S1[1], P) +
                 P->MLclosing;

            if (sc) {
              if (sc->f)
                fm += sc->f(i+1, length, 2, i, VRNA_DECOMP_PAIR_ML_EXT, sc->data);
            }

            if (fm<FcMd3) {
              FcMd3=fm; Md3i=i;
            }
          }

          eval = (hard_constraints[indx[length] + i + 1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) ? 1 : 0;
          eval = hc->up_ml[i] ? eval : 0;
          eval = hc->f(i + 1, length, 2, i - 1, VRNA_DECOMP_PAIR_ML_EXT, hc->data) ? eval : 0;
          if (eval) {
            fm = fM_d3[i-1] +
                 my_c[indx[length]+i+1]+
                 E_MLstem(type, S1[i], S1[1], P) +
                 P->MLclosing;

            if (sc) {
              if (sc->energy_up)
                fm += sc->energy_up[i][1];

              if (sc->f)
                fm += sc->f(i + 1, length, 2, i - 1, VRNA_DECOMP_PAIR_ML_EXT, sc->data);
            }

            if (fm<FcMd3) {
              FcMd3=fm; Md3i=-i;
            }
          }
        }
      }
    } else {
      if (hc->up_ml[1]) {
        for (i=turn+1; i<length-turn; i++) {
          for (u=2+turn; u<i-turn; u++) {
            fm = my_fML[indx[u]+2] +
                 my_fML[indx[i]+u+1];

            if (sc) {
              if (sc->energy_up)
                fm += sc->energy_up[1][1];

              if (sc->f)
                fm += sc->f(1, i, 2, i, VRNA_DECOMP_ML_ML, sc->data) +
                      sc->f(2, i, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
            }

            fM_d3[i] = MIN2(fM_d3[i], fm);
          }
        }

        for (i=2*turn+1; i<length-turn; i++) {
          type = ptype[indx[length]+i+1];
          if(hard_constraints[indx[length] + i + 1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP){
            fm = fM_d3[i]+
                  my_c[indx[length]+i+1]+
                  E_MLstem(type, -1, S1[1], P) +
                  P->MLclosing;

            if (sc) {
              if (sc->f)
                fm += sc->f(i + 1, length, 2, i, VRNA_DECOMP_PAIR_ML_EXT, sc->data);
            }

            if (fm<FcMd3) {
              FcMd3=fm; Md3i=i;
            }

            if(hc->up_ml[i]){
              fm = fM_d3[i-1]+
                    my_c[indx[length]+i+1]+
                    E_MLstem(type, S1[i], S1[1], P) +
                    P->MLclosing;

              if (sc) {
                if (sc->energy_up)
                  fm += sc->energy_up[i][1];

                if(sc->f)
                  fm += sc->f(i + 1, length, 2, i - 1, VRNA_DECOMP_PAIR_ML_EXT, sc->data);
              }

              if (fm<FcMd3) {
                FcMd3=fm; Md3i=-i;
              }
            }
          }
        }
      }
    }

    if (hc->f) {
      if (hc->up_ml[length]) {
        for (i=turn+1; i<length-turn; i++) {
          if (hc->f(i, length, i, length - 1, VRNA_DECOMP_ML_ML, hc->data)) {
            for (u=i+turn; u<length-turn; u++) {
              if (hc->f(i, length - 1, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
                fm = my_fML[indx[u]+i] +
                     my_fML[indx[length-1]+u+1];

                if (sc) {
                  if (sc->energy_up)
                    fm += sc->energy_up[length][1];

                  if (sc->f)
                    fm += sc->f(i, length, i, length - 1, VRNA_DECOMP_ML_ML, sc->data) +
                          sc->f(i, length - 1, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
                }

                fM_d5[i] = MIN2(fM_d5[i], fm);
              }
            }
          }
        }

        for (i=turn+1; i<length-2*turn; i++) {
          type = ptype[indx[i]+1];
          eval = (hard_constraints[indx[i]+1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) ? 1 : 0;
          eval = hc->f(1, i, i + 1, length - 1, VRNA_DECOMP_PAIR_ML_EXT, hc->data) ? eval : 0;
          if(eval){
            fm = E_MLstem(type, S1[length], -1, P) +
                  my_c[indx[i]+1] +
                  fM_d5[i+1] +
                  P->MLclosing;

            if (sc) {
              if (sc->f)
                fm += sc->f(1, i, i + 1, length - 1, VRNA_DECOMP_PAIR_ML_EXT, sc->data);
            }

            if (fm<FcMd5) {
              FcMd5=fm; Md5i=i;
            }
          }

          eval = (hard_constraints[indx[i]+1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) ? 1 : 0;
          eval = hc->f(1, i, i + 1, length - 1, VRNA_DECOMP_PAIR_ML_EXT, hc->data) ? eval : 0;
          eval = hc->up_ml[i+1] ? eval : 0;
          if (eval) {
            fm = E_MLstem(type, S1[length], S1[i+1], P) +
                  my_c[indx[i]+1] +
                  fM_d5[i+2] +
                  P->MLclosing;

            if (sc) {
              if (sc->energy_up)
                fm += sc->energy_up[i + 1][1];

              if (sc->f)
                fm += sc->f(i, length - 1, i + 1, length - 1, VRNA_DECOMP_ML_ML, sc->data) +
                      sc->f(1, i, i + 2, length - 1, VRNA_DECOMP_PAIR_ML_EXT, sc->data);
            }

            if (fm<FcMd5) {
              FcMd5=fm; Md5i=-i;
            }
          }
        }
      }
    } else {
      if (hc->up_ml[length]) {
        for (i=turn+1; i<length-turn; i++) {
          for (u=i+turn; u<length-turn; u++) {
            fm = my_fML[indx[u]+i] +
                  my_fML[indx[length-1]+u+1];

            if (sc) {
              if (sc->energy_up)
                fm += sc->energy_up[length][1];

              if (sc->f)
                fm += sc->f(i, length, i, length - 1, VRNA_DECOMP_ML_ML, sc->data) +
                      sc->f(i, length - 1, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
            }

            fM_d5[i] = MIN2(fM_d5[i], fm);
          }
        }

        for (i=turn+1; i<length-2*turn; i++) {
          if(hard_constraints[indx[i]+1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP){
            type = ptype[indx[i]+1];
            fm = E_MLstem(type, S1[length], -1, P) +
                  my_c[indx[i]+1] +
                  fM_d5[i+1] +
                  P->MLclosing;

            if (sc) {
              if (sc->f)
                fm += sc->f(1, i, i + 1, length - 1, VRNA_DECOMP_PAIR_ML_EXT, sc->data);
            }

            if (fm<FcMd5) {
              FcMd5=fm; Md5i=i;
            }

            if(hc->up_ml[i+1]){
              fm = E_MLstem(type, S1[length], S1[i+1], P) +
                    my_c[indx[i]+1] +
                    fM_d5[i+2] +
                    P->MLclosing;

              if (sc) {
                if (sc->energy_up)
                  fm += sc->energy_up[i + 1][1];

                if (sc->f)
                  fm += sc->f(i, length - 1, i + 1, length - 1, VRNA_DECOMP_ML_ML, sc->data) +
                        sc->f(1, i, i + 2, length - 1, VRNA_DECOMP_PAIR_ML_EXT, sc->data);
              }

              if (fm<FcMd5) {
                FcMd5=fm; Md5i=-i;
              }
            }
          }
        }
      }
    }

    if (FcMd5<MIN2(Fc,FcMd3)) {
      int real_i, sc_en = 0;

      /* looks like we have to do this ... */
      bt_stack[++(*bt)].i = 1;
      bt_stack[(*bt)].j = (Md5i>0)?Md5i:-Md5i;
      bt_stack[(*bt)].ml = 2;
      i = (Md5i>0)?Md5i+1 : -Md5i+2; /* let's backtrack fm_d5[Md5i+1] */
      real_i = (Md5i > 0) ? i : i - 1;

      if ((sc) && (sc->energy_up))
        sc_en += sc->energy_up[length][1];

      for (u=i+turn; u<length-turn; u++) {
        fm = my_fML[indx[u]+i] +
             my_fML[indx[length-1]+u+1] +
             sc_en;

        if (sc) {
          if (sc->energy_up)
            fm += sc->energy_up[real_i][i - real_i];

          if (sc->f)
            fm += sc->f(real_i, length, i, length - 1, VRNA_DECOMP_ML_ML, sc->data) +
                  sc->f(i, length - 1, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
        }

        if (fM_d5[i] == fm) {
          bt_stack[++(*bt)].i = i;
          bt_stack[(*bt)].j = u;
          bt_stack[(*bt)].ml = 1;
          bt_stack[++(*bt)].i =u+1;
          bt_stack[(*bt)].j = length-1;
          bt_stack[(*bt)].ml = 1;
          break;
        }
      }
      Fc = FcMd5;
    } else if (FcMd3<Fc) {
      int real_i, sc_en = 0;
      /* here we go again... */
      bt_stack[++(*bt)].i = (Md3i>0)?Md3i+1:-Md3i+1;
      bt_stack[(*bt)].j = length;
      bt_stack[(*bt)].ml = 2;
      i = (Md3i>0)? Md3i : -Md3i-1; /* let's backtrack fm_d3[Md3i] */
      real_i = (Md3i > 0) ? i : i + 1;

      if ((sc) && (sc->energy_up))
        sc_en += sc->energy_up[1][1];

      for (u=2+turn; u<i-turn; u++) {
        fm = my_fML[indx[u]+2] +
             my_fML[indx[i]+u+1] +
             sc_en;

        if (sc) {
          if (sc->energy_up)
            fm += sc->energy_up[real_i][real_i - i];

          if (sc->f)
            fm += sc->f(1, real_i, 2, i, VRNA_DECOMP_ML_ML, sc->data) +
                  sc->f(2, i, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
        }

        if (fM_d3[i] == fm) {
          bt_stack[++(*bt)].i = 2;
          bt_stack[(*bt)].j = u;
          bt_stack[(*bt)].ml = 1;
          bt_stack[++(*bt)].i =u+1;
          bt_stack[(*bt)].j = i;
          bt_stack[(*bt)].ml = 1;
          break;
        }
      }
      Fc = FcMd3;
    }
    free(fM_d3);
    free(fM_d5);
  }
  else if(Fc < 0){
    if (FcH==Fc) {
      bt_stack[++(*bt)].i = Hi;
      bt_stack[(*bt)].j = Hj;
      bt_stack[(*bt)].ml = 2;
    }
    else if (FcI==Fc) {
      bt_stack[++(*bt)].i = Ii;
      bt_stack[(*bt)].j = Ij;
      bt_stack[(*bt)].ml = 2;
      bt_stack[++(*bt)].i = Ip;
      bt_stack[(*bt)].j = Iq;
      bt_stack[(*bt)].ml = 2;
    }
    else if (FcM==Fc) { /* grumpf we found a Multiloop */
      int eee;
      /* backtrack in fM2 */
      fm = fM2[Mi+1];
      for (u=Mi+turn+1; u<length-turn; u++) {
        eee = my_fML[indx[u]+Mi+1] +
              my_fML[indx[length]+u+1];

        if (sc) {
          if (sc->f)
            eee += sc->f(Mi + 1, length, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
        }

        if (fm == eee) {
                bt_stack[++(*bt)].i=Mi+1;
                bt_stack[(*bt)].j=u;
                bt_stack[(*bt)].ml = 1;
                bt_stack[++(*bt)].i=u+1;
                bt_stack[(*bt)].j=length;
                bt_stack[(*bt)].ml = 1;
                break;
        }
      }
      bt_stack[++(*bt)].i = 1;
      bt_stack[(*bt)].j = Mi;
      bt_stack[(*bt)].ml = 1;
    }
  } else { /* unstructured */
    bt_stack[++(*bt)].i = 1;
    bt_stack[(*bt)].j = 1;
    bt_stack[(*bt)].ml = 0;
  }
  vc->matrices->FcH = FcH;
  vc->matrices->FcI = FcI;
  vc->matrices->FcM = FcM;
  vc->matrices->Fc  = Fc;
}

