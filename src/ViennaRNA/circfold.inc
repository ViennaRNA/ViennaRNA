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
  char              *ptype, *hard_constraints;
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

  if(hc->up_ext[1] >= length){
    Fc = 0;
    if(sc){
      if(sc->energy_up)
        Fc += sc->energy_up[1][length];
    }
  } else {
    Fc = INF;
  }

  for(i = 1; i < length; i++)
    for(j = i + turn + 1; j <= length; j++){
      u = length-j + i-1;
      if (u<TURN) continue;

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
  for (i=1; i<length-TURN; i++) {
    fM2[i] = INF;
    for (u=i+TURN; u<length-TURN; u++)
      fM2[i] = MIN2(fM2[i], my_fML[indx[u]+i] + my_fML[indx[length]+u+1]);
  }

  for (i=TURN+1; i<length-2*TURN; i++) {
    fm = my_fML[indx[i]+1]+fM2[i+1]+P->MLclosing;
    if (fm<FcM) {
      FcM=fm; Mi=i;
    }
  }
  Fc = MIN2(Fc, FcM);

  if ((dangle_model == 1) || (dangle_model == 3)){
    fM_d3 =  (int *) vrna_alloc(sizeof(int)*(length+2));
    fM_d5 =  (int *) vrna_alloc(sizeof(int)*(length+2));

    for (i=TURN+1; i<length-TURN; i++) {
      fM_d3[i] = INF;
      for (u=2+TURN; u<i-TURN; u++)
        fM_d3[i] = MIN2(fM_d3[i], my_fML[indx[u]+2] + my_fML[indx[i]+u+1]);
    }

    for (i=2*TURN+1; i<length-TURN; i++) {
      type = ptype[indx[length]+i+1];
      if(hard_constraints[indx[length] + i + 1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP){
        if(hc->up_ml[1]){
          fm = fM_d3[i]+my_c[indx[length]+i+1]+E_MLstem(type, -1, S1[1], P) + P->MLclosing;
          if (fm<FcMd3) {
            FcMd3=fm; Md3i=i;
          }

          if(hc->up_ml[i]){
            fm = fM_d3[i-1]+my_c[indx[length]+i+1]+E_MLstem(type, S1[i], S1[1], P) + P->MLclosing;
            if (fm<FcMd3) {
              FcMd3=fm; Md3i=-i;
            }
          }
        }
      }
    }

    for (i=TURN+1; i<length-TURN; i++) {
      fM_d5[i] = INF;
      for (u=i+TURN; u<length-TURN; u++)
        fM_d5[i] = MIN2(fM_d5[i], my_fML[indx[u]+i] + my_fML[indx[length-1]+u+1]);
    }

    for (i=TURN+1; i<length-2*TURN; i++) {
      if(hard_constraints[indx[i]+1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP){
        type = ptype[indx[i]+1];
        if(hc->up_ml[length]){
          fm = E_MLstem(type, S1[length], -1, P) + my_c[indx[i]+1] + fM_d5[i+1] + P->MLclosing;
          if (fm<FcMd5) {
            FcMd5=fm; Md5i=i;
          }
          if(hc->up_ml[i+1]){
            fm = E_MLstem(type, S1[length], S1[i+1], P) + my_c[indx[i]+1] + fM_d5[i+2] + P->MLclosing;
            if (fm<FcMd5) {
              FcMd5=fm; Md5i=-i;
            }
          }
        }
      }
    }

    if (FcMd5<MIN2(Fc,FcMd3)) {
      /* looks like we have to do this ... */
      bt_stack[++(*bt)].i = 1;
      bt_stack[(*bt)].j = (Md5i>0)?Md5i:-Md5i;
      bt_stack[(*bt)].ml = 2;
      i = (Md5i>0)?Md5i+1 : -Md5i+2; /* let's backtrack fm_d5[Md5i+1] */
      for (u=i+TURN; u<length-TURN; u++)
        if (fM_d5[i] == my_fML[indx[u]+i] + my_fML[indx[length-1]+u+1]) {
          bt_stack[++(*bt)].i = i;
          bt_stack[(*bt)].j = u;
          bt_stack[(*bt)].ml = 1;
          bt_stack[++(*bt)].i =u+1;
          bt_stack[(*bt)].j = length-1;
          bt_stack[(*bt)].ml = 1;
          break;
        }
      Fc = FcMd5;
    } else if (FcMd3<Fc) {
      /* here we go again... */
      bt_stack[++(*bt)].i = (Md3i>0)?Md3i+1:-Md3i+1;
      bt_stack[(*bt)].j = length;
      bt_stack[(*bt)].ml = 2;
      i = (Md3i>0)? Md3i : -Md3i-1; /* let's backtrack fm_d3[Md3i] */
      for (u=2+TURN; u<i-TURN; u++)
        if (fM_d3[i] == my_fML[indx[u]+2] + my_fML[indx[i]+u+1]) {
          bt_stack[++(*bt)].i = 2;
          bt_stack[(*bt)].j = u;
          bt_stack[(*bt)].ml = 1;
          bt_stack[++(*bt)].i =u+1;
          bt_stack[(*bt)].j = i;
          bt_stack[(*bt)].ml = 1;
          break;
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
      /* backtrack in fM2 */
      fm = fM2[Mi+1];
      for (u=Mi+TURN+1; u<length-TURN; u++)
        if (fm == my_fML[indx[u]+Mi+1] + my_fML[indx[length]+u+1]) {
                bt_stack[++(*bt)].i=Mi+1;
                bt_stack[(*bt)].j=u;
                bt_stack[(*bt)].ml = 1;
                bt_stack[++(*bt)].i=u+1;
                bt_stack[(*bt)].j=length;
                bt_stack[(*bt)].ml = 1;
                break;
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

