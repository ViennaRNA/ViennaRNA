/* -*-C-*- */
/* this file contains code for alifolding circular RNAs */
/* it's #include'd into alifold.c */

PRIVATE void
fill_arrays_comparative_circ( vrna_fold_compound_t *vc,
                              sect bt_stack[],
                              int *bt){

  /* variant of alifold() for circular RNAs */
  /* auxiliarry arrays:
     fM2 = multiloop region with exactly two stems, extending to 3' end
     for stupid dangles=1 case we also need:
     fM_d3 = multiloop region with >= 2 stems, starting at pos 2
             (a pair (k,n) will form 3' dangle with pos 1)
     fM_d5 = multiloop region with >= 2 stems, extending to pos n-1
             (a pair (1,k) will form a 5' dangle with pos n)
  */
  int               Fc, FcH, FcI, FcM, Hi, Hj, Ii, Ij, Ip, Iq, ip, iq, Mi, FcMd3, FcMd5;
  int               i, j, ij, u, length, new_c, fm, n_seq, s;
  int               *indx, *c, *fML, *fM2;
  char              *hard_constraints;
  unsigned short    **a2s;
  vrna_param_t      *P;
  vrna_hc_t         *hc;
  vrna_sc_t         **sc;

  a2s               = vc->a2s;
  P                 = vc->params;
  indx              = vc->jindx;     /* index for moving in the triangle matrices c[] and fMl[]*/
  c                 = vc->matrices->c;     /* energy array, given that i-j pair */
  fML               = vc->matrices->fML;     /* multi-loop auxiliary energy array */
  fM2               = vc->matrices->fM2;
  hc                = vc->hc;
  sc                = vc->scs;
  hard_constraints  = hc->matrix;
  n_seq             = vc->n_seq;
  length            = vc->length;

  /* extra arrays for circfold() */
  FcH = FcI= FcM = FcMd3= FcMd5= INF;

  if(hc->up_ext[1] >= length){
    Fc = 0;
    if(sc){
      for(s = 0; s < n_seq; s++){
        if(sc[s]->energy_up)
          Fc += sc[s]->energy_up[1][a2s[s][length] - 1];
      }
    }
  } else {
    Fc = INF;
  }

  for (i=1; i<length; i++)
    for (j=i+TURN+1; j <= length; j++) {
      u = length-j + i-1;
      if (u<TURN) continue;

      ij = indx[j]+i;

      if (!hard_constraints[ij])
        continue;

      /* exterior hairpin case */
      new_c =   vrna_E_hp_loop(vc, j, i)
              + c[ij];

      if (new_c<FcH) {
        FcH = new_c;
        Hi  = i;
        Hj  = j;
      }

      /* exterior interior loop case */
      ip = iq = 0;
      new_c =   vrna_E_ext_int_loop(vc, i, j, &ip, &iq)
              + c[ij];

      if(ip != 0){
        if(new_c < FcI){
          FcI = new_c;
          Ii  = i;
          Ij  = j;
          Ip  = ip;
          Iq  = iq;
        }
      }
    }
  Fc = MIN2(Fc, FcH);
  Fc = MIN2(Fc, FcI);

  /* compute the fM2 array (multi loops with exactly 2 helices) */
  /* to get a unique ML decomposition, just use fM1 instead of fML
     below. However, that will not work with dangles==1  */
  for (i=1; i<length-TURN; i++) {
    fM2[i] = INF;
    for (u=i+TURN; u<length-TURN; u++)
      fM2[i] = MIN2(fM2[i], fML[indx[u]+i] + fML[indx[length]+u+1]);
  }

  for (i=TURN+1; i<length-2*TURN; i++) {
    fm = fML[indx[i]+1]+fM2[i+1]+n_seq*P->MLclosing;
    if (fm<FcM) {
      FcM=fm; Mi=i;
    }
  }
  Fc = MIN2(Fc, FcM);

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
      if (fm == fML[indx[u]+Mi+1] + fML[indx[length]+u+1]) {
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
  } else { /* must be totally unstructured */
    bt_stack[++(*bt)].i = 1;
    bt_stack[(*bt)].j = 1;
    bt_stack[(*bt)].ml = 0;
  }
  vc->matrices->FcH = FcH;
  vc->matrices->FcI = FcI;
  vc->matrices->FcM = FcM;
  vc->matrices->Fc  = Fc;
}
