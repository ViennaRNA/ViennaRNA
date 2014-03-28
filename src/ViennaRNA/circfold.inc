/* -*-C-*- */
/* this file contains code for folding circular RNAs */
/* it's #include'd into fold.c */

PRIVATE void fill_arrays_circ(vrna_fold_compound *vc, sect bt_stack[], int *bt){
  /* variant of fold() for circular RNAs */
  /* auxiliarry arrays:
     fM2 = multiloop region with exactly two stems, extending to 3' end
     for stupid dangles=1 case we also need:
     fM_d3 = multiloop region with >= 2 stems, starting at pos 2
             (a pair (k,n) will form 3' dangle with pos 1)
     fM_d5 = multiloop region with >= 2 stems, extending to pos n-1
             (a pair (1,k) will form a 5' dangle with pos n)
  */
  int Hi, Hj, Ii, Ij, Ip, Iq, Mi;
  int *fM_d3, *fM_d5, Md3i, Md5i, FcMd3, FcMd5;
  int FcH, FcI, FcM, Fc, *fM2;
  int i,j, p,q,length, energy;

  /* the folding matrices */
  int   *my_c, *my_fML;

  paramT  *P          = vc->params;
  char *string        = vc->sequence;
  char *ptype         = vc->ptype;
  int   *indx         = vc->jindx;
  short *S1           = vc->sequence_encoding;

  int   dangle_model  = P->model_details.dangles;
  int   noGUclosure   = P->model_details.noGUclosure;
  int   *rtype        = &(P->model_details.rtype[0]);

  /* dereference of the folding matrices */
  my_c    = vc->matrices->c;
  my_fML  = vc->matrices->fML;
  fM2     = vc->matrices->fM2;

  length = (int) strlen(string);

  FcH = FcI= FcM = FcMd3= FcMd5= Fc = INF;
  for (i=1; i<length; i++)
    for (j=i+TURN+1; j <= length; j++) {
      int ij, type, u, new_c, no_close;
      u = length-j + i-1;
      if (u<TURN) continue;

      ij = indx[j]+i;
      type = ptype[ij];

      no_close = (((type==3)||(type==4))&&noGUclosure);

      /* if (j-i-1 > max_separation) type = 0; */  /* forces locality degree */

      type=rtype[type];
      if (!type) continue;
      if (no_close) new_c = FORBIDDEN;
      else {
        char loopseq[10];
        /*int si1, sj1;*/
        if (u<7) {
          strcpy(loopseq , string+j-1);
          strncat(loopseq, string, i);
        }
        new_c = E_Hairpin(u, type, S1[j+1], S1[i-1],  loopseq, P)+my_c[ij];
      }
      if (new_c<FcH) {
        FcH = new_c; Hi=i; Hj=j;
      }

      for (p = j+1; p < length ; p++) {
        int u1, qmin;
        u1 = p-j-1;
        if (u1+i-1>MAXLOOP) break;
        qmin = u1+i-1+length-MAXLOOP;
        if (qmin<p+TURN+1) qmin = p+TURN+1;
        for (q = qmin; q <=length; q++) {
          int u2, type_2/*, si1, sq1*/;
          type_2 = rtype[(unsigned char)ptype[indx[q]+p]];
          if (type_2==0) continue;
          u2 = i-1 + length-q;
          if (u1+u2>MAXLOOP) continue;
          energy = E_IntLoop(u1, u2, type, type_2, S1[j+1], S1[i-1], S1[p-1], S1[q+1], P);
          new_c = my_c[ij] + my_c[indx[q]+p] + energy;
          if (new_c<FcI) {
            FcI = new_c; Ii=i; Ij=j; Ip=p; Iq=q;
          }
        }
      }
    }
  Fc = MIN2(FcI, FcH);

  /* compute the fM2 array (multi loops with exactly 2 helices) */
  /* to get a unique ML decomposition, just use fM1 instead of fML
     below. However, that will not work with dangle_model==1  */
  for (i=1; i<length-TURN; i++) {
    int u;
    fM2[i] = INF;
    for (u=i+TURN; u<length-TURN; u++)
      fM2[i] = MIN2(fM2[i], my_fML[indx[u]+i] + my_fML[indx[length]+u+1]);
  }

  for (i=TURN+1; i<length-2*TURN; i++) {
    int fm;
    fm = my_fML[indx[i]+1]+fM2[i+1]+P->MLclosing;
    if (fm<FcM) {
      FcM=fm; Mi=i;
    }
  }
  Fc = MIN2(Fc, FcM);

  if (dangle_model==1) {
    int u;
    fM_d3 =  (int *) space(sizeof(int)*(length+2));
    fM_d5 =  (int *) space(sizeof(int)*(length+2));
    for (i=TURN+1; i<length-TURN; i++) {
      fM_d3[i] = INF;
      for (u=2+TURN; u<i-TURN; u++)
        fM_d3[i] = MIN2(fM_d3[i], my_fML[indx[u]+2] + my_fML[indx[i]+u+1]);
    }
    for (i=2*TURN+1; i<length-TURN; i++) {
      int fm, type;
      type = ptype[indx[length]+i+1];
      if (type==0) continue;
      fm = fM_d3[i]+my_c[indx[length]+i+1]+E_MLstem(type, -1, S1[1], P) + P->MLclosing;
      if (fm<FcMd3) {
        FcMd3=fm; Md3i=i;
      }
      fm = fM_d3[i-1]+my_c[indx[length]+i+1]+E_MLstem(type, S1[i], S1[1], P) + P->MLclosing;
      if (fm<FcMd3) {
        FcMd3=fm; Md3i=-i;
      }
    }

    for (i=TURN+1; i<length-TURN; i++) {
      fM_d5[i] = INF;
      for (u=i+TURN; u<length-TURN; u++)
        fM_d5[i] = MIN2(fM_d5[i], my_fML[indx[u]+i] + my_fML[indx[length-1]+u+1]);
    }
    for (i=TURN+1; i<length-2*TURN; i++) {
      int fm, type;
      type = ptype[indx[i]+1];
      if (type==0) continue;
      fm = E_MLstem(type, S1[length], -1, P) + my_c[indx[i]+1] + fM_d5[i+1] + P->MLclosing;
      if (fm<FcMd5) {
        FcMd5=fm; Md5i=i;
      }
      fm = E_MLstem(type, S1[length], S1[i+1], P) + my_c[indx[i]+1] + fM_d5[i+2] + P->MLclosing;
      if (fm<FcMd5) {
        FcMd5=fm; Md5i=-i;
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
  else if(Fc < INF){
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
      int fm, u;
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
  }
  vc->matrices->FcH = FcH;
  vc->matrices->FcI = FcI;
  vc->matrices->FcM = FcM;
  vc->matrices->Fc  = Fc;
}

