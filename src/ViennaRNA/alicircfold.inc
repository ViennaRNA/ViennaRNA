/* -*-C-*- */
/* this file contains code for alifolding circular RNAs */
/* it's #include'd into alifold.c */

PRIVATE void fill_arrays_circ(vrna_fold_compound *vc, sect bt_stack[], int *bt){
  /* variant of alifold() for circular RNAs */
  /* auxiliarry arrays:
     fM2 = multiloop region with exactly two stems, extending to 3' end
     for stupid dangles=1 case we also need:
     fM_d3 = multiloop region with >= 2 stems, starting at pos 2
             (a pair (k,n) will form 3' dangle with pos 1)
     fM_d5 = multiloop region with >= 2 stems, extending to pos n-1
             (a pair (1,k) will form a 5' dangle with pos n)
  */
  int Fc, FcH, FcI, FcM, Hi, Hj, Ii, Ij, Ip, Iq, Mi;
  int FcMd3, FcMd5;
  int i,j, p,q,length, energy;
  int n_seq, psc;
  int *type, s;

  short           **S     = vc->S;                                                                   
  short           **S5    = vc->S5;     /*S5[s][i] holds next base 5' of i in sequence s*/            
  short           **S3    = vc->S3;     /*Sl[s][i] holds next base 3' of i in sequence s*/            
  char            **Ss    = vc->Ss;                                                                   
  unsigned short  **a2s   = vc->a2s;                                                                   
  paramT          *P      = vc->params;                                                                   
  model_detailsT  *md     = &(P->model_details);
  int             *indx   = vc->jindx;     /* index for moving in the triangle matrices c[] and fMl[]*/  
  int             *c      = vc->matrices->c;     /* energy array, given that i-j pair */                       
  int             *fML    = vc->matrices->fML;     /* multi-loop auxiliary energy array */                       
  int             *fM2    = vc->matrices->fM2;
  int             *pscore = vc->pscore;     /* precomputed array of pair types */                         

  hard_constraintT *hc    = vc->hc;
  soft_constraintT **sc   = vc->scs;

  char **strings;
  strings = vc->sequences;
  n_seq   = vc->n_seq;
  length  = vc->length;

  int   dangle_model  = P->model_details.dangles;
  int   noGUclosure   = P->model_details.noGUclosure;
  int   *rtype        = &(P->model_details.rtype[0]);
  type                = (int *) space(n_seq*sizeof(int));

  oldAliEn=1; /* TODO allow new gap treatment in alicircfold() */


  /* extra arrays for circfold() */
  fM2 =  (int *) space(sizeof(int)*(length+2));
  FcH = FcI= FcM = FcMd3= FcMd5= Fc = INF;
  for (i=1; i<length; i++)
    for (j=i+TURN+1; j <= length; j++) {
      int ij, u, new_c;
      u = length-j + i-1;
      if (u<TURN) continue;

      ij = indx[j]+i;

      psc = pscore[indx[j]+i];

      if (psc<MINPSCORE) continue;

      new_c = 0;
      for (s=0; s<n_seq; s++) {
        char loopseq[10];
        int rt, si1, sj1;

        type[s] = md->pair[S[s][i]][S[s][j]];
        if (type[s]==0) type[s]=7;

        rt=rtype[type[s]];
        if (u<7) {
          strcpy(loopseq , strings[s]+j-1);
          strncat(loopseq, strings[s], i);
        }
        si1 = (i==1)?S[s][length] : S[s][i-1];
        sj1 = (j==length)?S[s][1] : S[s][j+1];
        new_c += E_Hairpin(u, rt, sj1, si1,  loopseq, P);
      }
      new_c += +c[ij];
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
          int u2, si1, sq1, type_2;

          if (pscore[indx[q]+p]<MINPSCORE) continue;

          u2 = i-1 + length-q;
          if (u1+u2>MAXLOOP) continue;
          for (energy = s=0; s<n_seq; s++) {
            int rt;
            rt=rtype[type[s]];
            type_2 = md->pair[S[s][q]][S[s][p]]; /* q,p not p,q! */
            if (type_2 ==0) type_2=7;
            si1 = (i==1)? S[s][length] : S[s][i-1];
            sq1 = (q==length)? S[s][1] : S[s][q+1];
            energy += E_IntLoop(u1, u2, rt, type_2,
                                 S[s][j+1], si1, S[s][p-1], sq1, P);
          }
          new_c = c[ij] + c[indx[q]+p] + energy;
          if (new_c<FcI) {
            FcI = new_c; Ii=i; Ij=j; Ip=p; Iq=q;
          }
        }
      }
    }
  Fc = MIN2(FcI, FcH);

  /* compute the fM2 array (multi loops with exactly 2 helices) */
  /* to get a unique ML decomposition, just use fM1 instead of fML
     below. However, that will not work with dangles==1  */
  for (i=1; i<length-TURN; i++) {
    int u;
    fM2[i] = INF;
    for (u=i+TURN; u<length-TURN; u++)
      fM2[i] = MIN2(fM2[i], fML[indx[u]+i] + fML[indx[length]+u+1]);
  }

  for (i=TURN+1; i<length-2*TURN; i++) {
    int fm;
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
    int fm, u;
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
  }
  vc->matrices->FcH = FcH;
  vc->matrices->FcI = FcI;
  vc->matrices->FcM = FcM;
  vc->matrices->Fc  = Fc;
}
