/* -*-C-*- */
/* this file contains code for folding circular RNAs */
/* it's #include'd into fold.c */

float circfold(const char *string, char *structure) {
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
  int i,j, p,q,length, energy, bonus=0, bonus_cnt=0;
  int s=0;
  /* if (!uniq_ML) {uniq_ML=1; init_length=-1;} */
  circular = 1;

  length = (int) strlen(string);

#ifdef USE_OPENMP
  /* always init everything since all global static variables are uninitialized when entering a thread */
  init_fold(length);
#else
  if (length>init_length) init_fold(length);
#endif

  if (fabs(P->temperature - temperature)>1e-6) update_fold_params();

  S   = encode_sequence(string, 0);
  S1  = encode_sequence(string, 1);

  BP = (int *)space(sizeof(int)*(length+2));
  make_ptypes(S, structure);

  /* compute all arrays used by fold() for linear RNA */
  energy = fill_arrays(string);

  FcH = FcI= FcM = FcMd3= FcMd5= Fc = INF;
  for (i=1; i<length; i++)
    for (j=i+TURN+1; j <= length; j++) {
      int ij, bonus=0, type, u, new_c, no_close;
      u = length-j + i-1;
      if (u<TURN) continue;

      ij = indx[j]+i;
      type = ptype[ij];

      /* enforcing structure constraints */
      if ((BP[i]==j)||(BP[i]==-1)||(BP[i]==-2)) bonus -= BONUS;
      if ((BP[j]==-1)||(BP[j]==-3)) bonus -= BONUS;
      if ((BP[i]==-4)||(BP[j]==-4)) type=0;

      no_close = (((type==3)||(type==4))&&no_closingGU&&(bonus==0));

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
        new_c = E_Hairpin(u, type, S1[j+1], S1[i-1],  loopseq, P)+bonus+c[ij];
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
          type_2 = rtype[ptype[indx[q]+p]];
          if (type_2==0) continue;
          u2 = i-1 + length-q;
          if (u1+u2>MAXLOOP) continue;
          energy = E_IntLoop(u1, u2, type, type_2, S1[j+1], S1[i-1], S1[p-1], S1[q+1], P);
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
    fm = fML[indx[i]+1]+fM2[i+1]+P->MLclosing;
    if (fm<FcM) {
      FcM=fm; Mi=i;
    }
  }
  Fc = MIN2(Fc, FcM);

  if (dangles==1) {
    int u;
    fM_d3 =  (int *) space(sizeof(int)*(length+2));
    fM_d5 =  (int *) space(sizeof(int)*(length+2));
    for (i=TURN+1; i<length-TURN; i++) {
      fM_d3[i] = INF;
      for (u=2+TURN; u<i-TURN; u++)
        fM_d3[i] = MIN2(fM_d3[i], fML[indx[u]+2] + fML[indx[i]+u+1]);
    }
    for (i=2*TURN+1; i<length-TURN; i++) {
      int fm, type;
      type = ptype[indx[length]+i+1];
      if (type==0) continue;
      fm = fM_d3[i]+c[indx[length]+i+1]+P->MLclosing+P->MLintern[type]+P->dangle3[type][S1[1]];
      if (fm<FcMd3) {
        FcMd3=fm; Md3i=i;
      }
      fm = fM_d3[i-1]+c[indx[length]+i+1]+P->MLclosing+P->MLintern[type]+P->dangle3[type][S1[1]]+P->dangle5[type][S1[i]];
      if (fm<FcMd3) {
        FcMd3=fm; Md3i=-i;
      }
    }

    for (i=TURN+1; i<length-TURN; i++) {
      fM_d5[i] = INF;
      for (u=i+TURN; u<length-TURN; u++)
        fM_d5[i] = MIN2(fM_d5[i], fML[indx[u]+i] + fML[indx[length-1]+u+1]);
    }
    for (i=TURN+1; i<length-2*TURN; i++) {
      int fm, type;
      type = ptype[indx[i]+1];
      if (type==0) continue;
      fm = P->dangle5[type][S1[length]]+c[indx[i]+1]+fM_d5[i+1]+P->MLclosing+P->MLintern[type];
      if (fm<FcMd5) {
        FcMd5=fm; Md5i=i;
      }
      fm = P->dangle5[type][S1[length]]+c[indx[i]+1]+P->dangle3[type][S1[i+1]]+
        fM_d5[i+2]+P->MLclosing+P->MLintern[type];
      if (fm<FcMd5) {
        FcMd5=fm; Md5i=-i;
      }
    }
    if (FcMd5<MIN2(Fc,FcMd3)) {
      /* looks like we have to do this ... */
      sector[++s].i = 1;
      sector[s].j = (Md5i>0)?Md5i:-Md5i;
      sector[s].ml = 2;
      i = (Md5i>0)?Md5i+1 : -Md5i+2; /* let's backtrack fm_d5[Md5i+1] */
      for (u=i+TURN; u<length-TURN; u++)
        if (fM_d5[i] == fML[indx[u]+i] + fML[indx[length-1]+u+1]) {
          sector[++s].i = i;
          sector[s].j = u;
          sector[s].ml = 1;
          sector[++s].i =u+1;
          sector[s].j = length-1;
          sector[s].ml = 1;
          break;
        }
      Fc = FcMd5;
    } else if (FcMd3<Fc) {
      /* here we go again... */
      sector[++s].i = (Md3i>0)?Md3i+1:-Md3i+1;
      sector[s].j = length;
      sector[s].ml = 2;
      i = (Md3i>0)? Md3i : -Md3i-1; /* let's backtrack fm_d3[Md3i] */
      for (u=2+TURN; u<i-TURN; u++)
        if (fM_d3[i] == fML[indx[u]+2] + fML[indx[i]+u+1]) {
          sector[++s].i = 2;
          sector[s].j = u;
          sector[s].ml = 1;
          sector[++s].i =u+1;
          sector[s].j = i;
          sector[s].ml = 1;
          break;
        }
      Fc = FcMd3;
    }
    free(fM_d3);
    free(fM_d5);
  }
  if(Fc < INF){
    if (FcH==Fc) {
      sector[++s].i = Hi;
      sector[s].j = Hj;
      sector[s].ml = 2;
    }
    else if (FcI==Fc) {
      sector[++s].i = Ii;
      sector[s].j = Ij;
      sector[s].ml = 2;
      sector[++s].i = Ip;
      sector[s].j = Iq;
      sector[s].ml = 2;
    }
    else if (FcM==Fc) { /* grumpf we found a Multiloop */
      int fm, u;
      /* backtrack in fM2 */
      fm = fM2[Mi+1];
      for (u=Mi+TURN+1; u<length-TURN; u++)
        if (fm == fML[indx[u]+Mi+1] + fML[indx[length]+u+1]) {
                sector[++s].i=Mi+1;
                sector[s].j=u;
                sector[s].ml = 1;
                sector[++s].i=u+1;
                sector[s].j=length;
                sector[s].ml = 1;
                break;
        }
      sector[++s].i = 1;
      sector[s].j = Mi;
      sector[s].ml = 1;
    }
  }
  backtrack(string, s);

#ifdef PAREN
  parenthesis_structure(structure, base_pair2, length);
#else
  letter_structure(structure, base_pair2, length);
#endif

  /*
  *  Backward compatibility:
  *  This block may be removed if deprecated functions
  *  relying on the global variable "base_pair" vanishs from within the package!
  */
  {
    if(base_pair) free(base_pair);
    base_pair = (bondT *)space(sizeof(bondT) * (1+length/2));
    memcpy(base_pair, base_pair2, sizeof(bondT) * (1+length/2));
  }

  /* check constraints */
  for(i=1;i<=length;i++) {
    if((BP[i]<0)&&(BP[i]>-4)) {
      bonus_cnt++;
      if((BP[i]==-3)&&(structure[i-1]==')')) bonus++;
      if((BP[i]==-2)&&(structure[i-1]=='(')) bonus++;
      if((BP[i]==-1)&&(structure[i-1]!='.')) bonus++;
    }

    if(BP[i]>i) {
      int l;
      bonus_cnt++;
      for(l=1; l<=base_pair2[0].i; l++)
        if((i==base_pair2[l].i)&&(BP[i]==base_pair2[l].j)) bonus++;
    }
  }

  if (bonus_cnt>bonus) fprintf(stderr,"\ncould not enforce all constraints\n");
  bonus*=BONUS;

  free(S); free(S1); free(BP);

  energy = Fc + bonus;      /*remove bonus energies from result */
  if(energy == INF)
    return (float) 0.0;
  else
    return (float) energy/100.;
}
