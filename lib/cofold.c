/* Last changed Time-stamp: <2008-12-03 17:44:38 ivo> */
/*
                  minimum free energy
                  RNA secondary structure prediction

                  c Ivo Hofacker, Chrisoph Flamm
                  original implementation by
                  Walter Fontana

                  Vienna RNA package
*/

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "utils.h"
#include "energy_par.h"
#include "fold_vars.h"
#include "pair_mat.h"
#include "params.h"
#include "subopt.h"
#include "fold.h"
#include "loop_energies.h"
#include "cofold.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/*@unused@*/
static char rcsid[] UNUSED = "$Id: cofold.c,v 1.12 2008/12/03 16:55:50 ivo Exp $";

#define PAREN

#define STACK_BULGE1      1       /* stacking energies for bulges of size 1 */
#define NEW_NINIO         1       /* new asymetry penalty */
#define MAXSECTORS        500     /* dimension for a backtrack array */
#define LOCALITY          0.      /* locality parameter for base-pairs */
#undef TURN
#define TURN              0       /* reset minimal base pair span for intermolecular pairings */
#define TURN2             3       /* used by zukersubopt */
#define SAME_STRAND(I,J)  (((I)>=cut_point)||((J)<cut_point))

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

PRIVATE float   mfe1, mfe2;       /* minimum free energies of the monomers */
PRIVATE int     *indx   = NULL;   /* index for moving in the triangle matrices c[] and fMl[]*/
PRIVATE int     *c      = NULL;   /* energy array, given that i-j pair */
PRIVATE int     *cc     = NULL;   /* linear array for calculating canonical structures */
PRIVATE int     *cc1    = NULL;   /*   "     "        */
PRIVATE int     *f5     = NULL;   /* energy of 5' end */
PRIVATE int     *fc     = NULL;   /* energy from i to cutpoint (and vice versa if i>cut) */
PRIVATE int     *fML    = NULL;   /* multi-loop auxiliary energy array */
PRIVATE int     *fM1    = NULL;   /* second ML array, only for subopt */
PRIVATE int     *Fmi    = NULL;   /* holds row i of fML (avoids jumps in memory) */
PRIVATE int     *DMLi   = NULL;   /* DMLi[j] holds MIN(fML[i,k]+fML[k+1,j])  */
PRIVATE int     *DMLi1  = NULL;   /*             MIN(fML[i+1,k]+fML[k+1,j])  */
PRIVATE int     *DMLi2  = NULL;   /*             MIN(fML[i+2,k]+fML[k+1,j])  */
PRIVATE char    *ptype  = NULL;   /* precomputed array of pair types */
PRIVATE short   *S = NULL, *S1 = NULL;
PRIVATE paramT  *P          = NULL;
PRIVATE int     init_length = -1;
PRIVATE int     zuker       = 0; /* Do Zuker style suboptimals? */
PRIVATE char    alpha[]     = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
PRIVATE sect    sector[MAXSECTORS];   /* stack for backtracking */
PRIVATE int     length;
PRIVATE bondT   *base_pair2 = NULL;
PRIVATE int     *BP; /* contains the structure constrainsts: BP[i]
                        -1: | = base must be paired
                        -2: < = base must be paired with j<i
                        -3: > = base must be paired with j>i
                        -4: x = base must not pair
                        positive int: base is paired with int      */
PRIVATE int     struct_constrained = 0;

#ifdef _OPENMP

/* NOTE: all variables are assumed to be uninitialized if they are declared as threadprivate
         thus we have to initialize them before usage by a seperate function!
         OR: use copyin in the PARALLEL directive!
         e.g.:
         #pragma omp parallel for copyin(P, init_length, ...)
*/
#pragma omp threadprivate(mfe1, mfe2, indx, c, cc, cc1, f5, fc, fML, fM1, Fmi, DMLi, DMLi1, DMLi2,\
                          ptype, S, S1, P, zuker, sector, length, base_pair2, BP, struct_constrained)

#endif

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/

PRIVATE void  init_cofold(int length, paramT *parameters);
PRIVATE void  get_arrays(unsigned int size);
/* PRIVATE void  scale_parameters(void); */
PRIVATE void  make_ptypes(const short *S, const char *structure);
PRIVATE void  backtrack(const char *sequence);
PRIVATE int   fill_arrays(const char *sequence);
PRIVATE void  free_end(int *array, int i, int start);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

/*--------------------------------------------------------------------------*/
PRIVATE void init_cofold(int length, paramT *parameters){

#ifdef _OPENMP
/* Explicitly turn off dynamic threads */
  omp_set_dynamic(0);
#endif

  if (length<1) nrerror("init_cofold: argument must be greater 0");
  free_co_arrays();
  get_arrays((unsigned) length);
  init_length=length;

  indx = get_indx((unsigned) length);

  update_cofold_params_par(parameters);
}

/*--------------------------------------------------------------------------*/

PRIVATE void get_arrays(unsigned int size){
  if(size >= (unsigned int)sqrt((double)INT_MAX))
    nrerror("get_arrays@cofold.c: sequence length exceeds addressable range");

  c     = (int *) space(sizeof(int)*((size*(size+1))/2+2));
  fML   = (int *) space(sizeof(int)*((size*(size+1))/2+2));
  if (uniq_ML)
    fM1    = (int *) space(sizeof(int)*((size*(size+1))/2+2));

  ptype = (char *) space(sizeof(char)*((size*(size+1))/2+2));
  f5    = (int *) space(sizeof(int)*(size+2));
  fc    = (int *) space(sizeof(int)*(size+2));
  cc    = (int *) space(sizeof(int)*(size+2));
  cc1   = (int *) space(sizeof(int)*(size+2));
  Fmi   = (int *) space(sizeof(int)*(size+1));
  DMLi  = (int *) space(sizeof(int)*(size+1));
  DMLi1 = (int *) space(sizeof(int)*(size+1));
  DMLi2 = (int *) space(sizeof(int)*(size+1));
  if (base_pair2) free(base_pair2);
  base_pair2 = (bondT *) space(sizeof(bondT)*(1+size/2));
}

/*--------------------------------------------------------------------------*/

PUBLIC void free_co_arrays(void){
  if(indx)        free(indx);
  if(c)           free(c);
  if(fML)         free(fML);
  if(f5)          free(f5);
  if(cc)          free(cc);
  if(cc1)         free(cc1);
  if(fc)          free(fc);
  if(ptype)       free(ptype);
  if(fM1)         free(fM1);
  if(base_pair2)  free(base_pair2);
  if(Fmi)         free(Fmi);
  if(DMLi)        free(DMLi);
  if(DMLi1)       free(DMLi1);
  if(DMLi2)       free(DMLi2);
  if(P)           free(P);
  indx = c = fML = f5 = cc = cc1 = fc = fM1 = Fmi = DMLi = DMLi1 = DMLi2 = NULL;
  ptype       = NULL;
  base_pair2  = NULL;
  P           = NULL;
  init_length = 0;
}


/*--------------------------------------------------------------------------*/

PUBLIC void export_cofold_arrays( int **f5_p,
                                  int **c_p,
                                  int **fML_p,
                                  int **fM1_p,
                                  int **fc_p,
                                  int **indx_p,
                                  char **ptype_p){

  /* make the DP arrays available to routines such as subopt() */
  *f5_p = f5; *c_p = c;
  *fML_p = fML; *fM1_p = fM1;
  *indx_p = indx; *ptype_p = ptype;
  *fc_p =fc;
}

/*--------------------------------------------------------------------------*/

PUBLIC float cofold(const char *string, char *structure) {
  return cofold_par(string, structure, NULL, fold_constrained);
}

PUBLIC float cofold_par(const char *string,
                        char *structure,
                        paramT *parameters,
                        int is_constrained){

  int i, length, energy, bonus=0, bonus_cnt=0;

  zuker = 0;

  struct_constrained  = is_constrained;
  length              = (int) strlen(string);

#ifdef _OPENMP
  /* always init everything since all global static variables are uninitialized when entering a thread */
  init_cofold(length, parameters);
#else
  if(parameters) init_cofold(length, parameters);
  else if (length>init_length) init_cofold(length, parameters);
  else if (fabs(P->temperature - temperature)>1e-6) update_cofold_params_par(parameters);
#endif

  S   = encode_sequence(string, 0);
  S1  = encode_sequence(string, 1);
  S1[0] = S[0]; /* store length at pos. 0 */

  BP = (int *)space(sizeof(int)*(length+2));
  make_ptypes(S, structure);

  energy = fill_arrays(string);

  backtrack(string);

#ifdef PAREN
  parenthesis_structure(structure, base_pair2, length);
#else
  letter_structure(structure, base_pair2, length);
#endif

  /*
  *  Backward compatibility:
  *  This block may be removed if deprecated functions
  *  relying on the global variable "base_pair" vanish from within the package!
  */
  base_pair = base_pair2;
  /*
  {
    if(base_pair) free(base_pair);
    base_pair = (bondT *)space(sizeof(bondT) * (1+length/2));
    memcpy(base_pair, base_pair2, sizeof(bondT) * (1+length/2));
  }
  */

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

  energy += bonus;      /*remove bonus energies from result */

  if (backtrack_type=='C')
    return (float) c[indx[length]+1]/100.;
  else if (backtrack_type=='M')
    return (float) fML[indx[length]+1]/100.;
  else
    return (float) energy/100.;
}

PRIVATE int fill_arrays(const char *string) {
  /* fill "c", "fML" and "f5" arrays and return  optimal energy */

  int   i, j, k, length, energy;
  int   decomp, new_fML, max_separation;
  int   no_close, type, type_2, tt, maxj;
  int   bonus=0;
  int   dangle_model  = P->model_details.dangles;
  int   noGUclosure   = P->model_details.noGUclosure;
  int   noLP          = P->model_details.noLP;

  length = (int) strlen(string);

  max_separation = (int) ((1.-LOCALITY)*(double)(length-2)); /* not in use */

  for (j=1; j<=length; j++) {
    Fmi[j]=DMLi[j]=DMLi1[j]=DMLi2[j]=INF;
    fc[j]=0;
  }

  for (j = 1; j<=length; j++)
    for (i=1; i<=j; i++) {
      c[indx[j]+i] = fML[indx[j]+i] = INF;
      if (uniq_ML) fM1[indx[j]+i] = INF;
    }

  for (i = length-TURN-1; i >= 1; i--) { /* i,j in [1..length] */

    maxj=(zuker)? (MIN2(i+cut_point-1,length)):length;
    for (j = i+TURN+1; j <= maxj; j++) {
      int p, q, ij;
      ij = indx[j]+i;
      bonus = 0;
      type = ptype[ij];

      /* enforcing structure constraints */
      if ((BP[i]==j)||(BP[i]==-1)||(BP[i]==-2)) bonus -= BONUS;
      if ((BP[j]==-1)||(BP[j]==-3)) bonus -= BONUS;
      if ((BP[i]==-4)||(BP[j]==-4)) type=0;

      no_close = (((type==3)||(type==4))&&noGUclosure&&(bonus==0));

      if (j-i-1 > max_separation) type = 0;  /* forces locality degree */

      if (type) {   /* we have a pair */
        int new_c=0, stackEnergy=INF;
        short si, sj;
        si  = SAME_STRAND(i, i+1) ? S1[i+1] : -1;
        sj  = SAME_STRAND(j-1, j) ? S1[j-1] : -1;
        /* hairpin ----------------------------------------------*/

        if (SAME_STRAND(i,j)) {
          if (no_close) new_c = FORBIDDEN;
          else
            new_c = E_Hairpin(j-i-1, type, si, sj, string+i-1, P);
        }
        else {
          if (dangle_model)
            new_c += E_ExtLoop(rtype[type], sj, si, P);
          else
            new_c += E_ExtLoop(rtype[type], -1, -1, P);
        }
        /*--------------------------------------------------------
          check for elementary structures involving more than one
          closing pair.
          --------------------------------------------------------*/

        for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1) ; p++) {
          int minq = j-i+p-MAXLOOP-2;
          if (minq<p+1+TURN) minq = p+1+TURN;
          for (q = minq; q < j; q++) {
            type_2 = ptype[indx[q]+p];

            if (type_2==0) continue;
            type_2 = rtype[type_2];

            if (noGUclosure)
              if (no_close||(type_2==3)||(type_2==4))
                if ((p>i+1)||(q<j-1)) continue;  /* continue unless stack */

            if (SAME_STRAND(i,p) && SAME_STRAND(q,j))
              energy = E_IntLoop(p-i-1, j-q-1, type, type_2, si, sj, S1[p-1], S1[q+1], P);
            else
              energy = E_IntLoop_Co(rtype[type], rtype[type_2],
                                    i, j, p, q,
                                    cut_point,
                                    si, sj,
                                    S1[p-1], S1[q+1],
                                    dangle_model,
                                    P);

            new_c = MIN2(energy+c[indx[q]+p], new_c);
            if ((p==i+1)&&(j==q+1)) stackEnergy = energy; /* remember stack energy */

          } /* end q-loop */
        } /* end p-loop */

        /* multi-loop decomposition ------------------------*/


        if (!no_close) {
          int MLenergy;

          if((si >= 0) && (sj >= 0)){
            decomp = DMLi1[j-1];
            tt = rtype[type];
            MLenergy = P->MLclosing;
            switch(dangle_model){
              case 0:   MLenergy += decomp + E_MLstem(tt, -1, -1, P);
                        break;
              case 2:   MLenergy += decomp + E_MLstem(tt, sj, si, P);
                        break;
              default:  decomp += E_MLstem(tt, -1, -1, P);
                        decomp = MIN2(decomp, DMLi1[j-2] + E_MLstem(tt, sj, -1, P) + P->MLbase);
                        decomp = MIN2(decomp, DMLi2[j-1] + E_MLstem(tt, -1, si, P) + P->MLbase);
                        decomp = MIN2(decomp, DMLi2[j-2] + E_MLstem(tt, sj, si, P) + 2*P->MLbase);
                        MLenergy += decomp;
                        break;
            }
            new_c = MIN2(new_c, MLenergy);
          }

          if (!SAME_STRAND(i,j)) { /* cut is somewhere in the multiloop*/
            decomp = fc[i+1] + fc[j-1];
            tt = rtype[type];
            switch(dangle_model){
              case 0:   decomp += E_ExtLoop(tt, -1, -1, P);
                        break;
              case 2:   decomp += E_ExtLoop(tt, sj, si, P);
                        break;
              default:  decomp += E_ExtLoop(tt, -1, -1, P);
                        decomp = MIN2(decomp, fc[i+2] + fc[j-2] + E_ExtLoop(tt, sj, si, P));
                        decomp = MIN2(decomp, fc[i+2] + fc[j-1] + E_ExtLoop(tt, -1, si, P));
                        decomp = MIN2(decomp, fc[i+1] + fc[j-2] + E_ExtLoop(tt, sj, -1, P));
                        break;
            }
            new_c = MIN2(new_c, decomp);
          }
        } /* end >> if (!no_close) << */

        /* coaxial stacking of (i.j) with (i+1.k) or (k+1.j-1) */

        if (dangle_model==3) {
          decomp = INF;
          for (k = i+2+TURN; k < j-2-TURN; k++) {
            type_2 = ptype[indx[k]+i+1]; type_2 = rtype[type_2];
            if (type_2)
              decomp = MIN2(decomp, c[indx[k]+i+1]+P->stack[type][type_2]+
                            fML[indx[j-1]+k+1]);
            type_2 = ptype[indx[j-1]+k+1]; type_2 = rtype[type_2];
            if (type_2)
              decomp = MIN2(decomp, c[indx[j-1]+k+1]+P->stack[type][type_2]+
                            fML[indx[k]+i+1]);
          }
          /* no TermAU penalty if coax stack */
          decomp += 2*P->MLintern[1] + P->MLclosing;
          new_c = MIN2(new_c, decomp);
        }

        new_c = MIN2(new_c, cc1[j-1]+stackEnergy);
        cc[j] = new_c + bonus;
        if (noLP){
          if (SAME_STRAND(i,i+1) && SAME_STRAND(j-1,j))
            c[ij] = cc1[j-1]+stackEnergy+bonus;
          else /* currently we don't allow stacking over the cut point */
            c[ij] = FORBIDDEN;
        }
        else
          c[ij] = cc[j];

      } /* end >> if (pair) << */

      else c[ij] = INF;


      /* done with c[i,j], now compute fML[i,j] */
      /* free ends ? -----------------------------------------*/
      new_fML=INF;
      if (SAME_STRAND(i-1,i)) {
        if (SAME_STRAND(i,i+1)) new_fML = fML[ij+1]+P->MLbase;
        if (SAME_STRAND(j-1,j)) new_fML = MIN2(fML[indx[j-1]+i]+P->MLbase, new_fML);
        if (SAME_STRAND(j,j+1)) {
          energy = c[ij];
          if(dangle_model == 2) energy += E_MLstem(type,(i>1) ? S1[i-1] : -1, (j<length) ? S1[j+1] : -1, P);
          else energy += E_MLstem(type, -1, -1, P);
          new_fML = MIN2(new_fML, energy);
          if(uniq_ML){
            fM1[ij] = energy;
            if(SAME_STRAND(j-1,j)) fM1[ij] = MIN2(energy, fM1[indx[j-1]+i] + P->MLbase);
          }
        }
        if (dangle_model%2==1) {  /* normal dangles */
          if (SAME_STRAND(i,i+1)) {
            tt = ptype[ij+1]; /* i+1,j */
            new_fML = MIN2(new_fML, c[ij+1] + P->MLbase + E_MLstem(tt, S1[i], -1, P));
          }
          if (SAME_STRAND(j-1,j)) {
            tt = ptype[indx[j-1]+i]; /* i,j-1 */
            new_fML = MIN2(new_fML, c[indx[j-1]+i] + P->MLbase + E_MLstem(tt, -1, S1[j], P));
          }
          if ((SAME_STRAND(j-1,j))&&(SAME_STRAND(i,i+1))) {
            tt = ptype[indx[j-1]+i+1]; /* i+1,j-1 */
            new_fML = MIN2(new_fML, c[indx[j-1]+i+1] + 2*P->MLbase + E_MLstem(tt, S1[i], S1[j], P));
          }
        }
      }

      /* modular decomposition -------------------------------*/

      {
        int stopp;     /*loop 1 up to cut, then loop 2*/
        stopp=(cut_point>0)? (cut_point):(j-2-TURN);
        for (decomp=INF, k = i+1+TURN; k<stopp; k++)
          decomp = MIN2(decomp, Fmi[k]+fML[indx[j]+k+1]);
        k++;
        for (;k <= j-2-TURN;k++)
          decomp = MIN2(decomp, Fmi[k]+fML[indx[j]+k+1]);
      }
      DMLi[j] = decomp;               /* store for use in ML decompositon */
      new_fML = MIN2(new_fML,decomp);

      /* coaxial stacking */
      if (dangle_model==3) {
        int stopp;
        stopp=(cut_point>0)? (cut_point):(j-2-TURN);
        /* additional ML decomposition as two coaxially stacked helices */
        for (decomp = INF, k = i+1+TURN; k<stopp; k++) {
          type = ptype[indx[k]+i]; type = rtype[type];
          type_2 = ptype[indx[j]+k+1]; type_2 = rtype[type_2];
          if (type && type_2)
            decomp = MIN2(decomp,
                          c[indx[k]+i]+c[indx[j]+k+1]+P->stack[type][type_2]);
        }
        k++;
        for (;k <= j-2-TURN; k++) {
          type = ptype[indx[k]+i]; type = rtype[type];
          type_2 = ptype[indx[j]+k+1]; type_2 = rtype[type_2];
          if (type && type_2)
            decomp = MIN2(decomp,
                          c[indx[k]+i]+c[indx[j]+k+1]+P->stack[type][type_2]);
        }

        decomp += 2*P->MLintern[1];

#if 0
        /* This is needed for Y shaped ML loops with coax stacking of
           interior pairs, but backtracking will fail if activated */
        DMLi[j] = MIN2(DMLi[j], decomp);
        if (SAME_STRAND(j-1,j)) DMLi[j] = MIN2(DMLi[j], DMLi[j-1]+P->MLbase);
        if (SAME_STRAND(i,i+1)) DMLi[j] = MIN2(DMLi[j], DMLi1[j]+P->MLbase);
        new_fML = MIN2(new_fML, DMLi[j]);
#endif
        new_fML = MIN2(new_fML, decomp);
      }

      fML[ij] = Fmi[j] = new_fML;     /* substring energy */

    }

    if (i==cut_point)
      for (j=i; j<=maxj; j++)
        free_end(fc, j, cut_point);
    if (i<cut_point)
      free_end(fc,i,cut_point-1);


    {
      int *FF; /* rotate the auxilliary arrays */
      FF = DMLi2; DMLi2 = DMLi1; DMLi1 = DMLi; DMLi = FF;
      FF = cc1; cc1=cc; cc=FF;
      for (j=1; j<=maxj; j++) {cc[j]=Fmi[j]=DMLi[j]=INF; }
    }
  }

  /* calculate energies of 5' and 3' fragments */

  for (i=1; i<=length; i++)
    free_end(f5, i, 1);

  if (cut_point>0) {
    mfe1=f5[cut_point-1];
    mfe2=fc[length];
    /* add DuplexInit, check whether duplex*/
    for (i=cut_point; i<=length; i++) {
      f5[i]=MIN2(f5[i]+P->DuplexInit, fc[i]+fc[1]);
    }
  }

  energy = f5[length];
  if (cut_point<1) mfe1=mfe2=energy;
  return energy;
}

PRIVATE void backtrack_co(const char *string, int s, int b /* b=0: start new structure, b \ne 0: add to existing structure */) {

  /*------------------------------------------------------------------
    trace back through the "c", "fc", "f5" and "fML" arrays to get the
    base pairing list. No search for equivalent structures is done.
    This is fast, since only few structure elements are recalculated.
    ------------------------------------------------------------------*/

  int   i, j, k, length, energy, new;
  int   no_close, type, type_2, tt;
  int   bonus;
  int   dangle_model  = P->model_details.dangles;
  int   noGUclosure   = P->model_details.noGUclosure;
  int   noLP          = P->model_details.noLP;

  /* int   b=0;*/

  length = strlen(string);
  if (s==0) {
    sector[++s].i = 1;
    sector[s].j   = length;
    sector[s].ml  = (backtrack_type=='M') ? 1 : ((backtrack_type=='C')?2:0);
  }
  while (s>0) {
    int ml, fij, fi, cij, traced, i1, j1, d3, d5, mm, p, q, jj=0;
    int canonical = 1;     /* (i,j) closes a canonical structure */
    i  = sector[s].i;
    j  = sector[s].j;
    ml = sector[s--].ml;   /* ml is a flag indicating if backtracking is to
                              occur in the fML- (1) or in the f-array (0) */
    if (ml==2) {
      base_pair2[++b].i = i;
      base_pair2[b].j   = j;
      goto repeat1;
    }

    if (j < i+TURN+1) continue; /* no more pairs in this interval */


    if (ml==0) {fij = f5[j]; fi = f5[j-1];}
    else if (ml==1) {fij = fML[indx[j]+i]; fi = fML[indx[j-1]+i]+P->MLbase;}
    else /* 3 or 4 */ {
      fij = fc[j];
      fi = (ml==3) ? INF : fc[j-1];
    }
    if (fij == fi) {  /* 3' end is unpaired */
      sector[++s].i = i;
      sector[s].j   = j-1;
      sector[s].ml  = ml;
      continue;
    }

    if (ml==0 || ml==4) { /* backtrack in f5 or fc[i=cut,j>cut] */
      int *ff;
      ff = (ml==4) ? fc : f5;
      switch(dangle_model){
        case 0:   /* j or j-1 is paired. Find pairing partner */
                  for (k=j-TURN-1,traced=0; k>=i; k--) {
                    int cc;
                    type = ptype[indx[j]+k];
                    if(type){
                      cc = c[indx[j]+k];
                      if(!SAME_STRAND(k,j)) cc += P->DuplexInit;
                      if(fij == ff[k-1] + cc + E_ExtLoop(type, -1, -1, P)){
                        traced = j; jj = k-1;
                      }
                    }
                    if(traced) break;
                  }

                  break;

        case 2:   /* j or j-1 is paired. Find pairing partner */
                  for (k=j-TURN-1,traced=0; k>=i; k--) {
                    int cc;
                    type = ptype[indx[j]+k];
                    if(type){
                      cc = c[indx[j]+k];
                      if(!SAME_STRAND(k,j)) cc += P->DuplexInit;
                      if(fij == ff[k-1] + cc + E_ExtLoop(type, (k>1) && SAME_STRAND(k-1,k) ? S1[k-1] : -1, (j<length) && SAME_STRAND(j,j+1) ? S1[j+1] : -1, P)){
                        traced = j; jj = k-1;
                      }
                    }
                    if(traced) break;
                  }
                  break;

        default:  for(k=j-TURN-1,traced=0; k>=i; k--){
                    int cc, en;
                    type = ptype[indx[j]+k];
                    if(type){
                      cc = c[indx[j]+k];
                      if(!SAME_STRAND(k,j)) cc += P->DuplexInit;
                      if(fij == ff[k-1] + cc + E_ExtLoop(type, -1, -1, P)){
                        traced = j; jj = k-1; break;
                      }
                      if((k>1) && SAME_STRAND(k-1,k))
                        if(fij == ff[k-2] + cc + E_ExtLoop(type, S1[k-1], -1, P)){
                              traced=j; jj=k-2; break;
                        }
                    }

                    type = ptype[indx[j-1]+k];
                    if(type && SAME_STRAND(j-1,j)){
                      cc = c[indx[j-1]+k];
                      if (!SAME_STRAND(k,j-1)) cc += P->DuplexInit; /*???*/
                      if (fij == cc + ff[k-1] + E_ExtLoop(type, -1, S1[j], P)){
                            traced=j-1; jj = k-1; break;
                      }
                      if(k>i){
                        if (fij == ff[k-2] + cc + E_ExtLoop(type, SAME_STRAND(k-1,k) ? S1[k-1] : -1, S1[j], P)){
                          traced=j-1; jj=k-2; break;
                        }
                      }
                    }
                  }

                  break;
      }
      if (!traced) nrerror("backtrack failed in f5 (or fc)");
      sector[++s].i = i;
      sector[s].j   = jj;
      sector[s].ml  = ml;

      i=k; j=traced;
      base_pair2[++b].i = i;
      base_pair2[b].j   = j;
      goto repeat1;
    }
    else if (ml==3) { /* backtrack in fc[i<cut,j=cut-1] */
      if (fc[i] == fc[i+1]) { /* 5' end is unpaired */
        sector[++s].i = i+1;
        sector[s].j   = j;
        sector[s].ml  = ml;
        continue;
      }
      /* i or i+1 is paired. Find pairing partner */
      switch(dangle_model){
        case 0:   for (k=i+TURN+1, traced=0; k<=j; k++){
                    jj=k+1;
                    type = ptype[indx[k]+i];
                    if (type) {
                      if(fc[i] == fc[k+1] + c[indx[k]+i] + E_ExtLoop(type, -1, -1, P)){
                        traced = i;
                      }
                    }
                    if (traced) break;
                  }
                  break;
        case 2:   for (k=i+TURN+1, traced=0; k<=j; k++){
                    jj=k+1;
                    type = ptype[indx[k]+i];
                    if(type){
                      if(fc[i] == fc[k+1] + c[indx[k]+i] + E_ExtLoop(type,(i>1 && SAME_STRAND(i-1,i)) ? S1[i-1] : -1,  SAME_STRAND(k,k+1) ? S1[k+1] : -1, P)){
                        traced = i;
                      }
                    }
                    if (traced) break;
                  }
                  break;
        default:  for(k=i+TURN+1, traced=0; k<=j; k++){
                    jj=k+1;
                    type = ptype[indx[k]+i];
                    if(type){
                      if(fc[i] == fc[k+1] + c[indx[k]+i] + E_ExtLoop(type, -1, -1, P)){
                        traced = i; break;
                      }
                      else if(fc[i] == fc[k+2] + c[indx[k]+i] + E_ExtLoop(type, -1, SAME_STRAND(k,k+1) ? S1[k+1] : -1, P)){
                        traced = i; jj=k+2; break;
                      }
                    }
                    type = ptype[indx[k]+i+1];
                    if(type){
                      if(fc[i] == fc[k+1] + c[indx[k]+i+1] + E_ExtLoop(type, SAME_STRAND(i, i+1) ? S1[i] : -1, -1, P)){
                        traced = i+1; break;
                      }
                      if(k<j){
                        if(fc[i] == fc[k+2] + c[indx[k]+i+1] + E_ExtLoop(type, SAME_STRAND(i, i+1) ? S1[i] : -1, SAME_STRAND(k, k+1) ? S1[k+1] : -1, P)){
                          traced = i+1; jj=k+2; break;
                        }
                      }
                    }
                  }
                  break;
      }

      if (!traced) nrerror("backtrack failed in fc[] 5' of cut");

      sector[++s].i = jj;
      sector[s].j   = j;
      sector[s].ml  = ml;

      j=k; i=traced;
      base_pair2[++b].i = i;
      base_pair2[b].j   = j;
      goto repeat1;
    }

    else { /* true multi-loop backtrack in fML */
      int cij1=INF, ci1j=INF, ci1j1=INF;
      if (fML[indx[j]+i+1]+P->MLbase == fij) { /* 5' end is unpaired */
        sector[++s].i = i+1;
        sector[s].j   = j;
        sector[s].ml  = ml;
        continue;
      }

      tt  = ptype[indx[j]+i];
      cij = c[indx[j]+i];
      switch(dangle_model){
        case 0:   if(fij == cij + E_MLstem(tt, -1, -1, P)){
                    base_pair2[++b].i  = i;
                    base_pair2[b].j    = j;
                    goto repeat1;
                  }
                  break;
        case 2:   if(fij == cij + E_MLstem(tt, (i>1) ? S1[i-1] : -1, (j<length) ? S1[j+1] : -1, P)){
                    base_pair2[++b].i  = i;
                    base_pair2[b].j    = j;
                    goto repeat1;
                  }
                  break;
        default:  if(fij == cij + E_MLstem(tt, -1, -1, P)){
                    base_pair2[++b].i  = i;
                    base_pair2[b].j    = j;
                    goto repeat1;
                  }
                  tt = ptype[indx[j]+i+1];
                  if(fij == c[indx[j]+i+1] + P->MLbase + E_MLstem(tt, S1[i], -1, P)){
                    i++;
                    base_pair2[++b].i  = i;
                    base_pair2[b].j    = j;
                    goto repeat1;
                  }
                  tt = ptype[indx[j-1]+i];
                  if(fij == c[indx[j-1]+i] + P->MLbase + E_MLstem(tt, -1, S1[j], P)){
                    j--;
                    base_pair2[++b].i  = i;
                    base_pair2[b].j    = j;
                    goto repeat1;
                  }
                  tt = ptype[indx[j-1]+i+1];
                  if(fij == c[indx[j-1]+i+1] + 2*P->MLbase + E_MLstem(tt, S1[i], S1[j], P)){
                    i++; j--;
                    base_pair2[++b].i  = i;
                    base_pair2[b].j    = j;
                    goto repeat1;
                  }
                  break;
      }

      /* find next component of multiloop */
      for (k = i+1+TURN; k <= j-2-TURN; k++)
        if (fij == (fML[indx[k]+i]+fML[indx[j]+k+1]))
          break;

      if ((dangle_model==3)&&(k>j-2-TURN)) { /* must be coax stack */
        ml = 2;
        for (k = i+1+TURN; k <= j-2-TURN; k++) {
          type = ptype[indx[k]+i];  type= rtype[type];
          type_2 = ptype[indx[j]+k+1]; type_2= rtype[type_2];
          if (type && type_2)
            if (fij == c[indx[k]+i]+c[indx[j]+k+1]+P->stack[type][type_2]+
                2*P->MLintern[1])
              break;
        }
      }

      sector[++s].i = i;
      sector[s].j   = k;
      sector[s].ml  = ml;
      sector[++s].i = k+1;
      sector[s].j   = j;
      sector[s].ml  = ml;

      if (k>j-2-TURN) nrerror("backtrack failed in fML");
      continue;
    }

  repeat1:

    /*----- begin of "repeat:" -----*/
    if (canonical)  cij = c[indx[j]+i];
    type = ptype[indx[j]+i];

    bonus = 0;

    if ((BP[i]==j)||(BP[i]==-1)||(BP[i]==-2)) bonus -= BONUS;
    if ((BP[j]==-1)||(BP[j]==-3)) bonus -= BONUS;

    if (noLP)
      if (cij == c[indx[j]+i]) {
        /* (i.j) closes canonical structures, thus
           (i+1.j-1) must be a pair                */
        type_2 = ptype[indx[j-1]+i+1]; type_2 = rtype[type_2];
        cij -= P->stack[type][type_2] + bonus;
        base_pair2[++b].i = i+1;
        base_pair2[b].j   = j-1;
        i++; j--;
        canonical=0;
        goto repeat1;
      }
    canonical = 1;


    no_close = (((type==3)||(type==4))&&noGUclosure&&(bonus==0));
    if (SAME_STRAND(i,j)) {
      if (no_close) {
        if (cij == FORBIDDEN) continue;
      } else
        if (cij == E_Hairpin(j-i-1, type, S1[i+1], S1[j-1],string+i-1, P)+bonus)
          continue;
    }
    else {
      int ee = 0;
      if(dangle_model){
        if(cij == E_ExtLoop(rtype[type], SAME_STRAND(j-1,j) ? S1[j-1] : -1, SAME_STRAND(i,i+1) ? S1[i+1] : -1, P)) continue;
      }
      else if(cij == E_ExtLoop(rtype[type], -1, -1, P)) continue;
    }

    for (p = i+1; p <= MIN2(j-2-TURN,i+MAXLOOP+1); p++) {
      int minq;
      minq = j-i+p-MAXLOOP-2;
      if (minq<p+1+TURN) minq = p+1+TURN;
      for (q = j-1; q >= minq; q--) {

        type_2 = ptype[indx[q]+p];
        if (type_2==0) continue;
        type_2 = rtype[type_2];
        if (noGUclosure)
          if (no_close||(type_2==3)||(type_2==4))
            if ((p>i+1)||(q<j-1)) continue;  /* continue unless stack */

        /* energy = oldLoopEnergy(i, j, p, q, type, type_2); */
        if (SAME_STRAND(i,p) && SAME_STRAND(q,j))
          energy = E_IntLoop(p-i-1, j-q-1, type, type_2,
                              S1[i+1], S1[j-1], S1[p-1], S1[q+1], P);
        else {
          energy = E_IntLoop_Co(rtype[type], rtype[type_2], i, j, p, q, cut_point, S1[i+1], S1[j-1], S1[p-1], S1[q+1], dangle_model, P);
        }

        new = energy+c[indx[q]+p]+bonus;
        traced = (cij == new);
        if (traced) {
          base_pair2[++b].i = p;
          base_pair2[b].j   = q;
          i = p, j = q;
          goto repeat1;
        }
      }
    }

    /* end of repeat: --------------------------------------------------*/

    /* (i.j) must close a fake or true multi-loop */
    tt = rtype[type];
    i1 = i+1;
    j1 = j-1;
    /* fake multi-loop */
    if(!SAME_STRAND(i,j)){
      int ii, jj, decomp;
      ii = jj = 0;
      decomp = fc[i1] + fc[j1];
      switch(dangle_model){
        case 0:   if(cij == decomp + E_ExtLoop(tt, -1, -1, P)){
                    ii=i1, jj=j1;
                  }
                  break;
        case 2:   if(cij == decomp + E_ExtLoop(tt, SAME_STRAND(j-1,j) ? S1[j-1] : -1, SAME_STRAND(i,i+1) ? S1[i+1] : -1, P)){
                    ii=i1, jj=j1;
                  }

                  break;
        default:  if(cij == decomp + E_ExtLoop(tt, -1, -1, P)){
                    ii=i1, jj=j1;
                  }
                  else if(cij == fc[i+2] + fc[j-1] + E_ExtLoop(tt, -1, SAME_STRAND(i,i+1) ? S1[i+1] : -1, P)){
                    ii = i+2; jj = j1;
                  }
                  else if(cij == fc[i+1] + fc[j-2] + E_ExtLoop(tt, SAME_STRAND(j-1,j) ? S1[j-1] : -1, -1, P)){
                    ii = i1; jj = j-2;
                  }
                  else if(cij == fc[i+2] + fc[j-2] + E_ExtLoop(tt, SAME_STRAND(j-1,j) ? S1[j-1] : -1, SAME_STRAND(i,i+1) ? S1[i+1] : -1, P)){
                    ii = i+2; jj = j-2;
                  }
                  break;
      }
      if(ii){
        sector[++s].i = ii;
        sector[s].j   = cut_point-1;
        sector[s].ml  = 3;
        sector[++s].i = cut_point;
        sector[s].j   = jj;
        sector[s].ml  = 4;
        continue;
      }
    }

    /* true multi-loop */
    mm = bonus + P->MLclosing;
    sector[s+1].ml  = sector[s+2].ml = 1;
    int ml0   = E_MLstem(tt, -1, -1, P);
    int ml5   = E_MLstem(tt, SAME_STRAND(j-1,j) ? S1[j-1] : -1, -1, P);
    int ml3   = E_MLstem(tt, -1, SAME_STRAND(i,i+1) ? S1[i+1] : -1, P);
    int ml53  = E_MLstem(tt, SAME_STRAND(j-1,j) ? S1[j-1] : -1, SAME_STRAND(i,i+1) ? S1[i+1] : -1, P);
    for (traced = 0, k = i+2+TURN; k < j-2-TURN; k++) {
      switch(dangle_model){
        case 0:   /* no dangles */
                  if(cij == mm + fML[indx[k]+i+1] + fML[indx[j-1]+k+1] + ml0)
                    traced = i+1;
                  break;
        case 2:   /*double dangles */
                  if(cij == mm + fML[indx[k]+i+1] + fML[indx[j-1]+k+1] + ml53)
                    traced = i+1;
                  break;
        default:  /* normal dangles */
                  if(cij == mm + fML[indx[k]+i+1] + fML[indx[j-1]+k+1] + ml0){
                    traced = i+1;
                    break;
                  }
                  else if (cij == fML[indx[k]+i+2] + fML[indx[j-1]+k+1] + ml3 + mm + P->MLbase){
                    traced = i1 = i+2;
                    break;
                  }
                  else if (cij == fML[indx[k]+i+1] + fML[indx[j-2]+k+1] + ml5 + mm + P->MLbase){
                    traced = i1 = i+1; j1 = j-2;
                    break;
                  }
                  else if (cij == fML[indx[k]+i+2] + fML[indx[j-2]+k+1] + ml53 + mm + 2*P->MLbase){
                    traced = i1 = i+2; j1 = j-2;
                    break;
                  }
                  break;
      }
      if(traced) break;
      /* coaxial stacking of (i.j) with (i+1.k) or (k.j-1) */
      /* use MLintern[1] since coax stacked pairs don't get TerminalAU */
      if (dangle_model==3) {
        int en;
        type_2 = ptype[indx[k]+i+1]; type_2 = rtype[type_2];
        if (type_2) {
          en = c[indx[k]+i+1]+P->stack[type][type_2]+fML[indx[j-1]+k+1];
          if (cij == en+2*P->MLintern[1]+P->MLclosing) {
            ml = 2;
            sector[s+1].ml  = 2;
            break;
          }
        }
        type_2 = ptype[indx[j-1]+k+1]; type_2 = rtype[type_2];
        if (type_2) {
          en = c[indx[j-1]+k+1]+P->stack[type][type_2]+fML[indx[k]+i+1];
          if (cij == en+2*P->MLintern[1]+P->MLclosing) {
            sector[s+2].ml = 2;
            break;
          }
        }
      }
    }
    if (k<=j-3-TURN) { /* found the decomposition */
      sector[++s].i = i1;
      sector[s].j   = k;
      sector[++s].i = k+1;
      sector[s].j   = j1;
    } else {
#if 0
      /* Y shaped ML loops don't work yet */
      if (dangle_model==3) {
        /* (i,j) must close a Y shaped ML loop with coax stacking */
        if (cij == fML[indx[j-2]+i+2] + mm + d3 + d5 + P->MLbase + P->MLbase) {
          i1 = i+2;
          j1 = j-2;
        } else if (cij ==  fML[indx[j-2]+i+1] + mm + d5 + P->MLbase)
          j1 = j-2;
        else if (cij ==  fML[indx[j-1]+i+2] + mm + d3 + P->MLbase)
          i1 = i+2;
        else /* last chance */
          if (cij != fML[indx[j-1]+i+1] + mm + P->MLbase)
            fprintf(stderr,  "backtracking failed in repeat");
        /* if we arrive here we can express cij via fML[i1,j1]+dangles */
        sector[++s].i = i1;
        sector[s].j   = j1;
      }
      else
#endif
        nrerror("backtracking failed in repeat");
    }

  } /* end >> while (s>0) << */

  base_pair2[0].i = b;    /* save the total number of base pairs */
}

PRIVATE void free_end(int *array, int i, int start) {
  int inc, type, energy, length, j, left, right;
  int dangle_model = P->model_details.dangles;

  inc = (i>start)? 1:-1;
  length = S[0];

  if (i==start) array[i]=0;
  else array[i] = array[i-inc];
  if (inc>0) {
    left = start; right=i;
  } else {
    left = i; right = start;
  }

  for (j=start; inc*(i-j)>TURN; j+=inc) {
    int d3, d5, ii, jj;
    short si, sj;
    if (i>j) { ii = j; jj = i;} /* inc>0 */
    else     { ii = i; jj = j;} /* inc<0 */
    type = ptype[indx[jj]+ii];
    if (type) {  /* i is paired with j */
      si = (ii>1)       && SAME_STRAND(ii-1,ii) ? S1[ii-1] : -1;
      sj = (jj<length)  && SAME_STRAND(jj,jj+1) ? S1[jj+1] : -1;
      energy = c[indx[jj]+ii];
      switch(dangle_model){
        case 0:   array[i] = MIN2(array[i], array[j-inc] + energy + E_ExtLoop(type, -1, -1, P));
                  break;
        case 2:   array[i] = MIN2(array[i], array[j-inc] + energy + E_ExtLoop(type, si, sj, P));
                  break;
        default:  array[i] = MIN2(array[i], array[j-inc] + energy + E_ExtLoop(type, -1, -1, P));
                  if(inc > 0){
                    if(j > left)
                    array[i] = MIN2(array[i], array[j-2] + energy + E_ExtLoop(type, si, -1, P));
                  }
                  else if(j < right)
                    array[i] = MIN2(array[i], array[j+2] + energy + E_ExtLoop(type, -1, sj, P));
                  break;
      }
    }
    if (dangle_model%2==1) {
      /* interval ends in a dangle (i.e. i-inc is paired) */
      if (i>j) { ii = j; jj = i-1;} /* inc>0 */
      else     { ii = i+1; jj = j;} /* inc<0 */
      type = ptype[indx[jj]+ii];
      if (!type) continue;

      si = (ii > left)  && SAME_STRAND(ii-1,ii) ? S1[ii-1] : -1;
      sj = (jj < right) && SAME_STRAND(jj,jj+1) ? S1[jj+1] : -1;
      energy = c[indx[jj]+ii];
      if(inc>0)
        array[i] = MIN2(array[i], array[j - inc] + energy + E_ExtLoop(type, -1, sj, P));
      else
        array[i] = MIN2(array[i], array[j - inc] + energy + E_ExtLoop(type, si, -1, P));
      if(j!= start){ /* dangle_model on both sides */
        array[i] = MIN2(array[i], array[j-2*inc] + energy + E_ExtLoop(type, si, sj, P));
      }
    }
  }
}

PUBLIC void update_cofold_params(void){
  update_cofold_params_par(NULL);
}

PUBLIC void update_cofold_params_par(paramT *parameters){
  if(P) free(P);

  if(parameters){
    P = get_parameter_copy(parameters);
  } else {
    model_detailsT md;
    set_model_details(&md);
    P = get_scaled_parameters(temperature, md);
  }
  make_pair_matrix();
  if (init_length < 0) init_length=0;
}

/*---------------------------------------------------------------------------*/

PRIVATE void make_ptypes(const short *S, const char *structure) {
  int n,i,j,k,l;
  int noLP = P->model_details.noLP;

  n=S[0];
  for (k=1; k<n-TURN; k++)
    for (l=1; l<=2; l++) {
      int type,ntype=0,otype=0;
      i=k; j = i+TURN+l; if (j>n) continue;
      type = pair[S[i]][S[j]];
      while ((i>=1)&&(j<=n)) {
        if ((i>1)&&(j<n)) ntype = pair[S[i-1]][S[j+1]];
        if (noLP && (!otype) && (!ntype))
          type = 0; /* i.j can only form isolated pairs */
        ptype[indx[j]+i] = (char) type;
        otype =  type;
        type  = ntype;
        i--; j++;
      }
    }

  if (struct_constrained && (structure != NULL))
    constrain_ptypes(structure, (unsigned int)n, ptype, BP, TURN, 0);
}

PUBLIC void get_monomere_mfes(float *e1, float *e2) {
  /*exports monomere free energies*/
  *e1 = mfe1;
  *e2 = mfe2;
}

PRIVATE void backtrack(const char *sequence) {
  /*routine to call backtrack_co from 1 to n, backtrack type??*/
  backtrack_co(sequence, 0,0);
}

PRIVATE comp_pair(const void *A, const void *B) {
  bondT *x,*y;
  int ex, ey;
  x = (bondT *) A;
  y = (bondT *) B;
  ex = c[indx[x->j]+x->i]+c[indx[x->i+length]+x->j];
  ey = c[indx[y->j]+y->i]+c[indx[y->i+length]+y->j];
  if (ex>ey) return 1;
  if (ex<ey) return -1;
  return (indx[x->j]+x->i - indx[y->j]+y->i);
}

PUBLIC SOLUTION *zukersubopt(const char *string) {
  return zukersubopt_par(string, NULL);
}

PUBLIC SOLUTION *zukersubopt_par(const char *string, paramT *parameters){

/* Compute zuker suboptimal. Here, we're abusing the cofold() code
   "double" sequence, compute dimerarray entries, track back every base pair.
   This is slightly wasteful compared to the normal solution */

  char      *doubleseq, *structure, *mfestructure, **todo;
  int       i, j, k, counter, num_pairs, psize, p;
  float     energy, mfenergy;
  SOLUTION  *zukresults;
  bondT     *pairlist;

  num_pairs       = counter = 0;
  zuker           = 1;
  length          = (int)strlen(string);
  doubleseq       = (char *)space((2*length+1)*sizeof(char));
  mfestructure    = (char *) space((unsigned) 2*length+1);
  structure       = (char *) space((unsigned) 2*length+1);
  zukresults      = (SOLUTION *)space(((length*(length-1))/2)*sizeof(SOLUTION));
  mfestructure[0] = '\0';
  BP              = (int *)space(sizeof(int)*(length+2));

  /* double the sequence */
  strcpy(doubleseq,string);
  strcat(doubleseq,string);
  cut_point = length + 1;

  /* get mfe and do forward recursion */
#ifdef _OPENMP
  /* always init everything since all global static variables are uninitialized when entering a thread */
  init_cofold(2 * length, parameters);
#else
  if(parameters) init_cofold(2 * length, parameters);
  else if ((2 * length) > init_length) init_cofold(2 * length, parameters);
  else if (fabs(P->temperature - temperature)>1e-6) update_cofold_params_par(parameters);
#endif


  S     = encode_sequence(doubleseq, 0);
  S1    = encode_sequence(doubleseq, 1);
  S1[0] = S[0]; /* store length at pos. 0 */
  make_ptypes(S, NULL); /* no constraint folding possible (yet?) with zukersubopt */

  (void)fill_arrays(doubleseq);
  mfenergy = f5[cut_point-1];

  psize     = length;
  pairlist  = (bondT *) space(sizeof(bondT)*(psize+1));
  todo      = (char **) space(sizeof(char *)*(length+1));
  for (i=1; i<length; i++) {
    todo[i] = (char *) space(sizeof(char)*(length+1));
  }

  /* Make a list of all base pairs */
  for (i=1; i<length; i++) {
    for (j=i+TURN2+1/*??*/; j<=length; j++) {
      if (ptype[indx[j]+i]==0) continue;
      if (num_pairs>=psize) {
        psize = 1.2*psize + 32;
        pairlist = xrealloc(pairlist, sizeof(bondT)*(psize+1));
      }
      pairlist[num_pairs].i   = i;
      pairlist[num_pairs++].j = j;
      todo[i][j]=1;
    }
  }
  qsort(pairlist, num_pairs, sizeof(bondT), comp_pair);

  for (p=0; p<num_pairs; p++) {
    i=pairlist[p].i;
    j=pairlist[p].j;
    if (todo[i][j]) {
      int k;
      sector[1].i   = i;
      sector[1].j   = j;
      sector[1].ml  = 2;
      backtrack_co(doubleseq, 1,0);
      sector[1].i   = j;
      sector[1].j   = i + length;
      sector[1].ml  = 2;
      backtrack_co(doubleseq, 1,base_pair2[0].i);
      energy = c[indx[j]+i]+c[indx[i+length]+j];
      parenthesis_zuker(structure, base_pair2, length);
      zukresults[counter].energy      = energy;
      zukresults[counter++].structure = strdup(structure);
      for (k = 1; k <= base_pair2[0].i; k++) { /* mark all pairs in structure as done */
        int x,y;
        x=base_pair2[k].i;
        y=base_pair2[k].j;
        if (x>length) x-=length;
        if (y>length) y-=length;
        if (x>y) {
          int temp;
          temp=x; x=y; y=temp;
        }
        todo[x][y] = 0;
      }
    }
  }
  /*free zeugs*/
  free(pairlist);
  for (i=1; i<length; i++)
    free(todo[i]);
  free(todo);
  free(structure);
  free(mfestructure);
  free(doubleseq);
  zuker=0;
  free(S); free(S1); free(BP);
  return zukresults;
}


/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

PUBLIC void initialize_cofold(int length){ /* DO NOTHING */ }
