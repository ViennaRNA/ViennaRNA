/** \file eval.c */


/*
                  Free energy evaluation

                  c Ivo Hofacker, Chrisoph Flamm
                  original implementation by
                  Walter Fontana

                  modification by Ronny Lorenz

                  Vienna RNA package
*/

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/utils.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/loop_energies.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/eval.h"

#define SAME_STRAND(I,J)  (((I)>=cut_point)||((J)<cut_point))

/*
#################################
# GLOBAL VARIABLES              #
#################################
*/
PUBLIC  int logML     = 0;  /* if nonzero use logarithmic ML energy in energy_of_struct */
PUBLIC  int cut_point = -1; /* set to first pos of second seq for cofolding */
PUBLIC  int eos_debug = 0;  /* verbose info from energy_of_struct */

/*
#################################
# PRIVATE VARIABLES             #
#################################
*/

PRIVATE short   *S = NULL, *S1 = NULL;
PRIVATE paramT  *P            = NULL;
PRIVATE short   *pair_table   = NULL;

#ifdef _OPENMP

/* NOTE: all variables are assumed to be uninitialized if they are declared as threadprivate
         thus we have to initialize them before usage by a seperate function!
         OR: use copyin in the PARALLEL directive!
         e.g.:
         #pragma omp parallel for copyin(P, init_length, min_hairpin)
*/
#pragma omp threadprivate(S, S1, P, pair_table)

#endif

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/
PRIVATE int   stack_energy(int i, const char *string, int verbostiy_level);
PRIVATE int   energy_of_extLoop_pt(int i, short *pair_table);
PRIVATE int   energy_of_ml_pt(int i, short *pt);
PRIVATE int   ML_Energy(int i, int is_extloop);
PRIVATE int   cut_in_loop(int i);


/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PRIVATE void  update_params(paramT *parameters){
  if(P) free(P);
  if(parameters){
    P = get_parameter_copy(parameters);
  } else {
    model_detailsT md;
    set_model_details(&md);
    P = get_scaled_parameters(temperature, md);
  }
  fill_pair_matrices(&(P->model_details));
/*
  make_pair_matrix();
*/
}

PUBLIC  float energy_of_structure(const char *string,
                                  const char *structure,
                                  int verbosity_level){

  return energy_of_struct_par(string, structure, NULL, verbosity_level);
}

PUBLIC float energy_of_struct_par(const char *string,
                                  const char *structure,
                                  paramT *parameters,
                                  int verbosity_level){
  int   energy;
  short *ss, *ss1;

  update_params(parameters);

  if (strlen(structure)!=strlen(string))
    nrerror("energy_of_struct: string and structure have unequal length");

  /* save the S and S1 pointers in case they were already in use */
  ss = S; ss1 = S1;
  S   = get_sequence_encoding(string, 0, &(P->model_details));
  S1  = get_sequence_encoding(string, 1, &(P->model_details));

  pair_table = vrna_pt_get(structure);

  energy = energy_of_struct_pt_par(string, pair_table, S, S1, parameters, verbosity_level);

  free(pair_table);
  free(S); free(S1);
  S=ss; S1=ss1;
  return  (float) energy/100.;
}

/*  returns a correction term that may be added to the energy retrieved
    from energy_of_struct_par() to correct misinterpreted loops. This
    correction is necessary since energy_of_struct_par() will forget 
    about the existance of gquadruplexes and just treat them as unpaired
    regions.

    recursive variant
*/
PRIVATE int en_corr_of_loop_gquad(int i,
                                  int j,
                                  const char *string,
                                  const char *structure,
                                  short *pt,
                                  int *loop_idx,
                                  const short *s1){

  int pos, energy, p, q, r, s, u, type, type2;
  int L, l[3];
  int *rtype;
  model_detailsT  *md;

  md    = &(P->model_details);
  rtype = &(md->rtype[0]);


  energy = 0;
  q = i;
  while((pos = parse_gquad(structure + q-1, &L, l)) > 0){
    q += pos-1;
    p = q - 4*L - l[0] - l[1] - l[2] + 1;
    if(q > j) break;
    /* we've found the first g-quadruplex at position [p,q] */
    energy += E_gquad(L, l, P);
    /* check if it's enclosed in a base pair */
    if(loop_idx[p] == 0){ q++; continue; /* g-quad in exterior loop */}
    else{
      energy += E_MLstem(0, -1, -1, P); /*  do not forget to remove this energy if
                                            the gquad is the only one surrounded by
                                            the enclosing pair
                                        */

      /*  find its enclosing pair */
      int num_elem, num_g, elem_i, elem_j, up_mis;
      num_elem  = 0;
      num_g     = 1;
      r         = p - 1;
      up_mis    = q - p + 1;

      /* seek for first pairing base located 5' of the g-quad */
      for(r = p - 1; !pt[r] && (r >= i); r--);
      if(r < i) nrerror("this should not happen");

      if(r < pt[r]){ /* found the enclosing pair */
        s = pt[r];
      } else {
        num_elem++;
        elem_i = pt[r];
        elem_j = r;
        r = pt[r]-1 ;
        /* seek for next pairing base 5' of r */
        for(; !pt[r] && (r >= i); r--);
        if(r < i) nrerror("so nich");
        if(r < pt[r]){ /* found the enclosing pair */
          s = pt[r];
        } else {
          /* hop over stems and unpaired nucleotides */
          while((r > pt[r]) && (r >= i)){
            if(pt[r]){ r = pt[r]; num_elem++;}
            r--;
          }
          if(r < i) nrerror("so nich");
          s = pt[r]; /* found the enclosing pair */
        }
      }
      /* now we have the enclosing pair (r,s) */

      u = q+1;
      /* we know everything about the 5' part of this loop so check the 3' part */
      while(u<s){
        if(structure[u-1] == '.') u++;
        else if (structure[u-1] == '+'){ /* found another gquad */
          pos = parse_gquad(structure + u - 1, &L, l);
          if(pos > 0){
            energy += E_gquad(L, l, P) + E_MLstem(0, -1, -1, P);
            up_mis += pos;
            u += pos;
            num_g++;
          }
        } else { /* we must have found a stem */
          if(!(u < pt[u])) nrerror("wtf!");
          num_elem++; elem_i = u; elem_j = pt[u];
          energy += en_corr_of_loop_gquad(u, pt[u], string, structure, pt, loop_idx, s1);
          u = pt[u] + 1;
        }
      }
      if(u!=s) nrerror("what the hell");
      else{ /* we are done since we've found no other 3' structure element */
        switch(num_elem){
          /* g-quad was misinterpreted as hairpin closed by (r,s) */
          case 0:   /* if(num_g == 1)
                      if((p-r-1 == 0) || (s-q-1 == 0))
                        nrerror("too few unpaired bases");
                    */
                    type = md->pair[s1[r]][s1[s]];
                    if(dangles == 2)
                      energy += P->mismatchI[type][s1[r+1]][s1[s-1]];
                    if(type > 2)
                      energy += P->TerminalAU;
                    energy += P->internal_loop[s - r - 1 - up_mis];
                    energy -= E_MLstem(0, -1, -1, P);
                    energy -= E_Hairpin(s - r - 1,
                                        type,
                                        s1[r + 1],
                                        s1[s - 1],
                                        string + r - 1,
                                        P);
                    break;
          /* g-quad was misinterpreted as interior loop closed by (r,s) with enclosed pair (elem_i, elem_j) */
          case 1:   type = md->pair[s1[r]][s1[s]];
                    type2 = md->pair[s1[elem_i]][s1[elem_j]];
                    energy += P->MLclosing
                              + E_MLstem(rtype[type], s1[s-1], s1[r+1], P)
                              + (elem_i - r - 1 + s - elem_j - 1 - up_mis) * P->MLbase
                              + E_MLstem(type2, s1[elem_i-1], s1[elem_j+1], P);
                    energy -= E_IntLoop(elem_i - r - 1,
                                        s - elem_j - 1,
                                        type,
                                        rtype[type2],
                                        s1[r + 1],
                                        s1[s - 1],
                                        s1[elem_i - 1],
                                        s1[elem_j + 1],
                                        P);
                    break;
          /* gquad was misinterpreted as unpaired nucleotides in a multiloop */
          default:  energy -= (up_mis) * P->MLbase;
                    break;
        }
      }
      q = s+1;
    }
  }
  return energy;
}

PUBLIC float
energy_of_gquad_structure(const char *string,
                          const char *structure,
                          int verbosity_level){

  return energy_of_gquad_struct_par(string, structure, NULL, verbosity_level);
}

PUBLIC float
energy_of_gquad_struct_par( const char *string,
                            const char *structure,
                            paramT *parameters,
                            int verbosity_level){

  int   energy, gge, *loop_idx;
  short *ss, *ss1;

  update_fold_params_par(parameters);

  if (strlen(structure)!=strlen(string))
    nrerror("energy_of_struct: string and structure have unequal length");

  /* save the S and S1 pointers in case they were already in use */
  ss = S; ss1 = S1;
  S   = get_sequence_encoding(string, 0, &(P->model_details));
  S1  = get_sequence_encoding(string, 1, &(P->model_details));

  /* the pair_table looses every information about the gquad position
     thus we have to find add the energy contributions for each loop
     that contains a gquad by ourself, substract all miscalculated
     contributions, i.e. loops that actually contain a gquad, from
     energy_of_structure_pt()
  */
  pair_table  = vrna_pt_get(structure);
  energy      = energy_of_structure_pt(string, pair_table, S, S1, verbosity_level);

  loop_idx    = vrna_get_loop_index(pair_table);
  gge         = en_corr_of_loop_gquad(1, S[0], string, structure, pair_table, loop_idx, S1);
  energy     += gge;

  free(pair_table);
  free(loop_idx);
  free(S); free(S1);
  S=ss; S1=ss1;
  return  (float) energy/100.;
}

PUBLIC int energy_of_structure_pt(const char *string,
                                  short *ptable,
                                  short *s,
                                  short *s1,
                                  int verbosity_level){
  return energy_of_struct_pt_par(string, ptable, s, s1, NULL, verbosity_level);
}

PUBLIC int energy_of_struct_pt_par( const char *string,
                                    short *ptable,
                                    short *s,
                                    short *s1,
                                    paramT *parameters,
                                    int verbosity_level){
  /* auxiliary function for kinfold,
     for most purposes call energy_of_struct instead */

  int   i, length, energy;
  short *ss, *ss1;

  update_params(parameters);

  pair_table = ptable;
  ss  = S;
  ss1 = S1;
  S = s;
  S1 = s1;

  length = S[0];
/*
  energy =  backtrack_type=='M' ? ML_Energy(0, 0) : ML_Energy(0, 1);
*/
  energy =  backtrack_type=='M' ? energy_of_ml_pt(0, ptable) : energy_of_extLoop_pt(0, ptable);
  if (verbosity_level>0)
    printf("External loop                           : %5d\n", energy);
  for (i=1; i<=length; i++) {
    if (pair_table[i]==0) continue;
    energy += stack_energy(i, string, verbosity_level);
    i=pair_table[i];
  }
  for (i=1; !SAME_STRAND(i,length); i++) {
    if (!SAME_STRAND(i,pair_table[i])) {
      energy+=P->DuplexInit;
      break;
    }
  }
  S   = ss;
  S1  = ss1;
  return energy;
}

PUBLIC float energy_of_circ_structure(const char *string,
                                      const char *structure,
                                      int verbosity_level){
  return energy_of_circ_struct_par(string, structure, NULL, verbosity_level);
}

PUBLIC  float energy_of_circ_struct_par(  const char *string,
                                          const char *structure,
                                          paramT *parameters,
                                          int verbosity_level){

  int   i, j, length, energy=0, en0, degree=0, type;
  short *ss, *ss1;

  update_params(parameters);

  int dangle_model = P->model_details.dangles;

  if (strlen(structure)!=strlen(string))
    nrerror("energy_of_struct: string and structure have unequal length");

  /* save the S and S1 pointers in case they were already in use */
  ss = S; ss1 = S1;
  S   = get_sequence_encoding(string, 0, &(P->model_details));
  S1  = get_sequence_encoding(string, 1, &(P->model_details));

  pair_table = vrna_pt_get(structure);

  length = S[0];

  for (i=1; i<=length; i++) {
    if (pair_table[i]==0) continue;
    degree++;
    energy += stack_energy(i, string, verbosity_level);
    i=pair_table[i];
  }

  if (degree==0) return 0.;
  for (i=1; pair_table[i]==0; i++);
  j = pair_table[i];
  type = P->model_details.pair[S[j]][S[i]];
  if (type==0) type=7;
  if (degree==1) {
    char loopseq[10];
    int u, si1, sj1;
    for (i=1; pair_table[i]==0; i++);
    u = length-j + i-1;
    if (u<7) {
      strcpy(loopseq , string+j-1);
      strncat(loopseq, string, i);
    }
    si1 = (i==1)?S1[length] : S1[i-1];
    sj1 = (j==length)?S1[1] : S1[j+1];
    en0 = E_Hairpin(u, type, sj1, si1, loopseq, P);
  } else
    if (degree==2) {
      int p,q, u1,u2, si1, sq1, type_2;
      for (p=j+1; pair_table[p]==0; p++);
      q=pair_table[p];
      u1 = p-j-1;
      u2 = i-1 + length-q;
      type_2 = P->model_details.pair[S[q]][S[p]];
      if (type_2==0) type_2=7;
      si1 = (i==1)? S1[length] : S1[i-1];
      sq1 = (q==length)? S1[1] : S1[q+1];
      en0 = E_IntLoop(u1, u2, type, type_2,
                       S1[j+1], si1, S1[p-1], sq1,P);
    } else { /* degree > 2 */
      en0 = ML_Energy(0, 0) - P->MLintern[0];
      if (dangle_model) {
        int d5, d3;
        if (pair_table[1]) {
          j = pair_table[1];
          type = P->model_details.pair[S[1]][S[j]];
          if (dangle_model==2)
            en0 += P->dangle5[type][S1[length]];
          else { /* dangle_model==1 */
            if (pair_table[length]==0) {
              d5 = P->dangle5[type][S1[length]];
              if (pair_table[length-1]!=0) {
                int tt;
                tt = P->model_details.pair[S[pair_table[length-1]]][S[length-1]];
                d3 = P->dangle3[tt][S1[length]];
                if (d3<d5) d5 = 0;
                else d5 -= d3;
              }
              en0 += d5;
            }
          }
        }
        if (pair_table[length]) {
          i = pair_table[length];
          type = P->model_details.pair[S[i]][S[length]];
          if (dangle_model==2)
            en0 += P->dangle3[type][S1[1]];
          else { /* dangle_model==1 */
            if (pair_table[1]==0) {
              d3 = P->dangle3[type][S1[1]];
              if (pair_table[2]) {
                int tt;
                tt = P->model_details.pair[S[2]][S[pair_table[2]]];
                d5 = P->dangle5[tt][1];
                if (d5<d3) d3=0;
                else d3 -= d5;
              }
              en0 += d3;
            }
          }
        }
      }
    }

  if (verbosity_level>0)
    printf("External loop                           : %5d\n", en0);
  energy += en0;
  /* fprintf(stderr, "ext loop degree %d tot %d\n", degree, energy); */
  free(S); free(S1);
  S=ss; S1=ss1;
  return  (float) energy/100.0;
}

/*---------------------------------------------------------------------------*/

PRIVATE int stack_energy( int i,
                          const char *string,
                          int verbosity_level){

  /* calculate energy of substructure enclosed by (i,j) */
  int ee, energy = 0;
  int j, p, q, type;
  int *rtype = &(P->model_details.rtype[0]);

  j=pair_table[i];
  type = P->model_details.pair[S[i]][S[j]];
  if (type==0) {
    type=7;
    if (verbosity_level>=0)
      fprintf(stderr,"WARNING: bases %d and %d (%c%c) can't pair!\n", i, j,
              string[i-1],string[j-1]);
  }

  p=i; q=j;
  while (p<q) { /* process all stacks and interior loops */
    int type_2;
    while (pair_table[++p]==0);
    while (pair_table[--q]==0);
    if ((pair_table[q]!=(short)p)||(p>q)) break;
    type_2 = P->model_details.pair[S[q]][S[p]];
    if (type_2==0) {
      type_2=7;
      if (verbosity_level>=0)
        fprintf(stderr,"WARNING: bases %d and %d (%c%c) can't pair!\n", p, q,
                string[p-1],string[q-1]);
    }
    /* energy += LoopEnergy(i, j, p, q, type, type_2); */
    if ( SAME_STRAND(i,p) && SAME_STRAND(q,j) )
      ee = E_IntLoop(p-i-1, j-q-1, type, type_2, S1[i+1], S1[j-1], S1[p-1], S1[q+1],P);
    else
      ee = energy_of_extLoop_pt(cut_in_loop(i), pair_table);
    if (verbosity_level>0)
      printf("Interior loop (%3d,%3d) %c%c; (%3d,%3d) %c%c: %5d\n",
             i,j,string[i-1],string[j-1],p,q,string[p-1],string[q-1], ee);
    energy += ee;
    i=p; j=q; type = rtype[type_2];
  } /* end while */

  /* p,q don't pair must have found hairpin or multiloop */

  if (p>q) {                       /* hair pin */
    if (SAME_STRAND(i,j))
      ee = E_Hairpin(j-i-1, type, S1[i+1], S1[j-1], string+i-1, P);
    else
      ee = energy_of_extLoop_pt(cut_in_loop(i), pair_table);
    energy += ee;
    if (verbosity_level>0)
      printf("Hairpin  loop (%3d,%3d) %c%c              : %5d\n",
             i, j, string[i-1],string[j-1], ee);

    return energy;
  }

  /* (i,j) is exterior pair of multiloop */
  while (p<j) {
    /* add up the contributions of the substructures of the ML */
    energy += stack_energy(p, string, verbosity_level);
    p = pair_table[p];
    /* search for next base pair in multiloop */
    while (pair_table[++p]==0);
  }
  {
    int ii;
    ii = cut_in_loop(i);
    ee = (ii==0) ? energy_of_ml_pt(i, pair_table) : energy_of_extLoop_pt(ii, pair_table);
  }
  energy += ee;
  if (verbosity_level>0)
    printf("Multi    loop (%3d,%3d) %c%c              : %5d\n",
           i,j,string[i-1],string[j-1],ee);

  return energy;
}

/*---------------------------------------------------------------------------*/



/**
*** Calculate the energy contribution of
*** stabilizing dangling-ends/mismatches
*** for all stems branching off the exterior
*** loop
**/
PRIVATE int energy_of_extLoop_pt(int i, short *pair_table) {
  int energy, mm5, mm3;
  int p, q, q_prev;
  int length = (int)pair_table[0];

  /* helper variables for dangles == 1 case */
  int E3_available;  /* energy of 5' part where 5' mismatch is available for current stem */
  int E3_occupied;   /* energy of 5' part where 5' mismatch is unavailable for current stem */

  int dangle_model = P->model_details.dangles;

  /* initialize vars */
  energy      = 0;
  p           = (i==0) ? 1 : i;
  q_prev      = -1;

  if(dangle_model%2 == 1){
    E3_available = INF;
    E3_occupied  = 0;
  }

  /* seek to opening base of first stem */
  while(p <= length && !pair_table[p]) p++;

  while(p < length){
    int tt;
    /* p must have a pairing partner */
    q  = (int)pair_table[p];
    /* get type of base pair (p,q) */
    tt = P->model_details.pair[S[p]][S[q]];
    if(tt==0) tt=7;

    switch(dangle_model){
      /* no dangles */
      case 0:   energy += E_ExtLoop(tt, -1, -1, P);
                break;

      /* the beloved double dangles */
      case 2:   mm5 = ((SAME_STRAND(p-1,p)) && (p>1))       ? S1[p-1] : -1;
                mm3 = ((SAME_STRAND(q,q+1)) && (q<length))  ? S1[q+1] : -1;
                energy += E_ExtLoop(tt, mm5, mm3, P);
                break;

      default:  {
                  int tmp;
                  if(q_prev + 2 < p){
                    E3_available = MIN2(E3_available, E3_occupied);
                    E3_occupied  = E3_available;
                  }
                  mm5 = ((SAME_STRAND(p-1,p)) && (p>1) && !pair_table[p-1])       ? S1[p-1] : -1;
                  mm3 = ((SAME_STRAND(q,q+1)) && (q<length) && !pair_table[q+1])  ? S1[q+1] : -1;
                  tmp = MIN2(
                                                E3_occupied  + E_ExtLoop(tt, -1, mm3, P),
                                                E3_available + E_ExtLoop(tt, mm5, mm3, P)
                                              );
                  E3_available =       MIN2(
                                                E3_occupied  + E_ExtLoop(tt, -1, -1, P),
                                                E3_available + E_ExtLoop(tt, mm5, -1, P)
                                              );
                  E3_occupied = tmp;
                }
                break;

    } /* end switch dangle_model */
    /* seek to the next stem */
    p = q + 1;
    q_prev = q;
    while (p <= length && !pair_table[p]) p++;
    if(p==i) break; /* cut was in loop */
  }

  if(dangle_model%2 == 1)
    energy = MIN2(E3_occupied, E3_available);

  return energy;
}

/**
*** i is the 5'-base of the closing pair
***
*** since each helix can coaxially stack with at most one of its
*** neighbors we need an auxiliarry variable  cx_energy
*** which contains the best energy given that the last two pairs stack.
*** energy  holds the best energy given the previous two pairs do not
*** stack (i.e. the two current helices may stack)
*** We don't allow the last helix to stack with the first, thus we have to
*** walk around the Loop twice with two starting points and take the minimum
***/
PRIVATE int energy_of_ml_pt(int i, short *pt){

  int energy, cx_energy, tmp, tmp2, best_energy=INF;
  int i1, j, p, q, q_prev, q_prev2, u, x, type, count, mm5, mm3, tt, ld5, new_cx, dang5, dang3, dang;
  int mlintern[NBPAIRS+1], mlclosing, mlbase;

  /* helper variables for dangles == 1|5 case */
  int E_mm5_available;  /* energy of 5' part where 5' mismatch of current stem is available */
  int E_mm5_occupied;   /* energy of 5' part where 5' mismatch of current stem is unavailable */
  int E2_mm5_available; /* energy of 5' part where 5' mismatch of current stem is available with possible 3' dangle for enclosing pair (i,j) */
  int E2_mm5_occupied;  /* energy of 5' part where 5' mismatch of current stem is unavailable with possible 3' dangle for enclosing pair (i,j) */
  int length = (int)pt[0];
  int dangle_model = P->model_details.dangles;
  int *rtype = &(P->model_details.rtype[0]);

  if(i >= pt[i])
    nrerror("energy_of_ml_pt: i is not 5' base of a closing pair!");

  j = (int)pt[i];

  /* init the variables */
  energy      = 0;
  p           = i+1;
  q_prev      = i-1;
  q_prev2     = i;

  for (x = 0; x <= NBPAIRS; x++) mlintern[x] = P->MLintern[x];
  mlclosing = P->MLclosing; mlbase = P->MLbase;

  /* seek to opening base of first stem */
  while(p <= j && !pair_table[p]) p++;
  u = p - i - 1;

  switch(dangle_model){
    case 0:   while(p < j){
                /* p must have a pairing partner */
                q  = (int)pair_table[p];
                /* get type of base pair (p,q) */
                tt = P->model_details.pair[S[p]][S[q]];
                if(tt==0) tt=7;
                energy += E_MLstem(tt, -1, -1, P);
                /* seek to the next stem */
                p = q + 1;
                q_prev = q_prev2 = q;
                while (p <= j && !pair_table[p]) p++;
                u += p - q - 1; /* add unpaired nucleotides */
              }
              /* now lets get the energy of the enclosing stem */
              type = P->model_details.pair[S[j]][S[i]]; if (type==0) type=7;
              energy += E_MLstem(type, -1, -1, P);
              break;

    case 2:   while(p < j){
                /* p must have a pairing partner */
                q  = (int)pair_table[p];
                /* get type of base pair (p,q) */
                tt = P->model_details.pair[S[p]][S[q]];
                if(tt==0) tt=7;
                mm5 = (SAME_STRAND(p-1,p))  ? S1[p-1] : -1;
                mm3 = (SAME_STRAND(q,q+1))  ? S1[q+1] : -1;
                energy += E_MLstem(tt, mm5, mm3, P);
                /* seek to the next stem */
                p = q + 1;
                q_prev = q_prev2 = q;
                while (p <= j && !pair_table[p]) p++;
                u += p - q - 1; /* add unpaired nucleotides */
              }
              type = P->model_details.pair[S[j]][S[i]]; if (type==0) type=7;
              mm5 = ((SAME_STRAND(j-1,j)) && !pair_table[j-1])  ? S1[j-1] : -1;
              mm3 = ((SAME_STRAND(i,i+1)) && !pair_table[i+1])  ? S1[i+1] : -1;
              energy += E_MLstem(type, S1[j-1], S1[i+1], P);
              break;

    case 3:   /* we treat helix stacking different */
              for (count=0; count<2; count++) { /* do it twice */
                ld5 = 0; /* 5' dangle energy on prev pair (type) */
                if ( i==0 ) {
                  j = (unsigned int)pair_table[0]+1;
                  type = 0;  /* no pair */
                }
                else {
                  j = (unsigned int)pair_table[i];
                  type = P->model_details.pair[S[j]][S[i]]; if (type==0) type=7;
                  /* prime the ld5 variable */
                  if (SAME_STRAND(j-1,j)) {
                    ld5 = P->dangle5[type][S1[j-1]];
                    if ((p=(unsigned int)pair_table[j-2]) && SAME_STRAND(j-2, j-1))
                    if (P->dangle3[P->model_details.pair[S[p]][S[j-2]]][S1[j-1]]<ld5) ld5 = 0;
                  }
                }
                i1=i; p = i+1; u=0;
                energy = 0; cx_energy=INF;
                do { /* walk around the multi-loop */
                  new_cx = INF;

                  /* hop over unpaired positions */
                  while (p <= (unsigned int)pair_table[0] && pair_table[p]==0) p++;

                  /* memorize number of unpaired positions */
                  u += p-i1-1;
                  /* get position of pairing partner */
                  if ( p == (unsigned int)pair_table[0]+1 ){
                    q = 0;tt = 0; /* virtual root pair */
                  } else {
                    q  = (unsigned int)pair_table[p];
                    /* get type of base pair P->q */
                    tt = P->model_details.pair[S[p]][S[q]]; if (tt==0) tt=7;
                  }

                  energy += mlintern[tt];
                  cx_energy += mlintern[tt];

                  dang5=dang3=0;
                  if ((SAME_STRAND(p-1,p))&&(p>1))
                    dang5=P->dangle5[tt][S1[p-1]];      /* 5'dangle of pq pair */
                  if ((SAME_STRAND(i1,i1+1))&&(i1<(unsigned int)S[0]))
                    dang3 = P->dangle3[type][S1[i1+1]]; /* 3'dangle of previous pair */

                  switch (p-i1-1) {
                    case 0:   /* adjacent helices */
                              if (i1!=0){
                                if (SAME_STRAND(i1,p)) {
                                  new_cx = energy + P->stack[rtype[type]][rtype[tt]];
                                  /* subtract 5'dangle and TerminalAU penalty */
                                  new_cx += -ld5 - mlintern[tt]-mlintern[type]+2*mlintern[1];
                                }
                                ld5=0;
                                energy = MIN2(energy, cx_energy);
                              }
                              break;
                    case 1:   /* 1 unpaired base between helices */
                              dang = MIN2(dang3, dang5);
                              energy = energy +dang; ld5 = dang - dang3;
                              /* may be problem here: Suppose
                                cx_energy>energy, cx_energy+dang5<energy
                                and the following helices are also stacked (i.e.
                                we'll subtract the dang5 again */
                              if (cx_energy+dang5 < energy) {
                                energy = cx_energy+dang5;
                                ld5 = dang5;
                              }
                              new_cx = INF;  /* no coax stacking with mismatch for now */
                              break;
                    default:  /* many unpaired base between helices */
                              energy += dang5 +dang3;
                              energy = MIN2(energy, cx_energy + dang5);
                              new_cx = INF;  /* no coax stacking possible */
                              ld5 = dang5;
                              break;
                  }
                  type = tt;
                  cx_energy = new_cx;
                  i1 = q; p=q+1;
                } while (q!=i);
                best_energy = MIN2(energy, best_energy); /* don't use cx_energy here */
                /* fprintf(stderr, "%6.2d\t", energy); */
                /* skip a helix and start again */
                while (pair_table[p]==0) p++;
                if (i == (unsigned int)pair_table[p]) break;
                i = (unsigned int)pair_table[p];
              } /* end doing it twice */
              energy = best_energy;
              break;

    default:  E_mm5_available = E2_mm5_available  = INF;
              E_mm5_occupied  = E2_mm5_occupied   = 0;
              while(p < j){
                /* p must have a pairing partner */
                q  = (int)pair_table[p];
                /* get type of base pair (p,q) */
                tt = P->model_details.pair[S[p]][S[q]];
                if(tt==0) tt=7;
                if(q_prev + 2 < p){
                  E_mm5_available = MIN2(E_mm5_available, E_mm5_occupied);
                  E_mm5_occupied  = E_mm5_available;
                }
                if(q_prev2 + 2 < p){
                  E2_mm5_available  = MIN2(E2_mm5_available, E2_mm5_occupied);
                  E2_mm5_occupied   = E2_mm5_available;
                }
                mm5 = ((SAME_STRAND(p-1,p)) && !pair_table[p-1])  ? S1[p-1] : -1;
                mm3 = ((SAME_STRAND(q,q+1)) && !pair_table[q+1])  ? S1[q+1] : -1;
                tmp =                   MIN2(
                                              E_mm5_occupied  + E_MLstem(tt, -1, mm3, P),
                                              E_mm5_available + E_MLstem(tt, mm5, mm3, P)
                                            );
                tmp   =                 MIN2(tmp, E_mm5_available + E_MLstem(tt, -1, mm3, P));
                tmp2  =                 MIN2(
                                              E_mm5_occupied  + E_MLstem(tt, -1, -1, P),
                                              E_mm5_available + E_MLstem(tt, mm5, -1, P)
                                            );
                E_mm5_available =       MIN2(tmp2, E_mm5_available  + E_MLstem(tt, -1, -1, P));
                E_mm5_occupied  = tmp;

                tmp =                  MIN2(
                                              E2_mm5_occupied  + E_MLstem(tt, -1, mm3, P),
                                              E2_mm5_available + E_MLstem(tt, mm5, mm3, P)
                                            );
                tmp =                   MIN2(tmp, E2_mm5_available + E_MLstem(tt, -1, mm3, P));
                tmp2 =                  MIN2(
                                              E2_mm5_occupied  + E_MLstem(tt, -1, -1, P),
                                              E2_mm5_available + E_MLstem(tt, mm5, -1, P)
                                            );
                E2_mm5_available =      MIN2(tmp2, E2_mm5_available + E_MLstem(tt, -1, -1, P));
                E2_mm5_occupied = tmp;
                /* seek to the next stem */
                p = q + 1;
                q_prev = q_prev2 = q;
                while (p <= j && !pair_table[p]) p++;
                u += p - q - 1; /* add unpaired nucleotides */
              }
              /* now lets see how we get the minimum including the enclosing stem */
              type = P->model_details.pair[S[j]][S[i]]; if (type==0) type=7;
              mm5 = ((SAME_STRAND(j-1,j)) && !pair_table[j-1])  ? S1[j-1] : -1;
              mm3 = ((SAME_STRAND(i,i+1)) && !pair_table[i+1])  ? S1[i+1] : -1;
              if(q_prev + 2 < p){
                E_mm5_available = MIN2(E_mm5_available, E_mm5_occupied);
                E_mm5_occupied  = E_mm5_available;
              }
              if(q_prev2 + 2 < p){
                E2_mm5_available  = MIN2(E2_mm5_available, E2_mm5_occupied);
                E2_mm5_occupied   = E2_mm5_available;
              }
              energy = MIN2(E_mm5_occupied  + E_MLstem(type, -1, -1, P),
                            E_mm5_available + E_MLstem(type, mm5, -1, P)
                          );
              energy = MIN2(energy, E_mm5_available   + E_MLstem(type, -1, -1, P));
              energy = MIN2(energy, E2_mm5_occupied   + E_MLstem(type, -1, mm3, P));
              energy = MIN2(energy, E2_mm5_occupied   + E_MLstem(type, -1, -1, P));
              energy = MIN2(energy, E2_mm5_available  + E_MLstem(type, mm5, mm3, P));
              energy = MIN2(energy, E2_mm5_available  + E_MLstem(type, -1, mm3, P));
              energy = MIN2(energy, E2_mm5_available  + E_MLstem(type, mm5, -1, P));
              energy = MIN2(energy, E2_mm5_available  + E_MLstem(type, -1, -1, P));
              break;
  }/* end switch dangle_model */

  energy += P->MLclosing;
  /* logarithmic ML loop energy if logML */
  if(logML && (u>6))
    energy += 6*P->MLbase+(int)(P->lxc*log((double)u/6.));
  else
    energy += (u*P->MLbase);

  return energy;
}


PUBLIC int loop_energy(short * ptable, short *s, short *s1, int i) {
  /* compute energy of a single loop closed by base pair (i,j) */
  int j, type, p,q, energy;
  short *Sold, *S1old, *ptold;

  ptold=pair_table;   Sold = S;   S1old = S1;
  pair_table = ptable;   S = s;   S1 = s1;

  if (i==0) { /* evaluate exterior loop */
    energy = energy_of_extLoop_pt(0,pair_table);
    pair_table=ptold; S=Sold; S1=S1old;
    return energy;
  }
  j = pair_table[i];
  if (j<i) nrerror("i is unpaired in loop_energy()");
  type = P->model_details.pair[S[i]][S[j]];
  if (type==0) {
    type=7;
    if (eos_debug>=0)
      fprintf(stderr,"WARNING: bases %d and %d (%c%c) can't pair!\n", i, j,
              get_encoded_char(S[i], &(P->model_details)), get_encoded_char(S[j], &(P->model_details)));
  }
  p=i; q=j;


  while (pair_table[++p]==0);
  while (pair_table[--q]==0);
  if (p>q) { /* Hairpin */
    char loopseq[8] = "";
    if (SAME_STRAND(i,j)) {
      if (j-i-1<7) {
        int u;
        for (u=0; i+u<=j; u++) loopseq[u] = get_encoded_char(S[i+u], &(P->model_details));
        loopseq[u] = '\0';
      }
      energy = E_Hairpin(j-i-1, type, S1[i+1], S1[j-1], loopseq, P);
    } else {
      energy = energy_of_extLoop_pt(cut_in_loop(i), pair_table);
    }
  }
  else if (pair_table[q]!=(short)p) { /* multi-loop */
    int ii;
    ii = cut_in_loop(i);
    energy = (ii==0) ? energy_of_ml_pt(i, pair_table) : energy_of_extLoop_pt(ii, pair_table);
  }
  else { /* found interior loop */
    int type_2;
    type_2 = P->model_details.pair[S[q]][S[p]];
    if (type_2==0) {
      type_2=7;
      if (eos_debug>=0)
        fprintf(stderr,"WARNING: bases %d and %d (%c%c) can't pair!\n", p, q,
              get_encoded_char(S[p], &(P->model_details)), get_encoded_char(S[q], &(P->model_details)));
    }
    /* energy += LoopEnergy(i, j, p, q, type, type_2); */
    if ( SAME_STRAND(i,p) && SAME_STRAND(q,j) )
      energy = E_IntLoop(p-i-1, j-q-1, type, type_2,
                          S1[i+1], S1[j-1], S1[p-1], S1[q+1], P);
    else
      energy = energy_of_extLoop_pt(cut_in_loop(i), pair_table);
  }

  pair_table=ptold; S=Sold; S1=S1old;
  return energy;
}

/*---------------------------------------------------------------------------*/


PUBLIC float energy_of_move(const char *string, const char *structure, int m1, int m2) {
  int   energy;
  short *ss, *ss1;

  if(P == NULL) update_fold_params();

  if (fabs(P->temperature - temperature)>1e-6) update_fold_params();

  if (strlen(structure)!=strlen(string))
    nrerror("energy_of_struct: string and structure have unequal length");

  /* save the S and S1 pointers in case they were already in use */
  ss = S; ss1 = S1;
  S   = get_sequence_encoding(string, 0, &(P->model_details));
  S1  = get_sequence_encoding(string, 1, &(P->model_details));

  pair_table = vrna_pt_get(structure);

  energy = energy_of_move_pt(pair_table, S, S1, m1, m2);

  free(pair_table);
  free(S); free(S1);
  S=ss; S1=ss1;
  return  (float) energy/100.;
}

/*---------------------------------------------------------------------------*/

PUBLIC int energy_of_move_pt(short *pt, short *s, short *s1, int m1, int m2) {
  /*compute change in energy given by move (m1,m2)*/
  int en_post, en_pre, i,j,k,l, len;

  len = pt[0];
  k = (m1>0)?m1:-m1;
  l = (m2>0)?m2:-m2;
  /* first find the enclosing pair i<k<l<j */
  for (j=l+1; j<=len; j++) {
    if (pt[j]<=0) continue; /* unpaired */
    if (pt[j]<k) break;   /* found it */
    if (pt[j]>j) j=pt[j]; /* skip substructure */
    else {
      fprintf(stderr, "%d %d %d %d ", m1, m2, j, pt[j]);
      nrerror("illegal move or broken pair table in energy_of_move()");
    }
  }
  i = (j<=len) ? pt[j] : 0;
  en_pre = loop_energy(pt, s, s1, i);
  en_post = 0;
  if (m1<0) { /*it's a delete move */
    en_pre += loop_energy(pt, s, s1, k);
    pt[k]=0;
    pt[l]=0;
  } else { /* insert move */
    pt[k]=l;
    pt[l]=k;
    en_post += loop_energy(pt, s, s1, k);
  }
  en_post += loop_energy(pt, s, s1, i);
  /*  restore pair table */
  if (m1<0) {
    pt[k]=l;
    pt[l]=k;
  } else {
    pt[k]=0;
    pt[l]=0;
  }
  return (en_post - en_pre);
}

PRIVATE int cut_in_loop(int i) {
  /* walk around the loop;  return j pos of pair after cut if
     cut_point in loop else 0 */
  int  p, j;
  p = j = pair_table[i];
  do {
    i  = pair_table[p];  p = i+1;
    while ( pair_table[p]==0 ) p++;
  } while (p!=j && SAME_STRAND(i,p));
  return SAME_STRAND(i,p) ? 0 : j;
}

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

PRIVATE int ML_Energy(int i, int is_extloop) {
  /* i is the 5'-base of the closing pair (or 0 for exterior loop)
     loop is scored as ML if extloop==0 else as exterior loop

     since each helix can coaxially stack with at most one of its
     neighbors we need an auxiliarry variable  cx_energy
     which contains the best energy given that the last two pairs stack.
     energy  holds the best energy given the previous two pairs do not
     stack (i.e. the two current helices may stack)
     We don't allow the last helix to stack with the first, thus we have to
     walk around the Loop twice with two starting points and take the minimum
  */

  int energy, cx_energy, best_energy=INF;
  int i1, j, p, q, u, x, type, count;
  int mlintern[NBPAIRS+1], mlclosing, mlbase;
  int dangle_model = P->model_details.dangles;
  int *rtype = &(P->model_details.rtype[0]);

  if (is_extloop) {
    for (x = 0; x <= NBPAIRS; x++)
      mlintern[x] = P->MLintern[x]-P->MLintern[1]; /* 0 or TerminalAU */
    mlclosing = mlbase = 0;
  } else {
    for (x = 0; x <= NBPAIRS; x++) mlintern[x] = P->MLintern[x];
    mlclosing = P->MLclosing; mlbase = P->MLbase;
  }

  /*  as we do not only have dangling end but also mismatch contributions,
  **  we do this a bit different to previous implementations
  */
  if(is_extloop){
    int mm5_prev, mm3_prev;
    int tt_prev = 0;
    energy = 0;
    i1  = i;
    p   = i+1;
    mm5_prev = mm3_prev = -1;

    int E_mm5_available, E_mm5_occupied;
    /* find out if we may have 5' mismatch for the next stem */
    while (p <= (int)pair_table[0] && pair_table[p]==0) p++;
    /* get position of pairing partner */
    if(p < (int)pair_table[0]){
        E_mm5_occupied  = (p - i - 1 > 0) ? INF : 0;
        E_mm5_available = (p - i - 1 > 0) ? 0 : INF;
    }

    if(p < (int)pair_table[0])
      do{
        int tt;
        /* p must have a pairing partner */
        q  = (int)pair_table[p];
        /* get type of base pair (p,q) */
        tt = P->model_details.pair[S[p]][S[q]];
        if(tt==0) tt=7;

        int mm5 = ((SAME_STRAND(p-1,p)) && (p>1)) ? S1[p-1]: -1;
        int mm3 = ((SAME_STRAND(q,q+1)) && (q<(unsigned int)pair_table[0])) ? S1[q+1]: -1;

        switch(dangle_model){
          /* dangle_model == 0 */
          case 0: energy += E_ExtLoop(tt, -1, -1, P);
                  break;
          /* dangle_model == 1 */
          case 1: {
                    /* check for unpaired nucleotide 3' to the current stem */
                    int u3 = ((q < pair_table[0]) && (pair_table[q+1] == 0)) ? 1 : 0;
                    if(pair_table[p-1] != 0) mm5 = -1;

                    if(!u3){
                      mm3 = -1;
                      E_mm5_occupied  = MIN2(
                                              E_mm5_occupied  + E_ExtLoop(tt, -1, -1, P),
                                              E_mm5_available + E_ExtLoop(tt, mm5, -1, P)
                                            );
                      E_mm5_available = E_mm5_occupied;
                    }
                    else{
                      E_mm5_occupied  = MIN2(
                                              E_mm5_occupied  + E_ExtLoop(tt, -1, mm3, P),
                                              E_mm5_available + E_ExtLoop(tt, mm5, mm3, P)
                                            );
                      E_mm5_available = MIN2(
                                              E_mm5_occupied  + E_ExtLoop(tt, -1, -1, P),
                                              E_mm5_available + E_ExtLoop(tt, mm5, -1, P)
                                            );
                    }
                  }
                  break;

          /* the beloved case dangle_model == 2 */
          case 2: energy += E_ExtLoop(tt, mm5, mm3, P);
                  break;

          /* dangle_model == 3 a.k.a. helix stacking */
          case 3: break;

        } /* end switch dangle_model */

        /* seek to the next stem */
        p = q + 1;
        while (p <= (int)pair_table[0] && pair_table[p]==0) p++;
        if(p == (int)pair_table[0] + 1){
          if(dangle_model == 1)
            energy = (p > q + 1) ? E_mm5_occupied : E_mm5_available;
          q = 0;
          break;
        }
      } while(q != i);
  }
  /* not exterior loop */
  else{
    for (count=0; count<2; count++) { /* do it twice */
      int ld5 = 0; /* 5' dangle energy on prev pair (type) */
      if ( i==0 ) {
        j = (unsigned int)pair_table[0]+1;
        type = 0;  /* no pair */
      }
      else {
        j = (unsigned int)pair_table[i];
        type = P->model_details.pair[S[j]][S[i]]; if (type==0) type=7;

        if (dangle_model==3) { /* prime the ld5 variable */
          if (SAME_STRAND(j-1,j)) {
            ld5 = P->dangle5[type][S1[j-1]];
            if ((p=(unsigned int)pair_table[j-2]) && SAME_STRAND(j-2, j-1))
                if (P->dangle3[P->model_details.pair[S[p]][S[j-2]]][S1[j-1]]<ld5) ld5 = 0;
          }
        }
      }
      i1=i; p = i+1; u=0;
      energy = 0; cx_energy=INF;
      do { /* walk around the multi-loop */
        int tt, new_cx = INF;

        /* hop over unpaired positions */
        while (p <= (unsigned int)pair_table[0] && pair_table[p]==0) p++;

        /* memorize number of unpaired positions */
        u += p-i1-1;
        /* get position of pairing partner */
        if ( p == (unsigned int)pair_table[0]+1 ){
          q = 0;tt = 0; /* virtual root pair */
        } else {
        q  = (unsigned int)pair_table[p];
          /* get type of base pair P->q */
        tt = P->model_details.pair[S[p]][S[q]]; if (tt==0) tt=7;
        }

        energy += mlintern[tt];
        cx_energy += mlintern[tt];

        if (dangle_model) {
          int dang5=0, dang3=0, dang;
          if ((SAME_STRAND(p-1,p))&&(p>1))
            dang5=P->dangle5[tt][S1[p-1]];      /* 5'dangle of pq pair */
          if ((SAME_STRAND(i1,i1+1))&&(i1<(unsigned int)S[0]))
            dang3 = P->dangle3[type][S1[i1+1]]; /* 3'dangle of previous pair */

          switch (p-i1-1) {
          case 0: /* adjacent helices */
            if (dangle_model==2)
              energy += dang3+dang5;
            else if (dangle_model==3 && i1!=0) {
              if (SAME_STRAND(i1,p)) {
                new_cx = energy + P->stack[rtype[type]][rtype[tt]];
                /* subtract 5'dangle and TerminalAU penalty */
                new_cx += -ld5 - mlintern[tt]-mlintern[type]+2*mlintern[1];
              }
              ld5=0;
              energy = MIN2(energy, cx_energy);
            }
            break;
          case 1: /* 1 unpaired base between helices */
            dang = (dangle_model==2)?(dang3+dang5):MIN2(dang3, dang5);
            if (dangle_model==3) {
              energy = energy +dang; ld5 = dang - dang3;
              /* may be problem here: Suppose
                 cx_energy>energy, cx_energy+dang5<energy
                 and the following helices are also stacked (i.e.
                 we'll subtract the dang5 again */
              if (cx_energy+dang5 < energy) {
                energy = cx_energy+dang5;
                ld5 = dang5;
              }
              new_cx = INF;  /* no coax stacking with mismatch for now */
            } else
              energy += dang;
            break;
          default: /* many unpaired base between helices */
            energy += dang5 +dang3;
            if (dangle_model==3) {
              energy = MIN2(energy, cx_energy + dang5);
              new_cx = INF;  /* no coax stacking possible */
              ld5 = dang5;
            }
          }
          type = tt;
        }
        if (dangle_model==3) cx_energy = new_cx;
        i1 = q; p=q+1;
      } while (q!=i);
      best_energy = MIN2(energy, best_energy); /* don't use cx_energy here */
      /* fprintf(stderr, "%6.2d\t", energy); */
      if (dangle_model!=3 || is_extloop) break;  /* may break cofold with co-ax */
      /* skip a helix and start again */
      while (pair_table[p]==0) p++;
      if (i == (unsigned int)pair_table[p]) break;
      i = (unsigned int)pair_table[p];
    }
    energy = best_energy;
    energy += mlclosing;
    /* logarithmic ML loop energy if logML */
    if ( (!is_extloop) && logML && (u>6) )
      energy += 6*mlbase+(int)(P->lxc*log((double)u/6.));
    else
      energy += mlbase*u;
    /* fprintf(stderr, "\n"); */
  }
  return energy;
}

PUBLIC float energy_of_struct(const char *string, const char *structure){
  return energy_of_structure(string, structure, eos_debug);
}

PUBLIC int energy_of_struct_pt(const char *string, short * ptable, short *s, short *s1){
  return energy_of_structure_pt(string, ptable, s, s1, eos_debug);
}

PUBLIC float energy_of_circ_struct(const char *string, const char *structure){
  return energy_of_circ_structure(string, structure, eos_debug);
}

