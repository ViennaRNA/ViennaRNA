/** \file eval.c */


/*
                  Free energy evaluation

                  c Ivo Hofacker, Chrisoph Flamm
                  original implementation by
                  Walter Fontana

                  ViennaRNA Package >= v2.0 by Ronny Lorenz

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
#include "ViennaRNA/model.h"
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
PUBLIC  int cut_point = -1; /* set to first pos of second seq for cofolding */
PUBLIC  int eos_debug = 0;  /* verbose info from energy_of_struct */

/*
#################################
# PRIVATE VARIABLES             #
#################################
*/

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/
PRIVATE int   stack_energy( int i,
                            const char *string,
                            const short *pt,
                            const short *s,
                            const short *s1,
                            FILE *file,
                            paramT *params,
                            int verbostiy_level);

PRIVATE int   energy_of_extLoop_pt( int i,
                                    const short *pt,
                                    const short *s,
                                    const short *s1,
                                    paramT *params);

PRIVATE int   energy_of_ml_pt(int i,
                              const short *pt,
                              const short *s,
                              const short *s1,
                              paramT *params);

PRIVATE int   cut_in_loop(int i, const short *pt);

PRIVATE int   eval_pt(const char *string,
                      const short *pt,
                      const short *s,
                      const short *s1,
                      FILE *file,
                      paramT *params,
                      int verbosity_level);

PRIVATE int   eval_circ_pt( const char *string,
                            const short *pt,
                            const short *s,
                            const short *s1,
                            paramT *P,
                            FILE *file,
                            int verbosity_level);

PRIVATE int   en_corr_of_loop_gquad(int i,
                                    int j,
                                    const char *string,
                                    const char *structure,
                                    const short *pt,
                                    const int *loop_idx,
                                    const short *s1,
                                    paramT *P);

PRIVATE paramT  *get_updated_params(paramT *parameters, int compat);

PRIVATE float   wrap_eval_structure(const char *string,
                                    const char *structure,
                                    const short *pt,
                                    const short *s,
                                    const short *s1,
                                    paramT *P,
                                    FILE *file,
                                    int verbosity);

PRIVATE int wrap_eval_loop_pt(const short *pt,
                              const short *s,
                              const short *s1,
                              int i,
                              paramT *params,
                              int verbosity);
/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PUBLIC  float
vrna_eval_structure(const char *string,
                    const char *structure,
                    paramT *parameters){

  paramT *P = get_updated_params(parameters, 0);
  short *s  = get_sequence_encoding(string, 0, &(P->model_details));
  short *s1 = get_sequence_encoding(string, 1, &(P->model_details));
  short *pt = vrna_pt_get(structure);
  float en  = wrap_eval_structure(string, structure, pt, s, s1, P, NULL, -1);

  free(P);
  free(pt);
  free(s);
  free(s1);
  return en;
}

PUBLIC float
vrna_eval_structure_verbose(const char *string,
                            const char *structure,
                            paramT *parameters,
                            FILE *file){

  paramT *P = get_updated_params(parameters, 0);
  short *s  = get_sequence_encoding(string, 0, &(P->model_details));
  short *s1 = get_sequence_encoding(string, 1, &(P->model_details));
  short *pt = vrna_pt_get(structure);
  float en  = wrap_eval_structure(string, structure, pt, s, s1, P, file, 1);

  free(P);
  free(pt);
  free(s);
  free(s1);
  return en;
}

PUBLIC int
vrna_eval_structure_pt( const char *string,
                        const short *pt,
                        paramT *parameters){

  if(pt && string){
    if(pt[0] != (short)strlen(string))
      nrerror("energy_of_struct: string and structure have unequal length");

    paramT *P = get_updated_params(parameters, 0); /* this ensures a proper P data structure */
    short *s  = get_sequence_encoding(string, 0, &(P->model_details));
    short *s1 = get_sequence_encoding(string, 1, &(P->model_details));

    int en = eval_pt(string, pt, s, s1, NULL, P, -1);
    free(s);
    free(s1);
    free(P);
    return en;
  } else
    return INF;
}

PUBLIC int
vrna_eval_structure_pt_verbose( const char *string,
                                const short *pt,
                                paramT *parameters,
                                FILE *file){

  if(pt && string){
    if(pt[0] != (short)strlen(string))
      nrerror("energy_of_struct: string and structure have unequal length");

    paramT *P = get_updated_params(parameters, 0); /* this ensures a proper P data structure */
    short *s  = get_sequence_encoding(string, 0, &(P->model_details));
    short *s1 = get_sequence_encoding(string, 1, &(P->model_details));

    int en = eval_pt(string, pt, s, s1, file, P, 1);
    free(s);
    free(s1);
    free(P);
    return en;
  } else
    return INF;
}

PUBLIC int
vrna_eval_structure_pt_fast(const char *string,
                            const short *pt,
                            const short *s,
                            const short *s1,
                            paramT *parameters){

  if(pt && string){
    if(pt[0] != (short)strlen(string))
      nrerror("energy_of_struct: string and structure have unequal length");

    paramT *P = get_updated_params(parameters, 0); /* this ensures a proper P data structure */

    int en = eval_pt(string, pt, s, s1, NULL, P, -1);
    free(P);
    return en;
  } else
    return INF;
}

PUBLIC int
vrna_eval_loop_pt(const short *pt,
                  const short *s,
                  const short *s1,
                  int i,
                  paramT *params){

  return wrap_eval_loop_pt(pt, s, s1, i, params, -1);
}

PRIVATE int
wrap_eval_loop_pt(const short *pt,
                  const short *s,
                  const short *s1,
                  int i,
                  paramT *params,
                  int verbosity){

  /* compute energy of a single loop closed by base pair (i,j) */
  int j, type, p,q, energy;

  paramT *P = get_updated_params(params, 0);

  if (i==0) { /* evaluate exterior loop */
    energy = energy_of_extLoop_pt(0,pt, s, s1, P);
    free(P);
    return energy;
  }
  j = pt[i];
  if (j<i) nrerror("i is unpaired in loop_energy()");
  type = P->model_details.pair[s[i]][s[j]];
  if (type==0) {
    type=7;
    if (verbosity>=0)
      fprintf(stderr,"WARNING: bases %d and %d (%c%c) can't pair!\n", i, j,
              get_encoded_char(s[i], &(P->model_details)), get_encoded_char(s[j], &(P->model_details)));
  }
  p=i; q=j;


  while (pt[++p]==0);
  while (pt[--q]==0);
  if (p>q) { /* Hairpin */
    char loopseq[8] = "";
    if (SAME_STRAND(i,j)) {
      if (j-i-1<7) {
        int u;
        for (u=0; i+u<=j; u++) loopseq[u] = get_encoded_char(s[i+u], &(P->model_details));
        loopseq[u] = '\0';
      }
      energy = E_Hairpin(j-i-1, type, s1[i+1], s1[j-1], loopseq, P);
    } else {
      energy = energy_of_extLoop_pt(cut_in_loop(i, (const short *)pt), (const short *)pt, s, s1, P);
    }
  }
  else if (pt[q]!=(short)p) { /* multi-loop */
    int ii;
    ii = cut_in_loop(i, (const short *)pt);
    energy = (ii==0) ? energy_of_ml_pt(i, (const short *)pt, s, s1, P) : energy_of_extLoop_pt(ii, (const short *)pt, s, s1, P);
  }
  else { /* found interior loop */
    int type_2;
    type_2 = P->model_details.pair[s[q]][s[p]];
    if (type_2==0) {
      type_2=7;
      if (verbosity>=0)
        fprintf(stderr,"WARNING: bases %d and %d (%c%c) can't pair!\n", p, q,
              get_encoded_char(s[p], &(P->model_details)), get_encoded_char(s[q], &(P->model_details)));
    }
    /* energy += LoopEnergy(i, j, p, q, type, type_2); */
    if ( SAME_STRAND(i,p) && SAME_STRAND(q,j) )
      energy = E_IntLoop(p-i-1, j-q-1, type, type_2,
                          s1[i+1], s1[j-1], s1[p-1], s1[q+1], P);
    else
      energy = energy_of_extLoop_pt(cut_in_loop(i, (const short *)pt), (const short *)pt, s, s1, P);
  }

  free(P);
  return energy;
}

PUBLIC float
vrna_eval_move( const char *string,
                const char *structure,
                int m1,
                int m2,
                paramT *params){

  if (strlen(structure)!=strlen(string))
    nrerror("vrna_eval_move: sequence and structure have unequal length");

  paramT *P = get_updated_params(params, 0);
  short *pt = vrna_pt_get(structure);
  short *s  = get_sequence_encoding(string, 0, &(P->model_details));
  short *s1 = get_sequence_encoding(string, 1, &(P->model_details));

  int en = vrna_eval_move_pt(pt, s, s1, m1, m2, P);

  free(s1);
  free(s);
  free(pt);
  free(P);
  return  (float)en/100.;
}

PUBLIC int
vrna_eval_move_pt(short *pt,
                  const short *s,
                  const short *s1,
                  int m1,
                  int m2,
                  paramT *params){

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
      nrerror("illegal move or broken pair table in vrna_eval_move_pt()");
    }
  }
  i = (j<=len) ? pt[j] : 0;
  en_pre = vrna_eval_loop_pt((const short *)pt, (const short *)s, (const short *)s1, i, params);
  en_post = 0;
  if (m1<0) { /*it's a delete move */
    en_pre += vrna_eval_loop_pt((const short *)pt, (const short *)s, (const short *)s1, k, params);
    pt[k]=0;
    pt[l]=0;
  } else { /* insert move */
    pt[k]=l;
    pt[l]=k;
    en_post += vrna_eval_loop_pt((const short *)pt, (const short *)s, (const short *)s1, k, params);
  }
  en_post += vrna_eval_loop_pt((const short *)pt, (const short *)s, (const short *)s1, i, params);
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

/*
#################################
# STATIC helper functions below #
#################################
*/

PRIVATE  paramT *
get_updated_params(paramT *parameters, int compat){
  paramT *P = NULL;
  if(parameters){
    P = get_parameter_copy(parameters);
  } else {
    model_detailsT md;
    if(compat)
      set_model_details(&md);
    else
      vrna_md_set_default(&md);
    P = get_scaled_parameters(temperature, md);
  }
  fill_pair_matrices(&(P->model_details));
  return P;
}

PRIVATE float
wrap_eval_structure(const char *string,
                    const char *structure,
                    const short *pt,
                    const short *s,
                    const short *s1,
                    paramT *P,
                    FILE *file,
                    int verbosity){

  int res;
  int gq;

  gq = P->model_details.gquad;
  P->model_details.gquad = 0;
  if(P->model_details.circ){
    res = eval_circ_pt(string, pt, s, s1, P, file, verbosity);
  } else {
    res = eval_pt(string, (const short *)pt, s, s1, file, P, verbosity);
  }
  P->model_details.gquad = gq;

  if(gq){
    int *loop_idx = vrna_get_loop_index(pt);
    int gge       = en_corr_of_loop_gquad(1, s[0], string, structure, pt, loop_idx, s1, P);
    res          += gge;
    free(loop_idx);
  }

  return (float)res/100.;
}

PRIVATE int
eval_pt(const char *string,
        const short *pt,
        const short *s,
        const short *s1,
        FILE *file,
        paramT *params,
        int verbosity_level){

  int   i, length, energy;
  short *ss, *ss1;

  FILE *out = (file) ? file : stdout;

  length = s[0];

  if(params->model_details.gquad)
    warn_user("vrna_eval_*_pt: No gquadruplex support!\nIgnoring potential gquads in structure!\nUse e.g. vrna_eval_structure() instead!");

  energy =  backtrack_type=='M' ? energy_of_ml_pt(0, pt, s, s1, params) : energy_of_extLoop_pt(0, pt, s, s1, params);
  if (verbosity_level>0)
    fprintf(out, "External loop                           : %5d\n", energy);
  for (i=1; i<=length; i++) {
    if (pt[i]==0) continue;
    energy += stack_energy(i, string, pt, s, s1, out, params, verbosity_level);
    i=pt[i];
  }
  for (i=1; !SAME_STRAND(i,length); i++) {
    if (!SAME_STRAND(i,pt[i])) {
      energy+=params->DuplexInit;
      break;
    }
  }

  return energy;
}

PRIVATE int
eval_circ_pt( const char *string,
              const short *pt,
              const short *s,
              const short *s1,
              paramT *P,
              FILE *file,
              int verbosity_level){

  int   i, j, length, energy=0, en0, degree=0, type;

  int dangle_model = P->model_details.dangles;

  FILE *out = (file) ? file : stdout;

  length = s[0];

  if(P->model_details.gquad)
    warn_user("vrna_eval_*_pt: No gquadruplex support!\nIgnoring potential gquads in structure!\nUse e.g. vrna_eval_structure() instead!");

  for (i=1; i<=length; i++) {
    if (pt[i]==0) continue;
    degree++;
    energy += stack_energy(i, string, (const short *)pt, s, s1, out, P, verbosity_level);
    i=pt[i];
  }

  if (degree==0){ return 0.;}
  for (i=1; pt[i]==0; i++);
  j = pt[i];
  type = P->model_details.pair[s[j]][s[i]];
  if (type==0) type=7;
  if (degree==1) {
    char loopseq[10];
    int u, si1, sj1;
    for (i=1; pt[i]==0; i++);
    u = length-j + i-1;
    if (u<7) {
      strcpy(loopseq , string+j-1);
      strncat(loopseq, string, i);
    }
    si1 = (i==1) ? s1[length] : s1[i-1];
    sj1 = (j==length) ? s1[1] : s1[j+1];
    en0 = E_Hairpin(u, type, sj1, si1, loopseq, P);
  } else
    if (degree==2) {
      int p,q, u1,u2, si1, sq1, type_2;
      for (p=j+1; pt[p]==0; p++);
      q=pt[p];
      u1 = p-j-1;
      u2 = i-1 + length-q;
      type_2 = P->model_details.pair[s[q]][s[p]];
      if (type_2==0) type_2=7;
      si1 = (i==1)? s1[length] : s1[i-1];
      sq1 = (q==length)? s1[1] : s1[q+1];
      en0 = E_IntLoop(u1, u2, type, type_2,
                       s1[j+1], si1, s1[p-1], sq1,P);
    } else { /* degree > 2 */
      en0 = energy_of_ml_pt(0, (const short *)pt, s, s1, P) - P->MLintern[0];
#if 0
      if (dangle_model) {
        int d5, d3;
        if (pt[1]) {
          j = pt[1];
          type = P->model_details.pair[s[1]][s[j]];
          if (dangle_model==2)
            en0 += P->dangle5[type][s1[length]];
          else { /* dangle_model==1 */
            if (pt[length]==0) {
              d5 = P->dangle5[type][s1[length]];
              if (pt[length-1]!=0) {
                int tt;
                tt = P->model_details.pair[s[pt[length-1]]][s[length-1]];
                d3 = P->dangle3[tt][s1[length]];
                if (d3<d5) d5 = 0;
                else d5 -= d3;
              }
              en0 += d5;
            }
          }
        }
        if (pt[length]) {
          i = pt[length];
          type = P->model_details.pair[s[i]][s[length]];
          if (dangle_model==2)
            en0 += P->dangle3[type][s1[1]];
          else { /* dangle_model==1 */
            if (pt[1]==0) {
              d3 = P->dangle3[type][s1[1]];
              if (pt[2]) {
                int tt;
                tt = P->model_details.pair[s[2]][s[pt[2]]];
                d5 = P->dangle5[tt][1];
                if (d5<d3) d3=0;
                else d3 -= d5;
              }
              en0 += d3;
            }
          }
        }
      }
#endif
    }

  if (verbosity_level>0)
    fprintf(out, "External loop                           : %5d\n", en0);
  energy += en0;
  /* fprintf(stderr, "ext loop degree %d tot %d\n", degree, energy); */

  return energy;
}



/*---------------------------------------------------------------------------*/
/*  returns a correction term that may be added to the energy retrieved
    from energy_of_struct_par() to correct misinterpreted loops. This
    correction is necessary since energy_of_struct_par() will forget 
    about the existance of gquadruplexes and just treat them as unpaired
    regions.

    recursive variant
*/
PRIVATE int
en_corr_of_loop_gquad(int i,
                      int j,
                      const char *string,
                      const char *structure,
                      const short *pt,
                      const int *loop_idx,
                      const short *s1,
                      paramT *P){

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
          energy += en_corr_of_loop_gquad(u, pt[u], string, structure, pt, loop_idx, s1, P);
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



PRIVATE int
stack_energy( int i,
              const char *string,
              const short *pt,
              const short *s,
              const short *s1,
              FILE *file,
              paramT *P,
              int verbosity_level){

  FILE *out = (file) ? file : stdout;

  /* calculate energy of substructure enclosed by (i,j) */
  int ee, energy = 0;
  int j, p, q, type;
  int *rtype = &(P->model_details.rtype[0]);

  j=pt[i];
  type = P->model_details.pair[s[i]][s[j]];
  if (type==0) {
    type=7;
    if (verbosity_level>=0)
      fprintf(stderr,"WARNING: bases %d and %d (%c%c) can't pair!\n", i, j,
              string[i-1],string[j-1]);
  }

  p=i; q=j;
  while (p<q) { /* process all stacks and interior loops */
    int type_2;
    while (pt[++p]==0);
    while (pt[--q]==0);
    if ((pt[q]!=(short)p)||(p>q)) break;
    type_2 = P->model_details.pair[s[q]][s[p]];
    if (type_2==0) {
      type_2=7;
      if (verbosity_level>=0)
        fprintf(stderr,"WARNING: bases %d and %d (%c%c) can't pair!\n", p, q,
                string[p-1],string[q-1]);
    }
    /* energy += LoopEnergy(i, j, p, q, type, type_2); */
    if ( SAME_STRAND(i,p) && SAME_STRAND(q,j) )
      ee = E_IntLoop(p-i-1, j-q-1, type, type_2, s1[i+1], s1[j-1], s1[p-1], s1[q+1],P);
    else
      ee = energy_of_extLoop_pt(cut_in_loop(i, pt), pt, s, s1, P);
    if (verbosity_level>0)
      fprintf(out, "Interior loop (%3d,%3d) %c%c; (%3d,%3d) %c%c: %5d\n",
             i,j,string[i-1],string[j-1],p,q,string[p-1],string[q-1], ee);
    energy += ee;
    i=p; j=q; type = rtype[type_2];
  } /* end while */

  /* p,q don't pair must have found hairpin or multiloop */

  if (p>q) {                       /* hair pin */
    if (SAME_STRAND(i,j))
      ee = E_Hairpin(j-i-1, type, s1[i+1], s1[j-1], string+i-1, P);
    else
      ee = energy_of_extLoop_pt(cut_in_loop(i, pt), pt, s, s1, P);
    energy += ee;
    if (verbosity_level>0)
      fprintf(out, "Hairpin  loop (%3d,%3d) %c%c              : %5d\n",
             i, j, string[i-1],string[j-1], ee);

    return energy;
  }

  /* (i,j) is exterior pair of multiloop */
  while (p<j) {
    /* add up the contributions of the substructures of the ML */
    energy += stack_energy(p, string, pt, s, s1, out, P, verbosity_level);
    p = pt[p];
    /* search for next base pair in multiloop */
    while (pt[++p]==0);
  }
  {
    int ii;
    ii = cut_in_loop(i, pt);
    ee = (ii==0) ? energy_of_ml_pt(i, pt, s, s1, P) : energy_of_extLoop_pt(ii, pt, s, s1, P);
  }
  energy += ee;
  if (verbosity_level>0)
    fprintf(out, "Multi    loop (%3d,%3d) %c%c              : %5d\n",
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
PRIVATE int
energy_of_extLoop_pt( int i,
                      const short *pt,
                      const short *s,
                      const short *s1,
                      paramT *P){

  int energy, mm5, mm3;
  int p, q, q_prev;
  int length = (int)pt[0];

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
  while(p <= length && !pt[p]) p++;

  while(p < length){
    int tt;
    /* p must have a pairing partner */
    q  = (int)pt[p];
    /* get type of base pair (p,q) */
    tt = P->model_details.pair[s[p]][s[q]];
    if(tt==0) tt=7;

    switch(dangle_model){
      /* no dangles */
      case 0:   energy += E_ExtLoop(tt, -1, -1, P);
                break;

      /* the beloved double dangles */
      case 2:   mm5 = ((SAME_STRAND(p-1,p)) && (p>1))       ? s1[p-1] : -1;
                mm3 = ((SAME_STRAND(q,q+1)) && (q<length))  ? s1[q+1] : -1;
                energy += E_ExtLoop(tt, mm5, mm3, P);
                break;

      default:  {
                  int tmp;
                  if(q_prev + 2 < p){
                    E3_available = MIN2(E3_available, E3_occupied);
                    E3_occupied  = E3_available;
                  }
                  mm5 = ((SAME_STRAND(p-1,p)) && (p>1) && !pt[p-1])       ? s1[p-1] : -1;
                  mm3 = ((SAME_STRAND(q,q+1)) && (q<length) && !pt[q+1])  ? s1[q+1] : -1;
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
    while (p <= length && !pt[p]) p++;
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
PRIVATE int
energy_of_ml_pt(int i,
                const short *pt,
                const short *s,
                const short *s1,
                paramT *P){

  int energy, cx_energy, tmp, tmp2, best_energy=INF;
  int i1, j, p, q, q_prev, q_prev2, u, x, type, count, mm5, mm3, tt, ld5, new_cx, dang5, dang3, dang;
  int e_stem, e_stem5, e_stem3, e_stem53;
  int mlintern[NBPAIRS+1];

  /* helper variables for dangles == 1|5 case */
  int E_mm5_available;  /* energy of 5' part where 5' mismatch of current stem is available */
  int E_mm5_occupied;   /* energy of 5' part where 5' mismatch of current stem is unavailable */
  int E2_mm5_available; /* energy of 5' part where 5' mismatch of current stem is available with possible 3' dangle for enclosing pair (i,j) */
  int E2_mm5_occupied;  /* energy of 5' part where 5' mismatch of current stem is unavailable with possible 3' dangle for enclosing pair (i,j) */
  int dangle_model = P->model_details.dangles;
  int *rtype = &(P->model_details.rtype[0]);

  if(i >= pt[i])
    nrerror("energy_of_ml_pt: i is not 5' base of a closing pair!");

  if(i == 0)
    j = (int)pt[0] + 1;
  else
    j = (int)pt[i];

  /* init the variables */
  energy      = 0;
  p           = i+1;
  q_prev      = i-1;
  q_prev2     = i;

  for (x = 0; x <= NBPAIRS; x++) mlintern[x] = P->MLintern[x];

  /* seek to opening base of first stem */
  while(p <= j && !pt[p]) p++;
  u = p - i - 1;

  switch(dangle_model){
    case 0:   while(p < j){
                /* p must have a pairing partner */
                q  = (int)pt[p];
                /* get type of base pair (p,q) */
                tt = P->model_details.pair[s[p]][s[q]];
                if(tt==0) tt=7;
                energy += E_MLstem(tt, -1, -1, P);
                /* seek to the next stem */
                p = q + 1;
                q_prev = q_prev2 = q;
                while (p <= j && !pt[p]) p++;
                u += p - q - 1; /* add unpaired nucleotides */
              }
              /* now lets get the energy of the enclosing stem */
              if(i > 0){  /* actual closing pair */
                type = P->model_details.pair[s[j]][s[i]]; if (type==0) type=7;
                energy += E_MLstem(type, -1, -1, P);
              } else {  /* virtual closing pair */
                energy += E_MLstem(0, -1, -1, P);
              }
              break;

    case 2:   while(p < j){
                /* p must have a pairing partner */
                q  = (int)pt[p];
                /* get type of base pair (p,q) */
                tt = P->model_details.pair[s[p]][s[q]];
                if(tt==0) tt=7;
                mm5 = (SAME_STRAND(p-1,p))  ? s1[p-1] : -1;
                mm3 = (SAME_STRAND(q,q+1))  ? s1[q+1] : -1;
                energy += E_MLstem(tt, mm5, mm3, P);
                /* seek to the next stem */
                p = q + 1;
                q_prev = q_prev2 = q;
                while (p <= j && !pt[p]) p++;
                u += p - q - 1; /* add unpaired nucleotides */
              }
              if(i > 0){  /* actual closing pair */
                type = P->model_details.pair[s[j]][s[i]]; if (type==0) type=7;
                mm5 = ((SAME_STRAND(j-1,j)) && !pt[j-1])  ? s1[j-1] : -1;
                mm3 = ((SAME_STRAND(i,i+1)) && !pt[i+1])  ? s1[i+1] : -1;
                energy += E_MLstem(type, s1[j-1], s1[i+1], P);
              } else {  /* virtual closing pair */
                energy += E_MLstem(0, -1, -1, P);
              }
              break;

    case 3:   /* we treat helix stacking different */
              for (count=0; count<2; count++) { /* do it twice */
                ld5 = 0; /* 5' dangle energy on prev pair (type) */
                if ( i==0 ) {
                  j = (unsigned int)pt[0]+1;
                  type = 0;  /* no pair */
                }
                else {
                  j = (unsigned int)pt[i];
                  type = P->model_details.pair[s[j]][s[i]]; if (type==0) type=7;
                  /* prime the ld5 variable */
                  if (SAME_STRAND(j-1,j)) {
                    ld5 = P->dangle5[type][s1[j-1]];
                    if ((p=(unsigned int)pt[j-2]) && SAME_STRAND(j-2, j-1))
                    if (P->dangle3[P->model_details.pair[s[p]][s[j-2]]][s1[j-1]]<ld5) ld5 = 0;
                  }
                }
                i1=i; p = i+1; u=0;
                energy = 0; cx_energy=INF;
                do { /* walk around the multi-loop */
                  new_cx = INF;

                  /* hop over unpaired positions */
                  while (p <= (unsigned int)pt[0] && pt[p]==0) p++;

                  /* memorize number of unpaired positions */
                  u += p-i1-1;
                  /* get position of pairing partner */
                  if ( p == (unsigned int)pt[0]+1 ){
                    q = 0;tt = 0; /* virtual root pair */
                  } else {
                    q  = (unsigned int)pt[p];
                    /* get type of base pair P->q */
                    tt = P->model_details.pair[s[p]][s[q]]; if (tt==0) tt=7;
                  }

                  energy += mlintern[tt];
                  cx_energy += mlintern[tt];

                  dang5=dang3=0;
                  if ((SAME_STRAND(p-1,p))&&(p>1))
                    dang5=P->dangle5[tt][s1[p-1]];      /* 5'dangle of pq pair */
                  if ((SAME_STRAND(i1,i1+1))&&(i1<(unsigned int)s[0]))
                    dang3 = P->dangle3[type][s1[i1+1]]; /* 3'dangle of previous pair */

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
                while (pt[p]==0) p++;
                if (i == (unsigned int)pt[p]) break;
                i = (unsigned int)pt[p];
              } /* end doing it twice */
              energy = best_energy;
              break;

    default:  E_mm5_available = E2_mm5_available  = INF;
              E_mm5_occupied  = E2_mm5_occupied   = 0;
              while(p < j){
                /* p must have a pairing partner */
                q  = (int)pt[p];
                /* get type of base pair (p,q) */
                tt = P->model_details.pair[s[p]][s[q]];
                if(tt==0) tt=7;
                if(q_prev + 2 < p){
                  E_mm5_available = MIN2(E_mm5_available, E_mm5_occupied);
                  E_mm5_occupied  = E_mm5_available;
                }
                if(q_prev2 + 2 < p){
                  E2_mm5_available  = MIN2(E2_mm5_available, E2_mm5_occupied);
                  E2_mm5_occupied   = E2_mm5_available;
                }
                mm5 = ((SAME_STRAND(p-1,p)) && !pt[p-1])  ? s1[p-1] : -1;
                mm3 = ((SAME_STRAND(q,q+1)) && !pt[q+1])  ? s1[q+1] : -1;
                e_stem    = E_MLstem(tt, -1, -1, P);
                e_stem5   = E_MLstem(tt, mm5, -1, P);
                e_stem3   = E_MLstem(tt, -1, mm3, P);
                e_stem53  = E_MLstem(tt, mm5, mm3, P);

                tmp   = E_mm5_occupied + e_stem3;
                tmp   = MIN2(tmp, E_mm5_available + e_stem53);
                tmp   = MIN2(tmp, E_mm5_available + e_stem3);
                tmp2  = E_mm5_occupied + e_stem;
                tmp2  = MIN2(tmp2, E_mm5_available + e_stem5);
                tmp2  = MIN2(tmp2, E_mm5_available + e_stem); 

                E_mm5_occupied  = tmp;
                E_mm5_available = tmp2;

                tmp   = E2_mm5_occupied + e_stem3;
                tmp   = MIN2(tmp, E2_mm5_available + e_stem53);
                tmp   = MIN2(tmp, E2_mm5_available + e_stem3);
                tmp2  = E2_mm5_occupied + e_stem;
                tmp2  = MIN2(tmp2, E2_mm5_available + e_stem5);
                tmp2  = MIN2(tmp2, E2_mm5_available + e_stem);
                
                E2_mm5_occupied   = tmp;
                E2_mm5_available  = tmp2;

                /* seek to the next stem */
                p = q + 1;
                q_prev = q_prev2 = q;
                while (p <= j && !pt[p]) p++;
                u += p - q - 1; /* add unpaired nucleotides */
              }
              if(i > 0){  /* actual closing pair */
                type = P->model_details.pair[s[j]][s[i]]; if (type==0) type=7;
                mm5 = ((SAME_STRAND(j-1,j)) && !pt[j-1])  ? s1[j-1] : -1;
                mm3 = ((SAME_STRAND(i,i+1)) && !pt[i+1])  ? s1[i+1] : -1;
                if(q_prev + 2 < p){
                  E_mm5_available = MIN2(E_mm5_available, E_mm5_occupied);
                  E_mm5_occupied  = E_mm5_available;
                }
                if(q_prev2 + 2 < p){
                  E2_mm5_available  = MIN2(E2_mm5_available, E2_mm5_occupied);
                  E2_mm5_occupied   = E2_mm5_available;
                }
                e_stem    = E_MLstem(type, -1, -1, P);
                e_stem5   = E_MLstem(type, mm5, -1, P);
                e_stem3   = E_MLstem(type, -1, mm3, P);
                e_stem53  = E_MLstem(type, mm5, mm3, P);
              } else {  /* virtual closing pair */
                e_stem = e_stem5 = e_stem3 = e_stem53 = E_MLstem(0, -1, -1, P);
              }
              /* now lets see how we get the minimum including the enclosing stem */
              energy = E_mm5_occupied  + e_stem;
              energy = MIN2(energy, E_mm5_available   + e_stem5);
              energy = MIN2(energy, E_mm5_available   + e_stem);
              energy = MIN2(energy, E2_mm5_occupied   + e_stem3);
              energy = MIN2(energy, E2_mm5_occupied   + e_stem);
              energy = MIN2(energy, E2_mm5_available  + e_stem53);
              energy = MIN2(energy, E2_mm5_available  + e_stem3);
              energy = MIN2(energy, E2_mm5_available  + e_stem5);
              energy = MIN2(energy, E2_mm5_available  + e_stem);
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



PRIVATE int
cut_in_loop(int i, const short *pt){

  /* walk around the loop;  return j pos of pair after cut if
     cut_point in loop else 0 */
  int  p, j;
  p = j = pt[i];
  do {
    i  = pt[p];  p = i+1;
    while ( pt[p]==0 ) p++;
  } while (p!=j && SAME_STRAND(i,p));
  return SAME_STRAND(i,p) ? 0 : j;
}

/*
#################################
# DEPRECATED functions below    #
#################################
*/

PUBLIC float
energy_of_struct( const char *string,
                  const char *structure){

  float en;
  paramT *P = get_updated_params(NULL, 1);

  if(eos_debug > 0)
    en = vrna_eval_structure_verbose(string, structure, P, NULL);
  else
    en = vrna_eval_structure(string, structure, P);

  free(P);
  return en;
}

PUBLIC int
energy_of_struct_pt(const char *string,
                    short *pt,
                    short *s,
                    short *s1){

  if(pt && string){
    if(pt[0] != (short)strlen(string))
      nrerror("energy_of_structure_pt: string and structure have unequal length");

    paramT *P = get_updated_params(NULL, 1);
    int en = eval_pt(string, pt, s, s1, NULL, P, eos_debug);
    free(P);
    return en;
  } else
    return INF;
}

PUBLIC float
energy_of_circ_struct(const char *string,
                      const char *structure){

  float en;
  paramT *P = get_updated_params(NULL, 1);
  P->model_details.circ = 1;

  if(eos_debug > 0)
    en = vrna_eval_structure_verbose(string, structure, P, NULL);
  else
    en = vrna_eval_structure(string, structure, P);

  free(P);
  return en;
}

PUBLIC  float
energy_of_structure(const char *string,
                    const char *structure,
                    int verbosity_level){

  float en;
  paramT *P = get_updated_params(NULL, 1);

  if(verbosity_level > 0)
    en = vrna_eval_structure_verbose(string, structure, P, NULL);
  else
    en = vrna_eval_structure(string, structure, P);

  free(P);
  return en;
}

PUBLIC float
energy_of_struct_par( const char *string,
                      const char *structure,
                      paramT *parameters,
                      int verbosity_level){

  float en;

  paramT *P = get_updated_params(parameters, 1);

  if(verbosity_level > 0)
    en = vrna_eval_structure_verbose(string, structure, P, NULL);
  else
    en = vrna_eval_structure(string, structure, P);

  free(P);

  return en;
}


PUBLIC float
energy_of_gquad_structure(const char *string,
                          const char *structure,
                          int verbosity_level){

  float en;
  paramT *P = get_updated_params(NULL, 1);
  P->model_details.gquad = 1;

  if(verbosity_level > 0)
    en = vrna_eval_structure_verbose(string, structure, P, NULL);
  else
    en = vrna_eval_structure(string, structure, P);

  free(P);

  return en;
}

PUBLIC float
energy_of_gquad_struct_par( const char *string,
                            const char *structure,
                            paramT *parameters,
                            int verbosity_level){

  float en;
  paramT *P = get_updated_params(parameters, 1);
  P->model_details.gquad = 1;

  if(verbosity_level > 0)
    en = vrna_eval_structure_verbose(string, structure, P, NULL);
  else
    en = vrna_eval_structure(string, structure, P);

  free(P);

  return en;
}

PUBLIC int
energy_of_structure_pt( const char *string,
                        short *pt,
                        short *s,
                        short *s1,
                        int verbosity_level){

  if(pt && string){
    if(pt[0] != (short)strlen(string))
      nrerror("energy_of_structure_pt: string and structure have unequal length");

    paramT *P = get_updated_params(NULL, 1);
    int en = eval_pt(string, pt, s, s1, NULL, P, verbosity_level);
    free(P);
    return en;
  } else
    return INF;
}

PUBLIC int
energy_of_struct_pt_par(const char *string,
                        short *pt,
                        short *s,
                        short *s1,
                        paramT *parameters,
                        int verbosity_level){

  if(pt && string){
    if(pt[0] != (short)strlen(string))
      nrerror("energy_of_structure_pt: string and structure have unequal length");

    paramT *P = get_updated_params(parameters, 1);
    int en = eval_pt(string, pt, s, s1, NULL, P, verbosity_level);
    free(P);
    return en;
  } else
    return INF;
}

PUBLIC float
energy_of_circ_structure( const char *string,
                          const char *structure,
                          int verbosity_level){

  float en;
  paramT *P = get_updated_params(NULL, 1);
  P->model_details.circ = 1;

  if(verbosity_level > 0)
    en = vrna_eval_structure_verbose(string, structure, P, NULL);
  else
    en = vrna_eval_structure(string, structure, P);

  free(P);

  return en;
}

PUBLIC float
energy_of_circ_struct_par(const char *string,
                          const char *structure,
                          paramT *parameters,
                          int verbosity_level){

  float en;
  paramT *P = get_updated_params(parameters, 1);
  P->model_details.circ = 1;

  if(verbosity_level > 0)
    en = vrna_eval_structure_verbose(string, structure, P, NULL);
  else
    en = vrna_eval_structure(string, structure, P);

  free(P);

  return en;
}

PUBLIC int
loop_energy(short *pt,
            short *s,
            short *s1,
            int i){

  paramT *P = get_updated_params(NULL,1);

  int en = wrap_eval_loop_pt(pt, s, s1, i, P, eos_debug);

  free(P);
  return en;
}


PUBLIC float
energy_of_move( const char *string,
                const char *structure,
                int m1,
                int m2){

  paramT *P = get_updated_params(NULL, 1);
  float en  = vrna_eval_move(string, structure, m1, m2, P);

  free(P);
  return en;
}

PUBLIC int
energy_of_move_pt(short *pt,
                  short *s,
                  short *s1,
                  int m1,
                  int m2){

  paramT *P = get_updated_params(NULL, 1);
  int en    = vrna_eval_move_pt(pt, s, s1, m1, m2, P);

  free(P);
  return en;
}

