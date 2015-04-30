#ifndef VIENNA_RNA_PACKAGE_MULTIBRANCH_LOOPS_H
#define VIENNA_RNA_PACKAGE_MULTIBRANCH_LOOPS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/energy_par.h>
#include <ViennaRNA/params.h>
#include <ViennaRNA/constraints.h>
#include <ViennaRNA/exterior_loops.h>
#include <ViennaRNA/gquad.h>

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#ifdef ON_SAME_STRAND
#undef ON_SAME_STRAND
#endif

#define ON_SAME_STRAND(I,J,C)  (((I)>=(C))||((J)<(C)))

/**
 *  @addtogroup   loops
 *
 *  @{
 *
 *  @file multibranch_loops.h
 *  @brief Energy evaluation of multibranch loops for MFE and partition function calculations
 */

/**
 *  @def E_MLstem(A,B,C,D)
 *  <H2>Compute the Energy contribution of a Multiloop stem</H2>
 *  This definition is a wrapper for the E_Stem() funtion.
 *  It is substituted by an E_Stem() funtion call with argument
 *  extLoop=0, so the energy contribution returned reflects a
 *  stem introduced in a multiloop.<BR>
 *  As for the parameters B (si1) and C (sj1) of the substituted
 *  E_Stem() function, you can inhibit to take 5'-, 3'-dangles
 *  or mismatch contributions to be taken into account by passing
 *  -1 to these parameters.
 * 
 *  @see    E_Stem()
 *  @param  A The pair type of the stem-closing pair
 *  @param  B The 5'-mismatching nucleotide
 *  @param  C The 3'-mismatching nucleotide
 *  @param  D The datastructure containing scaled energy parameters
 *  @return   The energy contribution of the introduced multiloop stem
 */
INLINE  PRIVATE int E_MLstem( int type,
                              int si1,
                              int sj1,
                              vrna_param_t *P);

/**
 *  @def exp_E_MLstem(A,B,C,D)
 *  This is the partition function variant of @ref E_MLstem()
 *  @see E_MLstem()
 *  @return The Boltzmann weighted energy contribution of the introduced multiloop stem
 */
INLINE  PRIVATE double exp_E_MLstem(int type,
                                    int si1,
                                    int sj1,
                                    vrna_exp_param_t *P);



/**
 *  @brief Evaluate energy of a multi branch helices stacking onto closing pair (i,j)
 *
 *  Computes total free energy for coaxial stacking of (i.j) with (i+1.k) or (k+1.j-1)
 */
INLINE PRIVATE int E_mb_loop_stack(int i, int j, vrna_fold_compound *vc);

/**
 *  @brief  Backtrack the decomposition of a multi branch loop closed by @f$ (i,j) @f$
 *
 *  @param    vc          The #vrna_fold_compound filled with all relevant data for backtracking
 *  @param    i           5' position of base pair closing the loop (will be set to 5' position
 *                        of leftmost decomposed block upon successful backtracking)
 *  @param    j           3' position of base pair closing the loop (will be set to 3' position
 *                        of rightmost decomposed block upon successful backtracking)
 *  @param    k           Split position that delimits leftmost from rightmost block, [i,k] and
 *                        [k+1, j], respectively. (Will be set upon successful backtracking)
 *  @param    en          The energy contribution of the substructure enclosed by @f$ (i,j) @f$
 *  @param    component1  Type of leftmost block (1 = ML, 2 = C)
 *  @param    component2  Type of rightmost block (1 = ML, 2 = C)
 *  @returns              1, if backtracking succeeded, 0 otherwise.
 */
INLINE PRIVATE int
vrna_BT_mb_loop(vrna_fold_compound *vc,
                int *i,
                int *j,
                int *k,
                int en,
                int *component1,
                int *component2);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/


INLINE PRIVATE int
E_mb_loop_fast( int i,
                int j,
                vrna_fold_compound *vc,
                int *dmli1,
                int *dmli2){

  int decomp, en;
  unsigned char type, tt;
  short S_i1, S_j1;

  int               cp      = vc->cutpoint;
  char              *ptype  = vc->ptype;
  short             *S      = vc->sequence_encoding;
  int               *indx   = vc->jindx;
  char              *hc     = vc->hc->matrix;
  int               *hc_up  = vc->hc->up_ml;
  vrna_sc_t         *sc     = vc->sc;
  int               *fc     = vc->matrices->fc;
  vrna_param_t      *P      = vc->params;

  int ij            = indx[j] + i;
  int hc_decompose  = hc[ij];
  int e             = INF;
  int dangle_model  = P->model_details.dangles;
  int *rtype        = &(P->model_details.rtype[0]);

  type              = (unsigned char)ptype[ij];

  if(cp < 0){
    S_i1    = S[i+1];
    S_j1    = S[j-1];
  } else {
    S_i1  = ON_SAME_STRAND(i, i + 1, cp) ? S[i+1] : -1;
    S_j1  = ON_SAME_STRAND(j - 1, j, cp) ? S[j-1] : -1;
  }

  if(hc_decompose & VRNA_CONSTRAINT_CONTEXT_MB_LOOP){
    if((S_i1 >= 0) && (S_j1 >= 0)){ /* regular multi branch loop */
      decomp = dmli1[j-1];
      tt = rtype[type];
      switch(dangle_model){
        /* no dangles */
        case 0:   if(decomp != INF){
                    decomp += E_MLstem(tt, -1, -1, P);
                    if(sc){
                      if(sc->en_basepair)
                        decomp += sc->en_basepair[ij];
                    }
                  }
                  break;

        /* double dangles */
        case 2:   decomp += E_MLstem(tt, S_j1, S_i1, P);
                  if(sc){
                    if(sc->en_basepair)
                      decomp += sc->en_basepair[ij];
                  }
                  break;

        /* normal dangles, aka dangles = 1 || 3 */
        default:  if(decomp != INF){
                    decomp += E_MLstem(tt, -1, -1, P);
                    if(sc){
                      if(sc->en_basepair)
                        decomp += sc->en_basepair[ij];
                    }
                  }
                  if(hc_up[i+1]){
                    if(dmli2[j-1] != INF){
                      en = dmli2[j-1] + E_MLstem(tt, -1, S_i1, P) + P->MLbase;
                      if(sc){
                        if(sc->free_energies)
                          en += sc->free_energies[i+1][1];

                        if(sc->en_basepair)
                          en += sc->en_basepair[ij];
                      }
                      decomp = MIN2(decomp, en);
                    }
                  }
                  if(hc_up[j-1] && hc_up[i+1]){
                    if(dmli2[j-2] != INF){
                      en = dmli2[j-2] + E_MLstem(tt, S_j1, S_i1, P) + 2*P->MLbase;
                      if(sc){
                        if(sc->free_energies)
                          en += sc->free_energies[i+1][1]
                                + sc->free_energies[j-1][1];

                        if(sc->en_basepair)
                          en += sc->en_basepair[ij];
                      }
                      decomp = MIN2(decomp, en);
                    }
                  }
                  if(hc_up[j-1]){
                    if(dmli1[j-2] != INF){
                      en = dmli1[j-2] + E_MLstem(tt, S_j1, -1, P) + P->MLbase;
                      if(sc){
                        if(sc->free_energies)
                          en += sc->free_energies[j-1][1];

                        if(sc->en_basepair)
                          en += sc->en_basepair[ij];
                      }
                      decomp = MIN2(decomp, en);
                    }
                  }
                  break;
      }
      if(decomp != INF)
        e = decomp + P->MLclosing;
    }

    if(!ON_SAME_STRAND(i, j, cp)){ /* multibrach like cofold structure with cut somewhere between i and j */
      decomp = fc[i+1] + fc[j-1];
      tt = rtype[type];
      switch(dangle_model){
        case 0:   decomp += E_ExtLoop(tt, -1, -1, P);
                  break;
        case 2:   decomp += E_ExtLoop(tt, S_j1, S_i1, P);
                  break;
        default:  decomp += E_ExtLoop(tt, -1, -1, P);
                  if((hc_up[i+1]) && (hc_up[j-1])){
                    en     = fc[i+2] + fc[j-2] + E_ExtLoop(tt, S_j1, S_i1, P);
                    decomp = MIN2(decomp, en);
                  }
                  if(hc_up[i+1]){
                    en     = fc[i+2] + fc[j-1] + E_ExtLoop(tt, -1, S_i1, P);
                    decomp = MIN2(decomp, en);
                  }
                  if(hc_up[j-1]){
                    en     = fc[i+1] + fc[j-2] + E_ExtLoop(tt, S_j1, -1, P);
                    decomp = MIN2(decomp, en);
                  }
                  break;
      }
      e = MIN2(e, decomp);
    }
  }
  return e;
}

INLINE PRIVATE int
E_mb_loop_stack(int i,
                int j,
                vrna_fold_compound *vc){

  int e, decomp, en, i1k, k1j1, ij, k;
  unsigned char type, type_2;

  int               *indx   = vc->jindx;
  char              *hc     = vc->hc->matrix;
  int               *c      = vc->matrices->c;
  int               *fML    = vc->matrices->fML;
  vrna_param_t      *P      = vc->params;
  vrna_md_t         *md     = &(P->model_details);
  int               turn    = md->min_loop_size;
  char              *ptype  = vc->ptype;
  int               *rtype  = &(md->rtype[0]);
  vrna_sc_t         *sc     = vc->sc;

  e     = INF;
  ij    = indx[j] + i;
  type  = ptype[ij];

  if(hc[ij] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP){
    decomp = INF;
    k1j1  = indx[j-1] + i + 2 + turn + 1;
    for (k = i+2+turn; k < j-2-turn; k++, k1j1++){
      i1k   = indx[k] + i + 1;
      if(hc[i1k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
        type_2  = rtype[(unsigned char)ptype[i1k]];
        en      = c[i1k]+P->stack[type][type_2]+fML[k1j1];
        decomp  = MIN2(decomp, en);
      }
      if(hc[k1j1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
        type_2  = rtype[(unsigned char)ptype[k1j1]];
        en      = c[k1j1]+P->stack[type][type_2]+fML[i1k];
        decomp  = MIN2(decomp, en);
      }
    }
    /* no TermAU penalty if coax stack */
    decomp += 2*P->MLintern[1] + P->MLclosing;
    if(sc){
      if(sc->en_basepair)
        decomp += sc->en_basepair[ij];
    }
    e = decomp;

  }
  return e;
}

INLINE PRIVATE int
E_ml_rightmost_stem(int i,
                    int j,
                    vrna_fold_compound *vc){

  int               en;
  vrna_param_t      *P            = vc->params;
  int               length        = vc->length;
  short             *S            = vc->sequence_encoding;
  int               *indx         = vc->jindx;
  char              *hc           = vc->hc->matrix;
  int               *hc_up        = vc->hc->up_ml;
  vrna_sc_t         *sc           = vc->sc;
  int               *c            = vc->matrices->c;
  int               *fm           = (P->model_details.uniq_ML) ? vc->matrices->fM1 : vc->matrices->fML;
  int               *ggg          = vc->matrices->ggg;
  int               ij            = indx[j] + i;
  int               type          = vc->ptype[ij];
  int               hc_decompose  = hc[ij];
  int               dangle_model  = P->model_details.dangles;
  int               with_gquad    = P->model_details.gquad;
  int               cp            = vc->cutpoint;
  int               e             = INF;

  if(ON_SAME_STRAND(i - 1, i, cp)){
    if(ON_SAME_STRAND(j, j + 1, cp)){
      if(hc_decompose & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
        e = c[ij];
        if(e != INF){
          switch(dangle_model){
            case 2:   e += E_MLstem(type, (i==1) ? S[length] : S[i-1], S[j+1], P);
                      break;
            default:  e += E_MLstem(type, -1, -1, P);
                      break;
          }
        }
      }

      if(with_gquad)
        if(ON_SAME_STRAND(i, j, cp)){
          en  = ggg[ij] + E_MLstem(0, -1, -1, P);
          e   = MIN2(e, en);
        }

    }

    if(ON_SAME_STRAND(j - 1, j, cp))
      if(hc_up[j]){
        if(fm[indx[j - 1] + i] != INF){
          en = fm[indx[j - 1] + i] + P->MLbase;
          if(sc)
            if(sc->free_energies)
              en += sc->free_energies[j][1];

          e = MIN2(e, en);
        }
      }
  }
  return e;
}

INLINE PRIVATE int
E_ml_stems_fast(int i,
                int j,
                vrna_fold_compound *vc,
                int *fmi,
                int *dmli){

  int k, en, decomp, mm5, mm3, type_2, k1j, stop;

  int               length        = (int)vc->length;
  char              *ptype        = vc->ptype;
  short             *S            = vc->sequence_encoding;
  int               *indx         = vc->jindx;
  char              *hc           = vc->hc->matrix;
  int               *hc_up        = vc->hc->up_ml;
  vrna_sc_t         *sc           = vc->sc;
  int               *c            = vc->matrices->c;
  int               *fm           = vc->matrices->fML;
  vrna_param_t      *P            = vc->params;
  int               ij            = indx[j] + i;
  int               dangle_model  = P->model_details.dangles;
  int               turn          = P->model_details.min_loop_size;
  int               type          = ptype[ij];
  int               *rtype        = &(P->model_details.rtype[0]);
  int               circular      = P->model_details.circ;
  int               cp            = vc->cutpoint;
  int               e             = INF;

  /*  extension with one unpaired nucleotide at the right (3' site)
      or full branch of (i,j)
  */
  e = E_ml_rightmost_stem(i,j,vc);

  /*  extension with one unpaired nucleotide at 5' site
      and all other variants which are needed for odd
      dangle models
  */
  if(ON_SAME_STRAND(i - 1, i, cp)){
    switch(dangle_model){
      /* no dangles */
      case 0:   /* fall through */

      /* double dangles */
      case 2:   if(ON_SAME_STRAND(i, i + 1, cp))
                  if(hc_up[i]){
                    if(fm[ij + 1] != INF){
                      en = fm[ij + 1] + P->MLbase;
                      if(sc)
                        if(sc->free_energies)
                          en += sc->free_energies[i][1];
                      e = MIN2(e, en);
                    }
                  }
                break;

      /* normal dangles, aka dangle_model = 1 || 3 */
      default:  mm5 = ((i>1) || circular) ? S[i] : -1;
                mm3 = ((j<length) || circular) ? S[j] : -1;
                if(ON_SAME_STRAND(i, i + 1, cp))
                  if(hc_up[i]){
                    if(fm[ij+1] != INF){
                      en = fm[ij+1] + P->MLbase;
                      if(sc)
                        if(sc->free_energies)
                          en += sc->free_energies[i][1];
                      e = MIN2(e, en);
                    }
                    if(hc[ij+1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
                      if(c[ij+1] != INF){
                        type = ptype[ij+1];
                        en = c[ij+1] + E_MLstem(type, mm5, -1, P) + P->MLbase;
                        if(sc)
                          if(sc->free_energies)
                            en += sc->free_energies[i][1];
                        e = MIN2(e, en);
                      }
                    }
                  }

                if(ON_SAME_STRAND(j - 1, j, cp))
                  if(hc_up[j]){
                    if(hc[indx[j-1]+i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
                      if(c[indx[j-1]+i] != INF){
                        type = ptype[indx[j-1]+i];
                        en = c[indx[j-1]+i] + E_MLstem(type, -1, mm3, P) + P->MLbase;
                        if(sc)
                          if(sc->free_energies)
                            en += sc->free_energies[j][1];
                        e = MIN2(e, en);
                      }
                    }
                  }

                if(ON_SAME_STRAND(j - 1, j, cp) && ON_SAME_STRAND(i, i + 1, cp))
                  if(hc[indx[j-1]+i+1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
                    if(hc_up[i] && hc_up[j]){
                      if(c[indx[j-1]+i+1] != INF){
                        type = ptype[indx[j-1]+i+1];
                        en = c[indx[j-1]+i+1] + E_MLstem(type, mm5, mm3, P) + 2*P->MLbase;
                        if(sc)
                          if(sc->free_energies)
                            en += sc->free_energies[j][1] + sc->free_energies[i][1];
                        e = MIN2(e, en);
                      }
                    }
                  }
                break;
    }
  }

  /* modular decomposition -------------------------------*/
  k1j   = indx[j] + i + turn + 2;
  stop  = (cp > 0) ? (cp - 1) : (j - 2 - turn);
  for (decomp = INF, k = i + 1 + turn; k <= stop; k++, k1j++){
    if((fmi[k] != INF ) && (fm[k1j] != INF)){
      en = fmi[k] + fm[k1j];
      decomp = MIN2(decomp, en);
    }
  }
  k++; k1j++;
  for (;k <= j - 2 - turn; k++, k1j++){
    if((fmi[k] != INF) && (fm[k1j] != INF)){
      en = fmi[k] + fm[k1j];
      decomp = MIN2(decomp, en);
    }
  }

  dmli[j] = decomp;               /* store for use in fast ML decompositon */
  e = MIN2(e, decomp);

  /* coaxial stacking */
  if (dangle_model==3) {
    /* additional ML decomposition as two coaxially stacked helices */
    int ik;
    k1j = indx[j]+i+turn+2;
    for (decomp = INF, k = i + 1 + turn; k <= stop; k++, k1j++){
      ik = indx[k]+i;
      if((hc[ik] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) && (hc[k1j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC)){
        type    = rtype[(unsigned char)ptype[ik]];
        type_2  = rtype[(unsigned char)ptype[k1j]];
        en      = c[ik] + c[k1j] + P->stack[type][type_2];
        decomp  = MIN2(decomp, en);
      }
    }
    k++; k1j++;
    for (; k <= j-2-turn; k++, k1j++){
      ik = indx[k]+i;
      if((hc[ik] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) && (hc[k1j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC)){
        type    = rtype[(unsigned char)ptype[ik]];
        type_2  = rtype[(unsigned char)ptype[k1j]];
        en      = c[ik] + c[k1j] + P->stack[type][type_2];
        decomp  = MIN2(decomp, en);
      }
    }

    decomp += 2*P->MLintern[1];        /* no TermAU penalty if coax stack */
#if 0
        /* This is needed for Y shaped ML loops with coax stacking of
           interior pairts, but backtracking will fail if activated */
        DMLi[j] = MIN2(DMLi[j], decomp);
        DMLi[j] = MIN2(DMLi[j], DMLi[j-1]+P->MLbase);
        DMLi[j] = MIN2(DMLi[j], DMLi1[j]+P->MLbase);
        new_fML = MIN2(new_fML, DMLi[j]);
#endif
    e = MIN2(e, decomp);
  }

  fmi[j] = e;

  return e;
}



INLINE  PRIVATE int E_MLstem(int type, int si1, int sj1, vrna_param_t *P){
  int energy = 0;
  if(si1 >= 0 && sj1 >= 0){
    energy += P->mismatchM[type][si1][sj1];
  }
  else if (si1 >= 0){
    energy += P->dangle5[type][si1];
  }
  else if (sj1 >= 0){
    energy += P->dangle3[type][sj1];
  }

  if(type > 2)
    energy += P->TerminalAU;

  energy += P->MLintern[type];

  return energy;
}



INLINE PRIVATE double
exp_E_MLstem( int type,
              int si1,
              int sj1,
              vrna_exp_param_t *P){

  double energy = 1.0;
  if(si1 >= 0 && sj1 >= 0){
    energy = P->expmismatchM[type][si1][sj1];
  }
  else if(si1 >= 0){
    energy = P->expdangle5[type][si1];
  }
  else if(sj1 >= 0){
    energy = P->expdangle3[type][sj1];
  }

  if(type > 2)
    energy *= P->expTermAU;

  energy *= P->expMLintern[type];
  return energy;
}

INLINE PRIVATE int
vrna_BT_mb_loop_split(vrna_fold_compound *vc,
                      int *i,
                      int *j,
                      int *k,
                      int *l,
                      bondT *bp_stack,
                      int *stack_count){

  int length, cp, ij, ii, jj, fij, fi, u, en, *my_c, *my_fML, *my_ggg, turn, *idx, with_gquad, dangle_model, *rtype;
  unsigned char type, type_2;
  char          *ptype;
  short         *S1;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;
  vrna_sc_t     *sc;

  length        = vc->length;
  cp            = vc->cutpoint;
  P             = vc->params;
  md            = &(P->model_details);
  hc            = vc->hc;
  sc            = vc->sc;
  idx           = vc->jindx;
  ptype         = vc->ptype;
  rtype         = &(md->rtype[0]);
  S1            = vc->sequence_encoding;

  my_c    = vc->matrices->c;
  my_fML  = vc->matrices->fML;
  my_ggg  = vc->matrices->ggg;
  turn    = md->min_loop_size;
  with_gquad    = md->gquad;
  dangle_model  = md->dangles;

  ii = *i;
  jj = *j;

  /* nibble off unpaired 3' bases */
  do{
    fij = my_fML[idx[jj] + ii];
    fi  = (hc->up_ml[jj]) ? my_fML[idx[jj - 1] + ii] + P->MLbase : INF;

    if(sc){
      if(sc->free_energies)
        fi += sc->free_energies[jj][1];
    }

    if(--jj == 0)
      break;

  } while (fij == fi);
  jj++;

  /* nibble off unpaired 5' bases */
  do{
    fij = my_fML[idx[jj] + ii];
    fi  = (hc->up_ml[ii]) ? my_fML[idx[jj] + ii + 1] + P->MLbase : INF;

    if(sc){
      if(sc->free_energies)
        fi += sc->free_energies[ii][1];
    }

    if(++ii == jj)
      break;

  } while (fij == fi);
  ii--;

  if (jj < ii + turn + 2){ /* no more pairs */
    return 0;
  }

  ij = idx[jj] + ii;

  /* 1. test for single component */

  if(with_gquad){
    if(fij == my_ggg[ij] + E_MLstem(0, -1, -1, P)){
      *i = *j = -1;
      *k = *l = -1;
      vrna_BT_gquad_mfe(vc, ii, jj, bp_stack, stack_count);
      return 1;
    }
  }

  type  = (unsigned char)ptype[ij];
  en    = my_c[ij];

  switch(dangle_model){
    case 0:   if(hc->matrix[ij] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
                if(fij == en + E_MLstem(type, -1, -1, P)){
                  bp_stack[++(*stack_count)].i = ii;
                  bp_stack[(*stack_count)].j   = jj;
                  *i = *j = -1;
                  *k = ii;
                  *l = jj;
                  return 1;
                }
              }
              break;

    case 2:   if(hc->matrix[ij] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
                if(fij == en + E_MLstem(type, S1[ii - 1], S1[jj + 1], P)){
                  bp_stack[++(*stack_count)].i = ii;
                  bp_stack[(*stack_count)].j   = jj;
                  *i = *j = -1;
                  *k = ii;
                  *l = jj;
                  return 1;
                }
              }
              break;

    default:  if(hc->matrix[ij] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
                if(fij == en + E_MLstem(type, -1, -1, P)){
                  bp_stack[++(*stack_count)].i = ii;
                  bp_stack[(*stack_count)].j   = jj;
                  *i = *j = -1;
                  *k = ii;
                  *l = jj;
                  return 1;
                }
              }
              if(hc->matrix[ij + 1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
                if(hc->up_ml[ii]){
                  int tmp_en = fij;
                  if(sc)
                    if(sc->free_energies)
                      tmp_en -= sc->free_energies[ii][1];

                  type = (unsigned char)ptype[ij + 1];
                  if(tmp_en == my_c[ij+1] + E_MLstem(type, S1[ii], -1, P) + P->MLbase){
                    bp_stack[++(*stack_count)].i = ii + 1;
                    bp_stack[(*stack_count)].j   = jj;
                    *i = *j = -1;
                    *k = ii + 1;
                    *l = jj;
                    return 1;
                  }
                }
              }
              if(hc->matrix[idx[jj - 1] + ii] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
                if(hc->up_ml[jj]){
                  int tmp_en = fij;
                  if(sc)
                    if(sc->free_energies)
                      tmp_en -= sc->free_energies[jj][1];

                  type = (unsigned char)ptype[idx[jj - 1] + ii];
                  if(tmp_en == my_c[idx[jj - 1] + ii] + E_MLstem(type, -1, S1[jj], P) + P->MLbase){
                    bp_stack[++(*stack_count)].i = ii;
                    bp_stack[(*stack_count)].j   = jj - 1;
                    *i = *j = -1;
                    *k = ii;
                    *l = jj - 1;
                    return 1;
                  }
                }
              }
              if(hc->matrix[idx[jj - 1] + ii + 1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
                if(hc->up_ml[ii] && hc->up_ml[jj]){
                  int tmp_en = fij;
                  if(sc)
                    if(sc->free_energies)
                      tmp_en -= sc->free_energies[ii][1] + sc->free_energies[jj][1];

                  type = (unsigned char)ptype[idx[jj - 1] + ii + 1];
                  if(tmp_en == my_c[idx[jj - 1] + ii + 1] + E_MLstem(type, S1[ii], S1[jj], P) + 2 * P->MLbase){
                    bp_stack[++(*stack_count)].i = ii + 1;
                    bp_stack[(*stack_count)].j   = jj - 1;
                    *i = *j = -1;
                    *k = ii + 1;
                    *l = jj - 1;
                    return 1;
                  }
                }
              }
              break;
  }

  /* 2. Test for possible split point */
  for(u = ii + 1 + turn; u <= jj - 2 - turn; u++){
    ij    = idx[jj] + u + 1;
    type  = (unsigned char)ptype[ij];
    en    = my_c[ij];

    switch(dangle_model){
      case 0:   if(hc->matrix[ij] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
                  if(fij == my_fML[idx[u] + ii] + en + E_MLstem(type, -1, -1, P)){
                    bp_stack[++(*stack_count)].i = u + 1;
                    bp_stack[(*stack_count)].j   = jj;
                    *i = ii;
                    *j = u;
                    *k = u + 1;
                    *l = jj;
                    return 1;
                  }
                }
                break;

      case 2:   if(hc->matrix[ij] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
                  if(fij == my_fML[idx[u] + ii] + en + E_MLstem(type, S1[u], S1[jj + 1], P)){
                    bp_stack[++(*stack_count)].i = u + 1;
                    bp_stack[(*stack_count)].j   = jj;
                    *i = ii;
                    *j = u;
                    *k = u + 1;
                    *l = jj;
                    return 1;
                  }
                }
                break;

      default:  if(hc->matrix[ij] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
                  if(fij == my_fML[idx[u] + ii] + en + E_MLstem(type, -1, -1, P)){
                    bp_stack[++(*stack_count)].i = u + 1;
                    bp_stack[(*stack_count)].j   = jj;
                    *i = ii;
                    *j = u;
                    *k = u + 1;
                    *l = jj;
                    return 1;
                  }
                }
                if(hc->matrix[ij + 1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
                  if(hc->up_ml[u + 1]){
                    en = my_c[ij + 1];
                    if(sc)
                      if(sc->free_energies)
                        en += sc->free_energies[u + 1][1];

                    type = (unsigned char)ptype[ij + 1];
                    if(fij == my_fML[idx[u] + ii] + en + E_MLstem(type, S1[u + 1], -1, P) + P->MLbase){
                      bp_stack[++(*stack_count)].i = u + 2;
                      bp_stack[(*stack_count)].j   = jj;
                      *i = ii;
                      *j = u;
                      *k = u + 2;
                      *l = jj;
                      return 1;
                    }
                  }
                }
                if(hc->matrix[idx[jj - 1] + u + 1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
                  if(hc->up_ml[jj]){
                    en = my_c[idx[jj - 1] + u + 1];
                    if(sc)
                      if(sc->free_energies)
                        en += sc->free_energies[jj][1];

                    type = (unsigned char)ptype[idx[jj - 1] + u + 1];
                    if(fij == my_fML[idx[u] + ii] + en + E_MLstem(type, -1, S1[jj], P) + P->MLbase){
                      bp_stack[++(*stack_count)].i = u + 1;
                      bp_stack[(*stack_count)].j   = jj - 1;
                      *i = ii;
                      *j = u;
                      *k = u + 1;
                      *l = jj - 1;
                      return 1;
                    }
                  }
                }
                if(hc->matrix[idx[jj - 1] + u + 2] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
                  if(hc->up_ml[u + 1] && hc->up_ml[jj]){
                    en = my_c[idx[jj - 1] + u + 2];
                    if(sc)
                      if(sc->free_energies)
                        en += sc->free_energies[u + 1][1] + sc->free_energies[jj][1];

                    type = (unsigned char)ptype[idx[jj - 1] + u + 2];
                    if(fij == my_fML[idx[u] + ii] + en + E_MLstem(type, S1[u + 1], S1[jj], P) + 2 * P->MLbase){
                      bp_stack[++(*stack_count)].i = u + 2;
                      bp_stack[(*stack_count)].j   = jj - 1;
                      *i = ii;
                      *j = u;
                      *k = u + 2;
                      *l = jj - 1;
                      return 1;
                    }
                  }
                }
                break;
    }
  }

  /* 3. last chance! Maybe coax stack */
  if(dangle_model==3){
    int ik, k1j, tmp_en;
    for (k1j = idx[jj] + ii + turn + 2, u = ii + 1 + turn; u <= jj - 2 - turn; u++, k1j++) {
      ik = idx[u] + ii;
      if((hc->matrix[ik] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) && (hc->matrix[k1j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC)){
        type    = rtype[(unsigned char)ptype[ik]];
        type_2  = rtype[(unsigned char)ptype[k1j]];
        tmp_en  = my_c[ik] + my_c[k1j] + P->stack[type][type_2] + 2*P->MLintern[1];
        if (fij == tmp_en){
          bp_stack[++(*stack_count)].i = u + 1;
          bp_stack[(*stack_count)].j   = jj;
          *i = ii;
          *j = u;
          *k = u + 1;
          *l = jj;
          return 1;
        }
      }
    }
  }

  return 0;
}


INLINE PRIVATE int
vrna_BT_mb_loop(vrna_fold_compound *vc,
                int *i,
                int *j,
                int *k,
                int en,
                int *component1,
                int *component2){

  int ij, p, q, r, e, cp, *idx, turn, dangle_model, *my_c, *my_fML, *my_fc, *rtype;
  unsigned char type, type_2;
  char          *ptype;
  short         s5, s3, *S1;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;
  vrna_sc_t     *sc;

  cp            = vc->cutpoint;
  idx           = vc->jindx;
  ij            = idx[*j] + *i;
  S1            = vc->sequence_encoding;
  P             = vc->params;
  md            = &(P->model_details);
  hc            = vc->hc;
  sc            = vc->sc;
  my_c          = vc->matrices->c;
  my_fML        = vc->matrices->fML;
  my_fc         = vc->matrices->fc;
  turn          = md->min_loop_size;
  ptype         = vc->ptype;
  rtype         = &(md->rtype[0]);
  type          = (unsigned char)ptype[ij];
  type          = rtype[type];
  dangle_model  = md->dangles;

  p = *i + 1;
  q = *j - 1;

  r = q - turn - 1;

  if(hc->matrix[ij] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP){

    /* is it a fake multi-loop? */
    if(!ON_SAME_STRAND(*i, *j, cp)){
      int ii, jj;
      ii = jj = 0;
      e = my_fc[p] + my_fc[q];
      if(sc)
        if(sc->en_basepair)
          e += sc->en_basepair[ij];

      s5 = ON_SAME_STRAND(q, *j, cp) ? S1[q] : -1;
      s3 = ON_SAME_STRAND(*i, p, cp) ? S1[p] : -1;

      switch(dangle_model){
        case 0:   if(en == e + E_ExtLoop(type, -1, -1, P)){
                    ii=p, jj=q;
                  }
                  break;
        case 2:   if(en == e + E_ExtLoop(type, s5, s3, P)){
                    ii=p, jj=q;
                  }
                  break;
        default:  {
                    if(en == e + E_ExtLoop(type, -1, -1, P)){
                      ii=p, jj=q;
                      break;
                    }
                    if(hc->up_ext[p]){
                      e = my_fc[p + 1] + my_fc[q];
                      if(sc){
                        if(sc->free_energies)
                          e += sc->free_energies[p][1];
                        if(sc->en_basepair)
                          e += sc->en_basepair[ij];
                      }
                      if(en == e + E_ExtLoop(type, -1, s3, P)){
                        ii = p + 1; jj = q;
                        break;
                      }
                    }
                    if(hc->up_ext[q]){
                      e = my_fc[p] + my_fc[q - 1];
                      if(sc){
                        if(sc->free_energies)
                          e += sc->free_energies[q][1];
                        if(sc->en_basepair)
                          e += sc->en_basepair[ij];
                      }
                      if(en == e + E_ExtLoop(type, s5, -1, P)){
                        ii = p; jj = q - 1;
                        break;
                      }
                    }
                    if((hc->up_ext[q]) && (hc->up_ext[p])){
                      e = my_fc[p + 1] + my_fc[q - 1];
                      if(sc){
                        if(sc->free_energies)
                          e += sc->free_energies[p][1] + sc->free_energies[q][1];
                        if(sc->en_basepair)
                          e += sc->en_basepair[ij];
                      }
                      if(en == e + E_ExtLoop(type, s5, s3, P)){
                        ii = p + 1; jj = q - 1;
                        break;
                      }
                    }
                  }
                  break;
      }

      if(ii){ /* found a decomposition */
        *component1 = 3;
        *i          = ii;
        *k          = cp - 1;
        *j          = jj;
        *component2 = 4;
        return 1;
      }
    }

    /* true multi loop? */
    *component1 = *component2 = 1;  /* both components are MB loop parts by default */

    s5 = ON_SAME_STRAND(q, *j, cp) ? S1[q] : -1;
    s3 = ON_SAME_STRAND(*i, p, cp) ? S1[p] : -1;
                
    switch(dangle_model){
      case 0:   e = en - E_MLstem(type, -1, -1, P) - P->MLclosing;
                if(sc){
                  if(sc->en_basepair)
                    e -= sc->en_basepair[ij];
                }
                for(r = *i + 2 + turn; r < *j - 2 - turn; ++r){
                  if(e == my_fML[idx[r] + p] + my_fML[idx[q] + r + 1])
                    break;
                }
                break;

      case 2:   e = en - E_MLstem(type, s5, s3, P) - P->MLclosing;
                if(sc){
                  if(sc->en_basepair)
                    e -= sc->en_basepair[ij];
                }
                for(r = p + turn + 1; r < q - turn - 1; ++r){
                    if(e == my_fML[idx[r] + p] + my_fML[idx[q] + r + 1])
                      break;
                }
                break;

      default:  e = en - P->MLclosing;
                if(sc){
                  if(sc->en_basepair)
                    e -= sc->en_basepair[ij];
                }
                for(r = p + turn + 1; r < q - turn - 1; ++r){
                  if(e == my_fML[idx[r] + p] + my_fML[idx[q] + r + 1] + E_MLstem(type, -1, -1, P)){
                    break;
                  }
                  if(hc->up_ml[p]){
                    int tmp_en = e;
                    if(sc)
                      if(sc->free_energies)
                        tmp_en -= sc->free_energies[p][1];

                    if(tmp_en == my_fML[idx[r] + p + 1] + my_fML[idx[q] + r + 1] + E_MLstem(type, -1, s3, P) + P->MLbase){
                      p += 1;
                      break;
                    }
                  }
                  if(hc->up_ml[q]){
                    int tmp_en = e;
                    if(sc)
                      if(sc->free_energies)
                        tmp_en -= sc->free_energies[q][1];

                    if(tmp_en == my_fML[idx[r] + p] + my_fML[idx[q - 1] + r + 1] + E_MLstem(type, s5, -1, P) + P->MLbase){
                      q -= 1;
                      break;
                    }
                  }
                  if((hc->up_ml[p]) && (hc->up_ml[q])){
                    int tmp_en = e;
                    if(sc)
                      if(sc->free_energies)
                        tmp_en -= sc->free_energies[p][1] + sc->free_energies[q][1];

                    if(tmp_en == my_fML[idx[r] + p + 1] + my_fML[idx[q - 1] + r + 1] + E_MLstem(type, s5, s3, P) + 2 * P->MLbase){
                      p += 1;
                      q -= 1;
                      break;
                    }
                  }
                  /* coaxial stacking of (i.j) with (i+1.r) or (r.j-1) */
                  /* use MLintern[1] since coax stacked pairs don't get TerminalAU */
                  if(dangle_model == 3){
                    int tmp_en = e;
                    type = rtype[type];

                    if(hc->matrix[idx[r] + p] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
                      type_2 = rtype[(unsigned char)ptype[idx[r] + p]];
                      tmp_en = my_c[idx[r] + p] + P->stack[type][type_2] + my_fML[idx[q] + r + 1];
                      if(sc){
                        if(sc->en_basepair)
                          tmp_en += sc->en_basepair[ij];
                      }
                      if(en == tmp_en + 2 * P->MLintern[1] + P->MLclosing){
                        *component1 = 2;
                        break;
                      }
                    }

                    if(hc->matrix[idx[q] + r + 1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC){
                      type_2 = rtype[(unsigned char)ptype[idx[q] + r + 1]];
                      tmp_en = my_c[idx[q] + r + 1] + P->stack[type][type_2] + my_fML[idx[r] + p];
                      if(sc){
                        if(sc->en_basepair)
                          tmp_en += sc->en_basepair[ij];
                      }
                      if (en == tmp_en + 2 * P->MLintern[1] + P->MLclosing){
                        *component2 = 2;
                        break;
                      }
                    }
                  }
                }
                break;
    }
  }

  if(r <= *j - turn - 1){
    *i = p;
    *k = r;
    *j = q;
    return 1;
  } else {
#if 0
    /* Y shaped ML loops fon't work yet */
    if (dangle_model==3) {
      d5 = P->dangle5[tt][S1[j-1]];
      d3 = P->dangle3[tt][S1[i+1]];
      /* (i,j) must close a Y shaped ML loop with coax stacking */
      if (cij ==  fML[indx[j-2]+i+2] + mm + d3 + d5 + P->MLbase + P->MLbase) {
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
      bt_stack[++s].i = i1;
      bt_stack[s].j   = j1;
    }
#endif
  }

  return 0;
}


/**
 * @}
 */

#endif
