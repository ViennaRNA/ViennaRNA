/** \file eval.c */


/*
 *                Free energy evaluation
 *
 *                c Ivo Hofacker, Chrisoph Flamm
 *                original implementation by
 *                Walter Fontana
 *
 *                ViennaRNA Package >= v2.0 by Ronny Lorenz
 *
 *                Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/eval/gquad.h"
#include "ViennaRNA/eval/structures.h"


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

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
/* consensus structure variants below */
PRIVATE int
covar_energy_of_struct_pt(vrna_fold_compound_t  *fc,
                          const short           *pt);


PRIVATE int
stack_energy_covar_pt(vrna_fold_compound_t  *fc,
                      int                   i,
                      const short           *ptable);


PRIVATE int
covar_en_corr_of_loop_gquad(vrna_fold_compound_t  *fc,
                            int                   i,
                            int                   j,
                            const char            *structure,
                            const short           *pt,
                            const int             *loop_idx);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC float
vrna_eval_covar_structure(vrna_fold_compound_t  *fc,
                          const char            *structure)
{
  short         *pt;
  unsigned int  n_seq;
  int           res, gq, *loop_idx;
  vrna_md_t     *md;

  res   = 0;
  n_seq = 1;

  if ((fc) &&
      (fc->type == VRNA_FC_TYPE_COMPARATIVE) &&
      (structure)) {
    n_seq     = fc->n_seq;
    pt        = vrna_ptable(structure);
    md        = &(fc->params->model_details);
    gq        = md->gquad;
    md->gquad = 0;
    res       = covar_energy_of_struct_pt(fc, pt);
    md->gquad = gq;

    if (gq) {
      loop_idx  = vrna_loopidx_from_ptable(pt);
      res       -= covar_en_corr_of_loop_gquad(fc,
                                               1,
                                               fc->length,
                                               structure,
                                               pt,
                                               (const int *)loop_idx);
      free(loop_idx);
    }

    free(pt);
  }

  return (float)res / (100. * (float)n_seq);
}


/*
 #################################
 # STATIC helper functions below #
 #################################
 */


/* below are the consensus structure evaluation functions */

PRIVATE int
covar_energy_of_struct_pt(vrna_fold_compound_t  *fc,
                          const short           *pt)
{
  unsigned int  i, length;
  int           e;

  e       = 0;
  length  = fc->length;

  for (i = 1; i <= length; i++) {
    if (pt[i] == 0)
      continue;

    e += stack_energy_covar_pt(fc, i, pt);
    i = pt[i];
  }

  return e;
}


PRIVATE int
covar_en_corr_of_loop_gquad(vrna_fold_compound_t  *fc,
                            int                   i,
                            int                   j,
                            const char            *structure,
                            const short           *pt,
                            const int             *loop_idx)
{
  short         **S;
  unsigned int  pos, L, l[3], n_seq;
  int           en_covar, p, q, r, s, u, gq_en[2], num_elem,
                num_g, up_mis;
  vrna_param_t  *P;

  en_covar  = 0;
  n_seq     = fc->n_seq;
  S         = fc->S;
  P         = fc->params;
  q         = i;

  while ((pos = vrna_gq_parse(structure + q - 1, &L, l)) > 0) {
    q += pos - 1;
    p = q - 4 * L - l[0] - l[1] - l[2] + 1;
    if (q > j)
      break;

    /* we've found the first g-quadruplex at position [p,q] */
    vrna_E_consensus_gquad(L,
                           l,
                           (unsigned int)p,
                           fc->length,
                           n_seq,
                           (const short **)S,
                           (const unsigned int **)fc->a2s,
                           P,
                           gq_en);
    en_covar += gq_en[1];
    /* check if it's enclosed in a base pair */
    if (loop_idx[p] == 0) {
      q++;
      continue;                          /* g-quad in exterior loop */
    } else {
      /*  find its enclosing pair */
      num_elem  = 0;
      num_g     = 1;
      r         = p - 1;
      up_mis    = q - p + 1;

      /* seek for first pairing base located 5' of the g-quad */
      for (r = p - 1; !pt[r] && (r >= i); r--);

      if (r < pt[r]) {
        /* found the enclosing pair */
        s = pt[r];
      } else {
        num_elem++;
        r = pt[r] - 1;
        /* seek for next pairing base 5' of r */
        for (; !pt[r] && (r >= i); r--);

        if (r < pt[r]) {
          /* found the enclosing pair */
          s = pt[r];
        } else {
          /* hop over stems and unpaired nucleotides */
          while ((r > pt[r]) && (r >= i)) {
            if (pt[r]) {
              r = pt[r];
              num_elem++;
            }

            r--;
          }

          s = pt[r]; /* found the enclosing pair */
        }
      }

      /* now we have the enclosing pair (r,s) */

      u = q + 1;
      /* we know everything about the 5' part of this loop so check the 3' part */
      while (u < s) {
        if (structure[u - 1] == '.') {
          u++;
        } else if (structure[u - 1] == '+') {
          /* found another gquad */
          pos = vrna_gq_parse(structure + u - 1, &L, l);
          if (pos > 0) {
            vrna_E_consensus_gquad(L,
                                   l,
                                   (unsigned int)u,
                                   fc->length,
                                   n_seq,
                                   (const short **)S,
                                   (const unsigned int **)fc->a2s,
                                   P,
                                   gq_en);
            en_covar  += gq_en[1];
            up_mis    += pos;
            u         += pos;
            num_g++;
          }
        } else {
          /* we must have found a stem */
          num_elem++;
          en_covar += covar_en_corr_of_loop_gquad(fc,
                                                  u,
                                                  pt[u],
                                                  structure,
                                                  pt,
                                                  loop_idx);
          u = pt[u] + 1;
        }
      }
      /* we are done since we've found no other 3' structure element */

      q = s + 1;
    }
  }
  return en_covar;
}


PRIVATE int
stack_energy_covar_pt(vrna_fold_compound_t  *fc,
                      int                   i,
                      const short           *pt)
{
  /* calculate energy of substructure enclosed by (i,j) */
  int j, p, q, energy, *indx, *pscore;

  energy  = 0;
  indx    = fc->jindx;  /* index for moving in the triangle matrices c[] and fMl[]*/
  pscore  = fc->pscore; /* precomputed array of pair types */
  j       = pt[i];
  p       = i;
  q       = j;

  while (p < q) {
    /* process all stacks and internal loops */
    while (pt[++p] == 0);
    while (pt[--q] == 0);
    if ((pt[q] != (short)p) || (p > q))
      break;

    energy  += pscore[indx[j] + i];
    i       = p;
    j       = q;
  }  /* end while */

  /* p,q don't pair must have found hairpin or multiloop */

  if (p > q) {
    /* hairpin case */
    energy += pscore[indx[j] + i];
    return energy;
  }

  /* (i,j) is exterior pair of multiloop */
  energy += pscore[indx[j] + i];
  while (p < j) {
    /* add up the contributions of the substructures of the ML */
    energy  += stack_energy_covar_pt(fc, p, pt);
    p       = pt[p];
    /* search for next base pair in multiloop */
    while (pt[++p] == 0);
  }

  return energy;
}
