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
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>

#include "ViennaRNA/utils.h"
#include "ViennaRNA/structure_utils.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/loop_energies.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/cofold.h"
#include "ViennaRNA/eval.h"

#include "ViennaRNA/color_output.inc"

#ifdef ON_SAME_STRAND
#undef ON_SAME_STRAND
#endif

#define ON_SAME_STRAND(I, J, C)  (((I) >= (C)) || ((J) < (C)))

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
PRIVATE int
stack_energy(vrna_fold_compound_t *vc,
             int                  i,
             const short          *pt,
             FILE                 *file,
             int                  verbostiy_level);


PRIVATE int
energy_of_extLoop_pt(vrna_fold_compound_t *vc,
                     int                  i,
                     const short          *pt);


PRIVATE int
energy_of_ml_pt(vrna_fold_compound_t  *vc,
                int                   i,
                const short           *pt);


PRIVATE int
cut_in_loop(int         i,
            const short *pt,
            int         cp);


PRIVATE int
eval_pt(vrna_fold_compound_t  *vc,
        const short           *pt,
        FILE                  *file,
        int                   verbosity_level);


PRIVATE int
eval_circ_pt(vrna_fold_compound_t *vc,
             const short          *pt,
             FILE                 *file,
             int                  verbosity_level);


PRIVATE int
en_corr_of_loop_gquad(vrna_fold_compound_t  *vc,
                      int                   i,
                      int                   j,
                      const char            *structure,
                      const short           *pt,
                      FILE                  *file,
                      int                   verbosity_level);


PRIVATE vrna_param_t *
get_updated_params(vrna_param_t *parameters,
                   int          compat);


PRIVATE float
wrap_eval_structure(vrna_fold_compound_t  *vc,
                    const char            *structure,
                    const short           *pt,
                    FILE                  *file,
                    int                   verbosity);


/* consensus structure variants below */
PRIVATE int
covar_energy_of_struct_pt(vrna_fold_compound_t  *vc,
                          const short           *pt);


PRIVATE int
stack_energy_covar_pt(vrna_fold_compound_t  *vc,
                      int                   i,
                      const short           *ptable);


PRIVATE int
en_corr_of_loop_gquad_ali(vrna_fold_compound_t  *vc,
                          int                   i,
                          int                   j,
                          const char            *structure,
                          const short           *pt,
                          const int             *loop_idx);


PRIVATE int
covar_en_corr_of_loop_gquad(vrna_fold_compound_t  *vc,
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
vrna_eval_structure_v(vrna_fold_compound_t  *vc,
                      const char            *structure,
                      int                   verbosity_level,
                      FILE                  *file)
{
  if (strlen(structure) != vc->length)
    vrna_message_error("vrna_eval_structure_*: string and structure have unequal length");

  short *pt = vrna_ptable(structure);
  float en  = wrap_eval_structure(vc, structure, pt, file, verbosity_level);

  free(pt);
  return en;
}


PUBLIC float
vrna_eval_covar_structure(vrna_fold_compound_t  *vc,
                          const char            *structure)
{
  int   res, gq, *loop_idx;
  short *pt;

  pt                              = vrna_ptable(structure);
  res                             = 0;
  gq                              = vc->params->model_details.gquad;
  vc->params->model_details.gquad = 0;

  if (vc->type == VRNA_FC_TYPE_COMPARATIVE) {
    res = (int)((float)covar_energy_of_struct_pt(vc, pt) / (float)vc->n_seq);

    vc->params->model_details.gquad = gq;
    if (gq) {
      loop_idx  = vrna_loopidx_from_ptable(pt);
      res       -=
        (int)((float)covar_en_corr_of_loop_gquad(vc, 1, vc->length, structure, pt,
                                                 (const int *)loop_idx) / (float)vc->n_seq);
      free(loop_idx);
    }
  }

  free(pt);

  return (float)res / 100.;
}


PUBLIC int
vrna_eval_structure_pt_v(vrna_fold_compound_t *vc,
                         const short          *pt,
                         int                  verbosity_level,
                         FILE                 *file)
{
  if (pt && vc) {
    if (pt[0] != (short)vc->length)
      vrna_message_error("vrna_eval_structure_*: string and structure have unequal length");

    return eval_pt(vc, pt, file, verbosity_level);
  } else {
    return INF;
  }
}


PUBLIC int
vrna_eval_loop_pt_v(vrna_fold_compound_t  *vc,
                    int                   i,
                    const short           *pt,
                    int                   verbosity_level)
{
  /* compute energy of a single loop closed by base pair (i,j) */
  int           j, type, p, q, energy, cp;
  short         *s;
  vrna_param_t  *P;

  energy = INF;

  if (pt && vc) {
    P   = vc->params;
    cp  = vc->cutpoint;
    s   = vc->sequence_encoding2;

    vrna_sc_prepare(vc, VRNA_OPTION_MFE);

    if (i == 0) {
      /* evaluate exterior loop */
      energy = energy_of_extLoop_pt(vc, 0, pt);
      return energy;
    }

    j = pt[i];
    if (j < i)
      vrna_message_error("i is unpaired in loop_energy()");

    type = P->model_details.pair[s[i]][s[j]];
    if (type == 0) {
      type = 7;
      if (verbosity_level > VRNA_VERBOSITY_QUIET) {
        vrna_message_warning("bases %d and %d (%c%c) can't pair!",
                             i, j,
                             vrna_nucleotide_decode(s[i], &(P->model_details)),
                             vrna_nucleotide_decode(s[j], &(P->model_details)));
      }
    }

    p = i;
    q = j;


    while (pt[++p] == 0);
    while (pt[--q] == 0);
    if (p > q) {
      /* Hairpin */
      energy = vrna_eval_hp_loop(vc, i, j);
    } else if (pt[q] != (short)p) {
      /* multi-loop */
      int ii;
      ii      = cut_in_loop(i, (const short *)pt, cp);
      energy  =
        (ii == 0) ? energy_of_ml_pt(vc, i, (const short *)pt) : energy_of_extLoop_pt(vc,
                                                                                     ii,
                                                                                     (const short *)pt);
    } else {
      /* found interior loop */
      int type_2;
      type_2 = P->model_details.pair[s[q]][s[p]];
      if (type_2 == 0) {
        type_2 = 7;
        if (verbosity_level > VRNA_VERBOSITY_QUIET) {
          vrna_message_warning("bases %d and %d (%c%c) can't pair!",
                               p, q,
                               vrna_nucleotide_decode(s[p], &(P->model_details)),
                               vrna_nucleotide_decode(s[q], &(P->model_details)));
        }
      }

      energy = vrna_eval_int_loop(vc, i, j, p, q);
    }
  }

  return energy;
}


PUBLIC int
vrna_eval_move_pt(vrna_fold_compound_t  *vc,
                  short                 *pt,
                  int                   m1,
                  int                   m2)
{
  /*compute change in energy given by move (m1,m2)*/
  int           en_post, en_pre, i, j, k, l, len, cp;
  vrna_param_t  *P;

  len = vc->length;
  cp  = vc->cutpoint;
  P   = vc->params;

  k = (m1 > 0) ? m1 : -m1;
  l = (m2 > 0) ? m2 : -m2;
  /* first find the enclosing pair i<k<l<j */
  for (j = l + 1; j <= len; j++) {
    if (pt[j] <= 0)
      continue;             /* unpaired */

    if (pt[j] < k)
      break;              /* found it */

    if (pt[j] > j)
      j = pt[j];          /* skip substructure */
    else
      vrna_message_error("illegal move or broken pair table in vrna_eval_move_pt()\n"
                         "%d %d %d %d ", m1, m2, j, pt[j]);
  }
  i       = (j <= len) ? pt[j] : 0;
  en_pre  = vrna_eval_loop_pt(vc, i, (const short *)pt);
  en_post = 0;
  if (m1 < 0) {
    /*it's a delete move */
    en_pre  += vrna_eval_loop_pt(vc, k, (const short *)pt);
    pt[k]   = 0;
    pt[l]   = 0;
  } else {
    /* insert move */
    pt[k]   = l;
    pt[l]   = k;
    en_post += vrna_eval_loop_pt(vc, k, (const short *)pt);
  }

  en_post += vrna_eval_loop_pt(vc, i, (const short *)pt);
  /*  restore pair table */
  if (m1 < 0) {
    pt[k] = l;
    pt[l] = k;
  } else {
    pt[k] = 0;
    pt[l] = 0;
  }

  /* Cofolding -- Check if move changes COFOLD-Penalty */
  if (!ON_SAME_STRAND(k, l, cp)) {
    int p, c;
    p = c = 0;
    for (p = 1; p < cp;) {
      /* Count basepairs between two strands */
      if (pt[p] != 0) {
        if (ON_SAME_STRAND(p, pt[p], cp)) /* Skip stuff */
          p = pt[p];
        else if (++c > 1)
          break;                 /* Count a basepair, break if we have more than one */
      }

      p++;
    }
    if (m1 < 0 && c == 1) /* First and only inserted basepair */
      return en_post - en_pre - P->DuplexInit;
    else
    if (c == 0) /* Must have been a delete move */
      return en_post - en_pre + P->DuplexInit;
  }

  return en_post - en_pre;
}


/*
 #################################
 # STATIC helper functions below #
 #################################
 */
PRIVATE INLINE int
eval_ext_int_loop(vrna_fold_compound_t  *vc,
                  int                   i,
                  int                   j,
                  int                   p,
                  int                   q)
{
  unsigned char type, type_2;
  short         **SS, **S5, **S3, *S, si, sj, sp, sq;
  unsigned int  s, n_seq, **a2s;
  int           e, length;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_sc_t     *sc, **scs;

  length  = vc->length;
  P       = vc->params;
  md      = &(P->model_details);
  S       = vc->sequence_encoding;
  e       = INF;

  switch (vc->type) {
    case VRNA_FC_TYPE_SINGLE:
      si      = S[j + 1];
      sj      = S[i - 1];
      sp      = S[p - 1];
      sq      = S[q + 1];
      type    = (unsigned char)md->pair[S[j]][S[i]];
      type_2  = (unsigned char)md->pair[S[q]][S[p]];
      sc      = vc->sc;

      if (type == 0)
        type = 7;

      if (type_2 == 0)
        type_2 = 7;

      e = ubf_eval_ext_int_loop(i, j, p, q,
                                i - 1, j + 1, p - 1, q + 1,
                                si, sj, sp, sq,
                                type, type_2,
                                length,
                                P, sc);
      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      n_seq = vc->n_seq;
      SS    = vc->S;
      S5    = vc->S5;
      S3    = vc->S3;
      a2s   = vc->a2s;
      n_seq = vc->n_seq;
      scs   = vc->scs;

      for (e = 0, s = 0; s < n_seq; s++) {
        type = (unsigned char)md->pair[SS[s][j]][SS[s][i]];
        if (type == 0)
          type = 7;

        type_2 = (unsigned char)md->pair[SS[s][q]][SS[s][p]];
        if (type_2 == 0)
          type_2 = 7;

        sc = (scs && scs[s]) ? scs[s] : NULL;

        e += ubf_eval_ext_int_loop(a2s[s][i], a2s[s][j], a2s[s][p], a2s[s][q],
                                   a2s[s][i - 1], a2s[s][j + 1], a2s[s][p - 1], a2s[s][q + 1],
                                   S3[s][j], S5[s][i], S5[s][p], S3[s][q],
                                   type, type_2,
                                   a2s[s][length],
                                   P, sc);
      }

      break;
  }

  return e;
}


PRIVATE float
wrap_eval_structure(vrna_fold_compound_t  *vc,
                    const char            *structure,
                    const short           *pt,
                    FILE                  *file,
                    int                   verbosity)
{
  int res, gq, L, l[3];

  res                             = INF;
  gq                              = vc->params->model_details.gquad;
  vc->params->model_details.gquad = 0;

  switch (vc->type) {
    case VRNA_FC_TYPE_SINGLE:
      if (vc->params->model_details.circ)
        res = eval_circ_pt(vc, pt, file, verbosity);
      else
        res = eval_pt(vc, pt, file, verbosity);

      vc->params->model_details.gquad = gq;

      if (gq && (parse_gquad(structure, &L, l) > 0)) {
        if (verbosity > 0)
          print_eval_sd_corr(file == NULL ? stdout : file);

        res += en_corr_of_loop_gquad(vc, 1, vc->length, structure, pt, file, verbosity);
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      if (vc->params->model_details.circ)
        res = (int)((float)eval_circ_pt(vc, pt, file, verbosity) / (float)vc->n_seq);
      else
        res = (int)((float)eval_pt(vc, pt, file, verbosity) / (float)vc->n_seq);

      vc->params->model_details.gquad = gq;

      if (gq) {
        int *loop_idx = vrna_loopidx_from_ptable(pt);
        res +=
          (int)((float)en_corr_of_loop_gquad_ali(vc, 1, vc->length, structure, pt,
                                                 (const int *)loop_idx) / (float)vc->n_seq);
        free(loop_idx);
      }

      break;

    default:                      /* do nothing */
      break;
  }
  return (float)res / 100.;
}


PRIVATE int
eval_pt(vrna_fold_compound_t  *vc,
        const short           *pt,
        FILE                  *file,
        int                   verbosity_level)
{
  int   i, length, energy, cp;
  FILE  *out;

  out     = (file) ? file : stdout;
  length  = vc->length;
  cp      = vc->cutpoint;

  if (vc->params->model_details.gquad)
    vrna_message_warning("vrna_eval_*_pt: No gquadruplex support!\n"
                         "Ignoring potential gquads in structure!\n"
                         "Use e.g. vrna_eval_structure() instead!");

  vrna_sc_prepare(vc, VRNA_OPTION_MFE);

  energy = vc->params->model_details.backtrack_type ==
           'M' ? energy_of_ml_pt(vc, 0, pt) : energy_of_extLoop_pt(vc, 0, pt);

  if (verbosity_level > 0)
    print_eval_ext_loop(out, energy);

  for (i = 1; i <= length; i++) {
    if (pt[i] == 0)
      continue;

    energy  += stack_energy(vc, i, pt, out, verbosity_level);
    i       = pt[i];
  }
  for (i = 1; !ON_SAME_STRAND(i, length, cp); i++) {
    if (!ON_SAME_STRAND(i, pt[i], cp)) {
      energy += vc->params->DuplexInit;
      break;
    }
  }

  return energy;
}


PRIVATE int
eval_circ_pt(vrna_fold_compound_t *vc,
             const short          *pt,
             FILE                 *file,
             int                  verbosity_level)
{
  unsigned int  s, n_seq;
  int           i, j, length, energy, en0, degree;
  unsigned int  **a2s;
  vrna_param_t  *P;
  vrna_sc_t     *sc, **scs;
  FILE          *out;

  energy  = 0;
  en0     = 0;
  degree  = 0;
  length  = vc->length;
  P       = vc->params;
  sc      = (vc->type == VRNA_FC_TYPE_SINGLE) ? vc->sc : NULL;
  scs     = (vc->type == VRNA_FC_TYPE_COMPARATIVE) ? vc->scs : NULL;
  out     = (file) ? file : stdout;

  if (P->model_details.gquad)
    vrna_message_warning("vrna_eval_*_pt: No gquadruplex support!\n"
                         "Ignoring potential gquads in structure!\n"
                         "Use e.g. vrna_eval_structure() instead!");

  vrna_sc_prepare(vc, VRNA_OPTION_MFE);

  /* evaluate all stems in exterior loop */
  for (i = 1; i <= length; i++) {
    if (pt[i] == 0)
      continue;

    degree++;
    energy  += stack_energy(vc, i, (const short *)pt, out, verbosity_level);
    i       = pt[i];
  }

  /* find first stem */
  for (i = 1; pt[i] == 0; i++);
  j = pt[i];

  /* evaluate exterior loop itself */
  switch (degree) {
    case 0:   /* unstructured */
      switch (vc->type) {
        case VRNA_FC_TYPE_SINGLE:
          if (sc)
            if (sc->energy_up)
              en0 += sc->energy_up[1][length];

          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          n_seq = vc->n_seq;
          a2s   = vc->a2s;
          if (scs) {
            for (s = 0; s < n_seq; s++)
              if (scs[s] && scs[s]->energy_up)
                en0 += scs[s]->energy_up[1][a2s[s][length]];
          }

          break;
      }
      break;
    case 1:   /* hairpin loop */
      en0 = vrna_eval_ext_hp_loop(vc, i, j);
      break;

    case 2:   /* interior loop */
    {
      int p, q;
      /* seek to next pair */
      for (p = j + 1; pt[p] == 0; p++);
      q = pt[p];

      en0 = eval_ext_int_loop(vc, i, j, p, q);
    }
    break;

    default:  /* multibranch loop */
      en0 = energy_of_ml_pt(vc, 0, (const short *)pt);

      if (vc->type == VRNA_FC_TYPE_SINGLE)
        en0 -= E_MLstem(0, -1, -1, P);         /* remove virtual closing pair */

      break;
  }

  if (verbosity_level > 0)
    print_eval_ext_loop(out, en0);

  energy += en0;

  return energy;
}


/*---------------------------------------------------------------------------*/
/*  returns a correction term that may be added to the energy retrieved
 *  from energy_of_struct_par() to correct misinterpreted loops. This
 *  correction is necessary since energy_of_struct_par() will forget
 *  about the existance of gquadruplexes and just treat them as unpaired
 *  regions.
 *
 *  recursive variant
 */
PRIVATE int
en_corr_of_loop_gquad(vrna_fold_compound_t  *vc,
                      int                   i,
                      int                   j,
                      const char            *structure,
                      const short           *pt,
                      FILE                  *out,
                      int                   verbosity_level)
{
  char          *sequence;
  int           pos, tmp_e, energy, p, q, r, s, u, type, type2, L, l[3], *rtype, *loop_idx;
  int           num_elem, num_g, elem_i, elem_j, up_mis;
  short         *s1, *s2;
  vrna_param_t  *P;
  vrna_md_t     *md;

  sequence  = vc->sequence;
  loop_idx  = vrna_loopidx_from_ptable(pt);
  s1        = vc->sequence_encoding;
  s2        = vc->sequence_encoding2;
  P         = vc->params;
  md        = &(P->model_details);
  rtype     = &(md->rtype[0]);

  energy  = 0;
  q       = i;

  if (out == NULL)
    out = stdout;

  while ((pos = parse_gquad(structure + q - 1, &L, l)) > 0) {
    q += pos - 1;
    p = q - 4 * L - l[0] - l[1] - l[2] + 1;
    if (q > j)
      break;

    /* we've found the first g-quadruplex at position [p,q] */
    tmp_e   = E_gquad(L, l, P);
    energy  += tmp_e;
    if (verbosity_level > 0)
      print_eval_gquad(out, p, L, l, tmp_e);

    /* check if it's enclosed in a base pair */
    if (loop_idx[p] == 0) {
      q++;
      continue;                         /* g-quad in exterior loop */
    } else {
      /*  find its enclosing pair */
      num_elem  = 0;
      num_g     = 1;
      r         = p - 1;
      up_mis    = q - p + 1;

      /* seek for first pairing base located 5' of the g-quad */
      for (r = p - 1; !pt[r] && (r >= i); r--);
      if (r < i)
        vrna_message_error("this should not happen");

      if (r < pt[r]) {
        /* found the enclosing pair */
        s = pt[r];
      } else {
        num_elem++;
        elem_i  = pt[r];
        elem_j  = r;
        r       = pt[r] - 1;
        /* seek for next pairing base 5' of r */
        for (; !pt[r] && (r >= i); r--);
        if (r < i)
          vrna_message_error("so nich");

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
          if (r < i)
            vrna_message_error("so nich");

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
          pos = parse_gquad(structure + u - 1, &L, l);
          if (pos > 0) {
            tmp_e = E_gquad(L, l, P);

            if (verbosity_level > 0)
              print_eval_gquad(out, pos, L, l, tmp_e);

            energy  += tmp_e;
            up_mis  += pos;
            u       += pos;
            num_g++;
          }
        } else {
          /* we must have found a stem */
          if (!(u < pt[u]))
            vrna_message_error("wtf!");

          num_elem++;
          elem_i  = u;
          elem_j  = pt[u];
          energy  += en_corr_of_loop_gquad(vc, u, pt[u], structure, pt, out, verbosity_level);
          u       = pt[u] + 1;
        }
      }

      if (u != s) {
        vrna_message_error("what the heck");
      } else {
        int e_minus, e_plus;

        e_plus = e_minus = 0;

        /* we are done since we've found no other 3' structure element */
        switch (num_elem) {
          /* g-quad was misinterpreted as hairpin closed by (r,s) */
          case 0:   /* if(num_g == 1)
                     * if((p-r-1 == 0) || (s-q-1 == 0))
                     *  vrna_message_error("too few unpaired bases");
                     */
            e_minus = vrna_eval_hp_loop(vc, r, s);
            if (verbosity_level > 0)
              print_eval_hp_loop_revert(out, r, s, sequence[r - 1], sequence[s - 1], e_minus);

            type = md->pair[s2[r]][s2[s]];

            /* if we consider the G-Quadruplex, we have */
            if (num_g == 1) {
              /* a) an interior loop like structure */
              if (dangles == 2)
                e_plus += P->mismatchI[type][s1[r + 1]][s1[s - 1]];

              if (type > 2)
                e_plus += P->TerminalAU;

              e_plus += P->internal_loop[s - r - 1 - up_mis];
              if (verbosity_level > 0)
                print_eval_int_loop(out, r, s, sequence[r - 1], sequence[s - 1], p, q,
                                    sequence[p - 1], sequence[q - 1], e_plus);
            } else {
              /* or b) a multibranch loop like structure */
              e_plus = P->MLclosing
                       + E_MLstem(rtype[type], s1[s - 1], s1[r + 1], P)
                       + num_g * E_MLstem(0, -1, -1, P)
                       + (s - r - 1 - up_mis) * P->MLbase;

              if (verbosity_level > 0)
                print_eval_mb_loop(out, r, s, sequence[r - 1], sequence[s - 1], e_plus);
            }

            energy += e_plus - e_minus;
            break;

          /* g-quad was misinterpreted as interior loop closed by (r,s) with enclosed pair (elem_i, elem_j) */
          case 1:
            type    = md->pair[s2[r]][s2[s]];
            type2   = md->pair[s2[elem_i]][s2[elem_j]];
            e_plus  = P->MLclosing
                      + E_MLstem(rtype[type], s1[s - 1], s1[r + 1], P)
                      + (elem_i - r - 1 + s - elem_j - 1 - up_mis) * P->MLbase
                      + E_MLstem(type2, s1[elem_i - 1], s1[elem_j + 1], P);

            e_plus += num_g * E_MLstem(0, -1, -1, P);

            e_minus = vrna_eval_int_loop(vc, r, s, elem_i, elem_j);

            energy += e_plus - e_minus;

            if (verbosity_level > 0) {
              print_eval_int_loop_revert(out,
                                         r,
                                         s,
                                         sequence[r - 1],
                                         sequence[j - 1],
                                         elem_i,
                                         elem_j,
                                         sequence[elem_i - 1],
                                         sequence[elem_j - 1],
                                         e_minus);
              print_eval_mb_loop(out, r, s, sequence[r - 1], sequence[s - 1], e_plus);
            }

            break;

          /* gquad was misinterpreted as unpaired nucleotides in a multiloop */
          default:
            e_minus = (up_mis) * P->MLbase;
            e_plus  = num_g * E_MLstem(0, -1, -1, P);
            energy  += e_plus - e_minus;

            if (verbosity_level > 0) {
              print_eval_mb_loop_revert(out, r, s, sequence[r - 1], sequence[s - 1], e_minus);
              print_eval_mb_loop(out, r, s, sequence[r - 1], sequence[s - 1], e_plus);
            }

            break;
        }
      }

      q = s + 1;
    }
  }

  free(loop_idx);
  return energy;
}


PRIVATE int
stack_energy(vrna_fold_compound_t *vc,
             int                  i,
             const short          *pt,
             FILE                 *file,
             int                  verbosity_level)
{
  /* recursively calculate energy of substructure enclosed by (i,j) */

  int           ee, energy, j, p, q, type, cp;
  char          *string;
  short         *s;
  FILE          *out;
  vrna_param_t  *P;
  vrna_md_t     *md;

  cp      = vc->cutpoint;
  s       = vc->sequence_encoding2;
  P       = vc->params;
  md      = &(P->model_details);
  energy  = 0;
  out     = (file) ? file : stdout;

  j = pt[i];

  if (vc->type == VRNA_FC_TYPE_SINGLE) {
    string = vc->sequence;
    if (md->pair[s[i]][s[j]] == 0) {
      if (verbosity_level > VRNA_VERBOSITY_QUIET) {
        vrna_message_warning("bases %d and %d (%c%c) can't pair!",
                             i, j,
                             string[i - 1],
                             string[j - 1]);
      }
    }
  } else if (vc->type == VRNA_FC_TYPE_COMPARATIVE) {
    string = vc->cons_seq;
  } else {
    return INF;
  }

  p = i;
  q = j;

  while (p < q) {
    /* process all stacks and interior loops */
    while (pt[++p] == 0);
    while (pt[--q] == 0);
    if ((pt[q] != (short)p) || (p > q))
      break;

    ee = 0;
    if (vc->type == VRNA_FC_TYPE_SINGLE) {
      if (md->pair[s[q]][s[p]] == 0) {
        if (verbosity_level > VRNA_VERBOSITY_QUIET) {
          vrna_message_warning("bases %d and %d (%c%c) can't pair!",
                               p, q,
                               string[p - 1],
                               string[q - 1]);
        }
      }
    }

    ee = vrna_eval_int_loop(vc, i, j, p, q);

    if (verbosity_level > 0)
      print_eval_int_loop(out, i, j, string[i - 1], string[j - 1],
                          p, q, string[p - 1], string[q - 1],
                          ee);

    energy  += ee;
    i       = p;
    j       = q;
  } /* end while */

  /* p,q don't pair must have found hairpin or multiloop */

  if (p > q) {
    /* hairpin */
    ee      = vrna_eval_hp_loop(vc, i, j);
    energy  += ee;

    if (verbosity_level > 0)
      print_eval_hp_loop(out, i, j, string[i - 1], string[j - 1], ee);

    return energy;
  }

  /* (i,j) is exterior pair of multiloop */
  while (p < j) {
    /* add up the contributions of the substructures of the ML */
    energy  += stack_energy(vc, p, pt, out, verbosity_level);
    p       = pt[p];
    /* search for next base pair in multiloop */
    while (pt[++p] == 0);
  }

  switch (vc->type) {
    case VRNA_FC_TYPE_SINGLE:
    {
      int ii = cut_in_loop(i, pt, cp);
      ee = (ii == 0) ? energy_of_ml_pt(vc, i, pt) : energy_of_extLoop_pt(vc, ii, pt);
    }
    break;

    case VRNA_FC_TYPE_COMPARATIVE:
      ee = energy_of_ml_pt(vc, i, pt);
      break;

    default:
      break;                             /* this should never happen */
  }

  energy += ee;
  if (verbosity_level > 0)
    print_eval_mb_loop(out, i, j, string[i - 1], string[j - 1], ee);

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
energy_of_extLoop_pt(vrna_fold_compound_t *vc,
                     int                  i,
                     const short          *pt)
{
  int           energy, mm5, mm3, bonus, p, q, q_prev, length, dangle_model, n_seq, cp, ss, u,
                start;
  short         *s, *s1, **S, **S5, **S3;
  unsigned int  **a2s;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_sc_t     *sc, **scs;


  /* helper variables for dangles == 1 case */
  int           E3_available;   /* energy of 5' part where 5' mismatch is available for current stem */
  int           E3_occupied;    /* energy of 5' part where 5' mismatch is unavailable for current stem */


  /* initialize vars */
  length        = vc->length;
  cp            = vc->cutpoint;
  P             = vc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;

  energy  = 0;
  bonus   = 0;
  p       = start = (i == 0) ? 1 : i;
  q_prev  = -1;

  if (dangle_model % 2 == 1) {
    E3_available  = INF;
    E3_occupied   = 0;
  }

  /* seek to opening base of first stem */
  while (p <= length && !pt[p])
    p++;

  switch (vc->type) {
    case VRNA_FC_TYPE_SINGLE:
      s   = vc->sequence_encoding2;
      s1  = vc->sequence_encoding;
      sc  = vc->sc;

      /* add soft constraints for first unpaired nucleotides */
      if (sc) {
        if (sc->energy_up)
          bonus += sc->energy_up[start][p - start];

        /* how do we handle generalized soft constraints here ? */
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      S     = vc->S;
      S5    = vc->S5;
      S3    = vc->S3;
      a2s   = vc->a2s;
      n_seq = vc->n_seq;
      scs   = vc->scs;

      /* add soft constraints for first unpaired nucleotides */
      if (scs) {
        for (ss = 0; ss < n_seq; ss++) {
          if (scs[ss]) {
            if (scs[ss]->energy_up) {
              u     = a2s[ss][p] - a2s[ss][start];
              bonus += scs[ss]->energy_up[a2s[ss][start]][u];
            }

            /* how do we handle generalized soft constraints here ? */
          }
        }
      }

      break;

    default:
      return INF;
      break;
  }

  while (p < length) {
    int tt;
    /* p must have a pairing partner */
    q = (int)pt[p];

    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:     /* get type of base pair (p,q) */
        tt = md->pair[s[p]][s[q]];
        if (tt == 0)
          tt = 7;

        switch (dangle_model) {
          /* no dangles */
          case 0:
            energy += E_ExtLoop(tt, -1, -1, P);
            break;

          /* the beloved double dangles */
          case 2:
            mm5     = ((ON_SAME_STRAND(p - 1, p, cp)) && (p > 1))       ? s1[p - 1] : -1;
            mm3     = ((ON_SAME_STRAND(q, q + 1, cp)) && (q < length))  ? s1[q + 1] : -1;
            energy  += E_ExtLoop(tt, mm5, mm3, P);
            break;

          default:
          {
            int tmp;
            if (q_prev + 2 < p) {
              E3_available  = MIN2(E3_available, E3_occupied);
              E3_occupied   = E3_available;
            }

            mm5 = ((ON_SAME_STRAND(p - 1, p, cp)) && (p > 1) && !pt[p - 1])       ? s1[p - 1] : -1;
            mm3 = ((ON_SAME_STRAND(q, q + 1, cp)) && (q < length) && !pt[q + 1])  ? s1[q + 1] : -1;
            tmp = MIN2(
              E3_occupied + E_ExtLoop(tt, -1, mm3, P),
              E3_available + E_ExtLoop(tt, mm5, mm3, P)
              );
            E3_available = MIN2(
              E3_occupied + E_ExtLoop(tt, -1, -1, P),
              E3_available + E_ExtLoop(tt, mm5, -1, P)
              );
            E3_occupied = tmp;
          }
          break;
        }                             /* end switch dangle_model */
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        for (ss = 0; ss < n_seq; ss++) {
          /* get type of base pair (p,q) */
          tt = md->pair[S[ss][p]][S[ss][q]];
          if (tt == 0)
            tt = 7;

          switch (dangle_model) {
            case 0:
              energy += E_ExtLoop(tt, -1, -1, P);
              break;

            case 2:
              mm5     = (a2s[ss][p] > 1) ? S5[ss][p] : -1;
              mm3     = (a2s[ss][q] < a2s[ss][S[0][0]]) ? S3[ss][q] : -1;  /* why S[0][0] ??? */
              energy  += E_ExtLoop(tt, mm5, mm3, P);
              break;

            default:
              break;                                     /* odd dangles not implemented yet */
          }
        }
        break;

      default:
        break;                             /* this should never happen */
    }

    /* seek to the next stem */
    p       = q + 1;
    q_prev  = q;
    while (p <= length && !pt[p])
      p++;

    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:     /* add soft constraints for unpaired region */
        if (sc && (q_prev + 1 <= length)) {
          if (sc->energy_up)
            bonus += sc->energy_up[q_prev + 1][p - q_prev - 1];

          /* how do we handle generalized soft constraints here ? */
        }

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        if (scs) {
          for (ss = 0; ss < n_seq; ss++) {
            if (scs[ss]) {
              if (scs[ss]->energy_up) {
                u     = a2s[ss][p] - a2s[ss][q_prev + 1];
                bonus += scs[ss]->energy_up[a2s[ss][q_prev + 1]][u];
              }
            }
          }
        }

        break;

      default:
        break;                             /* this should never happen */
    }

    if (p == i)
      break; /* cut was in loop */
  }

  if (dangle_model % 2 == 1)
    energy = MIN2(E3_occupied, E3_available);

  return energy + bonus;
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
energy_of_ml_pt(vrna_fold_compound_t  *vc,
                int                   i,
                const short           *pt)
{
  int           energy, cx_energy, tmp, tmp2, best_energy = INF, bonus, *idx, cp, dangle_model,
                logML, circular, *rtype, ss, n, n_seq;
  int           i1, j, p, q, q_prev, q_prev2, u, uu, x, type, count, mm5, mm3, tt, ld5, new_cx,
                dang5, dang3, dang;
  int           e_stem, e_stem5, e_stem3, e_stem53;
  int           mlintern[NBPAIRS + 1];
  short         *s, *s1, **S, **S5, **S3;
  unsigned int  **a2s;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_sc_t     *sc, **scs;

  /* helper variables for dangles == 1|5 case */
  int           E_mm5_available;    /* energy of 5' part where 5' mismatch of current stem is available */
  int           E_mm5_occupied;     /* energy of 5' part where 5' mismatch of current stem is unavailable */
  int           E2_mm5_available;   /* energy of 5' part where 5' mismatch of current stem is available with possible 3' dangle for enclosing pair (i,j) */
  int           E2_mm5_occupied;    /* energy of 5' part where 5' mismatch of current stem is unavailable with possible 3' dangle for enclosing pair (i,j) */

  n   = vc->length;
  cp  = vc->cutpoint;
  P   = vc->params;
  md  = &(P->model_details);
  idx = vc->jindx;

  circular      = md->circ;
  dangle_model  = md->dangles;
  logML         = md->logML;
  rtype         = &(md->rtype[0]);

  bonus = 0;

  if (i >= pt[i]) {
    vrna_message_warning("energy_of_ml_pt: i is not 5' base of a closing pair!");
    return INF;
  }

  j = (i == 0) ? n + 1 : (int)pt[i];

  switch (vc->type) {
    case VRNA_FC_TYPE_SINGLE:
      s   = vc->sequence_encoding2;
      s1  = vc->sequence_encoding;
      sc  = vc->sc;

      if (i != 0) {
        /* (i,j) is closing pair of multibranch loop, add soft constraints */
        if (sc)
          if (sc->energy_bp)
            bonus += sc->energy_bp[idx[j] + i];
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      S     = vc->S;
      S5    = vc->S5;
      S3    = vc->S3;
      a2s   = vc->a2s;
      n_seq = vc->n_seq;
      scs   = vc->scs;

      if ((dangle_model % 2) || (dangle_model > 2) || (dangle_model < 0)) {
        vrna_message_warning(
          "consensus structure evaluation for odd dangle models not implemented (yet)!");
        return INF;
      }

      if (i != 0) {
        /* (i,j) is closing pair of multibranch loop, add soft constraints */
        if (scs) {
          for (ss = 0; ss < n_seq; ss++)
            if (scs[ss] && scs[ss]->energy_bp)
              bonus += scs[ss]->energy_bp[idx[j] + i];
        }
      }

      break;

    default:
      return INF;
      break;
  }

  /* init the variables */
  energy  = 0;
  u       = 0;     /* the total number of unpaired nucleotides */
  p       = i + 1;
  q_prev  = i - 1;
  q_prev2 = i;


  for (x = 0; x <= NBPAIRS; x++)
    mlintern[x] = P->MLintern[x];

  /* seek to opening base of first stem */
  while (p <= j && !pt[p])
    p++;

  /* add bonus energies for first stretch of unpaired nucleotides */
  switch (vc->type) {
    case VRNA_FC_TYPE_SINGLE:
      u += p - i - 1;
      if (sc)
        if (sc->energy_up)
          bonus += sc->energy_up[i + 1][u];

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      if (scs) {
        for (ss = 0; ss < n_seq; ss++) {
          uu = a2s[ss][p] - a2s[ss][i + 1];
          if (scs[ss] && scs[ss]->energy_up)
            bonus += scs[ss]->energy_up[a2s[ss][i + 1]][uu];

          u += uu;
        }
      } else {
        for (ss = 0; ss < n_seq; ss++)
          u += a2s[ss][p] - a2s[ss][i + 1];
      }

      break;

    default:
      break;                             /* this should never happen */
  }

  switch (dangle_model) {
    case 0:
      switch (vc->type) {
        case VRNA_FC_TYPE_SINGLE:
          while (p < j) {
            /* p must have a pairing partner */
            q = (int)pt[p];
            /* get type of base pair (p,q) */
            tt = md->pair[s[p]][s[q]];
            if (tt == 0)
              tt = 7;

            energy += E_MLstem(tt, -1, -1, P);

            /* seek to the next stem */
            p       = q + 1;
            q_prev  = q_prev2 = q;
            while (p <= j && !pt[p])
              p++;
            u += p - q - 1;                                     /* add unpaired nucleotides */

            if (sc)
              if (sc->energy_up)
                bonus += sc->energy_up[q + 1][p - q - 1];
          }

          /* now lets get the energy of the enclosing stem */
          if (i > 0) {
            /* actual closing pair */
            tt = md->pair[s[j]][s[i]];
            if (tt == 0)
              tt = 7;

            energy += E_MLstem(tt, -1, -1, P);
          } else {
            /* virtual closing pair */
            energy += E_MLstem(0, -1, -1, P);
          }

          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          while (p < j) {
            /* p must have a pairing partner */
            q = (int)pt[p];
            for (ss = 0; ss < n_seq; ss++) {
              /* get type of base pair (p,q) */
              tt = md->pair[S[ss][p]][S[ss][q]];
              if (tt == 0)
                tt = 7;

              energy += E_MLstem(tt, -1, -1, P);
            }

            /* seek to the next stem */
            p       = q + 1;
            q_prev  = q_prev2 = q;
            while (p <= j && !pt[p])
              p++;

            /* add unpaired nucleotides and possible soft constraints */
            if (scs) {
              for (ss = 0; ss < n_seq; ss++) {
                uu = a2s[ss][p] - a2s[ss][q + 1];
                if (scs[ss] && scs[ss]->energy_up)
                  bonus += sc->energy_up[a2s[ss][q + 1]][uu];

                u += uu;
              }
            } else {
              for (ss = 0; ss < n_seq; ss++)
                u += a2s[ss][p] - a2s[ss][q + 1];
            }
          }

          /* now lets get the energy of the enclosing stem */
          if (i > 0) {
            /* actual closing pair */
            for (ss = 0; ss < n_seq; ss++) {
              tt = md->pair[S[ss][j]][S[ss][i]];
              if (tt == 0)
                tt = 7;

              energy += E_MLstem(tt, -1, -1, P);
            }
          }

          break;

        default:
          break;                                     /* this should never happen */
      }
      break;

    case 2:
      switch (vc->type) {
        case VRNA_FC_TYPE_SINGLE:
          while (p < j) {
            /* p must have a pairing partner */
            q = (int)pt[p];
            /* get type of base pair (p,q) */
            tt = md->pair[s[p]][s[q]];
            if (tt == 0)
              tt = 7;

            mm5     = ON_SAME_STRAND(p - 1, p, cp) ? s1[p - 1] : -1;
            mm3     = ON_SAME_STRAND(q, q + 1, cp) ? s1[q + 1] : -1;
            energy  += E_MLstem(tt, mm5, mm3, P);

            /* seek to the next stem */
            p       = q + 1;
            q_prev  = q_prev2 = q;
            while (p <= j && !pt[p])
              p++;
            u += p - q - 1;                                     /* add unpaired nucleotides */

            if (sc)
              if (sc->energy_up)
                bonus += sc->energy_up[q + 1][p - q - 1];
          }
          if (i > 0) {
            /* actual closing pair */
            tt = md->pair[s[j]][s[i]];
            if (tt == 0)
              tt = 7;

            mm5     = ON_SAME_STRAND(j - 1, j, cp) ? s1[j - 1] : -1;
            mm3     = ON_SAME_STRAND(i, i + 1, cp) ? s1[i + 1] : -1;
            energy  += E_MLstem(tt, mm5, mm3, P);
          } else {
            /* virtual closing pair */
            energy += E_MLstem(0, -1, -1, P);
          }

          break;

        case VRNA_FC_TYPE_COMPARATIVE:
          while (p < j) {
            /* p must have a pairing partner */
            q = (int)pt[p];

            for (ss = 0; ss < n_seq; ss++) {
              /* get type of base pair (p,q) */
              tt = md->pair[S[ss][p]][S[ss][q]];
              if (tt == 0)
                tt = 7;

              mm5     = ((a2s[ss][p] > 1) || circular) ? S5[ss][p] : -1;
              mm3     = ((a2s[ss][q] < a2s[ss][S[0][0]]) || circular) ? S3[ss][q] : -1;
              energy  += E_MLstem(tt, mm5, mm3, P);
            }

            /* seek to the next stem */
            p       = q + 1;
            q_prev  = q_prev2 = q;
            while (p <= j && !pt[p])
              p++;

            /* add unpaired nucleotides and possible soft constraints */
            if (scs) {
              for (ss = 0; ss < n_seq; ss++) {
                uu = a2s[ss][p] - a2s[ss][q + 1];
                if (scs[ss] && scs[ss]->energy_up)
                  bonus += sc->energy_up[a2s[ss][q + 1]][uu];

                u += uu;
              }
            } else {
              for (ss = 0; ss < n_seq; ss++)
                u += a2s[ss][p] - a2s[ss][q + 1];
            }
          }

          if (i > 0) {
            /* actual closing pair */
            for (ss = 0; ss < n_seq; ss++) {
              tt = md->pair[S[ss][j]][S[ss][i]];
              if (tt == 0)
                tt = 7;

              mm5     = S5[ss][j];
              mm3     = S3[ss][i];
              energy  += E_MLstem(tt, mm5, mm3, P);
            }
          }

          break;

        default:
          break;                                     /* this should never happen */
      }
      break;

    case 3:   /* we treat helix stacking different */
      for (count = 0; count < 2; count++) {
        /* do it twice */
        ld5 = 0;         /* 5' dangle energy on prev pair (type) */
        if (i == 0) {
          j     = (unsigned int)pt[0] + 1;
          type  = 0;         /* no pair */
        } else {
          j     = (unsigned int)pt[i];
          type  = P->model_details.pair[s[j]][s[i]];
          if (type == 0)
            type = 7;

          /* prime the ld5 variable */
          if (ON_SAME_STRAND(j - 1, j, cp)) {
            ld5 = P->dangle5[type][s1[j - 1]];
            if ((p = (unsigned int)pt[j - 2]) && ON_SAME_STRAND(j - 2, j - 1, cp))
              if (P->dangle3[P->model_details.pair[s[p]][s[j - 2]]][s1[j - 1]] < ld5)
                ld5 = 0;
          }
        }

        i1        = i;
        p         = i + 1;
        u         = 0;
        energy    = 0;
        cx_energy = INF;
        do {
          /* walk around the multi-loop */
          new_cx = INF;

          /* hop over unpaired positions */
          while (p <= (unsigned int)pt[0] && pt[p] == 0)
            p++;

          /* memorize number of unpaired positions */
          u += p - i1 - 1;

          if (sc)
            if (sc->energy_up)
              bonus += sc->energy_up[i1 + 1][p - i1 - 1];

          /* get position of pairing partner */
          if (p == (unsigned int)pt[0] + 1) {
            q   = 0;
            tt  = 0;              /* virtual root pair */
          } else {
            q = (unsigned int)pt[p];
            /* get type of base pair P->q */
            tt = P->model_details.pair[s[p]][s[q]];
            if (tt == 0)
              tt = 7;
          }

          energy    += mlintern[tt];
          cx_energy += mlintern[tt];

          dang5 = dang3 = 0;
          if ((ON_SAME_STRAND(p - 1, p, cp)) && (p > 1))
            dang5 = P->dangle5[tt][s1[p - 1]];          /* 5'dangle of pq pair */

          if ((ON_SAME_STRAND(i1, i1 + 1, cp)) && (i1 < (unsigned int)s[0]))
            dang3 = P->dangle3[type][s1[i1 + 1]];       /* 3'dangle of previous pair */

          switch (p - i1 - 1) {
            case 0:           /* adjacent helices */
              if (i1 != 0) {
                if (ON_SAME_STRAND(i1, p, cp)) {
                  new_cx = energy + P->stack[rtype[type]][rtype[tt]];
                  /* subtract 5'dangle and TerminalAU penalty */
                  new_cx += -ld5 - mlintern[tt] - mlintern[type] + 2 * mlintern[1];
                }

                ld5     = 0;
                energy  = MIN2(energy, cx_energy);
              }

              break;
            case 1:           /* 1 unpaired base between helices */
              dang    = MIN2(dang3, dang5);
              energy  = energy + dang;
              ld5     = dang - dang3;
              /* may be problem here: Suppose
               * cx_energy>energy, cx_energy+dang5<energy
               * and the following helices are also stacked (i.e.
               * we'll subtract the dang5 again */
              if (cx_energy + dang5 < energy) {
                energy  = cx_energy + dang5;
                ld5     = dang5;
              }

              new_cx = INF;   /* no coax stacking with mismatch for now */
              break;
            default:          /* many unpaired base between helices */
              energy  += dang5 + dang3;
              energy  = MIN2(energy, cx_energy + dang5);
              new_cx  = INF;                 /* no coax stacking possible */
              ld5     = dang5;
              break;
          }
          type      = tt;
          cx_energy = new_cx;
          i1        = q;
          p         = q + 1;
        } while (q != i);
        best_energy = MIN2(energy, best_energy);         /* don't use cx_energy here */
        /* skip a helix and start again */
        while (pt[p] == 0)
          p++;
        if (i == (unsigned int)pt[p])
          break;

        i = (unsigned int)pt[p];
      }         /* end doing it twice */
      energy = best_energy;
      break;

    default:
      E_mm5_available = E2_mm5_available = INF;
      E_mm5_occupied  = E2_mm5_occupied = 0;
      while (p < j) {
        /* p must have a pairing partner */
        q = (int)pt[p];
        /* get type of base pair (p,q) */
        tt = P->model_details.pair[s[p]][s[q]];
        if (tt == 0)
          tt = 7;

        if (q_prev + 2 < p) {
          E_mm5_available = MIN2(E_mm5_available, E_mm5_occupied);
          E_mm5_occupied  = E_mm5_available;
        }

        if (q_prev2 + 2 < p) {
          E2_mm5_available  = MIN2(E2_mm5_available, E2_mm5_occupied);
          E2_mm5_occupied   = E2_mm5_available;
        }

        mm5       = ((ON_SAME_STRAND(p - 1, p, cp)) && !pt[p - 1])  ? s1[p - 1] : -1;
        mm3       = ((ON_SAME_STRAND(q, q + 1, cp)) && !pt[q + 1])  ? s1[q + 1] : -1;
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
        p       = q + 1;
        q_prev  = q_prev2 = q;
        while (p <= j && !pt[p])
          p++;
        u += p - q - 1;         /* add unpaired nucleotides */

        if (sc)
          if (sc->energy_up)
            bonus += sc->energy_up[q + 1][p - q - 1];
      }
      if (i > 0) {
        /* actual closing pair */
        type = P->model_details.pair[s[j]][s[i]];
        if (type == 0)
          type = 7;

        mm5 = ((ON_SAME_STRAND(j - 1, j, cp)) && !pt[j - 1])  ? s1[j - 1] : -1;
        mm3 = ((ON_SAME_STRAND(i, i + 1, cp)) && !pt[i + 1])  ? s1[i + 1] : -1;
        if (q_prev + 2 < p) {
          E_mm5_available = MIN2(E_mm5_available, E_mm5_occupied);
          E_mm5_occupied  = E_mm5_available;
        }

        if (q_prev2 + 2 < p) {
          E2_mm5_available  = MIN2(E2_mm5_available, E2_mm5_occupied);
          E2_mm5_occupied   = E2_mm5_available;
        }

        e_stem    = E_MLstem(type, -1, -1, P);
        e_stem5   = E_MLstem(type, mm5, -1, P);
        e_stem3   = E_MLstem(type, -1, mm3, P);
        e_stem53  = E_MLstem(type, mm5, mm3, P);
      } else {
        /* virtual closing pair */
        e_stem = e_stem5 = e_stem3 = e_stem53 = E_MLstem(0, -1, -1, P);
      }

      /* now lets see how we get the minimum including the enclosing stem */
      energy  = E_mm5_occupied + e_stem;
      energy  = MIN2(energy, E_mm5_available + e_stem5);
      energy  = MIN2(energy, E_mm5_available + e_stem);
      energy  = MIN2(energy, E2_mm5_occupied + e_stem3);
      energy  = MIN2(energy, E2_mm5_occupied + e_stem);
      energy  = MIN2(energy, E2_mm5_available + e_stem53);
      energy  = MIN2(energy, E2_mm5_available + e_stem3);
      energy  = MIN2(energy, E2_mm5_available + e_stem5);
      energy  = MIN2(energy, E2_mm5_available + e_stem);
      break;
  }/* end switch dangle_model */

  switch (vc->type) {
    case VRNA_FC_TYPE_SINGLE:
      energy += P->MLclosing;
      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      energy += P->MLclosing * n_seq;
      break;

    default:
      break;
  }

  /* logarithmic ML loop energy if logML */
  /* does this work for comparative predictions as well? */
  if (logML && (u > 6))
    energy += 6 * P->MLbase + (int)(P->lxc * log((double)u / 6.));
  else
    energy += (u * P->MLbase);

  return energy + bonus;
}


PRIVATE int
cut_in_loop(int         i,
            const short *pt,
            int         cp)
{
  /* walk around the loop;  return j pos of pair after cut if
   * cut_point in loop else 0 */
  int p, j;

  p = j = pt[i];
  do {
    i = pt[p];
    p = i + 1;
    while (pt[p] == 0)
      p++;
  } while (p != j && ON_SAME_STRAND(i, p, cp));
  return ON_SAME_STRAND(i, p, cp) ? 0 : j;
}


/* below are the consensus structure evaluation functions */

PRIVATE int
covar_energy_of_struct_pt(vrna_fold_compound_t  *vc,
                          const short           *pt)
{
  int e       = 0;
  int length  = vc->length;
  int i;

  for (i = 1; i <= length; i++) {
    if (pt[i] == 0)
      continue;

    e += stack_energy_covar_pt(vc, i, pt);
    i = pt[i];
  }

  return e;
}


PRIVATE int
en_corr_of_loop_gquad_ali(vrna_fold_compound_t  *vc,
                          int                   i,
                          int                   j,
                          const char            *structure,
                          const short           *pt,
                          const int             *loop_idx)
{
  int           pos, cnt, tmp_e, energy, p, q, r, s, u, type, gq_en[2];
  int           num_elem, num_g, elem_i, elem_j, up_mis;
  int           L, l[3];

  short         **S           = vc->S;
  short         **S5          = vc->S5;
  short         **S3          = vc->S3;
  vrna_param_t  *P            = vc->params;
  vrna_md_t     *md           = &(P->model_details);
  int           n_seq         = vc->n_seq;
  int           dangle_model  = md->dangles;

  energy  = 0;
  q       = i;
  while ((pos = parse_gquad(structure + q - 1, &L, l)) > 0) {
    q += pos - 1;
    p = q - 4 * L - l[0] - l[1] - l[2] + 1;
    if (q > j)
      break;

    /* we've found the first g-quadruplex at position [p,q] */
    E_gquad_ali_en(p, L, l, (const short **)S, n_seq, gq_en, P);
    tmp_e   = gq_en[0];
    energy  += tmp_e;
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
      if (r < i)
        vrna_message_error("this should not happen");

      if (r < pt[r]) {
        /* found the enclosing pair */
        s = pt[r];
      } else {
        num_elem++;
        elem_i  = pt[r];
        elem_j  = r;
        r       = pt[r] - 1;
        /* seek for next pairing base 5' of r */
        for (; !pt[r] && (r >= i); r--);
        if (r < i)
          vrna_message_error("so nich");

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
          if (r < i)
            vrna_message_error("so nich");

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
          pos = parse_gquad(structure + u - 1, &L, l);
          if (pos > 0) {
            E_gquad_ali_en(u, L, l, (const short **)S, n_seq, gq_en, P);
            tmp_e   = gq_en[0];
            energy  += tmp_e;
            up_mis  += pos;
            u       += pos;
            num_g++;
          }
        } else {
          /* we must have found a stem */
          if (!(u < pt[u]))
            vrna_message_error("wtf!");

          num_elem++;
          elem_i  = u;
          elem_j  = pt[u];
          energy  += en_corr_of_loop_gquad_ali(vc,
                                               u,
                                               pt[u],
                                               structure,
                                               pt,
                                               loop_idx);
          u = pt[u] + 1;
        }
      }
      if (u != s) {
        vrna_message_error("what the ...");
      } else {
        int e_minus, e_plus, e_temp;

        e_plus = e_minus = 0;

        /* we are done since we've found no other 3' structure element */
        switch (num_elem) {
          /* g-quad was misinterpreted as hairpin closed by (r,s) */
          case 0:   /* if(num_g == 1)
                     * if((p-r-1 == 0) || (s-q-1 == 0))
                     *  vrna_message_error("too few unpaired bases");
                     */
            e_minus = vrna_eval_hp_loop(vc, r, s);
            /* if we consider the G-Quadruplex, we have */
            if (num_g == 1) {
              /* a) an interior loop like structure */
              for (cnt = 0; cnt < n_seq; cnt++) {
                type = md->pair[S[cnt][r]][S[cnt][s]];
                if (type == 0)
                  type = 7;

                if (dangle_model == 2)
                  e_plus += P->mismatchI[type][S3[cnt][r]][S5[cnt][s]];

                if (type > 2)
                  e_plus += P->TerminalAU;
              }

              e_plus += n_seq * P->internal_loop[s - r - 1 - up_mis];
            } else {
              /* or b) a multibranch loop like structure */
              for (cnt = 0; cnt < n_seq; cnt++) {
                type = md->pair[S[cnt][s]][S[cnt][r]];
                if (type == 0)
                  type = 7;

                e_plus += E_MLstem(type, S5[cnt][s], S3[cnt][r], P);
              }

              e_temp = num_g * E_MLstem(0, -1, -1, P) +
                       P->MLclosing +
                       (elem_i - r - 1 + s - elem_j - 1 - up_mis) * P->MLbase;

              e_plus += n_seq * e_temp;
            }

            energy += e_plus - e_minus;
            break;

          /* g-quad was misinterpreted as interior loop closed by (r,s) with enclosed pair (elem_i, elem_j) */
          case 1:
            e_minus = vrna_eval_int_loop(vc, r, s, elem_i, elem_j);

            for (cnt = 0; cnt < n_seq; cnt++) {
              type = md->pair[S[cnt][s]][S[cnt][r]];
              if (type == 0)
                type = 7;

              e_plus += E_MLstem(type, S5[cnt][s], S3[cnt][r], P);

              type = md->pair[S[cnt][elem_i]][S[cnt][elem_j]];
              if (type == 0)
                type = 7;

              e_plus += E_MLstem(type, S5[cnt][elem_i], S3[cnt][elem_j], P);
            }

            e_temp = num_g * E_MLstem(0, -1, -1, P) +
                     P->MLclosing +
                     (elem_i - r - 1 + s - elem_j - 1 - up_mis) * P->MLbase;

            e_plus += n_seq * e_temp;

            energy += e_plus - e_minus;

            break;
          /* gquad was misinterpreted as unpaired nucleotides in a multiloop */
          default:
            e_minus = (up_mis) * P->MLbase * n_seq;
            e_plus  = num_g * E_MLstem(0, -1, -1, P) * n_seq;
            energy  += e_plus - e_minus;

            break;
        }
      }

      q = s + 1;
    }
  }

  return energy;
}


PRIVATE int
covar_en_corr_of_loop_gquad(vrna_fold_compound_t  *vc,
                            int                   i,
                            int                   j,
                            const char            *structure,
                            const short           *pt,
                            const int             *loop_idx)
{
  int           pos, en_covar, p, q, r, s, u, gq_en[2];
  int           num_elem, num_g, up_mis;
  int           L, l[3];

  short         **S   = vc->S;
  vrna_param_t  *P    = vc->params;
  int           n_seq = vc->n_seq;

  en_covar  = 0;
  q         = i;
  while ((pos = parse_gquad(structure + q - 1, &L, l)) > 0) {
    q += pos - 1;
    p = q - 4 * L - l[0] - l[1] - l[2] + 1;
    if (q > j)
      break;

    /* we've found the first g-quadruplex at position [p,q] */
    E_gquad_ali_en(p, L, l, (const short **)S, n_seq, gq_en, P);
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
      if (r < i)
        vrna_message_error("this should not happen");

      if (r < pt[r]) {
        /* found the enclosing pair */
        s = pt[r];
      } else {
        num_elem++;
        r = pt[r] - 1;
        /* seek for next pairing base 5' of r */
        for (; !pt[r] && (r >= i); r--);
        if (r < i)
          vrna_message_error("so nich");

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
          if (r < i)
            vrna_message_error("so nich");

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
          pos = parse_gquad(structure + u - 1, &L, l);
          if (pos > 0) {
            E_gquad_ali_en(u, L, l, (const short **)S, n_seq, gq_en, P);
            en_covar  += gq_en[1];
            up_mis    += pos;
            u         += pos;
            num_g++;
          }
        } else {
          /* we must have found a stem */
          if (!(u < pt[u]))
            vrna_message_error("wtf!");

          num_elem++;
          en_covar += covar_en_corr_of_loop_gquad(vc,
                                                  u,
                                                  pt[u],
                                                  structure,
                                                  pt,
                                                  loop_idx);
          u = pt[u] + 1;
        }
      }
      if (u != s) {
        vrna_message_error("what the ...");
      } else {
        /* we are done since we've found no other 3' structure element */
      }

      q = s + 1;
    }
  }
  return en_covar;
}


PRIVATE int
stack_energy_covar_pt(vrna_fold_compound_t  *vc,
                      int                   i,
                      const short           *pt)
{
  /* calculate energy of substructure enclosed by (i,j) */
  int *indx   = vc->jindx;                      /* index for moving in the triangle matrices c[] and fMl[]*/
  int *pscore = vc->pscore;                     /* precomputed array of pair types */

  int energy = 0;
  int j, p, q;

  j = pt[i];
  p = i;
  q = j;
  while (p < q) {
    /* process all stacks and interior loops */
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
    energy  += stack_energy_covar_pt(vc, p, pt);
    p       = pt[p];
    /* search for next base pair in multiloop */
    while (pt[++p] == 0);
  }

  return energy;
}
