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

#ifdef _WIN32
#ifdef __MINGW32__
#include <unistd.h>
#else
#include "ViennaRNA/unistd_win.h"
#endif
#else
#include <unistd.h>
#endif

#include <limits.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/structures.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/cofold.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/datastructures/char_stream.h"
#include "ViennaRNA/eval.h"

#include "ViennaRNA/color_output.inc"

#define   ADD_OR_INF(a, b)     (((a) != INF) && ((b) != INF) ?  (a) + (b) : INF)

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
stack_energy(vrna_fold_compound_t *fc,
             int                  i,
             const short          *pt,
             vrna_cstr_t          output_stream,
             int                  verbostiy_level);


PRIVATE int
energy_of_ml_pt(vrna_fold_compound_t  *fc,
                int                   i,
                const short           *pt);


PRIVATE int
eval_pt(vrna_fold_compound_t  *fc,
        const short           *pt,
        vrna_cstr_t           output_stream,
        int                   verbosity_level);


PRIVATE int
eval_circ_pt(vrna_fold_compound_t *fc,
             const short          *pt,
             vrna_cstr_t          output_stream,
             int                  verbosity_level);


PRIVATE int
en_corr_of_loop_gquad(vrna_fold_compound_t  *fc,
                      int                   i,
                      int                   j,
                      const char            *structure,
                      const short           *pt,
                      const int             *loop_idx,
                      vrna_cstr_t           output_stream,
                      int                   verbosity_level);


PRIVATE float
wrap_eval_structure(vrna_fold_compound_t  *fc,
                    const char            *structure,
                    const short           *pt,
                    vrna_cstr_t           output_stream,
                    int                   verbosity);


PRIVATE int
energy_of_extLoop_pt(vrna_fold_compound_t *fc,
                     unsigned int         begin,
                     const short          *pt);


PRIVATE int
energy_of_ext_loop_components(vrna_fold_compound_t  *fc,
                              const short           *pt,
                              vrna_cstr_t           output_stream,
                              int                   verbosity_level);


PRIVATE int
first_pair_after_last_nick(unsigned int i,
                           unsigned int j,
                           const short  *pt,
                           unsigned int *sn);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC float
vrna_eval_structure_v(vrna_fold_compound_t  *fc,
                      const char            *structure,
                      int                   verbosity_level,
                      FILE                  *file)
{
  short       *pt;
  float       en;
  vrna_cstr_t output_stream;

  en = (float)INF / 100.;

  if ((fc) &&
      (structure)) {
    if (strlen(structure) != fc->length) {
      vrna_message_warning("vrna_eval_structure_*: "
                           "string and structure have unequal length (%d vs. %d)",
                           fc->length,
                           strlen(structure));
      return en;
    }

    output_stream = vrna_cstr(fc->length, (file) ? file : stdout);
    pt            = vrna_ptable(structure);
    en            = wrap_eval_structure(fc,
                                        structure,
                                        pt,
                                        output_stream,
                                        verbosity_level);

    vrna_cstr_fflush(output_stream);
    vrna_cstr_free(output_stream);

    free(pt);
  }

  return en;
}


PUBLIC float
vrna_eval_structure_cstr(vrna_fold_compound_t *fc,
                         const char           *structure,
                         int                  verbosity_level,
                         vrna_cstr_t          output_stream)
{
  short *pt;
  float en;

  en = (float)INF / 100.;

  if ((fc) &&
      (structure)) {
    if (strlen(structure) != fc->length) {
      vrna_message_warning("vrna_eval_structure_*: "
                           "string and structure have unequal length (%d vs. %d)",
                           fc->length,
                           strlen(structure));
      return en;
    }

    pt  = vrna_ptable(structure);
    en  = wrap_eval_structure(fc,
                              structure,
                              pt,
                              output_stream,
                              verbosity_level);

    free(pt);
  }

  return en;
}


PUBLIC int
vrna_eval_structure_pt_v(vrna_fold_compound_t *fc,
                         const short          *pt,
                         int                  verbosity_level,
                         FILE                 *file)
{
  int         e;
  vrna_cstr_t output_stream;

  e = INF;

  if ((fc) &&
      (pt)) {
    if (pt[0] != (short)fc->length) {
      vrna_message_warning("vrna_eval_structure_*: "
                           "string and structure have unequal length (%d vs. %d)",
                           fc->length,
                           pt[0]);
      return INF;
    }

    output_stream = vrna_cstr(fc->length, (file) ? file : stdout);
    e             = eval_pt(fc,
                            pt,
                            output_stream,
                            verbosity_level);

    vrna_cstr_fflush(output_stream);
    vrna_cstr_free(output_stream);
  }

  return e;
}


PUBLIC int
vrna_eval_loop_pt_v(vrna_fold_compound_t  *fc,
                    int                   i,
                    const short           *pt,
                    int                   verbosity_level)
{
  /* compute energy of a single loop closed by base pair (i,j) */
  short         *s;
  unsigned int  *sn, begin;
  int           j, p, q, energy;
  vrna_md_t     *md;

  energy = INF;

  if ((fc) &&
      (pt)) {
    md  = &(fc->params->model_details);
    sn  = fc->strand_number;
    s   = fc->sequence_encoding2;

    vrna_sc_prepare(fc, VRNA_OPTION_MFE);

    /* evaluate exterior loop ? */
    if (i == 0)
      return energy_of_extLoop_pt(fc, 0, pt);

    j = pt[i];
    if (j < i) {
      vrna_message_warning("vrna_eval_loop_pt*: "
                           "i = %d is unpaired in loop_energy()",
                           i);
      return INF;
    }

    if (md->pair[s[i]][s[j]] == 0) {
      if (verbosity_level > VRNA_VERBOSITY_QUIET) {
        vrna_message_warning("bases %d and %d (%c%c) can't pair!",
                             i, j,
                             vrna_nucleotide_decode(s[i], md),
                             vrna_nucleotide_decode(s[j], md));
      }
    }

    p = i;
    q = j;


    while (pt[++p] == 0);
    while (pt[--q] == 0);

    /* check, whether this is a base pair enclosing an external loop, i.e. strand-nick in loop */
    if ((fc->strands > 1) &&
        ((begin = first_pair_after_last_nick(p, q, pt, sn)) != 0))
      return energy_of_extLoop_pt(fc, begin, pt);

    if (p > q) {
      /* Hairpin */
      energy = vrna_eval_hp_loop(fc, i, j);
    } else if (pt[q] != (short)p) {
      /* multi-loop */
      energy = energy_of_ml_pt(fc, i, (const short *)pt);
    } else {
      /* found interior loop */
      if (md->pair[s[q]][s[p]] == 0) {
        if (verbosity_level > VRNA_VERBOSITY_QUIET) {
          vrna_message_warning("bases %d and %d (%c%c) can't pair!",
                               p, q,
                               vrna_nucleotide_decode(s[p], md),
                               vrna_nucleotide_decode(s[q], md));
        }
      }

      energy = vrna_eval_int_loop(fc, i, j, p, q);
    }
  }

  return energy;
}


PUBLIC int
vrna_eval_move_pt(vrna_fold_compound_t  *fc,
                  short                 *pt,
                  int                   m1,
                  int                   m2)
{
  /*compute change in energy given by move (m1,m2)*/
  int en, en_post, en_pre, i, j, k, l, len;

  en = INF;

  if ((fc) &&
      (pt)) {
    len = fc->length;
    k   = (m1 > 0) ? m1 : -m1;
    l   = (m2 > 0) ? m2 : -m2;

    /* first find the enclosing pair i<k<l<j */
    for (j = l + 1; j <= len; j++) {
      if (pt[j] <= 0)
        continue;             /* unpaired */

      if (pt[j] < k)
        break;              /* found it */

      if (pt[j] > j) {
        j = pt[j];          /* skip substructure */
      } else {
        vrna_message_warning("vrna_eval_move_pt: "
                             "illegal move or broken pair table in vrna_eval_move_pt()\n"
                             "%d %d %d %d ", m1, m2, j, pt[j]);
        return en;
      }
    }

    i       = (j <= len) ? pt[j] : 0;
    en_pre  = vrna_eval_loop_pt(fc, i, (const short *)pt);
    en_post = 0;

    if (m1 < 0) {
      /*it's a delete move */
      en_pre  += vrna_eval_loop_pt(fc, k, (const short *)pt);
      pt[k]   = 0;
      pt[l]   = 0;
    } else {
      /* insert move */
      pt[k]   = l;
      pt[l]   = k;
      en_post += vrna_eval_loop_pt(fc, k, (const short *)pt);
    }

    en_post += vrna_eval_loop_pt(fc, i, (const short *)pt);

    /*  restore pair table */
    if (m1 < 0) {
      pt[k] = l;
      pt[l] = k;
    } else {
      pt[k] = 0;
      pt[l] = 0;
    }

    en = en_post - en_pre;
  }

  return en;
}


/*
 #################################
 # STATIC helper functions below #
 #################################
 */
PRIVATE INLINE int
eval_ext_int_loop(vrna_fold_compound_t  *fc,
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

  length  = fc->length;
  P       = fc->params;
  md      = &(P->model_details);
  S       = fc->sequence_encoding;
  e       = INF;

  switch (fc->type) {
    case VRNA_FC_TYPE_COMPARATIVE:
      n_seq = fc->n_seq;
      SS    = fc->S;
      S5    = fc->S5;
      S3    = fc->S3;
      a2s   = fc->a2s;
      n_seq = fc->n_seq;
      scs   = fc->scs;

      for (e = 0, s = 0; s < n_seq; s++) {
        type    = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
        type_2  = vrna_get_ptype_md(SS[s][q], SS[s][p], md);

        sc = (scs && scs[s]) ? scs[s] : NULL;

        e += ubf_eval_ext_int_loop(a2s[s][i], a2s[s][j], a2s[s][p], a2s[s][q],
                                   a2s[s][i - 1], a2s[s][j + 1], a2s[s][p - 1], a2s[s][q + 1],
                                   S3[s][j], S5[s][i], S5[s][p], S3[s][q],
                                   type, type_2,
                                   a2s[s][length],
                                   P, sc);
      }

      break;

    default:
      si      = S[j + 1];
      sj      = S[i - 1];
      sp      = S[p - 1];
      sq      = S[q + 1];
      type    = vrna_get_ptype_md(S[j], S[i], md);
      type_2  = vrna_get_ptype_md(S[q], S[p], md);
      sc      = fc->sc;

      e = ubf_eval_ext_int_loop(i, j, p, q,
                                i - 1, j + 1, p - 1, q + 1,
                                si, sj, sp, sq,
                                type, type_2,
                                length,
                                P, sc);
      break;
  }

  return e;
}


PRIVATE float
wrap_eval_structure(vrna_fold_compound_t  *fc,
                    const char            *structure,
                    const short           *pt,
                    vrna_cstr_t           output_stream,
                    int                   verbosity)
{
  unsigned int  n_seq;
  int           res, gq, *loop_idx, L, l[3];
  float         energy;
  vrna_md_t     *md;

  energy    = (float)INF / 100.;
  n_seq     = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  md        = &(fc->params->model_details);
  gq        = md->gquad;
  md->gquad = 0;

  if (md->circ)
    res = eval_circ_pt(fc, pt, output_stream, verbosity);
  else
    res = eval_pt(fc, pt, output_stream, verbosity);

  md->gquad = gq;

  if (gq && (parse_gquad(structure, &L, l) > 0)) {
    if (verbosity > 0)
      vrna_cstr_print_eval_sd_corr(output_stream);

    loop_idx  = vrna_loopidx_from_ptable(pt);
    res       += en_corr_of_loop_gquad(fc,
                                       1,
                                       fc->length,
                                       structure,
                                       pt,
                                       (const int *)loop_idx,
                                       output_stream,
                                       verbosity);
    free(loop_idx);
  }

  energy = (float)res / (100. * (float)n_seq);

  return energy;
}


PRIVATE int
eval_pt(vrna_fold_compound_t  *fc,
        const short           *pt,
        vrna_cstr_t           output_stream,
        int                   verbosity_level)
{
  int ee, energy;

  if (fc->params->model_details.gquad)
    vrna_message_warning("vrna_eval_*_pt: No gquadruplex support!\n"
                         "Ignoring potential gquads in structure!\n"
                         "Use e.g. vrna_eval_structure() instead!");

  vrna_sc_prepare(fc, VRNA_OPTION_MFE);

  energy = energy_of_extLoop_pt(fc, 0, pt);

  if (verbosity_level > 0) {
    vrna_cstr_print_eval_ext_loop(output_stream,
                                  (fc->type == VRNA_FC_TYPE_COMPARATIVE) ?
                                  (int)energy / (int)fc->n_seq :
                                  energy);
  }

  ee      = energy_of_ext_loop_components(fc, pt, output_stream, verbosity_level);
  energy  = ADD_OR_INF(energy, ee);

  return energy;
}


PRIVATE int
energy_of_extLoop_pt(vrna_fold_compound_t *fc,
                     unsigned int         begin,
                     const short          *pt)
{
  short         *s, *s1, s5, s3, **S, **S5, **S3;
  unsigned int  a, n, u, tt, *so, *sn, *ss, sss, strand, last_strand, i, j, last_i,
                start, n_seq, **a2s, strand_start, strand_end;
  int           energy, dangle_model, bonus, e, e_mm3_occupied, e_mm3_available;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_sc_t     *sc, **scs;

  energy        = 0;
  n             = fc->length;
  so            = fc->strand_order;
  sn            = fc->strand_number;
  ss            = fc->strand_start;
  P             = fc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;

  switch (fc->type) {
    case VRNA_FC_TYPE_COMPARATIVE:
      n_seq = fc->n_seq;
      s     = NULL;
      s1    = NULL;
      S     = fc->S;
      S5    = fc->S5;
      S3    = fc->S3;
      a2s   = fc->a2s;
      sc    = NULL;
      scs   = fc->scs;
      break;

    default:
      n_seq = 1;
      s     = fc->sequence_encoding2;
      s1    = fc->sequence_encoding;
      S     = NULL;
      S5    = NULL;
      S3    = NULL;
      a2s   = NULL;
      sc    = fc->sc;
      scs   = NULL;
      break;
  }

  strand_start = strand_end = 0;

  /*
   *  if begin == 0, or there exists only a single strand,
   *  we evaluiate the entire external loop
   */
  if ((begin == 0) || (fc->strands < 2)) {
    strand_end = fc->strands - 1;
  } else {
    /*
     *  otherwise, we either evaluate an unconnected strand or the
     *  external loop part that connects the current strand (sn[begin])
     */
    for (a = strand_start; a < fc->strands; a++)
      if (so[a] == sn[begin]) {
        strand_start = strand_end = a;
        break;
      }
  }

  /* start scanning for enclosing pair from each strand nick */
  for (a = strand_start; a <= strand_end; a++) {
    e               = 0;
    bonus           = 0;
    e_mm3_available = INF;
    e_mm3_occupied  = 0;

    strand  = last_strand = so[a]; /* strand number in current permutation */
    i       = last_i = start = ss[strand];

    while (i <= n) {
      if (sn[i] != last_strand) {
        /* add energy of previous unpaired region */
        switch (fc->type) {
          case VRNA_FC_TYPE_COMPARATIVE:
            if (scs) {
              for (sss = 0; sss < n_seq; sss++) {
                if (scs[sss]) {
                  if (scs[sss]->energy_up) {
                    u     = a2s[sss][i] - a2s[sss][last_i];
                    bonus += scs[sss]->energy_up[a2s[sss][last_i]][u];
                  }
                }
              }
            }

            break;

          default:
            if (sc)
              if (sc->energy_up)
                bonus += sc->energy_up[last_i][i - last_i];

            break;
        }
        break;
      }

      if ((pt[i] != 0) && (pt[i] > i)) {
        /*
         * pairs down-stream
         * add energy of previous unpaired region
         */
        switch (fc->type) {
          case VRNA_FC_TYPE_COMPARATIVE:
            if (scs) {
              for (sss = 0; sss < n_seq; sss++) {
                if (scs[sss]) {
                  if (scs[sss]->energy_up) {
                    u     = a2s[sss][i] - a2s[sss][last_i];
                    bonus += scs[sss]->energy_up[a2s[sss][last_i]][u];
                  }
                }
              }
            }

            break;

          default:
            if (sc)
              if (sc->energy_up)
                bonus += sc->energy_up[last_i][i - last_i];

            break;
        }

        j = (unsigned int)pt[i];

        /* add energy of branch */
        switch (fc->type) {
          case VRNA_FC_TYPE_COMPARATIVE:
            for (sss = 0; sss < n_seq; sss++) {
              tt = vrna_get_ptype_md(S[sss][i], S[sss][j], md);

              switch (dangle_model) {
                case 0:
                  e += vrna_E_ext_stem(tt, -1, -1, P);
                  break;

                case 2:
                  s5  = ((sn[i - 1] == sn[i]) && (a2s[sss][i] > 1)) ? S5[sss][i] : -1;
                  s3  = ((sn[j] == sn[j + 1]) && (a2s[sss][j] < a2s[sss][n])) ? S3[sss][j] : -1;
                  e   += vrna_E_ext_stem(tt, s5, s3, P);
                  break;

                default:
                  /* odd dangles not implemented yet */
                  break;
              }
            }
            break;

          default:
            tt = vrna_get_ptype_md(s[i], s[j], md);

            switch (dangle_model) {
              case 0:
                e += vrna_E_ext_stem(tt, -1, -1, P);
                break;

              case 2:
                s5  = ((sn[i - 1] == sn[i]) && (i > 1)) ? s1[i - 1] : -1;
                s3  = ((sn[j] == sn[j + 1]) && (j < n)) ? s1[j + 1] : -1;
                e   += vrna_E_ext_stem(tt, s5, s3, P);
                break;

              default:
                s5  = ((sn[i - 1] == sn[i]) && (i > 1) && (!pt[i - 1])) ? s1[i - 1] : -1;
                s3  = ((sn[j] == sn[j + 1]) && (j < n) && (!pt[j + 1])) ? s1[j + 1] : -1;

                if ((last_i + 1 < i) || ((last_i == start) && (last_i < i))) {
                  e_mm3_available = MIN2(e_mm3_available, e_mm3_occupied);
                  e_mm3_occupied  = e_mm3_available;
                }

                e = MIN2(
                  e_mm3_occupied + vrna_E_ext_stem(tt, -1, s3, P),
                  e_mm3_available + vrna_E_ext_stem(tt, s5, s3, P)
                  );
                e_mm3_available = MIN2(
                  e_mm3_occupied + vrna_E_ext_stem(tt, -1, -1, P),
                  e_mm3_available + vrna_E_ext_stem(tt, s5, -1, P)
                  );
                e_mm3_occupied = e;
                break;
            }
            break;
        }

        i           = j;      /* skip the branch we'ce just evaluated */
        last_i      = i + 1;  /* update last unpaired nt */
        last_strand = sn[i];  /* update current strand number */
      } else if (pt[i] != 0) {
        /*
         * found 3' end of enclosing pair
         * add energy of previous unpaired region
         */
        switch (fc->type) {
          case VRNA_FC_TYPE_COMPARATIVE:
            if (scs) {
              for (sss = 0; sss < n_seq; sss++) {
                if (scs[sss]) {
                  if (scs[sss]->energy_up) {
                    u     = a2s[sss][i] - a2s[sss][last_i];
                    bonus += scs[sss]->energy_up[a2s[sss][last_i]][u];
                  }
                }
              }
            }

            break;

          default:
            if (sc)
              if (sc->energy_up)
                bonus += sc->energy_up[last_i][i - last_i];

            break;
        }

        j = i;
        i = (unsigned int)pt[i];

        /* add energy of enclosing base pair */
        switch (fc->type) {
          case VRNA_FC_TYPE_COMPARATIVE:
            for (sss = 0; sss < n_seq; sss++) {
              tt = vrna_get_ptype_md(S[sss][j], S[sss][i], md);

              switch (dangle_model) {
                case 0:
                  e += vrna_E_ext_stem(tt, -1, -1, P);
                  break;

                case 2:
                  s5  = (sn[j - 1] == sn[j]) ? S5[sss][j] : -1;
                  s3  = (sn[i] == sn[i + 1]) ? S3[sss][i] : -1;
                  e   += vrna_E_ext_stem(tt, s5, s3, P);
                  break;

                default:
                  /* odd dangles not implemented yet */
                  break;
              }
            }
            break;

          default:
            tt = vrna_get_ptype_md(s[j], s[i], md);

            switch (dangle_model) {
              case 0:
                e += vrna_E_ext_stem(tt, -1, -1, P);
                break;

              case 2:
                s5  = (sn[j - 1] == sn[j]) ? s1[j - 1] : -1;
                s3  = (sn[i] == sn[i + 1]) ? s1[i + 1] : -1;
                e   += vrna_E_ext_stem(tt, s5, s3, P);
                break;

              default:
                s5  = ((sn[j - 1] == sn[j]) && (!pt[j - 1])) ? s1[j - 1] : -1;
                s3  = ((sn[i] == sn[i + 1]) && (!pt[i + 1])) ? s1[i + 1] : -1;

                if ((last_i + 1 < j) || ((last_i == start) && (last_i < j))) {
                  e_mm3_available = MIN2(e_mm3_available, e_mm3_occupied);
                  e_mm3_occupied  = e_mm3_available;
                }

                e = MIN2(
                  e_mm3_occupied + vrna_E_ext_stem(tt, -1, s3, P),
                  e_mm3_available + vrna_E_ext_stem(tt, s5, s3, P)
                  );
                e_mm3_available = MIN2(
                  e_mm3_occupied + vrna_E_ext_stem(tt, -1, -1, P),
                  e_mm3_available + vrna_E_ext_stem(tt, s5, -1, P)
                  );
                e_mm3_occupied = e;
                break;
            }
            break;
        }

        /* add duplex initiation penalty */
        if (dangle_model % 2) {
          e_mm3_available += P->DuplexInit * n_seq;
          e_mm3_occupied  += P->DuplexInit * n_seq;
        } else {
          e += fc->params->DuplexInit * n_seq;
        }

        /* update index variables */
        last_i      = i + 1;
        last_strand = sn[i];
      } else if (i == n) {
        /*
         * end of outer-most external loop
         * add energy of unpaired region
         */
        switch (fc->type) {
          case VRNA_FC_TYPE_COMPARATIVE:
            if (scs) {
              for (sss = 0; sss < n_seq; sss++) {
                if (scs[sss]) {
                  if (scs[sss]->energy_up) {
                    u     = a2s[sss][i + 1] - a2s[sss][last_i];
                    bonus += scs[sss]->energy_up[a2s[sss][last_i]][u];
                  }
                }
              }
            }

            break;

          default:
            if (sc)
              if (sc->energy_up)
                bonus += sc->energy_up[last_i][i - last_i + 1];

            break;
        }
      }

      i++;
    }

    if (dangle_model % 2)
      e = MIN2(e_mm3_available, e_mm3_occupied);

    e += bonus;

    energy += e;
  }

  return energy;
}


PRIVATE int
energy_of_ext_loop_components(vrna_fold_compound_t  *fc,
                              const short           *pt,
                              vrna_cstr_t           output_stream,
                              int                   verbosity_level)
{
  unsigned int  last_s, s, i, n, a, *so, *sn, *ss;
  int           energy = 0;

  n       = fc->length;
  so      = fc->strand_order;
  sn      = fc->strand_number;
  ss      = fc->strand_start;
  energy  = 0;

  /* start scanning for enclosing pair from each strand start site in 5'->3' direction */
  for (a = 0; a < fc->strands; a++) {
    s = last_s = so[a]; /* strand number in current permutation */
    i = ss[s];

    while (i <= n) {
      if (sn[i] != last_s)
        break;

      if (pt[i] != 0) {
        if (pt[i] > i) {
          /*
           * pairs down-stream
           * add energy of enclosed substem
           */
          energy  += stack_energy(fc, i, pt, output_stream, verbosity_level);
          i       = (unsigned int)pt[i];
          last_s  = sn[i]; /* update current strand number */
        } else {
          /*
           * found 3' end of enclosing pair
           * update index variables
           */
          i       = (unsigned int)pt[i];
          last_s  = sn[i];
        }
      }

      i++;
    }
  }

  return energy;
}


PRIVATE int
eval_circ_pt(vrna_fold_compound_t *fc,
             const short          *pt,
             vrna_cstr_t          output_stream,
             int                  verbosity_level)
{
  unsigned int  s, n_seq, **a2s;
  int           i, j, length, energy, en0, degree;
  vrna_param_t  *P;
  vrna_sc_t     *sc, **scs;

  energy  = 0;
  en0     = 0;
  degree  = 0;
  length  = fc->length;
  P       = fc->params;

  switch (fc->type) {
    case VRNA_FC_TYPE_COMPARATIVE:
      n_seq = fc->n_seq;
      sc    = NULL;
      scs   = fc->scs;
      a2s   = fc->a2s;
      break;

    default:
      n_seq = 1;
      sc    = fc->sc;
      scs   = NULL;
      a2s   = NULL;
      break;
  }

  if (P->model_details.gquad)
    vrna_message_warning("vrna_eval_*_pt: No gquadruplex support!\n"
                         "Ignoring potential gquads in structure!\n"
                         "Use e.g. vrna_eval_structure() instead!");

  vrna_sc_prepare(fc, VRNA_OPTION_MFE);

  /* evaluate all stems in exterior loop */
  for (i = 1; i <= length; i++) {
    if (pt[i] == 0)
      continue;

    degree++;
    energy += stack_energy(fc,
                           i,
                           (const short *)pt,
                           output_stream,
                           verbosity_level);
    i = pt[i];
  }

  /* find first stem */
  for (i = 1; (i <= length) && (!pt[i]); i++);
  j = pt[i];

  /* evaluate exterior loop itself */
  switch (degree) {
    case 0:   /* unstructured */
      switch (fc->type) {
        case VRNA_FC_TYPE_COMPARATIVE:
          if (scs) {
            for (s = 0; s < n_seq; s++)
              if (scs[s] && scs[s]->energy_up)
                en0 += scs[s]->energy_up[1][a2s[s][length]];
          }

          break;

        default:
          if (sc)
            if (sc->energy_up)
              en0 += sc->energy_up[1][length];

          break;
      }
      break;

    case 1:   /* hairpin loop */
      en0 = vrna_eval_ext_hp_loop(fc, i, j);
      break;

    case 2:   /* interior loop */
    {
      int p, q;
      /* seek to next pair */
      for (p = j + 1; pt[p] == 0; p++);
      q = pt[p];

      en0 = eval_ext_int_loop(fc, i, j, p, q);
    }
    break;

    default:  /* multibranch loop */
      en0 = energy_of_ml_pt(fc, 0, (const short *)pt);

      if (fc->type == VRNA_FC_TYPE_SINGLE)
        en0 -= E_MLstem(0, -1, -1, P);         /* remove virtual closing pair */

      break;
  }

  if (verbosity_level > 0)
    vrna_cstr_print_eval_ext_loop(output_stream,
                                  (int)en0 / (int)n_seq);

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
en_corr_of_loop_gquad(vrna_fold_compound_t  *fc,
                      int                   i,
                      int                   j,
                      const char            *structure,
                      const short           *pt,
                      const int             *loop_idx,
                      vrna_cstr_t           output_stream,
                      int                   verbosity_level)
{
  char          *sequence;
  short         *s1, *s2, **S, **S5, **S3;
  unsigned int  cnt, n_seq;
  int           pos, tmp_e, energy, p, q, r, s, u, type, type2,
                L, l[3], num_elem, num_g, elem_i, elem_j,
                up_mis, gq_en[2], dangle_model;
  vrna_param_t  *P;
  vrna_md_t     *md;

  s1            = fc->sequence_encoding;
  s2            = fc->sequence_encoding2;
  S             = fc->S;
  S5            = fc->S5;
  S3            = fc->S3;
  P             = fc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;

  switch (fc->type) {
    case VRNA_FC_TYPE_COMPARATIVE:
      n_seq     = fc->n_seq;
      sequence  = fc->cons_seq;
      break;

    default:
      n_seq     = 1;
      sequence  = fc->sequence;
      break;
  }

  energy  = 0;
  q       = i;

  while ((pos = parse_gquad(structure + q - 1, &L, l)) > 0) {
    q += pos - 1;
    p = q - 4 * L - l[0] - l[1] - l[2] + 1;
    if (q > j)
      break;

    /* we've found the first g-quadruplex at position [p,q] */
    switch (fc->type) {
      case VRNA_FC_TYPE_COMPARATIVE:
        E_gquad_ali_en(p, L, l, (const short **)S, fc->a2s, n_seq, P, gq_en);
        tmp_e = gq_en[0];
        break;

      default:
        tmp_e = E_gquad(L, l, P);
        break;
    }

    energy += tmp_e;

    if (verbosity_level > 0) {
      vrna_cstr_print_eval_gquad(output_stream,
                                 p,
                                 L,
                                 l,
                                 (int)tmp_e / (int)n_seq);
    }

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
          pos = parse_gquad(structure + u - 1, &L, l);
          if (pos > 0) {
            switch (fc->type) {
              case VRNA_FC_TYPE_COMPARATIVE:
                E_gquad_ali_en(u, L, l, (const short **)S, fc->a2s, n_seq, P, gq_en);
                tmp_e = gq_en[0];
                break;

              default:
                tmp_e = E_gquad(L, l, P);
                break;
            }

            if (verbosity_level > 0) {
              vrna_cstr_print_eval_gquad(output_stream,
                                         pos,
                                         L,
                                         l,
                                         (int)tmp_e / (int)n_seq);
            }

            energy  += tmp_e;
            up_mis  += pos;
            u       += pos;
            num_g++;
          }
        } else {
          /* we must have found a stem */
          num_elem++;
          elem_i  = u;
          elem_j  = pt[u];
          energy  += en_corr_of_loop_gquad(fc,
                                           u,
                                           pt[u],
                                           structure,
                                           pt,
                                           loop_idx,
                                           output_stream,
                                           verbosity_level);
          u = pt[u] + 1;
        }
      }

      /* here, u == s */
      int e_minus, e_plus, e_temp;

      e_plus = e_minus = 0;

      /* we are done since we've found no other 3' structure element */
      switch (num_elem) {
        /* g-quad was misinterpreted as hairpin closed by (r,s) */
        case 0:
          e_minus = vrna_eval_hp_loop(fc, r, s);

          if (verbosity_level > 0) {
            vrna_cstr_print_eval_hp_loop_revert(output_stream,
                                                r,
                                                s,
                                                sequence[r - 1],
                                                sequence[s - 1],
                                                (int)e_minus / (int)n_seq);
          }

          /* if we consider the G-Quadruplex, we have */
          if (num_g == 1) {
            /* a) an interior loop like structure */
            switch (fc->type) {
              case VRNA_FC_TYPE_COMPARATIVE:
                for (cnt = 0; cnt < n_seq; cnt++) {
                  type = vrna_get_ptype_md(S[cnt][r], S[cnt][s], md);

                  if (dangle_model == 2)
                    e_plus += P->mismatchI[type][S3[cnt][r]][S5[cnt][s]];

                  if (type > 2)
                    e_plus += P->TerminalAU;
                }
                break;

              default:
                type = md->pair[s2[r]][s2[s]];
                if (dangle_model == 2)
                  e_plus += P->mismatchI[type][s1[r + 1]][s1[s - 1]];

                if (type > 2)
                  e_plus += P->TerminalAU;

                break;
            }

            e_plus += n_seq * P->internal_loop[s - r - 1 - up_mis];

            if (verbosity_level > 0) {
              vrna_cstr_print_eval_int_loop(output_stream,
                                            r,
                                            s,
                                            sequence[r - 1],
                                            sequence[s - 1],
                                            p,
                                            q,
                                            sequence[p - 1],
                                            sequence[q - 1],
                                            (int)e_plus / (int)n_seq);
            }
          } else {
            /* or b) a multibranch loop like structure */
            e_temp = num_g * E_MLstem(0, -1, -1, P) +
                     P->MLclosing +
                     (elem_i - r - 1 + s - elem_j - 1 - up_mis) * P->MLbase;

            e_plus = n_seq * e_temp;

            switch (fc->type) {
              case VRNA_FC_TYPE_COMPARATIVE:
                for (cnt = 0; cnt < n_seq; cnt++) {
                  type    = vrna_get_ptype_md(S[cnt][s], S[cnt][r], md);
                  e_plus  += E_MLstem(type, S5[cnt][s], S3[cnt][r], P);
                }
                break;

              default:
                type    = md->pair[s2[s]][s2[r]];
                e_plus  += E_MLstem(type, s1[s - 1], s1[r + 1], P);
                break;
            }

            if (verbosity_level > 0) {
              vrna_cstr_print_eval_mb_loop(output_stream,
                                           r,
                                           s,
                                           sequence[r - 1],
                                           sequence[s - 1],
                                           (int)e_plus / (int)n_seq);
            }
          }

          energy += e_plus - e_minus;
          break;

        /* g-quad was misinterpreted as interior loop closed by (r,s) with enclosed pair (elem_i, elem_j) */
        case 1:
          e_temp = num_g * E_MLstem(0, -1, -1, P) +
                   P->MLclosing +
                   (elem_i - r - 1 + s - elem_j - 1 - up_mis) * P->MLbase;
          e_plus = n_seq * e_temp;

          switch (fc->type) {
            case VRNA_FC_TYPE_COMPARATIVE:
              for (cnt = 0; cnt < n_seq; cnt++) {
                type    = vrna_get_ptype_md(S[cnt][s], S[cnt][r], md);
                type2   = vrna_get_ptype_md(S[cnt][elem_i], S[cnt][elem_j], md);
                e_plus  += E_MLstem(type, S5[cnt][s], S3[cnt][r], P) +
                           E_MLstem(type, S5[cnt][elem_i], S3[cnt][elem_j], P);
              }
              break;

            default:
              type    = md->pair[s2[s]][s2[r]];
              type2   = md->pair[s2[elem_i]][s2[elem_j]];
              e_plus  += E_MLstem(type, s1[s - 1], s1[r + 1], P) +
                         E_MLstem(type2, s1[elem_i - 1], s1[elem_j + 1], P);
              break;
          }

          e_minus = vrna_eval_int_loop(fc, r, s, elem_i, elem_j);
          energy  += e_plus - e_minus;

          if (verbosity_level > 0) {
            vrna_cstr_print_eval_int_loop_revert(output_stream,
                                                 r,
                                                 s,
                                                 sequence[r - 1],
                                                 sequence[j - 1],
                                                 elem_i,
                                                 elem_j,
                                                 sequence[elem_i - 1],
                                                 sequence[elem_j - 1],
                                                 (int)e_minus / (int)n_seq);

            vrna_cstr_print_eval_mb_loop(output_stream,
                                         r,
                                         s,
                                         sequence[r - 1],
                                         sequence[s - 1],
                                         (int)e_plus / (int)n_seq);
          }

          break;

        /* gquad was misinterpreted as unpaired nucleotides in a multiloop */
        default:
          e_minus = (up_mis) * P->MLbase * n_seq;
          e_plus  = num_g * E_MLstem(0, -1, -1, P) * n_seq;
          energy  += e_plus - e_minus;

          if (verbosity_level > 0) {
            vrna_cstr_print_eval_mb_loop_revert(output_stream,
                                                r,
                                                s,
                                                sequence[r - 1],
                                                sequence[s - 1],
                                                (int)e_minus / (int)n_seq);

            vrna_cstr_print_eval_mb_loop(output_stream,
                                         r,
                                         s,
                                         sequence[r - 1],
                                         sequence[s - 1],
                                         (int)e_plus / (int)n_seq);
          }

          break;
      }

      q = s + 1;
    }
  }

  return energy;
}


/*
 * returns first base pairing partner
 * after the last strand nick in the loop
 * enclosed by (i,j).
 */
PRIVATE int
first_pair_after_last_nick(unsigned int i,
                           unsigned int j,
                           const short  *pt,
                           unsigned int *sn)
{
  unsigned int p, r, first_strand, last_strand;

  first_strand  = sn[i];
  last_strand   = sn[j];
  p             = j;

  if (first_strand != last_strand) {
    /* go backwards from j to first strand nick */
    for (r = j - 1; r > i; r--) {
      if (sn[r] != last_strand)
        break;

      if (pt[r] != 0) {
        /* hop over base pair and store 5' pairing partner */
        r           = p = pt[r];
        last_strand = sn[p];
      }
    }
  }

  return (last_strand == first_strand) ? 0 : p;
}


PRIVATE int
stack_energy(vrna_fold_compound_t *fc,
             int                  i,
             const short          *pt,
             vrna_cstr_t          output_stream,
             int                  verbosity_level)
{
  /* recursively calculate energy of substructure enclosed by (i,j) */
  char          *string;
  short         *s;
  unsigned int  *sn, n_seq;
  int           ee, energy, j, p, q;
  vrna_param_t  *P;
  vrna_md_t     *md;

  sn      = fc->strand_number;
  s       = fc->sequence_encoding2;
  P       = fc->params;
  md      = &(P->model_details);
  energy  = 0;

  j = pt[i];

  switch (fc->type) {
    case VRNA_FC_TYPE_COMPARATIVE:
      string  = fc->cons_seq;
      n_seq   = fc->n_seq;
      break;

    default:
      string  = fc->sequence;
      n_seq   = 1;
      if (md->pair[s[i]][s[j]] == 0) {
        if (verbosity_level > VRNA_VERBOSITY_QUIET) {
          vrna_message_warning("bases %d and %d (%c%c) can't pair!",
                               i, j,
                               string[i - 1],
                               string[j - 1]);
        }
      }

      break;
  }

  p = i;
  q = j;

  while (p < q) {
    /* process all stacks and interior loops */
    while (pt[++p] == 0);
    while (pt[--q] == 0);
    if ((pt[q] != (short)p) || (p > q))
      break;

    if ((sn[i] == sn[p]) &&
        (sn[q] == sn[j])) {
      if (fc->type == VRNA_FC_TYPE_SINGLE) {
        if (md->pair[s[q]][s[p]] == 0) {
          if (verbosity_level > VRNA_VERBOSITY_QUIET) {
            vrna_message_warning("bases %d and %d (%c%c) can't pair!",
                                 p, q,
                                 string[p - 1],
                                 string[q - 1]);
          }
        }
      }

      ee = vrna_eval_int_loop(fc, i, j, p, q);

      if (verbosity_level > 0) {
        vrna_cstr_print_eval_int_loop(output_stream,
                                      i, j,
                                      string[i - 1], string[j - 1],
                                      p, q,
                                      string[p - 1], string[q - 1],
                                      (int)ee / (int)n_seq);
      }

      energy  += ee;
      i       = p;
      j       = q;
    } else {
      return energy;
    }
  } /* end while */

  /* p,q don't pair must have found hairpin or multiloop */

  if (p > q) {
    /* hairpin */
    if (sn[i] == sn[j]) {
      ee = vrna_eval_hp_loop(fc, i, j);
      if (verbosity_level > 0) {
        vrna_cstr_print_eval_hp_loop(output_stream,
                                     i, j,
                                     string[i - 1], string[j - 1],
                                     (int)ee / (int)n_seq);
      }

      energy += ee;
    }

    return energy;
  }

  /* (i,j) is exterior pair of multiloop or external loop */
  if (!first_pair_after_last_nick(i, j, pt, sn)) {
    while (p < j) {
      /* add up the contributions of the substructures of the ML */
      energy  += stack_energy(fc, p, pt, output_stream, verbosity_level);
      p       = pt[p];
      /* search for next base pair in multiloop */
      while (pt[++p] == 0);
    }

    ee = energy_of_ml_pt(fc, i, pt);

    if (verbosity_level > 0) {
      vrna_cstr_print_eval_mb_loop(output_stream,
                                   i, j,
                                   string[i - 1], string[j - 1],
                                   (int)ee / (int)n_seq);
    }

    energy += ee;
  } else {
    return energy;
  }

  return energy;
}


/*---------------------------------------------------------------------------*/


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
energy_of_ml_pt(vrna_fold_compound_t  *fc,
                int                   i,
                const short           *pt)
{
  short         *s, *s1, **S, **S5, **S3;
  unsigned int  *sn, **a2s, n_seq;
  int           energy, cx_energy, tmp, tmp2, best_energy = INF, bonus, *idx, dangle_model,
                logML, circular, *rtype, ss, n, i1, j, p, q, q_prev, q_prev2, u, uu,
                x, type, count, mm5, mm3, tt, ld5, new_cx, dang5, dang3, dang, e_stem,
                e_stem5, e_stem3, e_stem53, mlintern[NBPAIRS + 1];
  int           E_mm5_available;    /* energy of 5' part where 5' mismatch of current stem is available */
  int           E_mm5_occupied;     /* energy of 5' part where 5' mismatch of current stem is unavailable */
  int           E2_mm5_available;   /* energy of 5' part where 5' mismatch of current stem is available with possible 3' dangle for enclosing pair (i,j) */
  int           E2_mm5_occupied;    /* energy of 5' part where 5' mismatch of current stem is unavailable with possible 3' dangle for enclosing pair (i,j) */
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_sc_t     *sc, **scs;

  n   = fc->length;
  sn  = fc->strand_number;
  P   = fc->params;
  md  = &(P->model_details);
  idx = fc->jindx;

  circular      = md->circ;
  dangle_model  = md->dangles;
  logML         = md->logML;
  rtype         = &(md->rtype[0]);

  switch (fc->type) {
    case VRNA_FC_TYPE_COMPARATIVE:
      s     = NULL;
      s1    = NULL;
      sc    = NULL;
      S     = fc->S;
      S5    = fc->S5;
      S3    = fc->S3;
      a2s   = fc->a2s;
      n_seq = fc->n_seq;
      scs   = fc->scs;
      break;

    default:
      s     = fc->sequence_encoding2;
      s1    = fc->sequence_encoding;
      sc    = fc->sc;
      S     = NULL;
      S5    = NULL;
      S3    = NULL;
      a2s   = NULL;
      n_seq = 1;
      scs   = NULL;
      break;
  }

  bonus = 0;

  if (i >= pt[i]) {
    vrna_message_warning("energy_of_ml_pt: i is not 5' base of a closing pair!");
    return INF;
  }

  j = (i == 0) ? n + 1 : (int)pt[i];

  switch (fc->type) {
    case VRNA_FC_TYPE_COMPARATIVE:
      if ((dangle_model % 2) ||
          (dangle_model > 2) ||
          (dangle_model < 0)) {
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
      if (i != 0) {
        /* (i,j) is closing pair of multibranch loop, add soft constraints */
        if (sc)
          if (sc->energy_bp)
            bonus += sc->energy_bp[idx[j] + i];
      }

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
  switch (fc->type) {
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
      u += p - i - 1;
      if (sc)
        if (sc->energy_up)
          bonus += sc->energy_up[i + 1][u];

      break;
  }

  switch (dangle_model) {
    case 0:
      switch (fc->type) {
        case VRNA_FC_TYPE_COMPARATIVE:
          while (p < j) {
            /* p must have a pairing partner */
            q = (int)pt[p];
            for (ss = 0; ss < n_seq; ss++) {
              /* get type of base pair (p,q) */
              tt = vrna_get_ptype_md(S[ss][p], S[ss][q], md);

              energy += E_MLstem(tt, -1, -1, P);
            }

            /* seek to the next stem */
            p       = q + 1;
            q_prev  = q_prev2 = q;
            while (p < j && !pt[p])
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
              tt = vrna_get_ptype_md(S[ss][j], S[ss][i], md);

              energy += E_MLstem(tt, -1, -1, P);
            }
          }

          break;

        default:
          while (p < j) {
            /* p must have a pairing partner */
            q = (int)pt[p];
            /* get type of base pair (p,q) */
            tt = vrna_get_ptype_md(s[p], s[q], md);

            energy += E_MLstem(tt, -1, -1, P);

            /* seek to the next stem */
            p       = q + 1;
            q_prev  = q_prev2 = q;
            while (p < j && !pt[p])
              p++;
            u += p - q - 1;                                     /* add unpaired nucleotides */

            if (sc)
              if (sc->energy_up)
                bonus += sc->energy_up[q + 1][p - q - 1];
          }

          /* now lets get the energy of the enclosing stem */
          if (i > 0) {
            /* actual closing pair */
            tt = vrna_get_ptype_md(s[j], s[i], md);

            energy += E_MLstem(tt, -1, -1, P);
          } else {
            /* virtual closing pair */
            energy += E_MLstem(0, -1, -1, P);
          }

          break;
      }
      break;

    case 2:
      switch (fc->type) {
        case VRNA_FC_TYPE_COMPARATIVE:
          while (p < j) {
            /* p must have a pairing partner */
            q = (int)pt[p];

            for (ss = 0; ss < n_seq; ss++) {
              /* get type of base pair (p,q) */
              tt = vrna_get_ptype_md(S[ss][p], S[ss][q], md);

              mm5     = ((a2s[ss][p] > 1) || circular) ? S5[ss][p] : -1;
              mm3     = ((a2s[ss][q] < a2s[ss][n]) || circular) ? S3[ss][q] : -1;
              energy  += E_MLstem(tt, mm5, mm3, P);
            }

            /* seek to the next stem */
            p       = q + 1;
            q_prev  = q_prev2 = q;
            while (p < j && !pt[p])
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
              tt = vrna_get_ptype_md(S[ss][j], S[ss][i], md);

              mm5     = S5[ss][j];
              mm3     = S3[ss][i];
              energy  += E_MLstem(tt, mm5, mm3, P);
            }
          }

          break;

        default:
          while (p < j) {
            /* p must have a pairing partner */
            q = (int)pt[p];
            /* get type of base pair (p,q) */
            tt = vrna_get_ptype_md(s[p], s[q], md);

            mm5     = sn[p - 1] == sn[p] ? s1[p - 1] : -1;
            mm3     = sn[q] == sn[q + 1] ? s1[q + 1] : -1;
            energy  += E_MLstem(tt, mm5, mm3, P);

            /* seek to the next stem */
            p       = q + 1;
            q_prev  = q_prev2 = q;
            while (p < j && !pt[p])
              p++;
            u += p - q - 1;                                     /* add unpaired nucleotides */

            if (sc)
              if (sc->energy_up)
                bonus += sc->energy_up[q + 1][p - q - 1];
          }
          if (i > 0) {
            /* actual closing pair */
            tt = vrna_get_ptype_md(s[j], s[i], md);

            mm5     = sn[j - 1] == sn[j] ? s1[j - 1] : -1;
            mm3     = sn[i] == sn[i + 1] ? s1[i + 1] : -1;
            energy  += E_MLstem(tt, mm5, mm3, P);
          } else {
            /* virtual closing pair */
            energy += E_MLstem(0, -1, -1, P);
          }

          break;
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
          type  = vrna_get_ptype_md(s[j], s[i], md);

          /* prime the ld5 variable */
          if (sn[j - 1] == sn[j]) {
            ld5 = P->dangle5[type][s1[j - 1]];
            if ((p = (unsigned int)pt[j - 2]) && (sn[j - 2] == sn[j - 1]))
              if (P->dangle3[md->pair[s[p]][s[j - 2]]][s1[j - 1]] < ld5)
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
            tt = vrna_get_ptype_md(s[p], s[q], md);
          }

          energy    += mlintern[tt];
          cx_energy += mlintern[tt];

          dang5 = dang3 = 0;
          if ((sn[p - 1] == sn[p]) && (p > 1))
            dang5 = P->dangle5[tt][s1[p - 1]];          /* 5'dangle of pq pair */

          if ((sn[i1] == sn[i1 + 1]) && (i1 < (unsigned int)s[0]))
            dang3 = P->dangle3[type][s1[i1 + 1]];       /* 3'dangle of previous pair */

          switch (p - i1 - 1) {
            case 0:           /* adjacent helices */
              if (i1 != 0) {
                if (sn[i1] == sn[p]) {
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
        tt = vrna_get_ptype_md(s[p], s[q], md);

        if (q_prev + 2 < p) {
          E_mm5_available = MIN2(E_mm5_available, E_mm5_occupied);
          E_mm5_occupied  = E_mm5_available;
        }

        if (q_prev2 + 2 < p) {
          E2_mm5_available  = MIN2(E2_mm5_available, E2_mm5_occupied);
          E2_mm5_occupied   = E2_mm5_available;
        }

        mm5       = ((sn[p - 1] == sn[p]) && !pt[p - 1])  ? s1[p - 1] : -1;
        mm3       = ((sn[q] == sn[q + 1]) && !pt[q + 1])  ? s1[q + 1] : -1;
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
        while (p < j && !pt[p])
          p++;
        u += p - q - 1;         /* add unpaired nucleotides */

        if (sc)
          if (sc->energy_up)
            bonus += sc->energy_up[q + 1][p - q - 1];
      }
      if (i > 0) {
        /* actual closing pair */
        type = vrna_get_ptype_md(s[j], s[i], md);

        mm5 = ((sn[j - 1] == sn[j]) && !pt[j - 1])  ? s1[j - 1] : -1;
        mm3 = ((sn[i] == sn[i + 1]) && !pt[i + 1])  ? s1[i + 1] : -1;
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

  energy += P->MLclosing * n_seq;

  /*
   * logarithmic ML loop energy if logML
   * does this work for comparative predictions as well?
   */
  if (logML && (u > 6))
    energy += 6 * P->MLbase + (int)(P->lxc * log((double)u / 6.));
  else
    energy += (u * P->MLbase);

  return energy + bonus;
}
