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
#include "ViennaRNA/intern/unistd_win.h"
#endif
#else
#include <unistd.h>
#endif

#include <limits.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/structures/pairtable.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/eval/hairpin.h"
#include "ViennaRNA/eval/internal.h"
#include "ViennaRNA/eval/multibranch.h"
#include "ViennaRNA/eval/gquad.h"
#include "ViennaRNA/cofold.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/datastructures/char_stream.h"
#include "ViennaRNA/datastructures/array.h"
#include "ViennaRNA/eval/structures.h"

#include "ViennaRNA/intern/color_output.h"

#define   ADD_OR_INF(a, b)     (((a) != INF) && ((b) != INF) ?  (a) + (b) : INF)

typedef enum {
  VRNA_STRUCTURE_ELEM_EXT_LOOP  = 0,
  VRNA_STRUCTURE_ELEM_HP_LOOP   = 1,
  VRNA_STRUCTURE_ELEM_INT_LOOP  = 2,
  VRNA_STRUCTURE_ELEM_MB_LOOP   = 3,
  VRNA_STRUCTURE_ELEM_GQUAD     = 4,
  VRNA_STRUCTURE_ELEM_UD        = 5,
  VRNA_STRUCTURE_ELEM_STACK     = 6,
  VRNA_STRUCTURE_ELEM_BULGE     = 7
} vrna_struct_elem_type_e;


typedef struct {
  vrna_struct_elem_type_e type;
  unsigned int            i;
  unsigned int            j;
  unsigned int            p;
  unsigned int            q;
  unsigned int            L;
  unsigned int            l[3];
  int                     energy;
} vrna_struct_elem_t;


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
stack_energy(vrna_fold_compound_t           *fc,
             int                            i,
             const short                    *pt,
             vrna_array(vrna_struct_elem_t) *elements);


PRIVATE int
energy_of_ml_pt(vrna_fold_compound_t  *fc,
                unsigned int          i,
                const short           *pt);


PRIVATE int
eval_pt(vrna_fold_compound_t            *fc,
        const short                     *pt,
        vrna_array(vrna_struct_elem_t)  *elements);


PRIVATE int
eval_circ_pt(vrna_fold_compound_t           *fc,
             const short                    *pt,
             vrna_array(vrna_struct_elem_t) *elements);


PRIVATE int
en_corr_of_loop_gquad(vrna_fold_compound_t            *fc,
                      int                             i,
                      int                             j,
                      const char                      *structure,
                      const short                     *pt,
                      const int                       *loop_idx,
                      vrna_array(vrna_struct_elem_t)  *elements,
                      vrna_array(vrna_struct_elem_t)  *elements_rev);


PRIVATE int
en_corr_of_loop_gquad_circ(vrna_fold_compound_t           *fc,
                           unsigned int                   i,
                           unsigned int                   j,
                           const char                     *structure,
                           const short                    *pt,
                           const int                      *loop_idx,
                           vrna_array(vrna_struct_elem_t) *elements,
                           vrna_array(vrna_struct_elem_t) *elements_rev);


PRIVATE float
wrap_eval_structure(vrna_fold_compound_t            *fc,
                    const char                      *structure,
                    const short                     *pt,
                    vrna_array(vrna_struct_elem_t)  *elements);


PRIVATE int
energy_of_extLoop_pt(vrna_fold_compound_t *fc,
                     unsigned int         begin,
                     const short          *pt);


PRIVATE int
energy_of_ext_loop_components(vrna_fold_compound_t            *fc,
                              const short                     *pt,
                              vrna_array(vrna_struct_elem_t)  *elements);


PRIVATE int
first_pair_after_last_nick(unsigned int i,
                           unsigned int j,
                           const short  *pt,
                           unsigned int *sn);


PRIVATE void
print_structure_elements(const char                     *s,
                         vrna_array(vrna_struct_elem_t) elements,
                         vrna_cstr_t                    output_stream);


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
  char        *sequence;
  short       *pt;
  float       en;
  vrna_cstr_t output_stream;

  vrna_array(vrna_struct_elem_t)  elements = NULL;

  en = (float)INF / 100.;

  if ((fc) &&
      (structure)) {
    if (strlen(structure) != fc->length) {
      vrna_log_warning("Sequence and structure have unequal length (%d vs. %d)",
                       fc->length,
                       strlen(structure));
      return en;
    }

    if (verbosity_level > 0)
      vrna_array_init(elements);

    output_stream = vrna_cstr(fc->length, (file) ? file : stdout);
    pt            = vrna_ptable(structure);
    en            = wrap_eval_structure(fc,
                                        structure,
                                        pt,
                                        &elements);
    if (fc->type == VRNA_FC_TYPE_COMPARATIVE)
      sequence = fc->cons_seq;
    else
      sequence = fc->sequence;

    print_structure_elements(sequence,
                             elements,
                             output_stream);

    vrna_cstr_fflush(output_stream);
    vrna_cstr_free(output_stream);

    vrna_array_free(elements);
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
  char  *sequence;
  short *pt;
  float en;

  vrna_array(vrna_struct_elem_t) elements = NULL;

  en = (float)INF / 100.;

  if ((fc) &&
      (structure)) {
    if (strlen(structure) != fc->length) {
      vrna_log_warning("Sequence and structure have unequal length (%d vs. %d)",
                       fc->length,
                       strlen(structure));
      return en;
    }

    if (verbosity_level > 0)
      vrna_array_init(elements);

    pt  = vrna_ptable(structure);
    en  = wrap_eval_structure(fc,
                              structure,
                              pt,
                              &elements);

    if (fc->type == VRNA_FC_TYPE_COMPARATIVE)
      sequence = fc->cons_seq;
    else
      sequence = fc->sequence;

    print_structure_elements(sequence,
                             elements,
                             output_stream);


    vrna_array_free(elements);
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

  vrna_array(vrna_struct_elem_t) elements = NULL;


  e = INF;

  if ((fc) &&
      (pt)) {
    if (pt[0] != (short)fc->length) {
      vrna_log_warning("Sequence and structure have unequal length (%d vs. %d)",
                       fc->length,
                       pt[0]);
      return INF;
    }

    if (verbosity_level > 0)
      vrna_array_init(elements);

    output_stream = vrna_cstr(fc->length, (file) ? file : stdout);
    e             = eval_pt(fc,
                            pt,
                            &elements);

    vrna_cstr_fflush(output_stream);
    vrna_cstr_free(output_stream);
    vrna_array_free(elements);
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
      vrna_log_warning("i = %d is unpaired in loop_energy()", i);
      return INF;
    }

    if (md->pair[s[i]][s[j]] == 0) {
      if (verbosity_level > VRNA_VERBOSITY_QUIET) {
        vrna_log_warning("Bases %d and %d (%c%c) can't pair!",
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
      energy = vrna_eval_hairpin(fc, i, j, VRNA_EVAL_LOOP_NO_HC);
      if (energy == INF) {
        if (j - i - 1 < md->min_loop_size) {
          vrna_log_warning("Hairpin loop closed by %d and %d (%c,%c) too short",
                           i, j,
                           vrna_nucleotide_decode(s[i], md),
                           vrna_nucleotide_decode(s[j], md));
        } else {
          vrna_log_warning("Hairpin loop closed by %d and %d (%c,%c) forbidden",
                           i, j,
                           vrna_nucleotide_decode(s[i], md),
                           vrna_nucleotide_decode(s[j], md));
        }
      }
    } else if (pt[q] != (short)p) {
      /* multi-loop */
      energy = energy_of_ml_pt(fc, i, (const short *)pt);
    } else {
      /* found internal loop */
      if (md->pair[s[q]][s[p]] == 0) {
        if (verbosity_level > VRNA_VERBOSITY_QUIET) {
          vrna_log_warning("Bases %d and %d (%c%c) can't pair!",
                           p, q,
                           vrna_nucleotide_decode(s[p], md),
                           vrna_nucleotide_decode(s[q], md));
        }
      }

      energy = vrna_eval_internal(fc, i, j, p, q, VRNA_EVAL_LOOP_NO_HC);
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
        vrna_log_warning("illegal move or broken pair table\n"
                         "%d %d %d %d ",
                         m1, m2, j, pt[j]);
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
PRIVATE float
wrap_eval_structure(vrna_fold_compound_t            *fc,
                    const char                      *structure,
                    const short                     *pt,
                    vrna_array(vrna_struct_elem_t)  *elements)
{
  unsigned int  n_seq, L, l[3], gq;
  int           res, e, *loop_idx;
  float         energy;
  vrna_md_t     *md;

  energy  = (float)INF / 100.;
  n_seq   = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  md      = &(fc->params->model_details);
  gq      = md->gquad;
  /* temporarily turn-off gquad support to suppress warnings */
  md->gquad = 0;

  if (md->circ)
    res = eval_circ_pt(fc, pt, elements);
  else
    res = eval_pt(fc, pt, elements);

  /* re-set gquad support to previous state */
  md->gquad = gq;

  if ((gq) &&
      (vrna_gq_parse(structure, &L, l) > 0)) {
    vrna_array(vrna_struct_elem_t)  elements_rev = NULL;

    if (*elements)
      vrna_array_init(elements_rev);

    loop_idx  = vrna_loopidx_from_ptable(pt);
    e         = en_corr_of_loop_gquad(fc,
                                      1,
                                      fc->length,
                                      structure,
                                      pt,
                                      (const int *)loop_idx,
                                      elements,
                                      &elements_rev);
    res = ADD_OR_INF(res, e);

    /* revert loop energies for loops with gquads */
    if (elements_rev) {
      for (size_t i = 0; i < vrna_array_size(elements_rev); i++)
        /* remove reverted element in 'elements' */
        for (size_t j = 0; j < vrna_array_size(*elements); j++) {
          if (memcmp((*elements) + j, elements_rev + i, sizeof(vrna_struct_elem_t)) == 0) {
            if (j + 1 < vrna_array_size(*elements))
              memmove((*elements) + j, (*elements) + j + 1,
                      sizeof(vrna_struct_elem_t) * (vrna_array_size(*elements) - j));

            vrna_array_size(*elements) = vrna_array_size(*elements) - 1;
            break;
          }
        }
    }

    vrna_array_free(elements_rev);
    free(loop_idx);
  }

  energy = (float)res / (100. * (float)n_seq);

  return energy;
}


PRIVATE int
eval_pt(vrna_fold_compound_t            *fc,
        const short                     *pt,
        vrna_array(vrna_struct_elem_t)  *elements)
{
  int ee, energy;

  if (fc->params->model_details.gquad)
    vrna_log_warning("Missing G-Quadruplex support!\n"
                     "Ignoring potential gquads in structure!\n"
                     "Use e.g. vrna_eval_structure() instead!");

  vrna_sc_prepare(fc, VRNA_OPTION_MFE);

  energy = energy_of_extLoop_pt(fc, 0, pt);

  if (*elements)
    vrna_array_append(*elements,
                      ((vrna_struct_elem_t){
                        .type = VRNA_STRUCTURE_ELEM_EXT_LOOP,
                        .energy = (fc->type == VRNA_FC_TYPE_COMPARATIVE) ?
                                  (int)energy / (int)fc->n_seq :
                                  energy
                      }));

  ee      = energy_of_ext_loop_components(fc, pt, elements);
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

      if ((pt[i] != 0) && ((unsigned int)pt[i] > i)) {
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
                  e += vrna_E_exterior_stem(tt, -1, -1, P);
                  break;

                case 2:
                  s5  = ((sn[i - 1] == sn[i]) && (a2s[sss][i] > 1)) ? S5[sss][i] : -1;
                  s3  = ((sn[j] == sn[j + 1]) && (a2s[sss][j] < a2s[sss][n])) ? S3[sss][j] : -1;
                  e   += vrna_E_exterior_stem(tt, s5, s3, P);
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
                e += vrna_E_exterior_stem(tt, -1, -1, P);
                break;

              case 2:
                s5  = ((sn[i - 1] == sn[i]) && (i > 1)) ? s1[i - 1] : -1;
                s3  = ((sn[j] == sn[j + 1]) && (j < n)) ? s1[j + 1] : -1;
                e   += vrna_E_exterior_stem(tt, s5, s3, P);
                break;

              default:
                s5  = ((sn[i - 1] == sn[i]) && (i > 1) && (!pt[i - 1])) ? s1[i - 1] : -1;
                s3  = ((sn[j] == sn[j + 1]) && (j < n) && (!pt[j + 1])) ? s1[j + 1] : -1;

                if ((last_i + 1 < i) || ((last_i == start) && (last_i < i))) {
                  e_mm3_available = MIN2(e_mm3_available, e_mm3_occupied);
                  e_mm3_occupied  = e_mm3_available;
                }

                e = MIN2(
                  e_mm3_occupied + vrna_E_exterior_stem(tt, -1, s3, P),
                  e_mm3_available + vrna_E_exterior_stem(tt, s5, s3, P)
                  );
                e_mm3_available = MIN2(
                  e_mm3_occupied + vrna_E_exterior_stem(tt, -1, -1, P),
                  e_mm3_available + vrna_E_exterior_stem(tt, s5, -1, P)
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
                  e += vrna_E_exterior_stem(tt, -1, -1, P);
                  break;

                case 2:
                  s5  = (sn[j - 1] == sn[j]) ? S5[sss][j] : -1;
                  s3  = (sn[i] == sn[i + 1]) ? S3[sss][i] : -1;
                  e   += vrna_E_exterior_stem(tt, s5, s3, P);
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
                e += vrna_E_exterior_stem(tt, -1, -1, P);
                break;

              case 2:
                s5  = (sn[j - 1] == sn[j]) ? s1[j - 1] : -1;
                s3  = (sn[i] == sn[i + 1]) ? s1[i + 1] : -1;
                e   += vrna_E_exterior_stem(tt, s5, s3, P);
                break;

              default:
                s5  = ((sn[j - 1] == sn[j]) && (!pt[j - 1])) ? s1[j - 1] : -1;
                s3  = ((sn[i] == sn[i + 1]) && (!pt[i + 1])) ? s1[i + 1] : -1;

                if ((last_i + 1 < j) || ((last_i == start) && (last_i < j))) {
                  e_mm3_available = MIN2(e_mm3_available, e_mm3_occupied);
                  e_mm3_occupied  = e_mm3_available;
                }

                e = MIN2(
                  e_mm3_occupied + vrna_E_exterior_stem(tt, -1, s3, P),
                  e_mm3_available + vrna_E_exterior_stem(tt, s5, s3, P)
                  );
                e_mm3_available = MIN2(
                  e_mm3_occupied + vrna_E_exterior_stem(tt, -1, -1, P),
                  e_mm3_available + vrna_E_exterior_stem(tt, s5, -1, P)
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
energy_of_ext_loop_components(vrna_fold_compound_t            *fc,
                              const short                     *pt,
                              vrna_array(vrna_struct_elem_t)  *elements)
{
  unsigned int  last_s, s, i, n, a, *so, *sn, *ss;
  int           e, energy = 0;

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
        if ((unsigned int)pt[i] > i) {
          /*
           * pairs down-stream
           * add energy of enclosed substem
           */
          e       = stack_energy(fc, i, pt, elements);
          energy  = ADD_OR_INF(energy, e);
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
eval_circ_pt(vrna_fold_compound_t           *fc,
             const short                    *pt,
             vrna_array(vrna_struct_elem_t) *elements)
{
  unsigned int  s, n_seq, **a2s;
  int           i, j, length, energy, en0, degree;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_sc_t     *sc, **scs;

  energy  = 0;
  en0     = 0;
  degree  = 0;
  length  = fc->length;
  P       = fc->params;
  md      = &(P->model_details);

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

  if (md->gquad)
    vrna_log_warning("Missing G-Quadruplex support!\n"
                     "Ignoring potential gquads in structure!\n"
                     "Use e.g. vrna_eval_structure() instead!");

  vrna_sc_prepare(fc, VRNA_OPTION_MFE);

  /* evaluate all stems in exterior loop */
  for (i = 1; i <= length; i++) {
    if (pt[i] == 0)
      continue;

    degree++;
    en0 = stack_energy(fc,
                       i,
                       (const short *)pt,
                       elements);

    energy  = ADD_OR_INF(energy, en0);
    i       = pt[i];
  }

  /* find first stem */
  for (i = 1; (i <= length) && (!pt[i]); i++);
  j = pt[i];

  /* evaluate exterior loop itself */
  switch (degree) {
    case 0:   /* unstructured */
      switch (fc->type) {
        case VRNA_FC_TYPE_COMPARATIVE:
          for (s = 0; s < n_seq; s++)
            en0 += vrna_E_exterior_loop(a2s[s][length], md);
          break;

        default:
          en0 += vrna_E_exterior_loop(length, md);
          break;
      }

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
      en0 = vrna_eval_hairpin(fc, j, i, VRNA_EVAL_LOOP_NO_HC);
      break;

    case 2:   /* internal loop */
    {
      int p, q;
      /* seek to next pair */
      for (p = j + 1; pt[p] == 0; p++);
      q = pt[p];

      en0 = vrna_eval_internal(fc, i, j, p, q, VRNA_EVAL_LOOP_NO_HC);

    }
    break;

    default:  /* multibranch loop */
      en0 = energy_of_ml_pt(fc, 0, (const short *)pt);

      if (fc->type == VRNA_FC_TYPE_SINGLE)
        en0 -= vrna_E_multibranch_stem(0, -1, -1, P);         /* remove virtual closing pair */

      break;
  }

  if (*elements)
    vrna_array_append(*elements,
                      ((vrna_struct_elem_t){
                        .type = VRNA_STRUCTURE_ELEM_EXT_LOOP,
                        .energy = (int)en0 / (int)n_seq
                      }));

  energy = ADD_OR_INF(energy, en0);

  return energy;
}


PRIVATE int
en_corr_of_loop_gquad_circ(vrna_fold_compound_t           *fc,
                           unsigned int                   i,
                           unsigned int                   j,
                           const char                     *structure,
                           const short                    *pt,
                           const int                      *loop_idx,
                           vrna_array(vrna_struct_elem_t) *elements,
                           vrna_array(vrna_struct_elem_t) *elements_rev)
{
  short         *S, *S1, **SS, **S5, **S3;
  int           corr_en, tmp_e, e_minus, e_plus, gq_en[2];
  unsigned int  n_seq, n, p, num_elem, num_gq, up, up_mis, elem_i,
                elem_j, elem_p, elem_q, gq_p, gq_q, pos, u1, u2, u3,
                us1, us2, us3, type, type2, s, **a2s, L, l[3];
  vrna_param_t  *P;
  vrna_md_t     *md;

  /*
   *  we encountered a gquad [i,j] that is located in the exterior loop
   *  of a circRNA. Since it has not been present in the pair table used
   *  for energy evaluation before, we need to correct the exterior loop
   *  type and corresponding energy here.
   *
   *  Note, that the energy for gquad [i, j] has already been added in
   *  the calling function! So, we only need to add energies for further
   *  gquads in the exterior loop, as well as the loop-type dependent
   *  contributions here.
   */

  n     = fc->length;
  n_seq = (fc->type == VRNA_FC_TYPE_COMPARATIVE) ? fc->n_seq : 1;
  P     = fc->params;
  md    = &(P->model_details);
  S     = (fc->type == VRNA_FC_TYPE_COMPARATIVE) ? NULL : fc->sequence_encoding2;
  S1    = (fc->type == VRNA_FC_TYPE_COMPARATIVE) ? NULL : fc->sequence_encoding;
  SS    = (fc->type == VRNA_FC_TYPE_COMPARATIVE) ? fc->S : NULL;
  S5    = (fc->type == VRNA_FC_TYPE_COMPARATIVE) ? fc->S5 : NULL;
  S3    = (fc->type == VRNA_FC_TYPE_COMPARATIVE) ? fc->S3 : NULL;
  a2s   = (fc->type == VRNA_FC_TYPE_COMPARATIVE) ? fc->a2s : NULL;

  num_elem  = 0;                                        /* number of elements in exterior loop */
  up        = 0;                                        /* number of unpaired bases in exterior loop */
  num_gq    = 1;                                        /* total number of gquadruplexes in exterior loop */
  up_mis    = (i < j) ? j - i + 1 : (n - i + 1) + (j);  /* number of bases misinerpreted as unpaired */
  corr_en   = 0;
  e_minus   = 0;
  e_plus    = 0;

  elem_i = elem_j = elem_p = elem_q = gq_p = gq_q = 0;

  /*
   *  1. Let's determine, how many elements the exterior loop contains:
   */
  if (i < j) {
    /* 1.1 G-Quadruplex not spanning the n,1 junction */
    /*
     *  1.1.1 Seek to position 1. Note, that [i,j] already is the left-most
     *        gquad, so we only need to check for base pairs 5' of it...
     */
    for (p = i - 1; p > 0; p--) {
      if (pt[p] == 0) {
        up++;
        continue;
      } else if (p < (unsigned int)pt[p]) {
        vrna_log_error("Found base pair (%d,%d) enclosing gquad [%d,%d]",
                       p, pt[p], i, j);
        return INF;
      } else {
        tmp_e = en_corr_of_loop_gquad(fc,
                                      pt[p],
                                      p,
                                      structure,
                                      pt,
                                      loop_idx,
                                      elements,
                                      elements_rev);

        if (tmp_e == INF)
          return INF;

        corr_en += tmp_e;

        /* store base pair positions for any corrections we apply later on */
        if (num_elem == 0) {
          elem_i  = pt[p];
          elem_j  = p;
        } else if (num_elem == 1) {
          elem_p  = pt[p];
          elem_q  = p;
        }

        num_elem++;
        p = pt[p];
      }
    }

    /*
     *  1.1.2 Seek to position n. Here, we need to check for more gquads
     *        in the external loop, as well as for any stems that may need
     *        corrections as well.
     */
    for (p = j + 1; p <= n; p++) {
      if (p < (unsigned int)pt[p]) {
        /* base pair */
        tmp_e = en_corr_of_loop_gquad(fc,
                                      p,
                                      pt[p],
                                      structure,
                                      pt,
                                      loop_idx,
                                      elements,
                                      elements_rev);
        if (tmp_e == INF)
          return INF;

        corr_en += tmp_e;

        /* store base pair positions for any corrections we apply later on */
        if (num_elem == 0) {
          elem_i  = p;
          elem_j  = pt[p];
        } else if (num_elem == 1) {
          elem_p  = p;
          elem_q  = pt[p];
        }

        p = pt[p];
        num_elem++;
      } else if (structure[p - 1] == '.') {
        up++;
        continue;
      } else if (structure[p - 1] == '+') {
        /* found another gquad */
        pos = vrna_gq_parse(structure + p - 1, &L, l);

        if (pos > 0) {
          /* sanity check for malformed gquad (this one must not span the n,1 junction! */
          if (4 * L + l[0] + l[1] + l[2] > pos) {
            vrna_log_error("Malformed gquadruplex somewhere after position %u",
                           p - 1);
            return INF;
          }

          switch (fc->type) {
            case VRNA_FC_TYPE_COMPARATIVE:
              vrna_E_consensus_gquad(L,
                                     l,
                                     (unsigned int)p,
                                     n,
                                     n_seq,
                                     (const short **)SS,
                                     (const unsigned int **)a2s,
                                     P,
                                     gq_en);
              tmp_e = gq_en[0];
              break;

            default:
              tmp_e = vrna_E_gquad(L, l, P);
              break;
          }

          if (*elements)
            vrna_array_append(*elements,
                              ((vrna_struct_elem_t){
                                .type = VRNA_STRUCTURE_ELEM_GQUAD,
                                .i = p,
                                .j = p + pos - 1,
                                .L = L,
                                .l = { l[0], l[1], l[2] },
                                .energy = (int)tmp_e / (int)n_seq
                              }));

          gq_p    = p;
          gq_q    = p + pos - 1;
          corr_en += tmp_e;
          up_mis  += pos;
          p       += pos;
          num_gq++;
        }
      } else if (pt[p] != 0) {
        vrna_log_error("Found base pair (%d,%d) enclosing gquad [%d,%d]",
                       pt[p], p, i, j);
        return INF;
      }
    }
  } else {
    /* 1.2 G-Quadruplex spanning the n,1 junction */
    /*
     *  1.2.1 Seek to position i - 1. Here, we need to check for more gquads
     *        in the external loop, as well as for any stems that may need
     *        corrections as well.
     */
    for (p = j + 1; p < i; p++) {
      if (p < (unsigned int)pt[p]) {
        /* base pair */
        tmp_e = en_corr_of_loop_gquad(fc,
                                      p,
                                      pt[p],
                                      structure,
                                      pt,
                                      loop_idx,
                                      elements,
                                      elements_rev);
        if (tmp_e == INF)
          return INF;

        corr_en += tmp_e;

        /* store base pair positions for any corrections we apply later on */
        if (num_elem == 0) {
          elem_i  = p;
          elem_j  = pt[p];
        } else if (num_elem == 1) {
          elem_p  = p;
          elem_q  = pt[p];
        }

        p = pt[p];
        num_elem++;
      } else if (structure[p - 1] == '.') {
        up++;
        continue;
      } else if (structure[p - 1] == '+') {
        /* found another gquad */
        pos = vrna_gq_parse(structure + p - 1, &L, l);
        if (pos > 0) {
          /* sanity check for malformed gquad (this one must not span the n,1 junction! */
          if (4 * L + l[0] + l[1] + l[2] > pos) {
            vrna_log_error("Malformed gquadruplex somewhere after position %u",
                           p - 1);
            return INF;
          }

          switch (fc->type) {
            case VRNA_FC_TYPE_COMPARATIVE:
              vrna_E_consensus_gquad(L,
                                     l,
                                     (unsigned int)p,
                                     n,
                                     n_seq,
                                     (const short **)SS,
                                     (const unsigned int **)a2s,
                                     P,
                                     gq_en);
              tmp_e = gq_en[0];
              break;

            default:
              tmp_e = vrna_E_gquad(L, l, P);
              break;
          }

          if (*elements)
            vrna_array_append(*elements,
                              ((vrna_struct_elem_t){
                                .type = VRNA_STRUCTURE_ELEM_GQUAD,
                                .i = p,
                                .j = p + pos - 1,
                                .L = L,
                                .l = { l[0], l[1], l[2] },
                                .energy = (int)tmp_e / (int)n_seq
                              }));

          gq_p    = p;
          gq_q    = p + pos - 1;
          corr_en += tmp_e;
          up_mis  += pos;
          p       += pos;
          num_gq++;
        }
      } else if (pt[p] != 0) {
        vrna_log_error("Found base pair (%d,%d) enclosing gquad [%d,%d]",
                       pt[p], p, i, j);
        return INF;
      }
    }
  }

  /* now for the actual correction (if necessary) */
  switch (num_elem) {
    case 0: /* exterior loop has been falsly assumed to be entirely unpaired */
      e_minus += vrna_E_exterior_loop(n, md) * (int)n_seq;

      if (*elements_rev)
        vrna_array_append(*elements_rev,
                          ((vrna_struct_elem_t){
                            .type = VRNA_STRUCTURE_ELEM_EXT_LOOP,
                            .energy = (int)e_minus / (int)n_seq
                          }));

      tmp_e = 0;

      switch (num_gq) {
        case 1: /* actually a hairpin-like gquad structure */
          u1 = (i < j) ? (i - 1) + (n - j) : (i - j - 1);

          if (u1 < 3)
            return INF;

          tmp_e += vrna_E_hairpin(n - up_mis, 0, -1, -1, NULL, P) * (int)n_seq;

          break;

        case 2: /* actually an internal-loop-like structure */
          u1  = (i < j) ? i - 1 + n - gq_q : gq_p - j - 1;
          u2  = (i < j) ? gq_p - j - 1 : i - gq_q - 1;
          if (((u1 == 0) && (u2 < 3)) ||
              ((u1 < 3) && (u2 == 0)))
            return INF;

          tmp_e += P->internal_loop[n - up_mis];
          break;

        default:  /* actually a multibranch-loop-like structure */
          tmp_e += (P->MLclosing +
                    (up * P->MLbase) +
                    (num_gq * vrna_E_multibranch_stem(0, -1, -1, P))) *
                   (int)n_seq;
          break;
      }


      if (*elements)
        vrna_array_append(*elements,
                          ((vrna_struct_elem_t){
                            .type = VRNA_STRUCTURE_ELEM_EXT_LOOP,
                            .energy = (int)tmp_e / (int)n_seq
                          }));

      e_plus += tmp_e;

      break;

    case 1: /* exterior loop has been falsly assumed to be hairpin-like */
      e_minus += vrna_eval_hairpin(fc, elem_i, elem_j, VRNA_EVAL_LOOP_NO_HC);

      if (*elements_rev)
        vrna_array_append(*elements_rev,
                          ((vrna_struct_elem_t){
                            .type = VRNA_STRUCTURE_ELEM_EXT_LOOP,
                            .energy = (int)e_minus / (int)n_seq
                          }));

      tmp_e = 0;

      switch (num_gq) {
        case 1: /* actually an internal-loop-like structure */
          switch (fc->type) {
            case VRNA_FC_TYPE_SINGLE:
              type = vrna_get_ptype_md(S[elem_j], S[elem_i], md);
              if (dangles == 2)
                tmp_e += P->mismatchI[type][S1[elem_j + 1]][S1[elem_i - 1]];

              if (type > 2)
                tmp_e += P->TerminalAU;

              if (i < j) {
                if (i > elem_j) {
                  u1  = elem_i - 1;
                  u2  = i - elem_j - 1;
                  u3  = n - j;
                } else {
                  u1  = i - 1;
                  u2  = elem_i - j - 1;
                  u3  = n - elem_j;
                }
              } else {
                u1  = i - elem_j - 1;
                u2  = elem_i - j - 1;
                u3  = 0;
              }

              tmp_e += P->internal_loop[u1 + u2 + u3];
              break;

            case VRNA_FC_TYPE_COMPARATIVE:
              for (s = 0; s < n_seq; s++) {
                type = vrna_get_ptype_md(SS[s][elem_j], SS[s][elem_i], md);
                if (md->dangles == 2)
                  tmp_e += P->mismatchI[type][S3[s][elem_j]][S5[s][elem_i]];

                if (type > 2)
                  tmp_e += P->TerminalAU;

                if (i < j) {
                  if (i > elem_j) {
                    us1 = (elem_i > 1) ? a2s[s][elem_i - 1] - a2s[s][1]: 0;
                    us2 = a2s[s][i - 1] - a2s[s][elem_j];
                    us3 = a2s[s][n] - a2s[s][j];
                  } else {
                    us1 = (i > 1) ? a2s[s][i - 1] - a2s[s][1] : 0;
                    us2 = a2s[s][elem_i - 1] - a2s[s][j];
                    us3 = a2s[s][n] - a2s[s][elem_j];
                  }
                } else {
                  us1 = a2s[s][elem_i - 1] - a2s[s][j];
                  us2 = a2s[s][i - 1] - a2s[s][elem_j];
                  us3 = 0;
                }

                tmp_e += P->internal_loop[us1 + us2 + us3];
              }
              break;
          }

          break;

        default: /* actually a multibranch-loop-like structure */
          tmp_e = P->MLclosing +
                  (up * P->MLbase) +
                  (num_gq * vrna_E_multibranch_stem(0, -1, -1, P));

          tmp_e *= (int)n_seq;

          /* contribution of the base pair branching-off the mb-loop */
          switch (fc->type) {
            case VRNA_FC_TYPE_COMPARATIVE:
              if (dangles == 2) {
                for (s = 0; s < n_seq; s++) {
                  type  = vrna_get_ptype_md(SS[s][elem_i], SS[s][elem_j], md);
                  tmp_e += vrna_E_multibranch_stem(type, S5[s][elem_i], S3[s][elem_j], P);
                }
              } else {
                for (s = 0; s < n_seq; s++) {
                  type  = vrna_get_ptype_md(SS[s][elem_i], SS[s][elem_j], md);
                  tmp_e += vrna_E_multibranch_stem(type, -1, -1, P);
                }
              }

              break;
            default:
              type = vrna_get_ptype_md(S[elem_i], S[elem_j], md);
              if (dangles == 2)
                tmp_e += vrna_E_multibranch_stem(type, S1[elem_i - 1], S1[elem_j + 1], P);
              else
                tmp_e += vrna_E_multibranch_stem(type, -1, -1, P);

              break;
          }
          break;
      }

      if (*elements)
        vrna_array_append(*elements,
                          ((vrna_struct_elem_t){
                            .type = VRNA_STRUCTURE_ELEM_EXT_LOOP,
                            .energy = (int)tmp_e / (int)n_seq
                          }));

      e_plus += tmp_e;
      break;

    case 2: /* exterior loop hase been falsly interpreted as internal-loop-like */
      e_minus += vrna_eval_internal(fc, elem_i, elem_j, elem_p, elem_q, VRNA_EVAL_LOOP_NO_HC);

      if (*elements_rev)
        vrna_array_append(*elements_rev,
                          ((vrna_struct_elem_t){
                            .type = VRNA_STRUCTURE_ELEM_EXT_LOOP,
                            .energy = (int)e_minus / (int)n_seq
                          }));

      /* we actually have a multibranch-loop-like structure here */
      tmp_e = P->MLclosing +
              (up * P->MLbase) +
              (num_gq * vrna_E_multibranch_stem(0, -1, -1, P));

      tmp_e *= (int)n_seq;

      /* add contributions of the two base pairs branching-off the mb-loop */
      switch (fc->type) {
        case VRNA_FC_TYPE_COMPARATIVE:
          if (dangles == 2) {
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][elem_i], SS[s][elem_j], md);
              type2 = vrna_get_ptype_md(SS[s][elem_p], SS[s][elem_q], md);
              tmp_e += vrna_E_multibranch_stem(type, S5[s][elem_i], S3[s][elem_j], P) +
                       vrna_E_multibranch_stem(type2, S5[s][elem_p], S3[s][elem_q], P);
            }
          } else {
            for (s = 0; s < n_seq; s++) {
              type  = vrna_get_ptype_md(SS[s][elem_i], SS[s][elem_j], md);
              type2 = vrna_get_ptype_md(SS[s][elem_p], SS[s][elem_q], md);
              tmp_e += vrna_E_multibranch_stem(type, -1, -1, P) +
                       vrna_E_multibranch_stem(type2, -1, -1, P);
            }
          }

          break;
        default:
          type  = vrna_get_ptype_md(S[elem_i], S[elem_j], md);
          type2 = vrna_get_ptype_md(S[elem_p], S[elem_q], md);
          if (dangles == 2) {
            tmp_e += vrna_E_multibranch_stem(type, S1[elem_i - 1], S1[elem_j + 1], P) +
                     vrna_E_multibranch_stem(type2, S1[elem_p - 1], S1[elem_q + 1], P);
          } else {
            tmp_e += vrna_E_multibranch_stem(type, -1, -1, P) +
                     vrna_E_multibranch_stem(type2, -1, -1, P);
          }

          break;
      }

      if (*elements)
        vrna_array_append(*elements,
                          ((vrna_struct_elem_t){
                            .type = VRNA_STRUCTURE_ELEM_EXT_LOOP,
                            .energy = (int)tmp_e / (int)n_seq
                          }));

      e_plus += tmp_e;
      break;

    default:  /* exterior loop has been treated as multibranch-loop-like */
      e_minus += (int)up_mis *
                 P->MLbase *
                 (int)n_seq;

      if (*elements_rev)
        vrna_array_append(*elements_rev,
                          ((vrna_struct_elem_t){
                            .type = VRNA_STRUCTURE_ELEM_EXT_LOOP,
                            .energy = (int)e_minus / (int)n_seq
                          }));

      tmp_e = (int)num_gq *
              vrna_E_multibranch_stem(0, -1, -1, P) *
              (int)n_seq;

      if (*elements)
        vrna_array_append(*elements,
                          ((vrna_struct_elem_t){
                            .type = VRNA_STRUCTURE_ELEM_EXT_LOOP,
                            .energy = (int)tmp_e / (int)n_seq
                          }));

      e_plus += tmp_e;

      break;
  }

  corr_en += e_plus -
             e_minus;

  return corr_en;
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
en_corr_of_loop_gquad(vrna_fold_compound_t            *fc,
                      int                             i,
                      int                             j,
                      const char                      *structure,
                      const short                     *pt,
                      const int                       *loop_idx,
                      vrna_array(vrna_struct_elem_t)  *elements,
                      vrna_array(vrna_struct_elem_t)  *elements_rev)
{
  short         *s1, *s2, **S, **S5, **S3;
  unsigned int  cnt, n_seq, L, l[3], pos;
  int           n, tmp_e, energy, p, q, r, s, u, type, type2,
                num_elem, num_g, elem_i, elem_j,
                up_mis, gq_en[2], dangle_model;
  vrna_param_t  *P;
  vrna_md_t     *md;

  n             = fc->length;
  s1            = fc->sequence_encoding;
  s2            = fc->sequence_encoding2;
  S             = fc->S;
  S5            = fc->S5;
  S3            = fc->S3;
  P             = fc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;
  elem_i        = 0;
  elem_j        = 0;

  switch (fc->type) {
    case VRNA_FC_TYPE_COMPARATIVE:
      n_seq = fc->n_seq;
      break;

    default:
      n_seq = 1;
      break;
  }

  energy  = 0;
  q       = i;

  while ((pos = vrna_gq_parse(structure + q - 1, &L, l)) > 0) {
    q += pos - 1;

    if (q > j)
      break;

    num_g = 1;

    /* check whether gquad spans n,1 junction */
    if (4 * L + l[0] + l[1] + l[2] > pos) {
      if ((md->circ) &&
          (num_g == 1)) {
        /* allow for n,1 junction spanning G-Quadruplex only once */
        p = n + pos + 1 - 4 * L - l[0] - l[1] - l[2];
      } else {
        break; /* not allowed for linear sequences */
      }
    } else if (q > j) {
      break;
    } else {
      p = q - 4 * L - l[0] - l[1] - l[2] + 1;
    }

    /* we've found the first g-quadruplex at position [p,q] */
    switch (fc->type) {
      case VRNA_FC_TYPE_COMPARATIVE:
        vrna_E_consensus_gquad(L,
                               l,
                               (unsigned int)p,
                               n,
                               n_seq,
                               (const short **)fc->S,
                               (const unsigned int **)fc->a2s,
                               P,
                               gq_en);
        tmp_e = gq_en[0];
        break;

      default:
        tmp_e = vrna_E_gquad(L, l, P);
        break;
    }

    energy += tmp_e;

    if (*elements)
      vrna_array_append(*elements,
                        ((vrna_struct_elem_t){
                          .type = VRNA_STRUCTURE_ELEM_GQUAD,
                          .i = p,
                          .j = q,
                          .L = L,
                          .l = { l[0], l[1], l[2] },
                          .energy = (int)tmp_e / (int)n_seq
                        }));

    /* check if it's enclosed in a base pair */
    if (loop_idx[p] == 0) {
      if (md->circ) {
        /*  gquad is located in exterior loop of circRNA,
         *  so we need to correct the exterior loop energy
         *  accordingly.
         */
        tmp_e = en_corr_of_loop_gquad_circ(fc,
                                           (unsigned int)p,
                                           (unsigned int)q,
                                           structure,
                                           pt,
                                           loop_idx,
                                           elements,
                                           elements_rev);

        energy = ADD_OR_INF(energy, tmp_e);
        break;
      } else {
        q++;
        continue;                         /* g-quad in exterior loop */
      }
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
          pos = vrna_gq_parse(structure + u - 1, &L, l);
          if (pos > 0) {
            /* sanity check for malformed gquad (this one must not span the n,1 junction! */
            if (4 * L + l[0] + l[1] + l[2] > pos) {
              vrna_log_error("Malformed gquadruplex somewhere after position %u",
                             u - 1);
              return INF;
            }

            switch (fc->type) {
              case VRNA_FC_TYPE_COMPARATIVE:
                vrna_E_consensus_gquad(L,
                                       l,
                                       (unsigned int)u,
                                       n,
                                       n_seq,
                                       (const short **)fc->S,
                                       (const unsigned int **)fc->a2s,
                                       P,
                                       gq_en);
                tmp_e = gq_en[0];
                break;

              default:
                tmp_e = vrna_E_gquad(L, l, P);
                break;
            }

            if (*elements)
              vrna_array_append(*elements,
                                ((vrna_struct_elem_t){
                                  .type = VRNA_STRUCTURE_ELEM_GQUAD,
                                  .i = u,
                                  .j = u + pos - 1,
                                  .L = L,
                                  .l = { l[0], l[1], l[2] },
                                  .energy = (int)tmp_e / (int)n_seq
                                }));

            energy  += tmp_e;
            up_mis  += u + pos - 1;
            u       += pos;
            num_g++;
          }
        } else {
          /* we must have found a stem */
          num_elem++;
          elem_i  = u;
          elem_j  = pt[u];
          tmp_e   = en_corr_of_loop_gquad(fc,
                                          u,
                                          pt[u],
                                          structure,
                                          pt,
                                          loop_idx,
                                          elements,
                                          elements_rev);
          if (tmp_e == INF)
            return INF;

          energy  += tmp_e;
          u       = pt[u] + 1;
        }
      }

      /* here, u == s */
      int e_minus, e_plus, e_temp;

      e_plus = e_minus = 0;

      /* we are done since we've found no other 3' structure element */
      switch (num_elem) {
        /* g-quad was misinterpreted as hairpin closed by (r,s) */
        case 0:
          e_minus = vrna_eval_hairpin(fc, r, s, VRNA_EVAL_LOOP_NO_HC);

          if (*elements_rev)
            vrna_array_append(*elements_rev,
                              ((vrna_struct_elem_t){
                                .type = VRNA_STRUCTURE_ELEM_HP_LOOP,
                                .i = r,
                                .j = s,
                                .energy = (int)e_minus / (int)n_seq
                              }));

          /* if we consider the G-Quadruplex, we have */
          if (num_g == 1) {
            /* a) an internal loop like structure */
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

            if (*elements)
              vrna_array_append(*elements,
                                ((vrna_struct_elem_t){
                                  .type = VRNA_STRUCTURE_ELEM_INT_LOOP,
                                  .i = r,
                                  .j = s,
                                  .p = p,
                                  .q = q,
                                  .energy = (int)e_plus / (int)n_seq
                                }));
          } else {
            /* or b) a multibranch loop like structure */
            e_temp = num_g * vrna_E_multibranch_stem(0, -1, -1, P) +
                     P->MLclosing +
                     (s - r - 1 - up_mis) * P->MLbase;

            e_plus = n_seq * e_temp;

            switch (fc->type) {
              case VRNA_FC_TYPE_COMPARATIVE:
                for (cnt = 0; cnt < n_seq; cnt++) {
                  type    = vrna_get_ptype_md(S[cnt][s], S[cnt][r], md);
                  e_plus  += vrna_E_multibranch_stem(type, S5[cnt][s], S3[cnt][r], P);
                }
                break;

              default:
                type    = md->pair[s2[s]][s2[r]];
                e_plus  += vrna_E_multibranch_stem(type, s1[s - 1], s1[r + 1], P);
                break;
            }

            if (*elements)
              vrna_array_append(*elements,
                                ((vrna_struct_elem_t){
                                  .type = VRNA_STRUCTURE_ELEM_MB_LOOP,
                                  .i = r,
                                  .j = s,
                                  .energy = (int)e_plus / (int)n_seq
                                }));
          }

          energy += e_plus - e_minus;
          break;

        /* g-quad was misinterpreted as internal loop closed by (r,s) with enclosed pair (elem_i, elem_j) */
        case 1:
          e_temp = num_g * vrna_E_multibranch_stem(0, -1, -1, P) +
                   P->MLclosing +
                   (s + elem_i - r - 1 - up_mis - elem_j - 1) * P->MLbase;
          e_plus = n_seq * e_temp;

          switch (fc->type) {
            case VRNA_FC_TYPE_COMPARATIVE:
              for (cnt = 0; cnt < n_seq; cnt++) {
                type    = vrna_get_ptype_md(S[cnt][s], S[cnt][r], md);
                type2   = vrna_get_ptype_md(S[cnt][elem_i], S[cnt][elem_j], md);
                e_plus  += vrna_E_multibranch_stem(type, S5[cnt][s], S3[cnt][r], P) +
                           vrna_E_multibranch_stem(type, S5[cnt][elem_i], S3[cnt][elem_j], P);
              }
              break;

            default:
              type    = md->pair[s2[s]][s2[r]];
              type2   = md->pair[s2[elem_i]][s2[elem_j]];
              e_plus  += vrna_E_multibranch_stem(type, s1[s - 1], s1[r + 1], P) +
                         vrna_E_multibranch_stem(type2, s1[elem_i - 1], s1[elem_j + 1], P);
              break;
          }

          e_minus = vrna_eval_internal(fc, r, s, elem_i, elem_j, VRNA_EVAL_LOOP_NO_HC);
          energy  += e_plus - e_minus;

          if (*elements_rev)
            vrna_array_append(*elements_rev,
                              ((vrna_struct_elem_t){
                                .type = VRNA_STRUCTURE_ELEM_INT_LOOP,
                                .i = r,
                                .j = s,
                                .p = elem_i,
                                .q = elem_j,
                                .energy = (int)e_minus / (int)n_seq
                              }));

          if (*elements)
            vrna_array_append(*elements,
                              ((vrna_struct_elem_t){
                                .type = VRNA_STRUCTURE_ELEM_MB_LOOP,
                                .i = r,
                                .j = s,
                                .energy = (int)e_plus / (int)n_seq
                              }));

          break;

        /* gquad was misinterpreted as unpaired nucleotides in a multiloop */
        default:
          e_minus = (up_mis) * P->MLbase * n_seq;
          e_plus  = num_g * vrna_E_multibranch_stem(0, -1, -1, P) * n_seq;
          energy  += e_plus - e_minus;

          if (*elements_rev)
            vrna_array_append(*elements_rev,
                              ((vrna_struct_elem_t){
                                .type = VRNA_STRUCTURE_ELEM_MB_LOOP,
                                .i = r,
                                .j = s,
                                .energy = (int)e_minus / (int)n_seq
                              }));

          if (*elements)
            vrna_array_append(*elements,
                              ((vrna_struct_elem_t){
                                .type = VRNA_STRUCTURE_ELEM_MB_LOOP,
                                .i = r,
                                .j = s,
                                .energy = (int)e_minus / (int)n_seq
                              }));

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
stack_energy(vrna_fold_compound_t           *fc,
             int                            i,
             const short                    *pt,
             vrna_array(vrna_struct_elem_t) *elements)
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
        vrna_log_warning("Bases %d and %d (%c%c) can't pair!",
                         i, j,
                         string[i - 1],
                         string[j - 1]);
      }

      break;
  }

  p = i;
  q = j;

  while (p < q) {
    /* process all stacks and internal loops */
    while (pt[++p] == 0);
    while (pt[--q] == 0);
    if ((pt[q] != (short)p) || (p > q))
      break;

    if ((sn[i] == sn[p]) &&
        (sn[q] == sn[j])) {
      if (fc->type == VRNA_FC_TYPE_SINGLE) {
        if (md->pair[s[q]][s[p]] == 0) {
          vrna_log_warning("Bases %d and %d (%c%c) can't pair!",
                           p, q,
                           string[p - 1],
                           string[q - 1]);
        }
      }

      ee = vrna_eval_internal(fc, i, j, p, q, VRNA_EVAL_LOOP_NO_HC);

      if (*elements)
        vrna_array_append(*elements,
                          ((vrna_struct_elem_t){
                            .type = VRNA_STRUCTURE_ELEM_INT_LOOP,
                            .i = i,
                            .j = j,
                            .p = p,
                            .q = q,
                            .energy = (int)ee / (int)n_seq
                          }));

      energy  = ADD_OR_INF(energy, ee);
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
      ee = vrna_eval_hairpin(fc, i, j, VRNA_EVAL_LOOP_NO_HC);

      if (ee == INF) {
        if (j - i - 1 < md->min_loop_size) {
          vrna_log_warning("Hairpin loop closed by %d and %d (%c,%c) too short",
                           i, j,
                           vrna_nucleotide_decode(s[i], md),
                           vrna_nucleotide_decode(s[j], md));
        } else {
          vrna_log_warning("Hairpin loop closed by %d and %d (%c,%c) forbidden",
                           i, j,
                           vrna_nucleotide_decode(s[i], md),
                           vrna_nucleotide_decode(s[j], md));
        }
      }

      if (*elements)
        vrna_array_append(*elements,
                          ((vrna_struct_elem_t){
                            .type = VRNA_STRUCTURE_ELEM_HP_LOOP,
                            .i = i,
                            .j = j,
                            .energy = (int)ee / (int)n_seq
                          }));

      energy = ADD_OR_INF(energy, ee);
    }

    return energy;
  }

  /* (i,j) is exterior pair of multiloop or external loop */
  if (!first_pair_after_last_nick(i, j, pt, sn)) {
    while (p < j) {
      /* add up the contributions of the substructures of the ML */
      energy  += stack_energy(fc, p, pt, elements);
      p       = pt[p];
      /* search for next base pair in multiloop */
      while (pt[++p] == 0);
    }

    ee = energy_of_ml_pt(fc, i, pt);

    if (*elements)
      vrna_array_append(*elements,
                        ((vrna_struct_elem_t){
                          .type = VRNA_STRUCTURE_ELEM_MB_LOOP,
                          .i = i,
                          .j = j,
                          .energy = (int)ee / (int)n_seq
                        }));

    energy = ADD_OR_INF(energy, ee);
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
                unsigned int          i,
                const short           *pt)
{
  short         *s, *s1, **S, **S5, **S3;
  unsigned int  *sn, **a2s, n_seq, ss, i1, j, p, q, q_prev, q_prev2, n, dangle_model, logML,
                circular;
  int           energy, cx_energy, tmp, tmp2, best_energy, bonus, *idx, *rtype, u, uu,
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

  best_energy = INF;

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

  if (i >= (unsigned int)pt[i]) {
    vrna_log_warning("i = %d (pt[%d] = %d) is not 5' base of a closing pair!", i, i, pt[i]);
    return INF;
  }

  j = (i == 0) ? n + 1 : (unsigned int)pt[i];

  switch (fc->type) {
    case VRNA_FC_TYPE_COMPARATIVE:
      if ((dangle_model % 2) ||
          (dangle_model > 2)) {
        vrna_log_warning(
          "Consensus structure evaluation for odd dangle models not implemented (yet)!");
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

              energy += vrna_E_multibranch_stem(tt, -1, -1, P);
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

              energy += vrna_E_multibranch_stem(tt, -1, -1, P);
            }
          }

          break;

        default:
          while (p < j) {
            /* p must have a pairing partner */
            q = (int)pt[p];
            /* get type of base pair (p,q) */
            tt = vrna_get_ptype_md(s[p], s[q], md);

            energy += vrna_E_multibranch_stem(tt, -1, -1, P);

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

            energy += vrna_E_multibranch_stem(tt, -1, -1, P);
          } else {
            /* virtual closing pair */
            energy += vrna_E_multibranch_stem(0, -1, -1, P);
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
              energy  += vrna_E_multibranch_stem(tt, mm5, mm3, P);
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
              energy  += vrna_E_multibranch_stem(tt, mm5, mm3, P);
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
            energy  += vrna_E_multibranch_stem(tt, mm5, mm3, P);

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
            energy  += vrna_E_multibranch_stem(tt, mm5, mm3, P);
          } else {
            /* virtual closing pair */
            energy += vrna_E_multibranch_stem(0, -1, -1, P);
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
        e_stem    = vrna_E_multibranch_stem(tt, -1, -1, P);
        e_stem5   = vrna_E_multibranch_stem(tt, mm5, -1, P);
        e_stem3   = vrna_E_multibranch_stem(tt, -1, mm3, P);
        e_stem53  = vrna_E_multibranch_stem(tt, mm5, mm3, P);

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

        e_stem    = vrna_E_multibranch_stem(type, -1, -1, P);
        e_stem5   = vrna_E_multibranch_stem(type, mm5, -1, P);
        e_stem3   = vrna_E_multibranch_stem(type, -1, mm3, P);
        e_stem53  = vrna_E_multibranch_stem(type, mm5, mm3, P);
      } else {
        /* virtual closing pair */
        e_stem = e_stem5 = e_stem3 = e_stem53 = vrna_E_multibranch_stem(0, -1, -1, P);
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


PRIVATE int
struct_elem_compare(const void  *a,
                    const void  *b)
{
  int                 ret = 0;

  vrna_struct_elem_t  *aa = (vrna_struct_elem_t *)a;
  vrna_struct_elem_t  *bb = (vrna_struct_elem_t *)b;

  ret = (int)aa->i - (int)bb->i;

  if (ret == 0)
    ret = (int)aa->j - (int)bb->j;

  return ret;
}


PRIVATE void
print_structure_elements(const char                     *s,
                         vrna_array(vrna_struct_elem_t) elements,
                         vrna_cstr_t                    output_stream)
{
  if (elements) {
    /* sort elements from 5' to 3' position */
    qsort(elements,
          vrna_array_size(elements),
          sizeof(vrna_struct_elem_t),
          struct_elem_compare);

    /* print individual energy contributions */
    for (size_t i = 0; i < vrna_array_size(elements); i++) {
      switch (elements[i].type) {
        case VRNA_STRUCTURE_ELEM_EXT_LOOP:
          vrna_cstr_print_eval_ext_loop(output_stream,
                                        elements[i].energy);
          break;

        case VRNA_STRUCTURE_ELEM_HP_LOOP:
          vrna_cstr_print_eval_hp_loop(output_stream,
                                       elements[i].i,
                                       elements[i].j,
                                       s[elements[i].i - 1],
                                       s[elements[i].j - 1],
                                       elements[i].energy);
          break;

        case VRNA_STRUCTURE_ELEM_INT_LOOP:
          vrna_cstr_print_eval_int_loop(output_stream,
                                        elements[i].i,
                                        elements[i].j,
                                        s[elements[i].i - 1],
                                        s[elements[i].j - 1],
                                        elements[i].p,
                                        elements[i].q,
                                        s[elements[i].p - 1],
                                        s[elements[i].q - 1],
                                        elements[i].energy);
          break;

        case VRNA_STRUCTURE_ELEM_MB_LOOP:
          vrna_cstr_print_eval_mb_loop(output_stream,
                                       elements[i].i,
                                       elements[i].j,
                                       s[elements[i].i - 1],
                                       s[elements[i].j - 1],
                                       elements[i].energy);
          break;

        case VRNA_STRUCTURE_ELEM_GQUAD:
          vrna_cstr_print_eval_gquad(output_stream,
                                     elements[i].i,
                                     elements[i].j,
                                     elements[i].L,
                                     elements[i].l,
                                     elements[i].energy);
          break;

        case VRNA_STRUCTURE_ELEM_UD:
        /* fall-through */
        case VRNA_STRUCTURE_ELEM_STACK:
        /* fall-through */
        case VRNA_STRUCTURE_ELEM_BULGE:
        /* fall-through */
        default:
          break;
      }
    }
  }
}
