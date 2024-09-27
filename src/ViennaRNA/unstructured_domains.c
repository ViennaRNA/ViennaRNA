/** \file unstructured_domains.c **/

/*
 *                Unstructured domains
 *
 *                This file contains everything necessary to
 *                deal with the default implementation for unstructured
 *                domains in secondary structures. This feature enables,
 *                for instance, ligand binding to unpaired stretches of
 *                an RNA secondary structure.
 *
 *                c 2016 Ronny Lorenz
 *
 *                ViennaRNA package
 */

#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/eval/structures.h"
#include "ViennaRNA/unstructured_domains.h"

/*
 #################################
 # PRIVATE MACROS                #
 #################################
 */

/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */

/*
 #################################
 # PRIVATE VARIABLES/STRUCTS     #
 #################################
 */

struct binding_segment {
  unsigned int  start;
  unsigned int  end;
  unsigned int  type;
};


struct ud_bt_stack {
  unsigned int    from;
  vrna_ud_motif_t *motifs;
  unsigned int    motif_cnt;
  unsigned int    motif_size;
};


struct default_outside {
  int         motif_num;
  FLT_OR_DBL  exp_energy;
};

/*
 *  Default data structure for ligand binding to unpaired stretches
 */
struct ligands_up_data_default {
  /*
   **********************************
   * pre-computed position-wise
   * motif list
   **********************************
   */
  int         n;
  int         **motif_list_ext;
  int         **motif_list_hp;
  int         **motif_list_int;
  int         **motif_list_mb;

  int         *dG;
  FLT_OR_DBL  *exp_dG;
  int         *len;

  /*
   **********************************
   * below are DP matrices to store
   * the production rule results
   **********************************
   */
  int         *energies_ext;
  int         *energies_hp;
  int         *energies_int;
  int         *energies_mb;
  FLT_OR_DBL  *exp_energies_ext;
  FLT_OR_DBL  *exp_energies_hp;
  FLT_OR_DBL  *exp_energies_int;
  FLT_OR_DBL  *exp_energies_mb;

  /*
   **********************************
   * below are lists to store the
   * outside partition function for
   * each motif starting at each
   * position
   **********************************
   */
  unsigned int            *outside_ext_count;
  struct default_outside  **outside_ext;
  unsigned int            *outside_hp_count;
  struct default_outside  **outside_hp;
  unsigned int            *outside_int_count;
  struct default_outside  **outside_int;
  unsigned int            *outside_mb_count;
  struct default_outside  **outside_mb;

  FLT_OR_DBL              (*default_cb[32])(int,
                                            int,
                                            struct ligands_up_data_default *);
  FLT_OR_DBL              *exp_e_mx[32];
};

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE struct binding_segment *
extract_binding_segments(const char   *structure,
                         unsigned int *segments_num);


PRIVATE void
fill_MFE_matrix(vrna_fold_compound_t  *fc,
                int                   *mx,
                unsigned int          from,
                unsigned int          to,
                unsigned int          type);


PRIVATE vrna_ud_motif_t *
backtrack_MFE_matrix(vrna_fold_compound_t *fc,
                     int                  *mx,
                     unsigned int         from,
                     unsigned int         to,
                     unsigned int         type);


PRIVATE vrna_ud_motif_t **
backtrack_MFE_matrix_exhaustive(vrna_fold_compound_t  *fc,
                                int                   *mx,
                                unsigned int          from,
                                unsigned int          to,
                                unsigned int          type);


PRIVATE void
fill_MEA_matrix(vrna_fold_compound_t  *fc,
                float                 *mx,
                unsigned int          from,
                unsigned int          to,
                float                 *pu,
                unsigned int          t);


PRIVATE vrna_ud_motif_t *
backtrack_MEA_matrix(vrna_fold_compound_t *fc,
                     float                *mx,
                     unsigned int         from,
                     unsigned int         to,
                     float                *pu,
                     unsigned int         t);


PRIVATE void
remove_ligands_up(vrna_fold_compound_t *vc);


PRIVATE void
init_ligands_up(vrna_fold_compound_t *vc);


PRIVATE void
add_ligand_motif(vrna_fold_compound_t *vc,
                 const char           *motif,
                 double               motif_en,
                 const char           *motif_name,
                 unsigned int         loop_type);


PRIVATE void
remove_default_data(void *d);


/* default implementations for unstructured domains feature */
PRIVATE void
default_prod_rule(vrna_fold_compound_t  *vc,
                  void                  *d);


PRIVATE void
default_exp_prod_rule(vrna_fold_compound_t  *vc,
                      void                  *d);


PRIVATE int
default_energy(vrna_fold_compound_t *vc,
               int                  i,
               int                  j,
               unsigned int         loop_type,
               void                 *d);


PRIVATE FLT_OR_DBL
default_exp_energy(vrna_fold_compound_t *vc,
                   int                  i,
                   int                  j,
                   unsigned int         loop_type,
                   void                 *d);


PRIVATE void
default_probs_add(vrna_fold_compound_t  *vc,
                  int                   i,
                  int                   j,
                  unsigned int          loop_type,
                  FLT_OR_DBL            exp_energy,
                  void                  *data);


PRIVATE FLT_OR_DBL
default_probs_get(vrna_fold_compound_t  *vc,
                  int                   i,
                  int                   j,
                  unsigned int          loop_type,
                  int                   motif,
                  void                  *data);


/* helper functions for default implementatations of unstructured domains feature */
PRIVATE int
default_energy_ext_motif(int                            i,
                         int                            j,
                         struct ligands_up_data_default *data);


PRIVATE int
default_energy_hp_motif(int                             i,
                        int                             j,
                        struct ligands_up_data_default  *data);


PRIVATE int
default_energy_int_motif(int                            i,
                         int                            j,
                         struct ligands_up_data_default *data);


PRIVATE int
default_energy_mb_motif(int                             i,
                        int                             j,
                        struct ligands_up_data_default  *data);


PRIVATE FLT_OR_DBL
default_exp_energy_ext_motif(int                            i,
                             int                            j,
                             struct ligands_up_data_default *data);


PRIVATE FLT_OR_DBL
default_exp_energy_hp_motif(int                             i,
                            int                             j,
                            struct ligands_up_data_default  *data);


PRIVATE FLT_OR_DBL
default_exp_energy_int_motif(int                            i,
                             int                            j,
                             struct ligands_up_data_default *data);


PRIVATE FLT_OR_DBL
default_exp_energy_mb_motif(int                             i,
                            int                             j,
                            struct ligands_up_data_default  *data);


PRIVATE void
free_default_data_matrices(struct ligands_up_data_default *data);


PRIVATE void
free_default_data_exp_matrices(struct ligands_up_data_default *data);


PRIVATE void
prepare_matrices(vrna_fold_compound_t           *vc,
                 struct ligands_up_data_default *data);


PRIVATE void
prepare_exp_matrices(vrna_fold_compound_t           *vc,
                     struct ligands_up_data_default *data);


PRIVATE struct ligands_up_data_default *
get_default_data(void);


PRIVATE void
prepare_default_data(vrna_fold_compound_t           *vc,
                     struct ligands_up_data_default *data);


PRIVATE void
free_default_data(struct ligands_up_data_default *data);


PRIVATE int *
get_motifs(vrna_fold_compound_t *vc,
           int                  i,
           unsigned int         loop_type);


PRIVATE void
annotate_ud(vrna_fold_compound_t  *vc,
            int                   start,
            int                   end,
            char                  l,
            vrna_ud_motif_t       **list,
            int                   *list_size,
            int                   *list_pos);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_ud_motif_t *
vrna_ud_motifs_centroid(vrna_fold_compound_t  *fc,
                        const char            *structure)
{
  unsigned int    i, j, m, s, motifs_size, motifs_count, t, num_segments;
  vrna_ud_motif_t *motifs;
  vrna_ud_t       *domains_up;

  motifs = NULL;

  if ((fc) &&
      (fc->domains_up) &&
      (fc->domains_up->probs_get) &&
      (structure)) {
    domains_up = fc->domains_up;
    struct binding_segment *segments = extract_binding_segments(structure, &num_segments);

    motifs_size   = 10;
    motifs_count  = 0;
    motifs        = (vrna_ud_motif_t *)vrna_alloc(sizeof(vrna_ud_motif_t) * (motifs_size + 1));

    for (s = 0; s < num_segments; s++) {
      t = segments[s].type;

      for (i = segments[s].start; i <= segments[s].end; i++) {
        for (m = 0; m < (unsigned int)domains_up->motif_count; m++) {
          j = i + domains_up->motif_size[m] - 1;
          if (j <= segments[s].end) {
            /* actually, the condition below should check whether the motif probability makes up more than 50% of the unpaired probability */
            if (domains_up->probs_get(fc, i, j, t, m, domains_up->data) > 0.5) {
              motifs[motifs_count].start  = i;
              motifs[motifs_count].number = m;
              motifs_count++;

              if (motifs_count == motifs_size) {
                motifs_size *= 1.4;
                motifs      = (vrna_ud_motif_t *)vrna_realloc(motifs,
                                                              sizeof(vrna_ud_motif_t) *
                                                              (motifs_size + 1));
              }
            }
          }
        }
      }
    }

    free(segments);

    if (motifs_count > 0) {
      /* add end of list marker */
      motifs[motifs_count].start  = 0;
      motifs[motifs_count].number = -1;
      motifs                      = (vrna_ud_motif_t *)vrna_realloc(motifs,
                                                                    sizeof(vrna_ud_motif_t) *
                                                                    (motifs_count + 1));
    } else {
      free(motifs);
      motifs = NULL;
    }
  }

  return motifs;
}


PUBLIC vrna_ud_motif_t *
vrna_ud_motifs_MEA(vrna_fold_compound_t *fc,
                   const char           *structure,
                   vrna_ep_t            *probability_list)
{
  unsigned int            from, to, i, n, s, motifs_size, motifs_count, t, num_segments;
  float                   *pu, *mx;
  vrna_ep_t               *ptr;
  vrna_ud_motif_t         *motifs;
  struct binding_segment  *segments;

  motifs = NULL;

  if ((fc) &&
      (fc->domains_up) &&
      (fc->domains_up->probs_get) &&
      (structure) &&
      (probability_list)) {
    n         = fc->length;
    segments  = extract_binding_segments(structure, &num_segments);
    pu        = (float *)vrna_alloc(sizeof(float) * (n + 1));
    mx        = (float *)vrna_alloc(sizeof(float) * (n + 1));

    /* determine probabilities to be unpaired */
    for (i = 1; i <= n; i++)
      pu[i] = 1.;

    for (ptr = probability_list; ptr->i > 0; ptr++) {
      if (ptr->type == VRNA_PLIST_TYPE_BASEPAIR) {
        pu[ptr->i]  -= ptr->p;
        pu[ptr->j]  -= ptr->p;
      } else if (ptr->type == VRNA_PLIST_TYPE_UD_MOTIF) {
        for (i = (unsigned int)ptr->i; i <= (unsigned int)ptr->j; i++)
          pu[i] -= ptr->p;
      }
    }

    motifs_count  = 0;
    motifs_size   = 10;
    motifs        = (vrna_ud_motif_t *)vrna_alloc(sizeof(vrna_ud_motif_t) * (motifs_size + 1));

    for (s = 0; s < num_segments; s++) {
      vrna_ud_motif_t *m;

      from  = segments[s].start;
      to    = segments[s].end;
      t     = segments[s].type;

      fill_MEA_matrix(fc, mx, from, to, pu, t);
      m = backtrack_MEA_matrix(fc, mx, from, to, pu, t);

      if (m) {
        /* determine number of new motifs */
        for (i = 0; m[i].start != 0; i++);

        /* resize target memory if necessary */
        if (motifs_count + i >= motifs_size) {
          motifs_size += motifs_size / 2 + 1 + i;
          motifs      = (vrna_ud_motif_t *)vrna_realloc(motifs,
                                                        sizeof(vrna_ud_motif_t) *
                                                        (motifs_size + 1));
        }

        /* append data */
        memcpy(motifs + motifs_count, m, sizeof(vrna_ud_motif_t) * i);

        /* increase motif counter */
        motifs_count += i;

        /* release backtracked motif list */
        free(m);
      }
    }
    free(mx);
    free(pu);
    free(segments);

    if (motifs_count > 0) {
      /* add end of list marker */
      motifs[motifs_count].start  = 0;
      motifs[motifs_count].number = -1;
      motifs                      = (vrna_ud_motif_t *)vrna_realloc(motifs,
                                                                    sizeof(vrna_ud_motif_t) *
                                                                    (motifs_count + 1));
    } else {
      free(motifs);
      motifs = NULL;
    }
  }

  return motifs;
}


PUBLIC vrna_ud_motif_t *
vrna_ud_motifs_MFE(vrna_fold_compound_t *fc,
                   const char           *structure)
{
  unsigned int            from, to, i, n, s, motifs_size, motifs_count, t, num_segments;
  int                     *mx;
  vrna_ud_motif_t         *motifs;
  struct binding_segment  *segments;

  motifs = NULL;

  if ((fc) && (fc->domains_up) && (fc->domains_up->probs_get) && (structure)) {
    n         = fc->length;
    segments  = extract_binding_segments(structure, &num_segments);
    mx        = (int *)vrna_alloc(sizeof(int) * (n + 1));

    motifs_count  = 0;
    motifs_size   = 10;
    motifs        = (vrna_ud_motif_t *)vrna_alloc(sizeof(vrna_ud_motif_t) * (motifs_size + 1));

    for (s = 0; s < num_segments; s++) {
      vrna_ud_motif_t *m;

      from  = segments[s].start;
      to    = segments[s].end;
      t     = segments[s].type;

      fill_MFE_matrix(fc, mx, from, to, t);
      m = backtrack_MFE_matrix(fc, mx, from, to, t);

      if (m) {
        /* determine number of new motifs */
        for (i = 0; m[i].start != 0; i++);

        /* resize target memory if necessary */
        if (motifs_count + i >= motifs_size) {
          motifs_size += motifs_size / 2 + 1 + i;
          motifs      = (vrna_ud_motif_t *)vrna_realloc(motifs,
                                                        sizeof(vrna_ud_motif_t) *
                                                        (motifs_size + 1));
        }

        /* append data */
        memcpy(motifs + motifs_count, m, sizeof(vrna_ud_motif_t) * i);

        /* increase motif counter */
        motifs_count += i;

        /* release backtracked motif list */
        free(m);
      }
    }
    free(mx);
    free(segments);

    if (motifs_count > 0) {
      /* add end of list marker */
      motifs[motifs_count].start  = 0;
      motifs[motifs_count].number = -1;
      motifs                      = (vrna_ud_motif_t *)vrna_realloc(motifs,
                                                                    sizeof(vrna_ud_motif_t) *
                                                                    (motifs_count + 1));
    } else {
      free(motifs);
      motifs = NULL;
    }
  }

  return motifs;
}


PUBLIC void
vrna_ud_remove(vrna_fold_compound_t *vc)
{
  if (vc && vc->domains_up)
    remove_ligands_up(vc);
}


PUBLIC void
vrna_ud_set_data(vrna_fold_compound_t       *vc,
                 void                       *data,
                 vrna_auxdata_free_f free_cb)
{
  if (vc) {
    /* init if not already present */
    if (!vc->domains_up)
      init_ligands_up(vc);

    /* free previous data if 'free_data' function present */
    if (vc->domains_up->free_data)
      vc->domains_up->free_data(vc->domains_up->data);

    /* set new data and free callback */
    vc->domains_up->free_data = free_cb;
    vc->domains_up->data      = data;
  }
}


PUBLIC void
vrna_ud_set_prod_rule_cb(vrna_fold_compound_t         *vc,
                         vrna_ud_production_f  pre_cb,
                         vrna_ud_f      e_cb)
{
  if (vc) {
    /* init if not already present */
    if (!vc->domains_up)
      init_ligands_up(vc);

    /* set new callback */
    vc->domains_up->prod_cb   = pre_cb;
    vc->domains_up->energy_cb = e_cb;
  }
}


PUBLIC void
vrna_ud_set_exp_prod_rule_cb(vrna_fold_compound_t             *vc,
                             vrna_ud_exp_production_f  pre_cb,
                             vrna_ud_exp_f      exp_e_cb)
{
  if (vc) {
    /* init if not already present */
    if (!vc->domains_up)
      init_ligands_up(vc);

    /* set new callback */
    vc->domains_up->exp_prod_cb   = pre_cb;
    vc->domains_up->exp_energy_cb = exp_e_cb;
  }
}


PUBLIC void
vrna_ud_set_prob_cb(vrna_fold_compound_t        *vc,
                    vrna_ud_add_probs_f  setter,
                    vrna_ud_get_probs_f  getter)
{
  if (vc) {
    /* init if not already present */
    if (!vc->domains_up)
      init_ligands_up(vc);

    /* set new callback */
    vc->domains_up->probs_add = setter;
    vc->domains_up->probs_get = getter;
  }
}


PUBLIC void
vrna_ud_add_motif(vrna_fold_compound_t  *vc,
                  const char            *motif,
                  double                motif_en,
                  const char            *motif_name,
                  unsigned int          loop_type)
{
  if (vc) {
    if (!vc->domains_up) {
      /* set all default callbacks */
      vrna_ud_set_prod_rule_cb(vc, &default_prod_rule, &default_energy);
      vrna_ud_set_exp_prod_rule_cb(vc, &default_exp_prod_rule, &default_exp_energy);
      vrna_ud_set_data(vc, get_default_data(), &remove_default_data);
      vrna_ud_set_prob_cb(vc, &default_probs_add, &default_probs_get);
    }

    add_ligand_motif(vc, motif, motif_en, motif_name, loop_type);
  }
}


PUBLIC int *
vrna_ud_get_motif_size_at(vrna_fold_compound_t  *vc,
                          int                   i,
                          unsigned int          loop_type)
{
  if (vc && vc->domains_up) {
    int k, l, cnt, *ret, *ptr;

    ret = NULL;
    if ((i > 0) &&
        ((unsigned int)i <= vc->length)) {
      ptr = get_motifs(vc, i, loop_type);
      if (ptr) {
        for (k = 0; ptr[k] != -1; k++) /* replace motif number with its size */
          ptr[k] = vc->domains_up->motif_size[ptr[k]];
        /* make the list unique */
        ret     = (int *)vrna_alloc(sizeof(int) * (k + 1));
        ret[0]  = -1;
        cnt     = 0;
        for (k = 0; ptr[k] != -1; k++) {
          for (l = 0; l < cnt; l++)
            if (ptr[k] == ret[l])
              break;

          /* we've already seen this size */

          if (l == cnt) {
            /* we've not seen this size before */
            ret[cnt]      = ptr[k];
            ret[cnt + 1]  = -1;
            cnt++;
          }
        }
        /* resize ret array */
        ret = (int *)vrna_realloc(ret, sizeof(int) * (cnt + 1));
      }

      free(ptr);
    }

    return ret;
  }

  return NULL;
}


PUBLIC int *
vrna_ud_get_motifs_at(vrna_fold_compound_t  *vc,
                      int                   i,
                      unsigned int          loop_type)
{
  if (vc && vc->domains_up)
    if ((i > 0) &&
        ((unsigned int)i <= vc->length))
      return get_motifs(vc, i, loop_type);

  return NULL;
}


vrna_ud_motif_t *
vrna_ud_detect_motifs(vrna_fold_compound_t  *vc,
                      const char            *structure)
{
  int             list_size, list_pos;
  vrna_ud_motif_t *motif_list;

  motif_list = NULL;

  if (structure && vc->domains_up) {
    unsigned int  l, start, end;
    char          last, *loops;

    l           = 0;
    list_pos    = 0;
    list_size   = 15;
    motif_list  = (vrna_ud_motif_t *)vrna_alloc(sizeof(vrna_ud_motif_t) * list_size);
    loops       = vrna_db_to_element_string(structure);

    while (l < vc->length) {
      /* skip uppercase encodings */
      while (l < vc->length) {
        if (islower(loops[l]))
          break;

        l++;
      }

      if (l < vc->length) {
        start = 1 + l;
        last  = loops[l];
        while (loops[l++] == last)
          if (l == vc->length)
            break;

        end = l - 1;
        annotate_ud(vc, start, end, last, &motif_list, &list_size, &list_pos);
      }
    }

    motif_list = (vrna_ud_motif_t *)vrna_realloc(motif_list,
                                                 sizeof(vrna_ud_motif_t) *
                                                 (list_pos + 1));
    motif_list[list_pos].start  = 0;
    motif_list[list_pos].number = -1;
    free(loops);
  }

  return motif_list;
}


PRIVATE vrna_ud_motif_t **
ud_get_motifs_MFE(vrna_fold_compound_t    *fc,
                  struct binding_segment  *segments,
                  unsigned int            segments_num)
{
  vrna_ud_motif_t **motif_lists;

  motif_lists = NULL;

  if (segments) {
    unsigned int    n, alt_cnt, shift, l, d, m, *merger_cnt, num_combinations, start, end, t;
    int             *mx;
    vrna_ud_motif_t **ptr, *ptr2, ***alternatives;

    alternatives  = (vrna_ud_motif_t ***)vrna_alloc(sizeof(vrna_ud_motif_t * *) * segments_num);
    alt_cnt       = 0;

    /* collect optimal configurations for each segment */
    for (n = 0; n < segments_num; n++) {
      start = segments[n].start;
      end   = segments[n].end;
      t     = segments[n].type;

      mx  = (int *)vrna_alloc(sizeof(int) * (end - start + 2));
      mx  -= start;

      fill_MFE_matrix(fc, mx, start, end, t);
      if ((ptr = backtrack_MFE_matrix_exhaustive(fc, mx, start, end, t)))
        alternatives[alt_cnt++] = ptr;

      mx += start;
      free(mx);
    }

    /* prepare for merging process */

    num_combinations  = 1;
    merger_cnt        = (unsigned int *)vrna_alloc(sizeof(unsigned int) * alt_cnt);
    vrna_ud_motif_t **merger = (vrna_ud_motif_t **)vrna_alloc(sizeof(vrna_ud_motif_t *) * alt_cnt);

    for (n = 0; n < alt_cnt; n++) {
      /* count number of alternatives for current segment */
      for (d = 0; alternatives[n][d]; d++);

      if (d)
        num_combinations *= d;

      /* set merger pointers */
      merger[n] = alternatives[n][0];

      /* count number of elements at top of list */
      for (ptr2 = merger[n]; ptr2->start != 0; ++ptr2);

      /* store motif counter for current configuration */
      merger_cnt[n] = ptr2 - merger[n];
    }

    /* finally merge the segments */
    motif_lists = (vrna_ud_motif_t **)vrna_alloc(sizeof(vrna_ud_motif_t *) *
                                                 (num_combinations + 1));

    for (n = 0; n < num_combinations; n++) {
      /* determine length of list */
      for (l = m = 0; m < alt_cnt; m++)
        l += merger_cnt[m];

      motif_lists[n] = (vrna_ud_motif_t *)vrna_alloc(sizeof(vrna_ud_motif_t) * (l + 1));
      for (shift = m = 0; m < alt_cnt; m++) {
        memcpy(motif_lists[n] + shift, merger[m], sizeof(vrna_ud_motif_t) * merger_cnt[m]);
        shift += merger_cnt[m];
      }

      motif_lists[n][l].start   = 0;
      motif_lists[n][l].number  = -1;

      /* update (increase) merger pointers */
      for (m = alt_cnt; m > 0; m--) {
        ++(merger[m - 1]);
        if (merger[m - 1]) {
          for (ptr2 = merger[m - 1]; ptr2->start != 0; ++ptr2);

          merger_cnt[m - 1] = ptr2 - merger[m - 1];
          break;
        } else if (m == 1) {
          break;
        } else {
          merger[m - 1] = alternatives[m - 1][0];
          for (ptr2 = merger[m - 1]; ptr2->start != 0; ++ptr2);

          merger_cnt[m - 1] = ptr2 - merger[m - 1];
        }
      }
    }

    free(merger);
    free(merger_cnt);

    for (n = 0; n < alt_cnt; n++) {
      for (m = 0; alternatives[n][m]; m++)
        free(alternatives[n][m]);
      free(alternatives[n]);
    }
    free(alternatives);

    motif_lists[num_combinations] = NULL;
  }

  return motif_lists;
}


PRIVATE vrna_ud_motif_t **
ud_get_motifs_energy(vrna_fold_compound_t   *fc VRNA_UNUSED,
                     struct binding_segment *segments,
                     unsigned int           segments_num VRNA_UNUSED,
                     int                    e VRNA_UNUSED)
{
  vrna_ud_motif_t **motif_lists;

  vrna_log_error("not implemented!");
  motif_lists = NULL;

  if (segments) {
  }

  return motif_lists;
}


PUBLIC vrna_ud_motif_t **
vrna_ud_extract_motifs(vrna_fold_compound_t *fc,
                       const char           *structure,
                       float                *energy)
{
  vrna_ud_motif_t **motif_lists;

  motif_lists = NULL;

  if ((fc) && (fc->domains_up) && (structure)) {
    unsigned int            segments_num;
    struct binding_segment  *segments;

    segments = extract_binding_segments(structure, &segments_num);

    if (!energy) {
      /* get MFE arrangement(s) */
      motif_lists = ud_get_motifs_MFE(fc, segments, segments_num);
    } else {
      /* get arrangement(s) that result in provided free energy */
      float e   = vrna_eval_structure(fc, structure);
      int   de  = (int)roundf(*energy - e) * 100; /* binding free energy in deka cal/mol */
      motif_lists = ud_get_motifs_energy(fc, segments, segments_num, de);
    }

    free(segments);
  }

  return motif_lists;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE struct ligands_up_data_default *
get_default_data(void)
{
  struct ligands_up_data_default *data = vrna_alloc(sizeof(struct ligands_up_data_default));

  data->n                 = 0;
  data->motif_list_ext    = NULL;
  data->motif_list_hp     = NULL;
  data->motif_list_int    = NULL;
  data->motif_list_mb     = NULL;
  data->dG                = NULL;
  data->exp_dG            = NULL;
  data->energies_ext      = NULL;
  data->energies_hp       = NULL;
  data->energies_int      = NULL;
  data->energies_mb       = NULL;
  data->exp_energies_ext  = NULL;
  data->exp_energies_hp   = NULL;
  data->exp_energies_int  = NULL;
  data->exp_energies_mb   = NULL;
  data->outside_ext       = NULL;
  data->outside_hp        = NULL;
  data->outside_int       = NULL;
  data->outside_mb        = NULL;
  data->outside_ext_count = NULL;
  data->outside_hp_count  = NULL;
  data->outside_int_count = NULL;
  data->outside_mb_count  = NULL;
  return data;
}


PRIVATE void
remove_ligands_up(vrna_fold_compound_t *vc)
{
  int i;

  /* free auxiliary data */
  if (vc->domains_up->free_data)
    vc->domains_up->free_data(vc->domains_up->data);

  /* free motif sequences */
  for (i = 0; i < vc->domains_up->motif_count; i++)
    free(vc->domains_up->motif[i]);

  /* free motif names */
  for (i = 0; i < vc->domains_up->motif_count; i++)
    free(vc->domains_up->motif_name[i]);

  free(vc->domains_up->motif);
  free(vc->domains_up->motif_name);
  free(vc->domains_up->motif_size);
  free(vc->domains_up->motif_en);
  free(vc->domains_up->motif_type);

  free(vc->domains_up->uniq_motif_size);

  free(vc->domains_up);

  vc->domains_up = NULL;
}


PRIVATE void
init_ligands_up(vrna_fold_compound_t *vc)
{
  vc->domains_up = (vrna_ud_t *)vrna_alloc(sizeof(vrna_ud_t));

  vc->domains_up->uniq_motif_count  = 0;
  vc->domains_up->uniq_motif_size   = NULL;
  vc->domains_up->motif_count       = 0;
  vc->domains_up->motif             = NULL;
  vc->domains_up->motif_name        = NULL;
  vc->domains_up->motif_size        = NULL;
  vc->domains_up->motif_en          = NULL;
  vc->domains_up->motif_type        = NULL;
  vc->domains_up->prod_cb           = NULL;
  vc->domains_up->exp_prod_cb       = NULL;
  vc->domains_up->energy_cb         = NULL;
  vc->domains_up->exp_energy_cb     = NULL;
  vc->domains_up->data              = NULL;
  vc->domains_up->free_data         = NULL;
  vc->domains_up->probs_add         = NULL;
  vc->domains_up->probs_get         = NULL;
}


/*
 *  Given a secondary structure in dot-bracket notation,
 *  extract all consecutive unpaired nucleotides as individual
 *  segments.
 *  After successful execution, this function returns a list of
 *  segments and the number of list elements is stored in the
 *  variable passed as segments_num
 */
PRIVATE struct binding_segment *
extract_binding_segments(const char   *structure,
                         unsigned int *segments_num)
{
  struct binding_segment  *segments;
  char                    *loops;
  unsigned int            n, pos, start, segments_length;

  segments  = NULL;
  n         = strlen(structure);
  loops     = vrna_db_to_element_string(structure);

  *segments_num   = 0;
  segments_length = 15;
  segments        = (struct binding_segment *)vrna_alloc(sizeof(struct binding_segment) *
                                                         segments_length);

  pos = 1;

  /* extract segments that possibly harbor bound ligands */
  while (pos <= n) {
    /* skip uppercase encodings, i.e. paired nucleotides */
    for (; isupper(loops[pos - 1]) && (pos <= n); pos++);

    /* no more unpaired segments */
    if (pos > n)
      break;

    start = pos;

    /* find next uppercase encoding, i.e. paired nucleotides */
    for (; islower(loops[pos - 1]) && (pos <= n); pos++);

    segments[(*segments_num)].start = start;
    segments[(*segments_num)].end   = pos - 1;
    segments[(*segments_num)].type  = 0;

    if (loops[start - 1] == 'e')
      segments[(*segments_num)].type = VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP;
    else if (loops[start - 1] == 'h')
      segments[(*segments_num)].type = VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP;
    else if (loops[start - 1] == 'i')
      segments[(*segments_num)].type = VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP;
    else if (loops[start - 1] == 'm')
      segments[(*segments_num)].type = VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP;

    (*segments_num)++;

    /* resize memory if necessary */
    if ((*segments_num) == segments_length) {
      segments_length *= 1.4;
      segments        = (struct binding_segment *)vrna_realloc(segments,
                                                               sizeof(struct binding_segment) *
                                                               segments_length);
    }
  }

  segments = (struct binding_segment *)vrna_realloc(segments,
                                                    sizeof(struct binding_segment) *
                                                    (*segments_num));

  free(loops);

  return segments;
}


PRIVATE void
fill_MFE_matrix(vrna_fold_compound_t  *fc,
                int                   *mx,
                unsigned int          from,
                unsigned int          to,
                unsigned int          type)
{
  unsigned int  i, d, m, u;
  int           e, ee;
  vrna_ud_t     *domains_up;

  domains_up = fc->domains_up;

  e = 0;
  for (m = 0; m < (unsigned int)domains_up->uniq_motif_count; m++) {
    if (domains_up->uniq_motif_size[m] == 1) {
      ee = domains_up->energy_cb(fc,
                                 to,
                                 to,
                                 type | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                 domains_up->data);
      e = MIN2(e, ee);
    }
  }

  mx[to] = e;

  for (d = 2, i = to - 1; i >= from; i--, d++) {
    e = mx[i + 1];

    for (m = 0; m < (unsigned int)domains_up->uniq_motif_count; m++) {
      u = domains_up->uniq_motif_size[m];
      if (u <= d) {
        ee = domains_up->energy_cb(fc,
                                   i,
                                   i + u - 1,
                                   type | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                   domains_up->data);

        if (u < d)
          ee += mx[i + u];

        e = MIN2(e, ee);
      }
    }
    mx[i] = e;
  }
}


PRIVATE vrna_ud_motif_t *
backtrack_MFE_matrix(vrna_fold_compound_t *fc,
                     int                  *mx,
                     unsigned int         from,
                     unsigned int         to,
                     unsigned int         type)
{
  unsigned int    d, i, k, m, u, motif_cnt, motif_size;
  int             e, ee, eee;
  vrna_ud_t       *domains_up;
  vrna_ud_motif_t *motif_list;

  domains_up  = fc->domains_up;
  motif_cnt   = 0;
  motif_size  = 10;
  motif_list  = (vrna_ud_motif_t *)vrna_alloc(sizeof(vrna_ud_motif_t) * (motif_size + 1));

  for (d = to - from + 1, i = from; i < to;) {
    e   = mx[i];
    ee  = mx[i + 1];

    if (e == ee) {
      i++;
      d--;
      continue;
    }

    for (m = 0; m < (unsigned int)domains_up->uniq_motif_count; m++) {
      u = domains_up->uniq_motif_size[m];

      if (u <= d) {
        eee = ee = domains_up->energy_cb(fc,
                                         i,
                                         i + u - 1,
                                         type | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                         domains_up->data);

        if (ee != INF) {
          if (u < d)
            ee += mx[i + u];

          if (e == ee) {
            /* determine actual motif number and add this motif to list */
            for (k = 0; k < (unsigned int)domains_up->motif_count; k++)
              if ((domains_up->motif_type[k] & type) && (domains_up->motif_size[k] == u))
                if (eee == (int)roundf(domains_up->motif_en[k] * 100.))
                  break;

            motif_list[motif_cnt].start   = i;
            motif_list[motif_cnt].number  = k;
            motif_cnt++;

            if (motif_cnt == motif_size) {
              motif_size  *= 1.4;
              motif_list  = (vrna_ud_motif_t *)vrna_realloc(motif_list,
                                                            sizeof(vrna_ud_motif_t) *
                                                            (motif_size + 1));
            }

            i += u;
            d -= u;
            break;
          }
        }
      }
    }
  }

  if (i == to) {
    e = mx[i];

    if (e != 0) {
      for (m = 0; m < (unsigned int)domains_up->uniq_motif_count; m++) {
        if (domains_up->uniq_motif_size[m] == 1) {
          ee = domains_up->energy_cb(fc,
                                     i,
                                     i,
                                     type | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                     domains_up->data);

          if (e == ee) {
            /* determine actual motif number and add this motif to list */
            for (k = 0; k < (unsigned int)domains_up->motif_count; k++)
              if ((domains_up->motif_type[k] & type) && (domains_up->motif_size[k] == 1))
                if (ee == (int)roundf(domains_up->motif_en[k] * 100.))
                  break;

            motif_list[motif_cnt].start   = i;
            motif_list[motif_cnt].number  = k;
            motif_cnt++;

            if (motif_cnt == motif_size) {
              motif_size  *= 1.4;
              motif_list  = (vrna_ud_motif_t *)vrna_realloc(motif_list,
                                                            sizeof(vrna_ud_motif_t) *
                                                            (motif_size + 1));
            }

            break;
          }
        }
      }
    }
  }

  if (motif_cnt == 0) {
    free(motif_list);
    motif_list = NULL;
  } else {
    motif_list = (vrna_ud_motif_t *)vrna_realloc(motif_list,
                                                 sizeof(vrna_ud_motif_t) *
                                                 (motif_cnt + 1));
    motif_list[motif_cnt].start   = 0;
    motif_list[motif_cnt].number  = -1;
  }

  return motif_list;
}


PRIVATE vrna_ud_motif_t **
backtrack_MFE_matrix_exhaustive(vrna_fold_compound_t  *fc,
                                int                   *mx,
                                unsigned int          from,
                                unsigned int          to,
                                unsigned int          type)
{
  vrna_ud_motif_t     **motif_lists, *local_list;
  vrna_ud_t           *domains_up;

  domains_up = fc->domains_up;

  unsigned int        i, k, m, u, motif_list_size, motif_list_num, bt_stack_size, bt_stack_pos,
                      local_cnt, local_size;
  int                 e, ee;
  struct ud_bt_stack  *bt_stack;

  /* prepare motif list */
  motif_list_size = 10;
  motif_list_num  = 0;
  motif_lists     = (vrna_ud_motif_t **)vrna_alloc(sizeof(vrna_ud_motif_t *) *
                                                   (motif_list_size + 1));

  /* prepare backtrack stack */
  bt_stack_pos  = 0;
  bt_stack_size = 10;
  bt_stack      = (struct ud_bt_stack *)vrna_alloc(sizeof(struct ud_bt_stack) * bt_stack_size);

  /* push start condition to backtrack stack */
  bt_stack[bt_stack_pos].from       = from;
  bt_stack[bt_stack_pos].motif_size = 10;
  bt_stack[bt_stack_pos].motifs     = (vrna_ud_motif_t *)vrna_alloc(sizeof(vrna_ud_motif_t) * 10);
  bt_stack[bt_stack_pos].motif_cnt  = 0;
  bt_stack_pos++;

  /* process backtrack stack */
  while (bt_stack_pos > 0) {
    /* pop condition from stack */
    bt_stack_pos--;
    i           = bt_stack[bt_stack_pos].from;
    local_list  = bt_stack[bt_stack_pos].motifs;
    local_cnt   = bt_stack[bt_stack_pos].motif_cnt;
    local_size  = bt_stack[bt_stack_pos].motif_size;

    if (i > to) {
      /* store result */
      if (local_list) {
        local_list = (vrna_ud_motif_t *)vrna_realloc(local_list,
                                                     sizeof(vrna_ud_motif_t *) *
                                                     (local_cnt + 1));
        local_list[local_cnt].start   = 0;
        local_list[local_cnt].number  = -1;

        motif_lists[motif_list_num++] = local_list;

        if (motif_list_num == motif_list_size) {
          motif_list_size *= 1.4;
          motif_lists     = (vrna_ud_motif_t **)vrna_realloc(motif_lists,
                                                             sizeof(vrna_ud_motif_t *) *
                                                             (motif_list_size + 1));
        }
      }

      continue;
    }

    /* backtrack motifs */
    e = mx[i];

    /* nibble off unpaired nucleotides at 5' end */
    while ((i + 1 <= to) && (mx[i + 1] == e))
      i++;

    /* detect motif */
    for (k = 0; k < (unsigned int)domains_up->uniq_motif_count; k++) {
      u = domains_up->uniq_motif_size[k];
      if (i + u - 1 <= to) {
        ee = domains_up->energy_cb(fc,
                                   i,
                                   i + u - 1,
                                   type | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                   domains_up->data);
        if (e == ee) {
          /* clone current list */
          vrna_ud_motif_t *ptr = (vrna_ud_motif_t *)vrna_alloc(sizeof(vrna_ud_motif_t) *
                                                               (local_cnt + 2));
          memcpy(ptr, local_list, sizeof(vrna_ud_motif_t) * local_cnt);

          /* determine actual motif number and add this motif to list */
          for (m = 0; m < (unsigned int)domains_up->uniq_motif_count; m++)
            if ((domains_up->motif_type[m] & type) && (domains_up->motif_size[m] == u)) {
              if (ee == (int)roundf(domains_up->motif_en[m] * 100.)) {
                k = m;
                break;
              }
            }

          ptr[local_cnt].start  = i;
          ptr[local_cnt].number = m;

          /* push back to stack */
          bt_stack[bt_stack_pos].from       = to + 1; /* mark end of backtracking */
          bt_stack[bt_stack_pos].motifs     = ptr;
          bt_stack[bt_stack_pos].motif_cnt  = local_cnt + 1;
          bt_stack[bt_stack_pos].motif_size = local_cnt + 2;
          bt_stack_pos++;
        }

        if (i + u - 1 < to) {
          if (e == ee + mx[i + u]) {
            /* clone current list */
            vrna_ud_motif_t *ptr = (vrna_ud_motif_t *)vrna_alloc(sizeof(vrna_ud_motif_t) *
                                                                 (local_cnt + local_size));
            memcpy(ptr, local_list, sizeof(vrna_ud_motif_t) * local_cnt);

            /* determine actual motif number and add this motif to list */
            for (m = 0; m < (unsigned int)domains_up->uniq_motif_count; m++)
              if ((domains_up->motif_type[m] & type) && (domains_up->motif_size[m] == u)) {
                if (ee == (int)roundf(domains_up->motif_en[m] * 100.)) {
                  k = m;
                  break;
                }
              }

            ptr[local_cnt].start  = i;
            ptr[local_cnt].number = m;

            /* push back to stack */
            bt_stack[bt_stack_pos].from       = i + u;
            bt_stack[bt_stack_pos].motifs     = ptr;
            bt_stack[bt_stack_pos].motif_cnt  = local_cnt + 1;
            bt_stack[bt_stack_pos].motif_size = local_cnt + local_size;
            bt_stack_pos++;
          }
        }
      }
    }

    free(local_list);
  }

  if (motif_list_num == 0) {
    free(motif_lists);
    motif_lists = NULL;
  } else {
    motif_lists = (vrna_ud_motif_t **)vrna_realloc(motif_lists,
                                                   sizeof(vrna_ud_motif_t *) *
                                                   (motif_list_num + 1));
    motif_lists[motif_list_num] = NULL;
  }

  free(bt_stack);

  return motif_lists;
}


PRIVATE void
fill_MEA_matrix(vrna_fold_compound_t  *fc,
                float                 *mx,
                unsigned int          from,
                unsigned int          to,
                float                 *pu,
                unsigned int          type)
{
  unsigned int  i, d, m, u;
  float         ea, p;
  vrna_ud_t     *domains_up;

  domains_up = fc->domains_up;

  ea = pu[to];

  for (m = 0; m < (unsigned int)domains_up->motif_count; m++) {
    if (!(type & domains_up->motif_type[m]))
      continue;

    if (domains_up->motif_size[m] == 1) {
      p   = domains_up->probs_get(fc, to, to, type, m, domains_up->data);
      ea  = MAX2(ea, p);
    }
  }

  mx[to] = ea;

  for (d = 2, i = to - 1; i >= from; i--, d++) {
    ea = mx[i + 1] + pu[i];

    for (m = 0; m < (unsigned int)domains_up->motif_count; m++) {
      if (!(type & domains_up->motif_type[m]))
        continue;

      u = domains_up->motif_size[m];
      if (u <= d) {
        p = domains_up->probs_get(fc, i, i + u - 1, type, m, domains_up->data);
        if (p > 0) {
          p *= u;

          if (u < d)
            p += mx[i + u];

          ea = MAX2(ea, p);
        }
      }
    }
    mx[i] = ea;
  }
}


PRIVATE vrna_ud_motif_t *
backtrack_MEA_matrix(vrna_fold_compound_t *fc,
                     float                *mx,
                     unsigned int         from,
                     unsigned int         to,
                     float                *pu,
                     unsigned int         type)
{
  unsigned int    i, u, m, d, motif_cnt, motif_size, found;
  float           mea, p, prec;
  vrna_ud_t       *domains_up;
  vrna_ud_motif_t *motif_list;

  domains_up  = fc->domains_up;
  motif_cnt   = 0;
  motif_size  = 10;
  motif_list  = (vrna_ud_motif_t *)vrna_alloc(sizeof(vrna_ud_motif_t) * (motif_size + 1));

  for (d = to - from + 1, i = from; i <= to;) {
    prec  = FLT_EPSILON * mx[i];
    mea   = mx[i];
    p     = pu[i];
    found = 0;

    if (i < to)
      p += mx[i + 1];

    if (mea <= p + prec) {
      /* nibble-off unpaired nucleotides */
      i++;
      d--;
      continue;
    }

    for (m = 0; m < (unsigned int)domains_up->motif_count; m++) {
      if (!(type & domains_up->motif_type[m]))
        continue;

      u = domains_up->motif_size[m];
      if (u <= d) {
        p = domains_up->probs_get(fc, i, i + u - 1, type, m, domains_up->data);
        if (p > 0.) {
          p *= u;

          if (u < d)
            p += mx[i + u];

          if (mea <= p + prec) {
            motif_list[motif_cnt].start   = i;
            motif_list[motif_cnt].number  = m;
            motif_cnt++;

            if (motif_cnt == motif_size) {
              motif_size  *= 1.4;
              motif_list  = (vrna_ud_motif_t *)vrna_realloc(motif_list,
                                                            sizeof(vrna_ud_motif_t) *
                                                            (motif_size + 1));
            }

            i += u;
            d -= u;

            found = 1;
            break;
          }
        }
      }
    }

    if (!found) {
      vrna_log_warning("Backtracking failed in unstructured domains MEA\n");
      motif_cnt = 0;
      break;
    }
  }

  if (motif_cnt == 0) {
    free(motif_list);
    motif_list = NULL;
  } else {
    motif_list = (vrna_ud_motif_t *)vrna_realloc(motif_list,
                                                 sizeof(vrna_ud_motif_t) *
                                                 (motif_cnt + 1));
    motif_list[motif_cnt].start   = 0;
    motif_list[motif_cnt].number  = -1;
  }

  return motif_list;
}


/*
 **********************************
 * Default implementation for
 * ligand binding to unpaired
 * stretches follows below
 **********************************
 */
PRIVATE void
add_ligand_motif(vrna_fold_compound_t *vc,
                 const char           *motif,
                 double               motif_en,
                 const char           *motif_name,
                 unsigned int         loop_type)
{
  unsigned int  i, n, same_size;
  vrna_ud_t     *ud;

  n   = (unsigned int)strlen(motif);
  ud  = vc->domains_up;

  /* First, we update the list of unique motif lengths */
  for (same_size = i = 0; i < (unsigned int)ud->uniq_motif_count; i++) {
    if (ud->uniq_motif_size[i] == n) {
      same_size = 1;
      break;
    }
  }

  if (!same_size) {
    ud->uniq_motif_size = (unsigned int *)vrna_realloc(ud->uniq_motif_size,
                                                       sizeof(unsigned int *) *
                                                       (ud->uniq_motif_count + 1));
    ud->uniq_motif_size[ud->uniq_motif_count] = n;
    ud->uniq_motif_count++;
  }

  /* And finally, we store the motif */
  ud->motif = (char **)vrna_realloc(ud->motif,
                                    sizeof(char *) *
                                    (ud->motif_count + 1));
  ud->motif[ud->motif_count] = strdup(motif);

  ud->motif_name = (char **)vrna_realloc(ud->motif_name,
                                         sizeof(char *) *
                                         (ud->motif_count + 1));
  ud->motif_name[ud->motif_count] = (motif_name) ? strdup(motif) : NULL;

  ud->motif_size = (unsigned int *)vrna_realloc(ud->motif_size,
                                                sizeof(unsigned int *) *
                                                (ud->motif_count + 1));
  ud->motif_size[ud->motif_count] = n;

  ud->motif_en = (double *)vrna_realloc(ud->motif_en,
                                        sizeof(double) *
                                        (ud->motif_count + 1));
  ud->motif_en[ud->motif_count] = motif_en;

  ud->motif_type = (unsigned int *)vrna_realloc(ud->motif_type,
                                                sizeof(double) *
                                                (ud->motif_count + 1));
  ud->motif_type[ud->motif_count] = loop_type;

  ud->motif_count++;
}


PRIVATE void
remove_default_data(void *d)
{
  struct ligands_up_data_default *data;

  data = (struct ligands_up_data_default *)d;

  free_default_data_matrices(data);
  free_default_data_exp_matrices(data);
  free_default_data(data);

  free(data);
}


PRIVATE void
free_default_data(struct ligands_up_data_default *data)
{
  int i;

  if (data->motif_list_ext) {
    for (i = 0; i <= data->n; i++)
      free(data->motif_list_ext[i]);
    free(data->motif_list_ext);
  }

  if (data->motif_list_hp) {
    for (i = 0; i <= data->n; i++)
      free(data->motif_list_hp[i]);
    free(data->motif_list_hp);
  }

  if (data->motif_list_int) {
    for (i = 0; i <= data->n; i++)
      free(data->motif_list_int[i]);
    free(data->motif_list_int);
  }

  if (data->motif_list_mb) {
    for (i = 0; i <= data->n; i++)
      free(data->motif_list_mb[i]);
    free(data->motif_list_mb);
  }

  free(data->len);
  free(data->dG);
  free(data->exp_dG);
}


PRIVATE void
free_default_data_matrices(struct ligands_up_data_default *data)
{
  /* the following four pointers may point to the same memory */
  if (data->energies_ext) {
    /* check whether one of the other b* points to the same memory location */
    if (data->energies_ext == data->energies_hp)
      data->energies_hp = NULL;

    if (data->energies_ext == data->energies_int)
      data->energies_int = NULL;

    if (data->energies_ext == data->energies_mb)
      data->energies_mb = NULL;

    free(data->energies_ext);
    data->energies_ext = NULL;
  }

  if (data->energies_hp) {
    /* check whether one of the other b* points to the same memory location */
    if (data->energies_hp == data->energies_int)
      data->energies_int = NULL;

    if (data->energies_hp == data->energies_mb)
      data->energies_mb = NULL;

    free(data->energies_hp);
    data->energies_hp = NULL;
  }

  if (data->energies_int) {
    /* check whether one of the other b* points to the same memory location */
    if (data->energies_int == data->energies_mb)
      data->energies_mb = NULL;

    free(data->energies_int);
    data->energies_int = NULL;
  }

  free(data->energies_mb);
  data->energies_mb = NULL;
}


PRIVATE void
free_default_data_exp_matrices(struct ligands_up_data_default *data)
{
  int i;

  /* the following four pointers may point to the same memory */
  if (data->exp_energies_ext) {
    /* check whether one of the other b* points to the same memory location */
    if (data->exp_energies_ext == data->exp_energies_hp)
      data->exp_energies_hp = NULL;

    if (data->exp_energies_ext == data->exp_energies_int)
      data->exp_energies_int = NULL;

    if (data->exp_energies_ext == data->exp_energies_mb)
      data->exp_energies_mb = NULL;

    free(data->exp_energies_ext);
    data->exp_energies_ext = NULL;
  }

  if (data->exp_energies_hp) {
    /* check whether one of the other b* points to the same memory location */
    if (data->exp_energies_hp == data->exp_energies_int)
      data->exp_energies_int = NULL;

    if (data->exp_energies_hp == data->exp_energies_mb)
      data->exp_energies_mb = NULL;

    free(data->exp_energies_hp);
    data->exp_energies_hp = NULL;
  }

  if (data->exp_energies_int) {
    /* check whether one of the other b* points to the same memory location */
    if (data->exp_energies_int == data->exp_energies_mb)
      data->exp_energies_mb = NULL;

    free(data->exp_energies_int);
    data->exp_energies_int = NULL;
  }

  free(data->exp_energies_mb);
  data->exp_energies_mb = NULL;

  if (data->outside_ext)
    for (i = 0; i <= data->n; i++)
      if (data->outside_ext[i])
        free(data->outside_ext[i]);

  free(data->outside_ext);
  free(data->outside_ext_count);

  if (data->outside_hp)
    for (i = 0; i <= data->n; i++)
      if (data->outside_hp[i])
        free(data->outside_hp[i]);

  free(data->outside_hp);
  free(data->outside_hp_count);

  if (data->outside_int)
    for (i = 0; i <= data->n; i++)
      if (data->outside_int[i])
        free(data->outside_int[i]);

  free(data->outside_int);
  free(data->outside_int_count);

  if (data->outside_mb)
    for (i = 0; i <= data->n; i++)
      if (data->outside_mb[i])
        free(data->outside_mb[i]);

  free(data->outside_mb);
  free(data->outside_mb_count);
}


PRIVATE int *
get_motifs(vrna_fold_compound_t *vc,
           int                  i,
           unsigned int         loop_type)
{
  int       k, j, u, n, *motif_list, cnt, guess;
  char      *sequence;
  vrna_ud_t *domains_up;

  sequence    = vc->sequence;
  n           = (int)vc->length;
  domains_up  = vc->domains_up;

  cnt         = 0;
  guess       = domains_up->motif_count;
  motif_list  = (int *)vrna_alloc(sizeof(int) * (guess + 1));

  /* collect list of motif numbers we find that start at position i */
  for (k = 0; k < domains_up->motif_count; k++) {
    if (!(domains_up->motif_type[k] & loop_type))
      continue;

    j = i + domains_up->motif_size[k] - 1;
    if (j <= n) {
      /* only consider motif that does not exceed sequence length (does not work for circular RNAs!) */
      for (u = i; u <= j; u++)
        if (!vrna_nucleotide_IUPAC_identity(sequence[u - 1], domains_up->motif[k][u - i]))
          break;

      if (u > j) /* got a complete motif match */
        motif_list[cnt++] = k;
    }
  }

  if (cnt == 0) {
    free(motif_list);
    return NULL;
  }

  motif_list      = (int *)vrna_realloc(motif_list, sizeof(int) * (cnt + 1));
  motif_list[cnt] = -1; /* end of list marker */

  return motif_list;
}


static void
annotate_ud(vrna_fold_compound_t  *vc,
            int                   start,
            int                   end,
            char                  l,
            vrna_ud_motif_t       **list,
            int                   *list_size,
            int                   *list_pos)
{
  int i, j;

  /* get motifs in segment [start,end] */
  for (i = start; i <= end; i++) {
    unsigned int type = 0;
    switch (l) {
      case 'e':
        type = VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP;
        break;
      case 'h':
        type = VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP;
        break;
      case 'i':
        type = VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP;
        break;
      case 'm':
        type = VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP;
        break;
    }

    int *m = vrna_ud_get_motifs_at(vc, i, type);
    if (m) {
      for (j = 0; m[j] != -1; j++) {
        int size = vc->domains_up->motif_size[m[j]];

        if (i + size - 1 <= end) {
          if (*list_pos == *list_size) {
            *list_size  *= 1.2;
            *list       = (vrna_ud_motif_t *)vrna_realloc(*list,
                                                          sizeof(vrna_ud_motif_t) *
                                                          (*list_size));
          }

          (*list)[*list_pos].start  = i;
          (*list)[*list_pos].number = m[j];
          (*list_pos)++;
        }
      }
    }

    free(m);
  }
}


PRIVATE void
prepare_matrices(vrna_fold_compound_t           *vc,
                 struct ligands_up_data_default *data)
{
  int       i, j, k, n, size;
  vrna_ud_t *domains_up;

  n           = (int)vc->length;
  size        = ((n + 1) * (n + 2)) / 2 + 1;
  domains_up  = vc->domains_up;

  free_default_data_matrices(data);

  /* here we save memory by re-using DP matrices */
  unsigned int  lt[4] = {
    VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP,
    VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
    VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
    VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP
  };
  int           **m[4], *mx;

  m[0]  = &data->energies_ext;
  m[1]  = &data->energies_hp;
  m[2]  = &data->energies_int;
  m[3]  = &data->energies_mb;

  for (i = 0; i < 4; i++) {
    unsigned int *col, *col2;
    if (*(m[i]))
      continue;

    mx      = (int *)vrna_alloc(sizeof(int) * size);
    col     = (unsigned int *)vrna_alloc(sizeof(unsigned int) * domains_up->motif_count);
    col2    = (unsigned int *)vrna_alloc(sizeof(unsigned int) * domains_up->motif_count);
    *(m[i]) = mx;

    for (k = 0; k < domains_up->motif_count; k++)
      col[k] = domains_up->motif_type[k] & lt[i];

    /* check if any of the remaining DP matrices can point to the same location */
    for (j = i + 1; j < 4; j++) {
      for (k = 0; k < domains_up->motif_count; k++) {
        col2[k] = domains_up->motif_type[k] & lt[j];
        if (col[k] != col2[k])
          break;
      }
      if (k == domains_up->motif_count)
        *(m[j]) = mx;
    }

    free(col);
    free(col2);
  }
}


PRIVATE void
prepare_exp_matrices(vrna_fold_compound_t           *vc,
                     struct ligands_up_data_default *data)
{
  int       i, j, k, n, size;
  vrna_ud_t *domains_up;

  n           = (int)vc->length;
  size        = ((n + 1) * (n + 2)) / 2 + 1;
  domains_up  = vc->domains_up;

  free_default_data_exp_matrices(data);

  /* here we save memory by re-using DP matrices */
  unsigned int  lt[4] = {
    VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP,
    VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
    VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
    VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP
  };
  FLT_OR_DBL    **m[4], *mx;

  m[0]  = &data->exp_energies_ext;
  m[1]  = &data->exp_energies_hp;
  m[2]  = &data->exp_energies_int;
  m[3]  = &data->exp_energies_mb;

  for (i = 0; i < 4; i++) {
    unsigned int *col, *col2;
    if (*(m[i]))
      continue;

    mx      = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);
    col     = (unsigned int *)vrna_alloc(sizeof(unsigned int) * domains_up->motif_count);
    col2    = (unsigned int *)vrna_alloc(sizeof(unsigned int) * domains_up->motif_count);
    *(m[i]) = mx;

    for (k = 0; k < domains_up->motif_count; k++)
      col[k] = domains_up->motif_type[k] & lt[i];

    /* check if any of the remaining DP matrices can point to the same location */
    for (j = i + 1; j < 4; j++) {
      for (k = 0; k < domains_up->motif_count; k++) {
        col2[k] = domains_up->motif_type[k] & lt[j];
        if (col[k] != col2[k])
          break;
      }
      if (k == domains_up->motif_count)
        *(m[j]) = mx;
    }

    free(col);
    free(col2);
  }

  /* now prepate memory for outside partition function */
  data->outside_ext = (struct default_outside **)vrna_alloc(
    sizeof(struct default_outside *) * (n + 2));
  data->outside_hp = (struct default_outside **)vrna_alloc(
    sizeof(struct default_outside *) * (n + 2));
  data->outside_int = (struct default_outside **)vrna_alloc(
    sizeof(struct default_outside *) * (n + 2));
  data->outside_mb = (struct default_outside **)vrna_alloc(
    sizeof(struct default_outside *) * (n + 2));
  data->outside_ext_count = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n + 2));
  data->outside_hp_count  = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n + 2));
  data->outside_int_count = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n + 2));
  data->outside_mb_count  = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n + 2));
}


PRIVATE void
prepare_default_data(vrna_fold_compound_t           *vc,
                     struct ligands_up_data_default *data)
{
  int       i, n;
  vrna_ud_t *domains_up;

  n           = (int)vc->length;
  domains_up  = vc->domains_up;

  data->n = n;
  free_default_data(data);

  /*
   *  create motif_list for associating a nucleotide position with all
   *  motifs that start there
   */
  data->motif_list_ext    = (int **)vrna_alloc(sizeof(int *) * (n + 1));
  data->motif_list_hp     = (int **)vrna_alloc(sizeof(int *) * (n + 1));
  data->motif_list_int    = (int **)vrna_alloc(sizeof(int *) * (n + 1));
  data->motif_list_mb     = (int **)vrna_alloc(sizeof(int *) * (n + 1));
  data->motif_list_ext[0] = NULL;
  data->motif_list_hp[0]  = NULL;
  data->motif_list_int[0] = NULL;
  data->motif_list_mb[0]  = NULL;
  for (i = 1; i <= n; i++) {
    data->motif_list_ext[i] = get_motifs(vc, i, VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP);
    data->motif_list_hp[i]  = get_motifs(vc, i, VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP);
    data->motif_list_int[i] = get_motifs(vc, i, VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP);
    data->motif_list_mb[i]  = get_motifs(vc, i, VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP);
  }

  data->default_cb[VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP] = default_exp_energy_ext_motif;
  data->default_cb[VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP]  = default_exp_energy_hp_motif;
  data->default_cb[VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP] = default_exp_energy_int_motif;
  data->default_cb[VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP]  = default_exp_energy_mb_motif;

  /*  store length of motifs in 'data' */
  data->len = (int *)vrna_alloc(sizeof(int) * domains_up->motif_count);
  for (i = 0; i < domains_up->motif_count; i++)
    data->len[i] = domains_up->motif_size[i];

  /*  precompute energy contributions of the motifs */
  data->dG = (int *)vrna_alloc(sizeof(int) * domains_up->motif_count);
  for (i = 0; i < domains_up->motif_count; i++)
    data->dG[i] = (int)roundf(domains_up->motif_en[i] * 100.);
}


PRIVATE void
default_prod_rule(vrna_fold_compound_t  *vc,
                  void                  *d)
{
  int                             i, j, k, l, u, n, e_ext, e_hp, e_int, e_mb, en, en2, *idx;
  struct ligands_up_data_default  *data;

  int                             *energies_ext;
  int                             *energies_hp;
  int                             *energies_int;
  int                             *energies_mb;

  n     = (int)vc->length;
  idx   = vc->jindx;
  data  = (struct ligands_up_data_default *)d;

  prepare_default_data(vc, data);
  prepare_matrices(vc, data);

  energies_ext  = data->energies_ext;
  energies_hp   = data->energies_hp;
  energies_int  = data->energies_int;
  energies_mb   = data->energies_mb;

  /* now we can start to fill the DP matrices */
  for (i = n; i > 0; i--) {
    int *list_ext = data->motif_list_ext[i];
    int *list_hp  = data->motif_list_hp[i];
    int *list_int = data->motif_list_int[i];
    int *list_mb  = data->motif_list_mb[i];
    for (j = i; j <= n; j++) {
      if (i < j) {
        e_ext = energies_ext[idx[j] + i + 1];
        e_hp  = energies_hp[idx[j] + i + 1];
        e_int = energies_int[idx[j] + i + 1];
        e_mb  = energies_mb[idx[j] + i + 1];
      } else {
        e_ext = INF;
        e_hp  = INF;
        e_int = INF;
        e_mb  = INF;
      }

      if (list_ext) {
        for (k = 0; -1 != (l = list_ext[k]); k++) {
          u   = i + data->len[l] - 1;
          en  = data->dG[l];
          if (u <= j) {
            e_ext = MIN2(e_ext, en);
            if (u < j) {
              en2   = en + energies_ext[idx[j] + u + 1];
              e_ext = MIN2(e_ext, en2);
            }
          }
        }
      }

      if (list_hp) {
        for (k = 0; -1 != (l = list_hp[k]); k++) {
          u   = i + data->len[l] - 1;
          en  = data->dG[l];
          if (u <= j) {
            e_hp = MIN2(e_hp, en);
            if (u < j) {
              en2   = en + energies_hp[idx[j] + u + 1];
              e_hp  = MIN2(e_hp, en2);
            }
          }
        }
      }

      if (list_int) {
        for (k = 0; -1 != (l = list_int[k]); k++) {
          u   = i + data->len[l] - 1;
          en  = data->dG[l];
          if (u <= j) {
            e_int = MIN2(e_int, en);
            if (u < j) {
              en2   = en + energies_int[idx[j] + u + 1];
              e_int = MIN2(e_int, en2);
            }
          }
        }
      }

      if (list_mb) {
        for (k = 0; -1 != (l = list_mb[k]); k++) {
          u   = i + data->len[l] - 1;
          en  = data->dG[l];
          if (u <= j) {
            e_mb = MIN2(e_mb, en);
            if (u < j) {
              en2   = en + energies_mb[idx[j] + u + 1];
              e_mb  = MIN2(e_mb, en2);
            }
          }
        }
      }

      energies_ext[idx[j] + i]  = e_ext;
      energies_hp[idx[j] + i]   = e_hp;
      energies_int[idx[j] + i]  = e_int;
      energies_mb[idx[j] + i]   = e_mb;
    }
  }
}


PRIVATE void
default_exp_prod_rule(vrna_fold_compound_t  *vc,
                      void                  *d)
{
  int                             i, j, k, l, u, n, *idx;
  FLT_OR_DBL                      q_ext, q_hp, q_int, q_mb, q;
  vrna_ud_t                       *domains_up;
  struct ligands_up_data_default  *data;

  FLT_OR_DBL                      *exp_energies_ext;
  FLT_OR_DBL                      *exp_energies_hp;
  FLT_OR_DBL                      *exp_energies_int;
  FLT_OR_DBL                      *exp_energies_mb;
  double                          kT;

  n           = (int)vc->length;
  idx         = vc->iindx;
  domains_up  = vc->domains_up;
  data        = (struct ligands_up_data_default *)d;
  kT          = vc->exp_params->kT;

  prepare_default_data(vc, data);
  prepare_exp_matrices(vc, data);

  exp_energies_ext  = data->exp_energies_ext;
  exp_energies_hp   = data->exp_energies_hp;
  exp_energies_int  = data->exp_energies_int;
  exp_energies_mb   = data->exp_energies_mb;

  data->exp_e_mx[VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP] = data->exp_energies_ext;
  data->exp_e_mx[VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP]  = data->exp_energies_hp;
  data->exp_e_mx[VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP] = data->exp_energies_int;
  data->exp_e_mx[VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP]  = data->exp_energies_mb;

  /*  precompute energy contributions of the motifs */
  data->exp_dG = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * domains_up->motif_count);
  for (i = 0; i < domains_up->motif_count; i++) {
    double GT = domains_up->motif_en[i] * 1000.; /* in cal/mol */
    data->exp_dG[i] = (FLT_OR_DBL)exp(-GT / kT);
  }

  /* now we can start to fill the DP matrices */
  for (i = n; i > 0; i--) {
    int *list_ext = data->motif_list_ext[i];
    int *list_hp  = data->motif_list_hp[i];
    int *list_int = data->motif_list_int[i];
    int *list_mb  = data->motif_list_mb[i];
    for (j = i; j <= n; j++) {
      if (i < j) {
        q_ext = exp_energies_ext[idx[i + 1] - j];
        q_hp  = exp_energies_hp[idx[i + 1] - j];
        q_int = exp_energies_int[idx[i + 1] - j];
        q_mb  = exp_energies_mb[idx[i + 1] - j];
      } else {
        q_ext = 0;
        q_hp  = 0;
        q_int = 0;
        q_mb  = 0;
      }

      if (list_ext) {
        for (k = 0; -1 != (l = list_ext[k]); k++) {
          u = i + data->len[l] - 1;
          q = data->exp_dG[l];
          if (u <= j) {
            q_ext += q;
            if (u < j)
              q_ext += q * exp_energies_ext[idx[u + 1] - j];
          }
        }
      }

      if (list_hp) {
        for (k = 0; -1 != (l = list_hp[k]); k++) {
          u = i + data->len[l] - 1;
          q = data->exp_dG[l];
          if (u <= j) {
            q_hp += q;
            if (u < j)
              q_hp += q * exp_energies_hp[idx[u + 1] - j];
          }
        }
      }

      if (list_int) {
        for (k = 0; -1 != (l = list_int[k]); k++) {
          u = i + data->len[l] - 1;
          q = data->exp_dG[l];
          if (u <= j) {
            q_int += q;
            if (u < j)
              q_int += q * exp_energies_int[idx[u + 1] - j];
          }
        }
      }

      if (list_mb) {
        for (k = 0; -1 != (l = list_mb[k]); k++) {
          u = i + data->len[l] - 1;
          q = data->exp_dG[l];
          if (u <= j) {
            q_mb += q;
            if (u < j)
              q_mb += q * exp_energies_mb[idx[u + 1] - j];
          }
        }
      }

      exp_energies_ext[idx[i] - j]  = q_ext;
      exp_energies_hp[idx[i] - j]   = q_hp;
      exp_energies_int[idx[i] - j]  = q_int;
      exp_energies_mb[idx[i] - j]   = q_mb;
    }
  }
}


PRIVATE int
default_energy(vrna_fold_compound_t *vc,
               int                  i,
               int                  j,
               unsigned int         loop_type,
               void                 *d)
{
  int                             en, ij, *idx = vc->jindx;
  struct ligands_up_data_default  *data = (struct ligands_up_data_default *)d;

  en  = INF;
  ij  = idx[j] + i;

  if (j < i)
    return INF;

  if (loop_type & VRNA_UNSTRUCTURED_DOMAIN_MOTIF) {
    if (loop_type & VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP)
      en = default_energy_ext_motif(i, j, data);
    else if (loop_type & VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP)
      en = default_energy_hp_motif(i, j, data);
    else if (loop_type & VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP)
      en = default_energy_int_motif(i, j, data);
    else if (loop_type & VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP)
      en = default_energy_mb_motif(i, j, data);
  } else {
    if (loop_type & VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP) {
      if (data->energies_ext)
        en = data->energies_ext[ij];
    } else if (loop_type & VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP) {
      if (data->energies_hp)
        en = data->energies_hp[ij];
    } else if (loop_type & VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP) {
      if (data->energies_int)
        en = data->energies_int[ij];
    } else if (loop_type & VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP) {
      if (data->energies_mb)
        en = data->energies_mb[ij];
    }
  }

  return en;
}


PRIVATE FLT_OR_DBL
default_exp_energy(vrna_fold_compound_t *vc,
                   int                  i,
                   int                  j,
                   unsigned int         loop_type,
                   void                 *d)
{
  FLT_OR_DBL                      q;
  int                             ij, *idx;
  struct ligands_up_data_default  *data;

  q     = 0;
  data  = (struct ligands_up_data_default *)d;

  if (j < i)
    return 0.;

  if (loop_type & VRNA_UNSTRUCTURED_DOMAIN_MOTIF) {
    q = data->default_cb[loop_type & ~(VRNA_UNSTRUCTURED_DOMAIN_MOTIF)](i, j, data);
  } else {
    idx = vc->iindx;
    ij  = idx[i] - j;
    q   = data->exp_e_mx[loop_type][ij];
  }

  return q;
}


PRIVATE int
default_energy_ext_motif(int                            i,
                         int                            j,
                         struct ligands_up_data_default *data)
{
  int k, m;
  int e = INF;

  if (data->motif_list_ext[i]) {
    k = 0;
    while (-1 != (m = data->motif_list_ext[i][k])) {
      if ((i + data->len[m] - 1) == j)
        e = MIN2(e, data->dG[m]);

      k++;
    }
  }

  return e;
}


PRIVATE int
default_energy_hp_motif(int                             i,
                        int                             j,
                        struct ligands_up_data_default  *data)
{
  int k, m;
  int e = INF;

  if (data->motif_list_hp[i]) {
    k = 0;
    while (-1 != (m = data->motif_list_hp[i][k])) {
      if ((i + data->len[m] - 1) == j)
        e = MIN2(e, data->dG[m]);

      k++;
    }
  }

  return e;
}


PRIVATE int
default_energy_int_motif(int                            i,
                         int                            j,
                         struct ligands_up_data_default *data)
{
  int k, m;
  int e = INF;

  if (data->motif_list_int[i]) {
    k = 0;
    while (-1 != (m = data->motif_list_int[i][k])) {
      if ((i + data->len[m] - 1) == j)
        e = MIN2(e, data->dG[m]);

      k++;
    }
  }

  return e;
}


PRIVATE int
default_energy_mb_motif(int                             i,
                        int                             j,
                        struct ligands_up_data_default  *data)
{
  int k, m;
  int e = INF;

  if (data->motif_list_mb[i]) {
    k = 0;
    while (-1 != (m = data->motif_list_mb[i][k])) {
      if ((i + data->len[m] - 1) == j)
        e = MIN2(2, data->dG[m]);

      k++;
    }
  }

  return e;
}


PRIVATE FLT_OR_DBL
default_exp_energy_ext_motif(int                            i,
                             int                            j,
                             struct ligands_up_data_default *data)
{
  int         k, m;
  FLT_OR_DBL  q = 0;

  if (data->motif_list_ext[i]) {
    k = 0;
    while (-1 != (m = data->motif_list_ext[i][k])) {
      if ((i + data->len[m] - 1) == j)
        q += data->exp_dG[m];

      k++;
    }
  }

  return q;
}


PRIVATE FLT_OR_DBL
default_exp_energy_hp_motif(int                             i,
                            int                             j,
                            struct ligands_up_data_default  *data)
{
  int         k, m;
  FLT_OR_DBL  q = 0;

  if (data->motif_list_hp[i]) {
    k = 0;
    while (-1 != (m = data->motif_list_hp[i][k])) {
      if ((i + data->len[m] - 1) == j)
        q += data->exp_dG[m];

      k++;
    }
  }

  return q;
}


PRIVATE FLT_OR_DBL
default_exp_energy_int_motif(int                            i,
                             int                            j,
                             struct ligands_up_data_default *data)
{
  int         k, m;
  FLT_OR_DBL  q = 0;

  if (data->motif_list_int[i]) {
    k = 0;
    while (-1 != (m = data->motif_list_int[i][k])) {
      if ((i + data->len[m] - 1) == j)
        q += data->exp_dG[m];

      k++;
    }
  }

  return q;
}


PRIVATE FLT_OR_DBL
default_exp_energy_mb_motif(int                             i,
                            int                             j,
                            struct ligands_up_data_default  *data)
{
  int         k, m;
  FLT_OR_DBL  q = 0;

  if (data->motif_list_mb[i]) {
    k = 0;
    while (-1 != (m = data->motif_list_mb[i][k])) {
      if ((i + data->len[m] - 1) == j)
        q += data->exp_dG[m];

      k++;
    }
  }

  return q;
}


PRIVATE void
default_probs_add(vrna_fold_compound_t  *vc,
                  int                   i,
                  int                   j,
                  unsigned int          loop_type,
                  FLT_OR_DBL            exp_energy,
                  void                  *data)
{
  int                             **motif_list, k, l, m;
  unsigned int                    *size, *cnt, o;
  struct ligands_up_data_default  *d;
  struct default_outside          **storage, **st;

  d = (struct ligands_up_data_default *)data;

  if (loop_type & VRNA_UNSTRUCTURED_DOMAIN_MOTIF) {
    if (j < i)
      return;

    if (loop_type & VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP) {
      motif_list  = d->motif_list_ext;
      storage     = &(d->outside_ext[i]);
      size        = &(d->outside_ext_count[i]);
    } else if (loop_type & VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP) {
      motif_list  = d->motif_list_hp;
      storage     = &(d->outside_hp[i]);
      size        = &(d->outside_hp_count[i]);
    } else if (loop_type & VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP) {
      motif_list  = d->motif_list_int;
      storage     = &(d->outside_int[i]);
      size        = &(d->outside_int_count[i]);
    } else if (loop_type & VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP) {
      motif_list  = d->motif_list_mb;
      storage     = &(d->outside_mb[i]);
      size        = &(d->outside_mb_count[i]);
    } else {
      vrna_log_warning("Unknown unstructured domain loop type");
      return;
    }

    k = 0;
    while (-1 != (m = motif_list[i][k])) {
      if ((i + d->len[m] - 1) == j) {
        /* check for addition first */
        for (o = 0; o < *size; o++)
          if ((*storage)[o].motif_num == m) {
            /* found previously added motif constribution */
            (*storage)[o].exp_energy += exp_energy;
            break;
          }

        /* if we haven't added yet, create new list entry */
        if (o == *size) {
          *storage =
            (struct default_outside *)vrna_realloc(*storage,
                                                   sizeof(struct default_outside) * (*size + 1));
          (*storage)[*size].motif_num   = m;
          (*storage)[*size].exp_energy  = exp_energy;
          (*size)++;
        }
      }

      k++;
    }
  } else {
    if (j < i)
      return;

    FLT_OR_DBL pf, exp_e;
    pf = default_exp_energy(vc, i, j, loop_type, data);

    if (loop_type & VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP) {
      motif_list  = d->motif_list_ext;
      storage     = d->outside_ext;
      size        = d->outside_ext_count;
    } else if (loop_type & VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP) {
      motif_list  = d->motif_list_hp;
      storage     = d->outside_hp;
      size        = d->outside_hp_count;
    } else if (loop_type & VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP) {
      motif_list  = d->motif_list_int;
      storage     = d->outside_int;
      size        = d->outside_int_count;
    } else if (loop_type & VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP) {
      motif_list  = d->motif_list_mb;
      storage     = d->outside_mb;
      size        = d->outside_mb_count;
    } else {
      vrna_log_warning("Unknown unstructured domain loop type");
      return;
    }

    /* check for each motif starting at any k with i <= k <= j */
    for (k = i; k <= j; k++) {
      if (motif_list[k]) {
        st  = &(storage[k]);
        cnt = &(size[k]);
        for (l = 0; motif_list[k][l] != -1; l++) {
          m = motif_list[k][l];

          if (j < k + d->len[m] - 1) /* motifs must be sorted be length */
            continue;

          exp_e = d->exp_dG[m];
          FLT_OR_DBL p = exp_e / pf;

          /* add/insert contribution */

          /* check for addition first */
          for (o = 0; o < *cnt; o++)
            if ((*st)[o].motif_num == m) {
              /* found previously added motif constribution */
              (*st)[o].exp_energy += p * exp_energy;
              break;
            }

          /* if we haven't added yet, create new list entry */
          if (o == *cnt) {
            *st =
              (struct default_outside *)vrna_realloc(*st,
                                                     sizeof(struct default_outside) * (*cnt + 1));
            (*st)[*cnt].motif_num   = m;
            (*st)[*cnt].exp_energy  = p * exp_energy;
            (*cnt)++;
          }
        }
      }
    }
  }
}


PRIVATE FLT_OR_DBL
default_probs_get(vrna_fold_compound_t  *vc VRNA_UNUSED,
                  int                   i,
                  int                   j,
                  unsigned int          loop_type,
                  int                   motif,
                  void                  *data)
{
  FLT_OR_DBL                      outside = 0.;
  unsigned int                    *size, k;
  struct ligands_up_data_default  *d;
  struct default_outside          **storage;

  d = (struct ligands_up_data_default *)data;

  if (j < i)
    return 0.;

  if (loop_type & VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP) {
    if (d->outside_ext) {
      storage = &(d->outside_ext[i]);
      size    = &(d->outside_ext_count[i]);
      if ((storage) && (*storage)) {
        for (k = 0; k < *size; k++) {
          /* check for motif number match */
          if ((*storage)[k].motif_num == motif) {
            /* check for length match */
            if (i + d->len[motif] - 1 == j)
              outside += (*storage)[k].exp_energy;
          }
        }
      }
    }
  }

  if (loop_type & VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP) {
    if (d->outside_hp) {
      storage = &(d->outside_hp[i]);
      size    = &(d->outside_hp_count[i]);
      if ((storage) && (*storage)) {
        for (k = 0; k < *size; k++) {
          /* check for motif number match */
          if ((*storage)[k].motif_num == motif) {
            /* check for length match */
            if (i + d->len[motif] - 1 == j)
              outside += (*storage)[k].exp_energy;
          }
        }
      }
    }
  }

  if (loop_type & VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP) {
    if (d->outside_int) {
      storage = &(d->outside_int[i]);
      size    = &(d->outside_int_count[i]);
      if ((storage) && (*storage)) {
        for (k = 0; k < *size; k++) {
          /* check for motif number match */
          if ((*storage)[k].motif_num == motif) {
            /* check for length match */
            if (i + d->len[motif] - 1 == j)
              outside += (*storage)[k].exp_energy;
          }
        }
      }
    }
  }

  if (loop_type & VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP) {
    if (d->outside_mb) {
      storage = &(d->outside_mb[i]);
      size    = &(d->outside_mb_count[i]);
      if ((storage) && (*storage)) {
        for (k = 0; k < *size; k++) {
          /* check for motif number match */
          if ((*storage)[k].motif_num == motif) {
            /* check for length match */
            if (i + d->len[motif] - 1 == j)
              outside += (*storage)[k].exp_energy;
          }
        }
      }
    }
  }

  return outside;
}
