/*
 *      Utitlities for plot functions
 *
 *      c  Ronny Lorenz
 *      ViennaRNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>    
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/datastructures/string.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/plotting/utils.h"

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
 # PRIVATE VARIABLES             #
 #################################
 */

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC char **
vrna_annotate_covar_db(const char **alignment,
                       const char *structure,
                       vrna_md_t  *md_p)
{
  return vrna_annotate_covar_db_extended(alignment,
                                         structure,
                                         md_p,
                                         VRNA_BRACKETS_RND);
}


PUBLIC char **
vrna_annotate_covar_db_extended(const char   **alignment,
                                const char   *structure,
                                vrna_md_t    *md_p,
                                unsigned int options)
{
  short int *pt;
  char **result;
  vrna_string_t *res;
  
  pt  = vrna_ptable_from_string(structure, options);
  res = vrna_annotate_covar_pt(alignment, pt, md_p, 2, 0.2);

  result = (char **)vrna_alloc(sizeof(char) * 2);
  result[0] = strdup(res[0]);
  result[1] = strdup(res[1]);

  vrna_string_free(res[0]);
  vrna_string_free(res[1]);
  free(res);
  free(pt);

  return result;
}

PUBLIC vrna_string_t *
vrna_annotate_covar_pt(const char       **alignment,
                       const short int  *pt,
                       vrna_md_t        *md_p,
                       double           threshold,
                       double           min)
{
  /* produce annotation for colored drawings from vrna_file_PS_rnaplot_a() */
  size_t maxl;
  char  pps[64];
  int       i, N, n, s, pairings, a, b;
  double  tabs;
  double  sat_min = 0.2;
  double  sat_pairs[6] = {
    0.0,  /* red    */
    0.16, /* ochre  */
    0.32, /* turquoise */
    0.48, /* green  */
    0.65, /* blue   */
    0.81  /* violet */
  };

  vrna_md_t md;
  vrna_string_t *annotation, ps, colorps;

  if ((!alignment) || (!pt))
    return NULL;

  if (md_p)
    vrna_md_copy(&md, md_p);
  else
    vrna_md_set_default(&md);

  n     = strlen(alignment[0]);
  maxl  = 1024;

  annotation = (vrna_string_t *)vrna_alloc(sizeof(vrna_string_t) * 2);

  ps      = vrna_string_make(NULL);
  colorps = vrna_string_make(NULL);
  
  ps = vrna_string_make_space_for(ps, maxl);
  colorps = vrna_string_make_space_for(colorps, maxl);

  /* count number of sequences in alignment */
  for (N = 0; alignment[N] != NULL; N++);

  if ((min >= 0) && (min < 1))
    sat_min = min;

  if (threshold > (double)N)
    threshold = (double)N;

  if (threshold < 0) {
    /* default for any negative threshold input */
    tabs = 2;
  } else if (trunc(threshold) != threshold) {
    /*
       threshold is expressed as frequency, so
       count number of sequences to convert to
       absolute threshold
    */
    if (threshold < 1) { /* actual frequency */
      tabs = threshold * (double)N;
    } else { /* floating point number > 1. we will truncate this to integer */
      tabs = trunc(threshold);
    }
  } else {
    tabs = threshold;
  }

  snprintf(pps, 64, "0.8 -0.1 %f %f ConsLegend\n",
               threshold, 1. - min);
  colorps = vrna_string_append_cstring(colorps, pps);


  for (i = 1; i <= n; i++) {
    char  ci = '\0', cj = '\0';
    int   j, type, pfreq[8] = {
      0, 0, 0, 0, 0, 0, 0, 0
    }, vi = 0, vj = 0;

    if ((j = pt[i]) < i)
      continue;

    for (s = 0; alignment[s] != NULL; s++) {
      a     = vrna_nucleotide_encode(alignment[s][i - 1], &md);
      b     = vrna_nucleotide_encode(alignment[s][j - 1], &md);
      type  = md.pair[a][b];
      pfreq[type]++;
      if (type) {
        if (alignment[s][i - 1] != ci) {
          ci = alignment[s][i - 1];
          vi++;
        }

        if (alignment[s][j - 1] != cj) {
          cj = alignment[s][j - 1];
          vj++;
        }
      }
    }
    for (pairings = 0, s = 1; s <= 7; s++)
      if (pfreq[s])
        pairings++;

    if ((vrna_string_available_space(ps) < 192) || ((vrna_string_available_space(colorps)) < 64)) {
      maxl    *= 2;
      ps      = vrna_string_make_space_for(ps, maxl);
      colorps = vrna_string_make_space_for(colorps, maxl);
      if ((ps == NULL) || (colorps == NULL)) {
        vrna_log_error("out of memory in realloc");
        free(annotation);
        vrna_string_free(ps);
        vrna_string_free(colorps);
        return NULL;
      }
    }

    if ((pairings > 0) &&
        (pfreq[0] <= (int)tabs)) {
      double intensity = 1.;
      if (pfreq[0] > 0)
        intensity -= ((double)pfreq[0] / tabs) * (1. - sat_min);
      snprintf(pps, 64, "%d %d %.2f %.6f colorpair\n",
               i, j, sat_pairs[pairings - 1], intensity);
      colorps = vrna_string_append_cstring(colorps, pps);
    }

    if (pfreq[0] > 0) {
      snprintf(pps, 64, "%d %d %d gmark\n", i, j, pfreq[0]);
      ps = vrna_string_append_cstring(ps, pps);
    }

    if (vi > 1) {
      snprintf(pps, 64, "%d cmark\n", i);
      ps = vrna_string_append_cstring(ps, pps);
    }

    if (vj > 1) {
      snprintf(pps, 64, "%d cmark\n", j);
      ps = vrna_string_append_cstring(ps, pps);
    }
  }

  annotation[0]  = colorps;
  annotation[1]  = ps;
  return annotation;
}


/* produce info for PS_color_dot_plot */
PUBLIC vrna_cpair_t *
vrna_annotate_covar_pairs(const char  **alignment,
                          vrna_ep_t   *pl,
                          vrna_ep_t   *mfel,
                          double      threshold,
                          vrna_md_t   *md_p)
{
  unsigned int  n_seq;
  int           i, n, s, a, b, z, j, c, pfreq[7];
  vrna_cpair_t  *cp;
  vrna_ep_t     *ptr;
  vrna_md_t     md;

  if ((!alignment) || (!pl))
    return NULL;

  if (md_p)
    vrna_md_copy(&md, md_p);
  else
    vrna_md_set_default(&md);

  /* count number of sequences */
  for (n_seq = 0; alignment[n_seq] != NULL; n_seq++);

  /* count number of entries in pl */
  for (n = 0; pl[n].i > 0; n++);

  c   = 0;
  cp  = (vrna_cpair_t *)vrna_alloc(sizeof(vrna_cpair_t) * (n + 1));

  for (i = 0; i < n; i++) {
    int ncomp = 0;
    if (pl[i].p > threshold) {
      cp[c].i     = pl[i].i;
      cp[c].j     = pl[i].j;
      cp[c].p     = pl[i].p;
      cp[c].type  = pl[i].type;
      for (z = 0; z < 7; z++)
        pfreq[z] = 0;
      for (s = 0; s < n_seq; s++) {
        a = vrna_nucleotide_encode(alignment[s][cp[c].i - 1], &md);
        b = vrna_nucleotide_encode(alignment[s][cp[c].j - 1], &md);
        if ((alignment[s][cp[c].j - 1] == '~') || (alignment[s][cp[c].i - 1] == '~'))
          continue;

        if ((md.gquad) && (a == 3) && (b == 3))
          continue;

        pfreq[md.pair[a][b]]++;
      }
      for (z = 1; z < 7; z++)
        if (pfreq[z] > 0)
          ncomp++;

      cp[c].hue = MAX2(0.0, (ncomp - 1.0) / 6.2);   /* hue<6/6.9 (hue=1 ==  hue=0) */
      cp[c].sat = 1 - MIN2(1.0, (float)(pfreq[0] * 2. /*pi[i].bp[0]*/ / (n_seq)));
      c++;
    }
  }

  /*
   *  check whether all MFE base pairs are present in above list, and
   *  insert them if they are missing
   */
  if (mfel) {
    for (ptr = mfel; ptr->i > 0; ptr++) {
      int nofound = 1;
      for (j = 0; j < c; j++) {
        if ((cp[j].i == ptr->i) && (cp[j].j == ptr->j)) {
          cp[j].mfe = 1;
          nofound   = 0;
          break;
        }
      }
      if (nofound) {
        vrna_log_warning("mfe base pair with very low prob in pf: %d %d",
                             ptr->i,
                             ptr->j);

        cp          = (vrna_cpair_t *)vrna_realloc(cp, sizeof(vrna_cpair_t) * (c + 2));
        cp[c].i     = ptr->i;
        cp[c].j     = ptr->j;
        cp[c].p     = 0.;
        cp[c].type  = VRNA_PLIST_TYPE_BASEPAIR;
        cp[c].hue   = 0;
        cp[c].sat   = 0;
        cp[c].mfe   = 1;
        c++;
        cp[c].i = cp[c].j = 0;
      }
    }
  }

  return cp;
}
