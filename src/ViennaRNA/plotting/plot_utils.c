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
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/structures.h"
#include "ViennaRNA/alphabet.h"
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
  /* produce annotation for colored drawings from vrna_file_PS_rnaplot_a() */
  char      *ps, *colorps, **A;
  int       i, n, s, pairings, maxl, a, b;
  short     *ptable;
  char      *colorMatrix[6][3] = {
    { "0.0 1",  "0.0 0.6",  "0.0 0.2"  },   /* red    */
    { "0.16 1", "0.16 0.6", "0.16 0.2" },   /* ochre  */
    { "0.32 1", "0.32 0.6", "0.32 0.2" },   /* turquoise */
    { "0.48 1", "0.48 0.6", "0.48 0.2" },   /* green  */
    { "0.65 1", "0.65 0.6", "0.65 0.2" },   /* blue   */
    { "0.81 1", "0.81 0.6", "0.81 0.2" }    /* violet */
  };
  vrna_md_t md;

  if ((!alignment) || (!structure))
    return NULL;

  if (md_p)
    vrna_md_copy(&md, md_p);
  else
    vrna_md_set_default(&md);

  n     = strlen(alignment[0]);
  maxl  = 1024;

  A       = (char **)vrna_alloc(sizeof(char *) * 2);
  ps      = (char *)vrna_alloc(maxl);
  colorps = (char *)vrna_alloc(maxl);

  ptable  = vrna_ptable_from_string(structure, options);

  for (i = 1; i <= n; i++) {
    char  pps[64], ci = '\0', cj = '\0';
    int   j, type, pfreq[8] = {
      0, 0, 0, 0, 0, 0, 0, 0
    }, vi = 0, vj = 0;

    if ((j = ptable[i]) < i)
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

    if ((maxl - strlen(ps) < 192) || ((maxl - strlen(colorps)) < 64)) {
      maxl    *= 2;
      ps      = (char *)vrna_realloc(ps, sizeof(char) * maxl);
      colorps = (char *)vrna_realloc(colorps, sizeof(char) * maxl);
      if ((ps == NULL) || (colorps == NULL))
        vrna_message_error("out of memory in realloc");
    }

    if (pfreq[0] <= 2 && pairings > 0) {
      snprintf(pps, 64, "%d %d %s colorpair\n",
               i, j, colorMatrix[pairings - 1][pfreq[0]]);
      strcat(colorps, pps);
    }

    if (pfreq[0] > 0) {
      snprintf(pps, 64, "%d %d %d gmark\n", i, j, pfreq[0]);
      strcat(ps, pps);
    }

    if (vi > 1) {
      snprintf(pps, 64, "%d cmark\n", i);
      strcat(ps, pps);
    }

    if (vj > 1) {
      snprintf(pps, 64, "%d cmark\n", j);
      strcat(ps, pps);
    }
  }
  free(ptable);
  A[0]  = colorps;
  A[1]  = ps;
  return A;
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
        vrna_message_warning("mfe base pair with very low prob in pf: %d %d",
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
