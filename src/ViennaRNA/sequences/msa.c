/*
 *    ViennaRNA/utils/alignments.c
 *
 *    Helper functions related to alignments
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/params/ribosum.h"
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/io/utils.h"
#include "ViennaRNA/sequences/alignments.h"

#define NONE -10000 /* score for forbidden pairs */
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
 * IUP nucleotide classes indexed by a bit string of the present bases
 * A C AC G AG CG ACG U AU CU ACU GU AGU CGU ACGU
 */
static char IUP[17] = "-ACMGRSVUWYHKDBN";

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE char **
copy_alignment(const char   **alignment,
               unsigned int options);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_aln_mpi(const char **alignment)
{
  int   i, j, k, s, n_seq, n, pairnum = 0, sumident = 0;
  float ident = 0;

  if (alignment) {
    n = (int)strlen(alignment[0]);
    for (s = 0; alignment[s]; s++);
    n_seq = s;

    for (j = 0; j < n_seq - 1; j++)
      for (k = j + 1; k < n_seq; k++) {
        ident = 0;
        for (i = 1; i <= n; i++) {
          if (alignment[k][i] == alignment[j][i])
            ident++;

          pairnum++;
        }
        sumident += ident;
      }

    if (pairnum > 0)
      return (int)(sumident * 100 / pairnum);
  }

  return 0;
}


/*---------------------------------------------------------------------------*/
PRIVATE int
compare_pinfo(const void  *pi1,
              const void  *pi2)
{
  vrna_pinfo_t  *p1, *p2;
  int           i, nc1, nc2;

  p1  = (vrna_pinfo_t *)pi1;
  p2  = (vrna_pinfo_t *)pi2;
  for (nc1 = nc2 = 0, i = 1; i <= 6; i++) {
    if (p1->bp[i] > 0)
      nc1++;

    if (p2->bp[i] > 0)
      nc2++;
  }
  /* sort mostly by probability, add
   * epsilon * comp_mutations/(non-compatible+1) to break ties */
  return (p1->p + 0.01 * nc1 / (p1->bp[0] + 1.)) <
         (p2->p + 0.01 * nc2 / (p2->bp[0] + 1.)) ? 1 : -1;
}


PUBLIC vrna_pinfo_t *
vrna_aln_pinfo(vrna_fold_compound_t *vc,
               const char           *structure,
               double               threshold)
{
  int           i, j, num_p = 0, max_p = 64;
  vrna_pinfo_t  *pi;
  double        *duck, p;
  short         *ptable = NULL;

  short         **S       = vc->S;
  char          **AS      = vc->sequences;
  int           n_seq     = vc->n_seq;
  int           n         = vc->length;
  int           *my_iindx = vc->iindx;
  FLT_OR_DBL    *probs    = vc->exp_matrices->probs;
  vrna_md_t     *md       = &(vc->exp_params->model_details);

  max_p = 64;
  pi    = vrna_alloc(max_p * sizeof(vrna_pinfo_t));
  duck  = (double *)vrna_alloc((n + 1) * sizeof(double));
  if (structure)
    ptable = vrna_ptable(structure);

  for (i = 1; i < n; i++)
    for (j = i + 1; j <= n; j++) {
      if ((p = probs[my_iindx[i] - j]) >= threshold) {
        duck[i] -= p * log(p);
        duck[j] -= p * log(p);

        int type, s;
        pi[num_p].i   = i;
        pi[num_p].j   = j;
        pi[num_p].p   = p;
        pi[num_p].ent = duck[i] + duck[j] - p * log(p);

        for (type = 0; type < 8; type++)
          pi[num_p].bp[type] = 0;
        for (s = 0; s < n_seq; s++) {
          type = md->pair[S[s][i]][S[s][j]];
          if (S[s][i] == 0 && S[s][j] == 0)
            type = 7;                            /* gap-gap  */

          if ((AS[s][i - 1] == '-') || (AS[s][j - 1] == '-'))
            type = 7;

          if ((AS[s][i - 1] == '~') || (AS[s][j - 1] == '~'))
            type = 7;

          pi[num_p].bp[type]++;
        }
        if (ptable)
          pi[num_p].comp = (ptable[i] == j) ? 1 : 0;

        num_p++;
        if (num_p >= max_p) {
          max_p *= 2;
          pi    = vrna_realloc(pi, max_p * sizeof(vrna_pinfo_t));
        }
      }
    }
  free(duck);
  pi          = vrna_realloc(pi, (num_p + 1) * sizeof(vrna_pinfo_t));
  pi[num_p].i = 0;
  qsort(pi, num_p, sizeof(vrna_pinfo_t), compare_pinfo);

  free(ptable);
  return pi;
}


PUBLIC int
vrna_pscore(vrna_fold_compound_t  *fc,
            unsigned int          i,
            unsigned int          j)
{
  char          **AS;
  short         **S;
  unsigned int  n, tmp, u, s, n_seq;
  int           max_span;
  vrna_md_t     *md;

  if (i > j) {
    tmp = i;
    i   = j;
    j   = tmp;
  }

  if ((fc) &&
      (fc->type == VRNA_FC_TYPE_COMPARATIVE) &&
      (fc->length >= j)) {
    n         = fc->length;
    n_seq     = fc->n_seq;
    u         = j - i - 1;
    md        = &(fc->params->model_details);
    max_span  = md->max_bp_span;
    AS        = fc->sequences;
    S         = fc->S;

    if ((max_span < 2) ||
        (max_span > (int)n))
      max_span = (int)n;

    if ((u + 2) <= max_span) {
      unsigned int pfreq[8] = {
        0, 0, 0, 0, 0, 0, 0, 0
      };

      /* collect base pair frequencies */
      for (s = 0; s < n_seq; s++) {
        unsigned int type;

        if (S[s][i] == 0 && S[s][j] == 0) {
          type = 7;                             /* gap-gap  */
        } else {
          if ((AS[s][i] == '~') || (AS[s][j] == '~'))
            type = 7;
          else
            type = md->pair[S[s][i]][S[s][j]];
        }

        pfreq[type]++;
      }

      return vrna_pscore_freq(fc, (const unsigned int *)pfreq, 6);
    }
  }

  return NONE;
}


PUBLIC int
vrna_pscore_freq(vrna_fold_compound_t *fc,
                 const unsigned int   *frequencies,
                 unsigned int         pairs)
{
  if ((fc) &&
      (frequencies)) {
    unsigned int  k, l, n_seq;
    double        score;
    vrna_md_t     *md;

    n_seq = fc->n_seq;
    md    = &(fc->params->model_details);

    /*
     *  assume first (0) and last (pairs + 1) frequency to reflect
     *  counter-examples and gap-pairs, respectively
     */
    if (2 * frequencies[0] +
        frequencies[pairs + 1] < n_seq) {
      score = 0.;

      for (k = 1; k <= pairs; k++) /* ignore pairtype 7 (gap-gap) */
        for (l = k; l <= pairs; l++) {
          score += (double)frequencies[k] *
                   (double)frequencies[l] *
                   md->pair_dist[k][l];
        }

      /* counter examples score -1, gap-gap scores -0.25   */
      return md->cv_fact * ((UNIT * score) / n_seq -
                            md->nc_fact * UNIT * (frequencies[0] + frequencies[pairs + 1] * 0.25));
    }
  }

  return NONE;
}


PUBLIC int *
vrna_aln_pscore(const char  **alignment,
                vrna_md_t   *md)
{
  /*
   * calculate co-variance bonus for each pair depending on
   * compensatory/consistent mutations and incompatible seqs
   * should be 0 for conserved pairs, >0 for good pairs
   */


  int       i, j, k, l, s, n, n_seq, *indx, max_span;
  float     **dm;
  vrna_md_t md_default;
  int       *pscore;
  short     **S;

  int       olddm[7][7] = { { 0, 0, 0, 0, 0, 0, 0 },  /* hamming distance between pairs */
                            { 0, 0, 2, 2, 1, 2, 2 },  /* CG */
                            { 0, 2, 0, 1, 2, 2, 2 },  /* GC */
                            { 0, 2, 1, 0, 2, 1, 2 },  /* GU */
                            { 0, 1, 2, 2, 0, 2, 1 },  /* UG */
                            { 0, 2, 2, 1, 2, 0, 2 },  /* AU */
                            { 0, 2, 2, 2, 1, 2, 0 } /* UA */ };

  pscore = NULL;

  if (!md) {
    vrna_md_set_default(&md_default);
    md = &md_default;
  }

  if (alignment) {
    /* length of alignment */
    n = (int)strlen(alignment[0]);

    /* count number of sequences */
    for (s = 0; alignment[s]; s++);
    n_seq = s;

    /* make numeric encoding of sequences */
    S = (short **)vrna_alloc(sizeof(short *) * (n_seq + 1));
    for (s = 0; s < n_seq; s++)
      S[s] = vrna_seq_encode_simple(alignment[s], md);

    indx = vrna_idx_col_wise(n);

    pscore = (int *)vrna_alloc(sizeof(int) * ((n + 1) * (n + 2) / 2 + 2));

    if (md->ribo) {
      if (RibosumFile != NULL)
        dm = readribosum(RibosumFile);
      else
        dm = get_ribosum(alignment, n_seq, n);
    } else {
      /*use usual matrix*/
      dm = vrna_alloc(7 * sizeof(float *));
      for (i = 0; i < 7; i++) {
        dm[i] = vrna_alloc(7 * sizeof(float));
        for (j = 0; j < 7; j++)
          dm[i][j] = (float)olddm[i][j];
      }
    }

    max_span = md->max_bp_span;
    if ((max_span < 2) || (max_span > n))
      max_span = n;

    for (i = 1; i < n; i++) {
      for (j = i + 1; j <= n; j++) {
        int     pfreq[8] = {
          0, 0, 0, 0, 0, 0, 0, 0
        };
        double  score;

        if ((j - i + 1) > max_span) {
          pscore[indx[j] + i] = NONE;
          continue;
        }

        for (s = 0; s < n_seq; s++) {
          int type;
          if (S[s][i] == 0 && S[s][j] == 0) {
            type = 7;                             /* gap-gap  */
          } else {
            if ((alignment[s][i] == '~') || (alignment[s][j] == '~'))
              type = 7;
            else
              type = md->pair[S[s][i]][S[s][j]];
          }

          pfreq[type]++;
        }

        if (pfreq[0] * 2 + pfreq[7] >= n_seq) {
          pscore[indx[j] + i] = NONE;
          continue;
        }

        for (k = 1, score = 0; k <= 6; k++) /* ignore pairtype 7 (gap-gap) */
          for (l = k; l <= 6; l++)
            score += pfreq[k] * pfreq[l] * dm[k][l];
        /* counter examples score -1, gap-gap scores -0.25   */
        pscore[indx[j] + i] = md->cv_fact * ((UNIT * score) / n_seq -
                                             md->nc_fact * UNIT * (pfreq[0] + pfreq[7] * 0.25));
      }
    }

    if (md->noLP) {
      /* remove unwanted pairs */
      for (k = 1; k < n - 1; k++)
        for (l = 1; l <= 2; l++) {
          int type, ntype = 0, otype = 0;
          i     = k;
          j     = i + l;
          type  = pscore[indx[j] + i];
          while ((i >= 1) && (j <= n)) {
            if ((i > 1) && (j < n))
              ntype = pscore[indx[j + 1] + i - 1];

            if ((otype < md->cv_fact * MINPSCORE) && (ntype < md->cv_fact * MINPSCORE)) /* too many counterexamples */
              pscore[indx[j] + i] = NONE;                                               /* i.j can only form isolated pairs */

            otype = type;
            type  = ntype;
            i--;
            j++;
          }
        }
    }

    /*free dm */
    for (i = 0; i < 7; i++)
      free(dm[i]);
    free(dm);

    for (s = 0; s < n_seq; s++)
      free(S[s]);
    free(S);

    free(indx);
  }

  return pscore;
}


PUBLIC char **
vrna_aln_slice(const char   **alignment,
               unsigned int i,
               unsigned int j)
{
  char          **sub;
  unsigned int  n;
  int           n_seq, s;

  sub = NULL;

  if (alignment) {
    /* check if [i,j] is within range */
    n = strlen(alignment[0]);

    if ((i < j) && (j <= n)) {
      /* get number of sequences in alignment */
      for (n_seq = 0; alignment[n_seq] != NULL; n_seq++);

      sub = (char **)vrna_alloc(sizeof(char *) * (n_seq + 1));
      for (s = 0; s < n_seq; s++)
        sub[s] = vrna_alloc(sizeof(char) * (j - i + 2));
      sub[s] = NULL;

      /* copy subalignment */
      for (s = 0; s < n_seq; s++) {
        sub[s]              = memcpy(sub[s], alignment[s] + i - 1, sizeof(char) * (j - i + 1));
        sub[s][(j - i + 1)] = '\0';
      }
    }
  }

  return sub;
}


PUBLIC void
vrna_aln_free(char **alignment)
{
  int n_seq;

  if (alignment) {
    for (n_seq = 0; alignment[n_seq] != NULL; n_seq++)
      free(alignment[n_seq]);
    free(alignment);
  }
}


PUBLIC char **
vrna_aln_uppercase(const char **alignment)
{
  if (alignment)
    return copy_alignment(alignment, VRNA_ALN_UPPERCASE);

  return NULL;
}


PUBLIC char **
vrna_aln_toRNA(const char **alignment)
{
  if (alignment)
    return copy_alignment(alignment, VRNA_ALN_RNA);

  return NULL;
}


PUBLIC char **
vrna_aln_copy(const char    **alignment,
              unsigned int  options)
{
  if (alignment)
    return copy_alignment(alignment, options);

  return NULL;
}


PUBLIC float *
vrna_aln_conservation_struct(const char       **alignment,
                             const char       *structure,
                             const vrna_md_t  *md_p)
{
  unsigned int  i, j, n, s, n_seq;
  int           a, b;
  float         *conservation = NULL;
  short         *pt;
  vrna_md_t     md;

  if ((alignment) && (structure)) {
    n = strlen(structure);
    if (n > 0) {
      /* check alignment for consistency */
      for (s = 0; alignment[s]; s++) {
        if (strlen(alignment[s]) != n) {
          vrna_log_warning("vrna_aln_bpcons: Length of aligned sequence %u does not match consensus structure length\n"
                               "%s\n%s\n",
                               s + 1,
                               alignment[s],
                               structure);
          return NULL;
        }
      }

      n_seq = s;

      if (md_p)
        vrna_md_copy(&md, md_p);
      else
        vrna_md_set_default(&md);

      pt            = vrna_ptable(structure);
      conservation  = (float *)vrna_alloc(sizeof(float) * (n + 1));

      for (i = 1; i < n; i++) {
        if (i < pt[i]) {
          j = pt[i];
          for (s = 0; s < n_seq; s++) {
            a = vrna_nucleotide_encode(alignment[s][i - 1], &md);
            b = vrna_nucleotide_encode(alignment[s][j - 1], &md);
            if (md.pair[a][b]) {
              conservation[i] += 1.;
              conservation[j] += 1.;
            }
          }
          conservation[i] /= (float)n_seq;
          conservation[j] /= (float)n_seq;
        }
      }

      free(pt);
    } else {
      vrna_log_warning("vrna_aln_bpcons: Structure length is 0!");
    }
  }

  return conservation;
}


PUBLIC float *
vrna_aln_conservation_col(const char      **alignment,
                          const vrna_md_t *md_p,
                          unsigned int    options)
{
  unsigned int  i, n, s, n_seq;
  float         *conservation = NULL;
  vrna_md_t     md;

  if (alignment) {
    n = strlen(alignment[0]);
    if (n > 0) {
      /* check alignment for consistency */
      for (s = 1; alignment[s]; s++) {
        if (strlen(alignment[s]) != n) {
          vrna_log_warning("vrna_aln_conservation: Length of aligned sequence #%d does not match length of first sequence\n"
                               "%s\n\n",
                               s + 1,
                               alignment[s]);
          return NULL;
        }
      }

      n_seq = s;

      if (md_p)
        vrna_md_copy(&md, md_p);
      else
        vrna_md_set_default(&md);

      conservation = (float *)vrna_alloc(sizeof(float) * (n + 1));
      for (i = 1; i <= n; i++) {
        /*
         *  Here, we differentiate between 32 different nucleotide types.
         *  Change this if we ever cope with larger alphabets!
         */
        unsigned int nfreq[32] = {
          0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0
        };

        /* count frequencies of individual nucleotides */
        for (s = 0; s < n_seq; s++)
          nfreq[vrna_nucleotide_encode(alignment[s][i - 1], &md)]++;

        if (options & VRNA_MEASURE_SHANNON_ENTROPY) {
          double sum = 0;
          for (s = 0; s < 32; s++) {
            if (nfreq[s] > 0) {
              double p = (double)nfreq[s] / (double)n_seq;
              sum += p * log(p) / log(2.0);
            }
          }
          conservation[i] = (float)(-sum);
        }
      }
    } else {
      vrna_log_warning("vrna_aln_conservation: Length of first sequence in alignment is 0!");
    }
  }

  return conservation;
}


PUBLIC char *
vrna_aln_consensus_sequence(const char      **alignment,
                            const vrna_md_t *md_p)
{
  /* simple consensus sequence (most frequent character) */
  char          *consensus;
  unsigned int  i, n, s, n_seq;
  vrna_md_t     md;

  consensus = NULL;

  if (alignment) {
    n = strlen(alignment[0]);
    if (n > 0) {
      /* check alignment for consistency */
      for (s = 1; alignment[s]; s++) {
        if (strlen(alignment[s]) != n) {
          vrna_log_warning("vrna_aln_consensus_sequence: "
                               "Length of aligned sequence #%d does not match length of first sequence\n"
                               "%s\n\n",
                               s + 1,
                               alignment[s]);
          return NULL;
        }
      }

      n_seq = s;

      if (md_p)
        vrna_md_copy(&md, md_p);
      else
        vrna_md_set_default(&md);

      consensus = (char *)vrna_alloc((n + 1) * sizeof(char));

      for (i = 0; i < n; i++) {
        int c, fm;
        int freq[8] = {
          0, 0, 0, 0, 0, 0, 0, 0
        };

        for (s = 0; s < n_seq; s++)
          freq[vrna_nucleotide_encode(alignment[s][i], &md)]++;

        for (s = c = fm = 0; s < 8; s++) /* find the most frequent char */
          if (freq[s] > fm) {
            c   = s;
            fm  = freq[c];
          }

        if (s > 4)
          s++;        /* skip T */

        consensus[i] = vrna_nucleotide_decode(c, &md);
      }
    }
  }

  return consensus;
}


PUBLIC char *
vrna_aln_consensus_mis(const char       **alignment,
                       const vrna_md_t  *md_p)
{
  /*
   *  MIS displays the 'most informative sequence' (Freyhult et al 2004),
   *  elements in columns with frequency greater than the background
   *  frequency are projected into iupac notation. Columns where gaps are
   *  over-represented are in lower case.
   */
  char          *mis;
  unsigned char c;
  unsigned int  i, n, s, n_seq;
  unsigned int  bgfreq[8] = {
    0, 0, 0, 0, 0, 0, 0, 0
  };
  vrna_md_t     md;

  mis = NULL;

  if (alignment) {
    n = strlen(alignment[0]);
    if (n > 0) {
      /* check alignment for consistency */
      for (s = 1; alignment[s]; s++) {
        if (strlen(alignment[s]) != n) {
          vrna_log_warning("vrna_aln_consensus_mis: "
                               "Length of aligned sequence #%d does not match length of first sequence\n"
                               "%s\n\n",
                               s + 1,
                               alignment[s]);
          return NULL;
        }
      }

      n_seq = s;

      if (md_p)
        vrna_md_copy(&md, md_p);
      else
        vrna_md_set_default(&md);

      mis = (char *)vrna_alloc((n + 1) * sizeof(char));

      /* determien backgroud frequencies */
      for (i = 0; i < n; i++)
        for (s = 0; s < n_seq; s++) {
          c = (unsigned char)vrna_nucleotide_encode(alignment[s][i], &md);
          if (c > 4)
            c = 5;

          bgfreq[c]++;
        }

      /* generate MIS */
      for (i = 0; i < n; i++) {
        int           code    = 0;
        unsigned int  freq[8] = {
          0, 0, 0, 0, 0, 0, 0, 0
        };

        for (s = 0; s < n_seq; s++) {
          c = (unsigned char)vrna_nucleotide_encode(alignment[s][i], &md);
          if (c > 4)
            c = 5;

          freq[c]++;
        }
        for (c = 4; c > 0; c--) {
          code <<= 1;
          if (freq[c] * n >= bgfreq[c])
            code++;
        }
        mis[i] = IUP[code];
        if (freq[0] * n > bgfreq[0])
          mis[i] = tolower(IUP[code]);
      }
    }
  }

  return mis;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE char **
copy_alignment(const char   **alignment,
               unsigned int options)
{
  char          **output = NULL;
  unsigned int  s;

  /* determine number of sequence in the alignment */
  for (s = 0; alignment[s]; s++);

  /* allocate memory for output alignment */
  output = (char **)vrna_alloc(sizeof(char *) * (s + 1));

  /* copy alignment and convert to uppercase */
  for (s = 0; alignment[s]; s++) {
    output[s] = strdup(alignment[s]);

    if (options & VRNA_ALN_UPPERCASE)
      /* convert to uppercase letters */
      vrna_seq_toupper(output[s]);

    if (options & VRNA_ALN_RNA)
      /* convert DNA to RNA alphabet */
      vrna_seq_toRNA(output[s]);
  }

  /* mark end of alignment */
  output[s] = NULL;

  return output;
}


/*
 *###########################################
 *# deprecated functions below              #
 *###########################################
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

#define MAX_NUM_NAMES    500

PRIVATE void
encode_ali_sequence_old(const char      *sequence,
                        short           *S,
                        short           *s5,
                        short           *s3,
                        char            *ss,
                        unsigned short  *as,
                        int             circular);


int
read_clustal(FILE *clust,
             char *AlignedSeqs[],
             char *names[])
{
  char  *line, name[100] = "", *seq;
  int   n, nn = 0, num_seq = 0, i;

  if ((line = vrna_read_line(clust)) == NULL) {
    vrna_log_warning("Empty CLUSTAL file");
    return 0;
  }

  if ((strncmp(line, "CLUSTAL", 7) != 0) && (!strstr(line, "STOCKHOLM"))) {
    vrna_log_warning("This doesn't look like a CLUSTAL/STOCKHOLM file, sorry");
    free(line);
    return 0;
  }

  free(line);
  line = vrna_read_line(clust);

  while (line != NULL) {
    if (strncmp(line, "//", 2) == 0) {
      free(line);
      break;
    }

    if (((n = strlen(line)) < 4) || isspace((int)line[0])) {
      /* skip non-sequence line */
      free(line);
      line  = vrna_read_line(clust);
      nn    = 0; /* reset seqence number */
      continue;
    }

    /* skip comments */
    if (line[0] == '#') {
      free(line);
      line = vrna_read_line(clust);
      continue;
    }

    seq = (char *)vrna_alloc((n + 1) * sizeof(char));
    sscanf(line, "%99s %s", name, seq);

    for (i = 0; i < strlen(seq); i++) {
      if (seq[i] == '.')
        seq[i] = '-';                 /* replace '.' gaps by '-' */

      /* comment the next line and think about something more difficult to deal with
       * lowercase sequence letters if you really want to */
      seq[i] = toupper(seq[i]);
    }

    if (nn == num_seq) {
      /* first time */
      names[nn]       = strdup(name);
      AlignedSeqs[nn] = strdup(seq);
    } else {
      if (strcmp(name, names[nn]) != 0) {
        /* name doesn't match */
        vrna_log_warning("Sorry, your file is messed up (inconsitent seq-names)");
        free(line);
        free(seq);
        return 0;
      }

      AlignedSeqs[nn] = (char *)
                        vrna_realloc(AlignedSeqs[nn], strlen(seq) + strlen(AlignedSeqs[nn]) + 1);
      strcat(AlignedSeqs[nn], seq);
    }

    nn++;
    if (nn > num_seq)
      num_seq = nn;

    free(seq);
    free(line);
    if (num_seq >= MAX_NUM_NAMES) {
      vrna_log_warning("Too many sequences in CLUSTAL/STOCKHOLM file");
      return 0;
    }

    line = vrna_read_line(clust);
  }

  AlignedSeqs[num_seq]  = NULL;
  names[num_seq]        = NULL;
  if (num_seq == 0) {
    vrna_log_warning("No sequences found in CLUSTAL/STOCKHOLM file");
    return 0;
  }

  n = strlen(AlignedSeqs[0]);
  for (nn = 1; nn < num_seq; nn++) {
    if (strlen(AlignedSeqs[nn]) != n) {
      vrna_log_warning("Sorry, your file is messed up.\n"
                           "Unequal lengths!");
      return 0;
    }
  }

  vrna_log_info("%d sequences; length of alignment %d.", nn, n);
  return num_seq;
}


char *
consensus(const char *AS[])
{
  /* simple consensus sequence (most frequent character) */
  char  *string;
  int   i, n;

  string = NULL;

  if (AS) {
    n       = strlen(AS[0]);
    string  = (char *)vrna_alloc((n + 1) * sizeof(char));
    for (i = 0; i < n; i++) {
      int s, c, fm, freq[8] = {
        0, 0, 0, 0, 0, 0, 0, 0
      };
      for (s = 0; AS[s] != NULL; s++)
        freq[encode_char(AS[s][i])]++;
      for (s = c = fm = 0; s < 8; s++) /* find the most frequent char */
        if (freq[s] > fm)
          c = s, fm = freq[c];

      if (s > 4)
        s++;        /* skip T */

      string[i] = Law_and_Order[c];
    }
  }

  return string;
}


char *
consens_mis(const char *AS[])
{
  /* MIS displays the 'most informative sequence' (Freyhult et al 2004),
   * elements in columns with frequency greater than the background
   * frequency are projected into iupac notation. Columns where gaps are
   * over-represented are in lower case. */

  char  *cons;
  int   i, s, n, N, c;
  int   bgfreq[8] = {
    0, 0, 0, 0, 0, 0, 0, 0
  };

  cons = NULL;

  if (AS) {
    n = strlen(AS[0]);
    for (N = 0; AS[N] != NULL; N++);
    cons = (char *)vrna_alloc((n + 1) * sizeof(char));

    for (i = 0; i < n; i++)
      for (s = 0; s < N; s++) {
        c = encode_char(AS[s][i]);
        if (c > 4)
          c = 5;

        bgfreq[c]++;
      }

    for (i = 0; i < n; i++) {
      int freq[8] = {
        0, 0, 0, 0, 0, 0, 0, 0
      };
      int code = 0;
      for (s = 0; s < N; s++) {
        c = encode_char(AS[s][i]);
        if (c > 4)
          c = 5;

        freq[c]++;
      }
      for (c = 4; c > 0; c--) {
        code <<= 1;
        if (freq[c] * n >= bgfreq[c])
          code++;
      }
      cons[i] = IUP[code];
      if (freq[0] * n > bgfreq[0])
        cons[i] = tolower(IUP[code]);
    }
  }

  return cons;
}


PUBLIC char *
get_ungapped_sequence(const char *seq)
{
  char  *tmp_sequence, *b;
  int   i;

  tmp_sequence = strdup(seq);

  b = tmp_sequence;
  i = 0;
  do {
    if ((*b == '-') || (*b == '_') || (*b == '~') || (*b == '.'))
      continue;

    tmp_sequence[i] = *b;
    i++;
  } while (*(++b));

  tmp_sequence    = (char *)vrna_realloc(tmp_sequence, (i + 1) * sizeof(char));
  tmp_sequence[i] = '\0';

  return tmp_sequence;
}


PUBLIC int
get_mpi(char  *Alseq[],
        int   n_seq,
        int   length,
        int   *mini)
{
  int   i, j, k, pairnum = 0, sumident = 0;
  float ident = 0, minimum = 1.;

  for (j = 0; j < n_seq - 1; j++)
    for (k = j + 1; k < n_seq; k++) {
      ident = 0;
      for (i = 1; i <= length; i++) {
        if (Alseq[k][i] == Alseq[j][i])
          ident++;

        pairnum++;
      }
      if ((ident / length) < minimum)
        minimum = ident / (float)length;

      sumident += ident;
    }
  mini[0] = (int)(minimum * 100.);
  if (pairnum > 0)
    return (int)(sumident * 100 / pairnum);
  else
    return 0;
}


PUBLIC void
alloc_sequence_arrays(const char      **sequences,
                      short           ***S,
                      short           ***S5,
                      short           ***S3,
                      unsigned short  ***a2s,
                      char            ***Ss,
                      int             circ)
{
  unsigned int s, n_seq, length;

  if (sequences[0] != NULL) {
    length = strlen(sequences[0]);
    for (s = 0; sequences[s] != NULL; s++);
    n_seq = s;
    *S    = (short **)vrna_alloc((n_seq + 1) * sizeof(short *));
    *S5   = (short **)vrna_alloc((n_seq + 1) * sizeof(short *));
    *S3   = (short **)vrna_alloc((n_seq + 1) * sizeof(short *));
    *a2s  = (unsigned short **)vrna_alloc((n_seq + 1) * sizeof(unsigned short *));
    *Ss   = (char **)vrna_alloc((n_seq + 1) * sizeof(char *));
    for (s = 0; s < n_seq; s++) {
      if (strlen(sequences[s]) != length) {
        vrna_log_error("uneqal seqence lengths");
        *S = NULL;
        *S5 = NULL;
        *S3 = NULL;
        *a2s = NULL;
        *Ss = NULL;
        return;
      }

      (*S5)[s]  = (short *)vrna_alloc((length + 2) * sizeof(short));
      (*S3)[s]  = (short *)vrna_alloc((length + 2) * sizeof(short));
      (*a2s)[s] = (unsigned short *)vrna_alloc((length + 2) * sizeof(unsigned short));
      (*Ss)[s]  = (char *)vrna_alloc((length + 2) * sizeof(char));
      (*S)[s]   = (short *)vrna_alloc((length + 2) * sizeof(short));
      encode_ali_sequence_old(sequences[s], (*S)[s], (*S5)[s], (*S3)[s], (*Ss)[s], (*a2s)[s], circ);
    }
    (*S5)[n_seq]  = NULL;
    (*S3)[n_seq]  = NULL;
    (*a2s)[n_seq] = NULL;
    (*Ss)[n_seq]  = NULL;
    (*S)[n_seq]   = NULL;
  } else {
    vrna_log_error("no sequences in the alignment!");
    *S = NULL;
    *S5 = NULL;
    *S3 = NULL;
    *a2s = NULL;
    *Ss = NULL;
    return;
  }
}


PUBLIC void
free_sequence_arrays(unsigned int   n_seq,
                     short          ***S,
                     short          ***S5,
                     short          ***S3,
                     unsigned short ***a2s,
                     char           ***Ss)
{
  unsigned int s;

  for (s = 0; s < n_seq; s++) {
    free((*S)[s]);
    free((*S5)[s]);
    free((*S3)[s]);
    free((*a2s)[s]);
    free((*Ss)[s]);
  }
  free(*S);
  *S = NULL;
  free(*S5);
  *S5 = NULL;
  free(*S3);
  *S3 = NULL;
  free(*a2s);
  *a2s = NULL;
  free(*Ss);
  *Ss = NULL;
}


PUBLIC void
encode_ali_sequence(const char      *sequence,
                    short           *S,
                    short           *s5,
                    short           *s3,
                    char            *ss,
                    unsigned short  *as,
                    int             circular)
{
  if ((sequence) &&
      (S) &&
      (s5) &&
      (s3) &&
      (ss) &&
      (as))
    encode_ali_sequence_old(sequence, S, s5, s3, ss, as, circular);
}


PRIVATE void
encode_ali_sequence_old(const char      *sequence,
                        short           *S,
                        short           *s5,
                        short           *s3,
                        char            *ss,
                        unsigned short  *as,
                        int             circular)
{
  unsigned int    i, l;
  unsigned short  p;

  l     = strlen(sequence);
  S[0]  = (short)l;
  s5[0] = s5[1] = 0;

  /* make numerical encoding of sequence */
  for (i = 1; i <= l; i++) {
    short ctemp;
    ctemp = (short)encode_char(toupper(sequence[i - 1]));
    S[i]  = ctemp;
  }

  if (oldAliEn) {
    /* use alignment sequences in all energy evaluations */
    ss[0] = sequence[0];
    for (i = 1; i < l; i++) {
      s5[i] = S[i - 1];
      s3[i] = S[i + 1];
      ss[i] = sequence[i];
      as[i] = i;
    }
    ss[l]     = sequence[l];
    as[l]     = l;
    s5[l]     = S[l - 1];
    s3[l]     = 0;
    S[l + 1]  = S[1];
    s5[1]     = 0;
    if (circular) {
      s5[1]     = S[l];
      s3[l]     = S[1];
      ss[l + 1] = S[1];
    }
  } else {
    if (circular) {
      for (i = l; i > 0; i--) {
        char c5;
        c5 = sequence[i - 1];
        if ((c5 == '-') || (c5 == '_') || (c5 == '~') || (c5 == '.'))
          continue;

        s5[1] = S[i];
        break;
      }
      for (i = 1; i <= l; i++) {
        char c3;
        c3 = sequence[i - 1];
        if ((c3 == '-') || (c3 == '_') || (c3 == '~') || (c3 == '.'))
          continue;

        s3[l] = S[i];
        break;
      }
    } else {
      s5[1] = s3[l] = 0;
    }

    for (i = 1, p = 0; i <= l; i++) {
      char c5;
      c5 = sequence[i - 1];
      if ((c5 == '-') || (c5 == '_') || (c5 == '~') || (c5 == '.')) {
        s5[i + 1] = s5[i];
      } else {
        /* no gap */
        ss[p++]   = sequence[i - 1]; /*start at 0!!*/
        s5[i + 1] = S[i];
      }

      as[i] = p;
    }
    for (i = l; i >= 1; i--) {
      char c3;
      c3 = sequence[i - 1];
      if ((c3 == '-') || (c3 == '_') || (c3 == '~') || (c3 == '.'))
        s3[i - 1] = s3[i];
      else
        s3[i - 1] = S[i];
    }
  }
}


#endif
