/*
   Fast, but crude, pairwise structural Alignments of RNA sequences

   Possible structures of each RNA are encoded in a linear
   "probability profile", by computing for each base the probability
   of being unpaired, or paired upstream or downstream. These profiles
   can be aligned using standard string alignment.

   The is an extension of the old method in ProfileDist.c with the
   following changes:
   - use sequence as well as structure profile for scoring
   - use similarity alignment instead of distance (maybe add local alinment)
   - use affine gap costs

	  C Ivo L Hofacker, Vienna RNA Package
*/
/* Last changed Time-stamp: <2004-08-02 10:11:13 ivo> */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include "dist_vars.h"
#include "fold_vars.h"
#include "part_func.h"
#include "utils.h"
/*@unused@*/
static char rcsid[] = "$Id: ProfileAln.c,v 1.5 2006/01/18 13:00:30 ivo Exp $";

#define PUBLIC
#define PRIVATE        static
#define MIN(x,y)       (((x)<(y)) ? (x) : (y))
#define MAX(x,y)       (((x)>(y)) ? (x) : (y))
#define MAX3(x,y,z)    (MAX(  (MAX((x),(y))) ,(z)))
#define EQUAL(x,y)     (fabs((x)-(y)) <= fabs(x)*2*FLT_EPSILON)

PRIVATE int *alignment[2];

PUBLIC float    profile_aln(const float *T1, const char *seq1,
			    const float *T2, const char *seq2);

PUBLIC int set_paln_params(double gap_open, double gap_ext,
			   double seqweight, int free_ends);

PRIVATE void    sprint_aligned_bppm(const float *T1, const char *seq1,
				    const float *T2, const char *seq2);
PRIVATE double  PrfEditScore(const float *p1, const float *p2,
			     char c1, char c2);
PRIVATE double  average(double x, double y);

PRIVATE double  open=-1.5, ext=-0.666;  /* defaults from clustalw */
PRIVATE double  seqw=0.5;
PRIVATE int     free_ends=1;            /* whether to use free end gaps */
extern float    *Make_bp_profile(int length);

/*---------------------------------------------------------------------------*/

PRIVATE float **newmat(int l1, int l2) {
  float **a;
  int i;
  a = (float **) space((l1+1)*sizeof(float *));
  for (i=0; i<=l1; i++) a[i] = (float *) space((l2+1)*sizeof(float));
  return a;
}

PUBLIC float profile_aln(const float *T1, const char *seq1,
			 const float *T2, const char *seq2)
{
  /* align the 2 probability profiles T1, T2 */
  /* This is like a Needleman-Wunsch alignment, with affine gap-costs
     ala Gotoh. The score looks at both seq and pair profile */

  float  **S, **E, **F, tot_score;
  int    i, j, length1, length2;

  length1 = strlen(seq1);
  length2 = strlen(seq2);
  S = newmat(length1, length2);
  E = newmat(length1, length2);
  F = newmat(length1, length2);

  E[0][0] = F[0][0] = open - ext;
  S[0][0] = 0;
  for (i=1; i<=length1; i++) F[i][0] = -9999; /* impossible */
  for (j=1; j<=length2; j++) E[0][j] = -9999; /* impossible */
  if (!free_ends) {
    for (i=1; i<=length1; i++) S[i][0] = E[i][0] = E[i-1][0] +ext;
    for (j=1; j<=length2; j++) S[0][j] = F[0][j] = F[0][j-1] +ext;
  }

  for (i=1; i<=length1; i++) {
    for (j=1; j<=length2; j++) {
      float M;
      E[i][j] = MAX(E[i-1][j]+ext, S[i-1][j]+open);
      F[i][j] = MAX(F[i][j-1]+ext, S[i][j-1]+open);
      M = S[i-1][j-1] + PrfEditScore(T1+3*i,T2+3*j, seq1[i-1], seq2[j-1]);
      S[i][j] = MAX3(M, E[i][j], F[i][j]);
    }
  }

  if (edit_backtrack) {
    double score=0;
    char state = 'S';
    int pos, i,j;
    alignment[0] = (int *) space((length1+length2+1)*sizeof(int));
    alignment[1] = (int *) space((length1+length2+1)*sizeof(int));

    pos = length1+length2;
    i   = length1;
    j   = length2;

    tot_score = S[length1][length2];

    if (free_ends) {
      /* find starting point for backtracking,
	 search for highest entry in last row or column */
      int imax=0;
      for (i=1; i<=length1; i++) {
	if (S[i][length2]>score) {
	  score=S[i][length2];
	  imax=i;
	}
      }
      for (j=1; j<=length2; j++) {
	if (S[length1][j]>score) {
	  score=S[length1][j];
	  imax=-j;
	}
      }
      if (imax<0) {
	for (j=length2; j> -imax; j--) {
	  alignment[0][pos] = 0;
	  alignment[1][pos--] = j;
	}
	i=length1;
      } else {
	for (i=length1; i>imax; i--) {
	  alignment[0][pos] = i;
	  alignment[1][pos--] = 0;
	}
	j=length2;
      }
      tot_score=score;
    }

    while (i>0 && j>0) {
      switch (state) {
      case 'E':
	score = E[i][j];
	alignment[0][pos] = i;
	alignment[1][pos--] = 0;
	if (EQUAL(score, S[i-1][j] + open)) state = 'S';
	i--;
	break;
      case 'F':
	score = F[i][j];
	alignment[0][pos] = 0;
	alignment[1][pos--] = j;
	if (EQUAL(score, S[i][j-1] + open)) state = 'S';
	j--;
	break;
      case 'S':
	score = S[i][j];
	if (EQUAL(score, E[i][j])) state = 'E';
	else if (EQUAL(score, F[i][j])) state = 'F';
	else if (EQUAL(score, S[i-1][j-1] +
		       PrfEditScore(T1+3*i,T2+3*j, seq1[i-1], seq2[j-1]))) {
	  alignment[0][pos] = i;
	  alignment[1][pos--] = j;
	  i--; j--;
	}
	else nrerror("backtrack of alignment failed");
	break;
      }
    }

    for (; j>0; j--) {
      alignment[0][pos] = 0;
      alignment[1][pos--] = j;
    }
    for (; i>0; i--) {
      alignment[0][pos] = i;
      alignment[1][pos--] = 0;
    }

    for(i=pos+1; i<=length1+length2; i++){
      alignment[0][i-pos] = alignment[0][i];
      alignment[1][i-pos] = alignment[1][i];
    }
    alignment[0][0] = length1+length2-pos;   /* length of alignment */

    sprint_aligned_bppm(T1,seq1, T2,seq2);
    free(alignment[0]);
    free(alignment[1]);
  }
  for (i=0; i<=length1; i++) {
    free(S[i]); free(E[i]); free(F[i]);
  }
  free(S); free(E); free(F);

  return tot_score;
}


/*---------------------------------------------------------------------------*/
PRIVATE inline double average(double x, double y) {
  /*
     As in Bonhoeffer et al (1993) 'RNA Multi Structure Landscapes',
     Eur. Biophys. J. 22: 13-24 we have chosen  the geometric mean.
  */
  return (float) sqrt(x*y);
}

PRIVATE double PrfEditScore(const float *p1, const float *p2, char c1, char c2)
{
  double  score;
  int    k;

  for(score=0.,k=0; k<3; k++)
    score += average(p1[k],p2[k]);

  score *= (1- seqw);
  if (c1==c2) score +=  seqw;
  else if (((c1=='A') && (c2=='G')) ||
	   ((c1=='G') && (c2=='A')) ||
	   ((c1=='C') && (c2=='U')) ||
	   ((c1=='U') && (c2=='C')))
    score += 0.5*seqw;
  else score -= 0.9*seqw;
  return score;
}

/*---------------------------------------------------------------------------*/

PRIVATE void sprint_aligned_bppm(const float *T1, const char *seq1,
				 const float *T2, const char *seq2) {
   int     i, length;
   length = alignment[0][0];
   for (i=0; i<4; i++) {
     if (aligned_line[i] != NULL) free(aligned_line[i]);
     aligned_line[i] = (char *) space((length+1)*sizeof(char));
   }
   for(i=1; i<=length; i++){
      if (alignment[0][i]==0)
	aligned_line[0][i-1] = aligned_line[2][i-1] = '_';
      else {
	aligned_line[0][i-1] = bppm_symbol(T1+alignment[0][i]*3);
	aligned_line[2][i-1] = seq1[alignment[0][i]-1];
      }
      if (alignment[1][i]==0)
	aligned_line[1][i-1] = aligned_line[3][i-1] = '_';
      else {
	aligned_line[1][i-1] = bppm_symbol(T2+alignment[1][i]*3);
	aligned_line[3][i-1] = seq2[alignment[1][i]-1];
      }
   }
}

PUBLIC int set_paln_params(double gap_open, double gap_ext,
			   double seq_weight, int freeends) {
  open = (gap_open>0) ? -gap_open : gap_open;
  ext = (gap_ext>0) ? -gap_ext : gap_ext;
  if (open > ext) fprintf(stderr, "Gap extension penalty is smaller than "
			  "gap open. Do you realy want this?\n");
  seqw = seq_weight;
  if (seqw<0) {
    seqw = 0;
    fprintf(stderr, "Sequence weight set to 0 (must be in [0..1])\n");
  } else
  if (seqw>1) {
    seqw = 1;
    fprintf(stderr, "Sequence weight set to 1 (must be in [0..1])\n");
  }
  free_ends = (freeends) ? 1 : 0;
  return 0;
}

/*---------------------------------------------------------------------------*/
