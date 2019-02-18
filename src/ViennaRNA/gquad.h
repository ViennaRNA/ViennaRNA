#ifndef VIENNA_RNA_PACKAGE_GQUAD_H
#define VIENNA_RNA_PACKAGE_GQUAD_H

#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/params/basic.h>

#ifndef INLINE
#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif
#endif

/**
 *  @file       gquad.h
 *  @ingroup    paired_modules
 *  @brief      G-quadruplexes
 */

/**
 *  @addtogroup gquads
 *  @{
 *
 *  @brief Various functions related to G-quadruplex computations
 */


int         E_gquad(int           L,
                    int           l[3],
                    vrna_param_t  *P);


FLT_OR_DBL  exp_E_gquad(int               L,
                        int               l[3],
                        vrna_exp_param_t  *pf);


void E_gquad_ali_en(int           i,
                    int           L,
                    int           l[3],
                    const short   **S,
                    unsigned int  **a2s,
                    unsigned int  n_seq,
                    vrna_param_t  *P,
                    int           en[2]);


/**
 *  @brief Get a triangular matrix prefilled with minimum free energy
 *  contributions of G-quadruplexes.
 *
 *  At each position ij in the matrix, the minimum free energy of any
 *  G-quadruplex delimited by i and j is stored. If no G-quadruplex formation
 *  is possible, the matrix element is set to INF.
 *  Access the elements in the matrix via matrix[indx[j]+i]. To get
 *  the integer array indx see get_jindx().
 *
 *  @see get_jindx(), encode_sequence()
 *
 *  @param S  The encoded sequence
 *  @param P  A pointer to the data structure containing the precomputed energy contributions
 *  @return   A pointer to the G-quadruplex contribution matrix
 */
int *get_gquad_matrix(short         *S,
                      vrna_param_t  *P);


int *get_gquad_ali_matrix(short         *S_cons,
                          short         **S,
                          unsigned int  **a2s,
                          int           n_seq,
                          vrna_param_t  *P);


FLT_OR_DBL *get_gquad_pf_matrix(short             *S,
                                FLT_OR_DBL        *scale,
                                vrna_exp_param_t  *pf);


FLT_OR_DBL *get_gquad_pf_matrix_comparative(short             *S_cons,
                                            short             **S,
                                            unsigned int      **a2s,
                                            FLT_OR_DBL        *scale,
                                            unsigned int      n_seq,
                                            vrna_exp_param_t  *pf);


int **get_gquad_L_matrix(short        *S,
                         int          start,
                         int          maxdist,
                         int          n,
                         int          **g,
                         vrna_param_t *P);


void        vrna_gquad_mx_local_update(vrna_fold_compound_t *vc,
                                       int                  start);


void get_gquad_pattern_mfe(short        *S,
                           int          i,
                           int          j,
                           vrna_param_t *P,
                           int          *L,
                           int          l[3]);


void
get_gquad_pattern_exhaustive(short        *S,
                             int          i,
                             int          j,
                             vrna_param_t *P,
                             int          *L,
                             int          *l,
                             int          threshold);


void get_gquad_pattern_pf(short             *S,
                          int               i,
                          int               j,
                          vrna_exp_param_t  *pf,
                          int               *L,
                          int               l[3]);


plist *get_plist_gquad_from_pr(short            *S,
                               int              gi,
                               int              gj,
                               FLT_OR_DBL       *G,
                               FLT_OR_DBL       *probs,
                               FLT_OR_DBL       *scale,
                               vrna_exp_param_t *pf);


plist *get_plist_gquad_from_pr_max(short            *S,
                                   int              gi,
                                   int              gj,
                                   FLT_OR_DBL       *G,
                                   FLT_OR_DBL       *probs,
                                   FLT_OR_DBL       *scale,
                                   int              *L,
                                   int              l[3],
                                   vrna_exp_param_t *pf);


plist *get_plist_gquad_from_db(const char *structure,
                               float      pr);


plist *
vrna_get_plist_gquad_from_pr(vrna_fold_compound_t *fc,
                             int                  gi,
                             int                  gj);


plist *
vrna_get_plist_gquad_from_pr_max(vrna_fold_compound_t *fc,
                                 int                  gi,
                                 int                  gj,
                                 int                  *Lmax,
                                 int                  lmax[3]);


int         get_gquad_count(short *S,
                            int   i,
                            int   j);


int         get_gquad_layer_count(short *S,
                                  int   i,
                                  int   j);


void get_gquad_pattern_mfe_ali(short        **S,
                               unsigned int **a2s,
                               short        *S_cons,
                               int          n_seq,
                               int          i,
                               int          j,
                               vrna_param_t *P,
                               int          *L,
                               int          l[3]);


/**
 *  given a dot-bracket structure (possibly) containing gquads encoded
 *  by '+' signs, find first gquad, return end position or 0 if none found
 *  Upon return L and l[] contain the number of stacked layers, as well as
 *  the lengths of the linker regions.
 *  To parse a string with many gquads, call parse_gquad repeatedly e.g.
 *  end1 = parse_gquad(struc, &L, l); ... ;
 *  end2 = parse_gquad(struc+end1, &L, l); end2+=end1; ... ;
 *  end3 = parse_gquad(struc+end2, &L, l); end3+=end2; ... ;
 */
int parse_gquad(const char  *struc,
                int         *L,
                int         l[3]);


INLINE PRIVATE int backtrack_GQuad_IntLoop(int          c,
                                           int          i,
                                           int          j,
                                           int          type,
                                           short        *S,
                                           int          *ggg,
                                           int          *index,
                                           int          *p,
                                           int          *q,
                                           vrna_param_t *P);


INLINE PRIVATE int backtrack_GQuad_IntLoop_comparative(int          c,
                                                       int          i,
                                                       int          j,
                                                       unsigned int *type,
                                                       short        *S_cons,
                                                       short        **S5,
                                                       short        **S3,
                                                       unsigned int **a2s,
                                                       int          *ggg,
                                                       int          *index,
                                                       int          *p,
                                                       int          *q,
                                                       int          n_seq,
                                                       vrna_param_t *P);


INLINE PRIVATE int backtrack_GQuad_IntLoop_L(int          c,
                                             int          i,
                                             int          j,
                                             int          type,
                                             short        *S,
                                             int          **ggg,
                                             int          maxdist,
                                             int          *p,
                                             int          *q,
                                             vrna_param_t *P);


PRIVATE INLINE int
vrna_BT_gquad_int(vrna_fold_compound_t  *vc,
                  int                   i,
                  int                   j,
                  int                   en,
                  vrna_bp_stack_t       *bp_stack,
                  int                   *stack_count);


PRIVATE INLINE int
vrna_BT_gquad_mfe(vrna_fold_compound_t  *vc,
                  int                   i,
                  int                   j,
                  vrna_bp_stack_t       *bp_stack,
                  int                   *stack_count)
{
  /*
   * here we do some fancy stuff to backtrace the stacksize and linker lengths
   * of the g-quadruplex that should reside within position i,j
   */
  short         *S;
  int           l[3], L, a, n_seq;
  vrna_param_t  *P;

  if (vc) {
    P = vc->params;
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        S = vc->sequence_encoding2;
        L = -1;

        get_gquad_pattern_mfe(S, i, j, P, &L, l);
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        n_seq = vc->n_seq;
        L     = -1;
        get_gquad_pattern_mfe_ali(vc->S, vc->a2s, vc->S_cons, n_seq, i, j, P, &L, l);
        break;
    }

    if (L != -1) {
      /* fill the G's of the quadruplex into base_pair2 */
      for (a = 0; a < L; a++) {
        bp_stack[++(*stack_count)].i  = i + a;
        bp_stack[(*stack_count)].j    = i + a;
        bp_stack[++(*stack_count)].i  = i + L + l[0] + a;
        bp_stack[(*stack_count)].j    = i + L + l[0] + a;
        bp_stack[++(*stack_count)].i  = i + L + l[0] + L + l[1] + a;
        bp_stack[(*stack_count)].j    = i + L + l[0] + L + l[1] + a;
        bp_stack[++(*stack_count)].i  = i + L + l[0] + L + l[1] + L + l[2] + a;
        bp_stack[(*stack_count)].j    = i + L + l[0] + L + l[1] + L + l[2] + a;
      }
      return 1;
    } else {
      return 0;
    }
  }

  return 0;
}


PRIVATE INLINE int
vrna_BT_gquad_int(vrna_fold_compound_t  *vc,
                  int                   i,
                  int                   j,
                  int                   en,
                  vrna_bp_stack_t       *bp_stack,
                  int                   *stack_count)
{
  int           energy, dangles, *idx, ij, p, q, maxl, minl, c0, l1, *ggg;
  unsigned char type;
  char          *ptype;
  short         si, sj, *S, *S1;

  vrna_param_t  *P;
  vrna_md_t     *md;

  idx     = vc->jindx;
  ij      = idx[j] + i;
  P       = vc->params;
  md      = &(P->model_details);
  ptype   = vc->ptype;
  type    = (unsigned char)ptype[ij];
  S1      = vc->sequence_encoding;
  S       = vc->sequence_encoding2;
  dangles = md->dangles;
  si      = S1[i + 1];
  sj      = S1[j - 1];
  ggg     = vc->matrices->ggg;
  energy  = 0;

  if (dangles == 2)
    energy += P->mismatchI[type][si][sj];

  if (type > 2)
    energy += P->TerminalAU;

  p = i + 1;
  if (S1[p] == 3) {
    if (p < j - VRNA_GQUAD_MIN_BOX_SIZE) {
      minl  = j - i + p - MAXLOOP - 2;
      c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      minl  = MAX2(c0, minl);
      c0    = j - 3;
      maxl  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      maxl  = MIN2(c0, maxl);
      for (q = minl; q < maxl; q++) {
        if (S[q] != 3)
          continue;

        if (en == energy + ggg[idx[q] + p] + P->internal_loop[j - q - 1])
          return vrna_BT_gquad_mfe(vc, p, q, bp_stack, stack_count);
      }
    }
  }

  for (p = i + 2;
       p < j - VRNA_GQUAD_MIN_BOX_SIZE;
       p++) {
    l1 = p - i - 1;
    if (l1 > MAXLOOP)
      break;

    if (S1[p] != 3)
      continue;

    minl  = j - i + p - MAXLOOP - 2;
    c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
    minl  = MAX2(c0, minl);
    c0    = j - 1;
    maxl  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
    maxl  = MIN2(c0, maxl);
    for (q = minl; q < maxl; q++) {
      if (S1[q] != 3)
        continue;

      if (en == energy + ggg[idx[q] + p] + P->internal_loop[l1 + j - q - 1])
        return vrna_BT_gquad_mfe(vc, p, q, bp_stack, stack_count);
    }
  }

  q = j - 1;
  if (S1[q] == 3)
    for (p = i + 4;
         p < j - VRNA_GQUAD_MIN_BOX_SIZE;
         p++) {
      l1 = p - i - 1;
      if (l1 > MAXLOOP)
        break;

      if (S1[p] != 3)
        continue;

      if (en == energy + ggg[idx[q] + p] + P->internal_loop[l1])
        return vrna_BT_gquad_mfe(vc, p, q, bp_stack, stack_count);
    }

  return 0;
}


/**
 *  backtrack an interior loop like enclosed g-quadruplex
 *  with closing pair (i,j)
 *
 *  @param c      The total contribution the loop should resemble
 *  @param i      position i of enclosing pair
 *  @param j      position j of enclosing pair
 *  @param type   base pair type of enclosing pair (must be reverse type)
 *  @param S      integer encoded sequence
 *  @param ggg    triangular matrix containing g-quadruplex contributions
 *  @param index  the index for accessing the triangular matrix
 *  @param p      here the 5' position of the gquad is stored
 *  @param q      here the 3' position of the gquad is stored
 *  @param P      the datastructure containing the precalculated contibutions
 *
 *  @return       1 on success, 0 if no gquad found
 */
INLINE PRIVATE int
backtrack_GQuad_IntLoop(int           c,
                        int           i,
                        int           j,
                        int           type,
                        short         *S,
                        int           *ggg,
                        int           *index,
                        int           *p,
                        int           *q,
                        vrna_param_t  *P)
{
  int   energy, dangles, k, l, maxl, minl, c0, l1;
  short si, sj;

  dangles = P->model_details.dangles;
  si      = S[i + 1];
  sj      = S[j - 1];
  energy  = 0;

  if (dangles == 2)
    energy += P->mismatchI[type][si][sj];

  if (type > 2)
    energy += P->TerminalAU;

  k = i + 1;
  if (S[k] == 3) {
    if (k < j - VRNA_GQUAD_MIN_BOX_SIZE) {
      minl  = j - i + k - MAXLOOP - 2;
      c0    = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      minl  = MAX2(c0, minl);
      c0    = j - 3;
      maxl  = k + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      maxl  = MIN2(c0, maxl);
      for (l = minl; l < maxl; l++) {
        if (S[l] != 3)
          continue;

        if (c == energy + ggg[index[l] + k] + P->internal_loop[j - l - 1]) {
          *p  = k;
          *q  = l;
          return 1;
        }
      }
    }
  }

  for (k = i + 2;
       k < j - VRNA_GQUAD_MIN_BOX_SIZE;
       k++) {
    l1 = k - i - 1;
    if (l1 > MAXLOOP)
      break;

    if (S[k] != 3)
      continue;

    minl  = j - i + k - MAXLOOP - 2;
    c0    = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
    minl  = MAX2(c0, minl);
    c0    = j - 1;
    maxl  = k + VRNA_GQUAD_MAX_BOX_SIZE + 1;
    maxl  = MIN2(c0, maxl);
    for (l = minl; l < maxl; l++) {
      if (S[l] != 3)
        continue;

      if (c == energy + ggg[index[l] + k] + P->internal_loop[l1 + j - l - 1]) {
        *p  = k;
        *q  = l;
        return 1;
      }
    }
  }

  l = j - 1;
  if (S[l] == 3)
    for (k = i + 4;
         k < j - VRNA_GQUAD_MIN_BOX_SIZE;
         k++) {
      l1 = k - i - 1;
      if (l1 > MAXLOOP)
        break;

      if (S[k] != 3)
        continue;

      if (c == energy + ggg[index[l] + k] + P->internal_loop[l1]) {
        *p  = k;
        *q  = l;
        return 1;
      }
    }

  return 0;
}


INLINE PRIVATE int
backtrack_GQuad_IntLoop_comparative(int           c,
                                    int           i,
                                    int           j,
                                    unsigned int  *type,
                                    short         *S_cons,
                                    short         **S5,
                                    short         **S3,
                                    unsigned int  **a2s,
                                    int           *ggg,
                                    int           *index,
                                    int           *p,
                                    int           *q,
                                    int           n_seq,
                                    vrna_param_t  *P)
{
  int energy, dangles, k, l, maxl, minl, c0, l1, ss, tt, u1, u2, eee;

  dangles = P->model_details.dangles;
  energy  = 0;

  for (ss = 0; ss < n_seq; ss++) {
    tt = type[ss];
    if (tt == 0)
      tt = 7;

    if (dangles == 2)
      energy += P->mismatchI[tt][S3[ss][i]][S5[ss][j]];

    if (tt > 2)
      energy += P->TerminalAU;
  }

  k = i + 1;
  if (S_cons[k] == 3) {
    if (k < j - VRNA_GQUAD_MIN_BOX_SIZE) {
      minl  = j - i + k - MAXLOOP - 2;
      c0    = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      minl  = MAX2(c0, minl);
      c0    = j - 3;
      maxl  = k + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      maxl  = MIN2(c0, maxl);
      for (l = minl; l < maxl; l++) {
        if (S_cons[l] != 3)
          continue;

        eee = 0;

        for (ss = 0; ss < n_seq; ss++) {
          u1  = a2s[ss][j - 1] - a2s[ss][l];
          eee += P->internal_loop[u1];
        }

        if (c == energy + ggg[index[l] + k] + eee) {
          *p  = k;
          *q  = l;
          return 1;
        }
      }
    }
  }

  for (k = i + 2;
       k < j - VRNA_GQUAD_MIN_BOX_SIZE;
       k++) {
    l1 = k - i - 1;
    if (l1 > MAXLOOP)
      break;

    if (S_cons[k] != 3)
      continue;

    minl  = j - i + k - MAXLOOP - 2;
    c0    = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
    minl  = MAX2(c0, minl);
    c0    = j - 1;
    maxl  = k + VRNA_GQUAD_MAX_BOX_SIZE + 1;
    maxl  = MIN2(c0, maxl);
    for (l = minl; l < maxl; l++) {
      if (S_cons[l] != 3)
        continue;

      eee = 0;

      for (ss = 0; ss < n_seq; ss++) {
        u1  = a2s[ss][k - 1] - a2s[ss][i];
        u2  = a2s[ss][j - 1] - a2s[ss][l];
        eee += P->internal_loop[u1 + u2];
      }

      if (c == energy + ggg[index[l] + k] + eee) {
        *p  = k;
        *q  = l;
        return 1;
      }
    }
  }

  l = j - 1;
  if (S_cons[l] == 3)
    for (k = i + 4;
         k < j - VRNA_GQUAD_MIN_BOX_SIZE;
         k++) {
      l1 = k - i - 1;
      if (l1 > MAXLOOP)
        break;

      if (S_cons[k] != 3)
        continue;

      eee = 0;

      for (ss = 0; ss < n_seq; ss++) {
        u1  = a2s[ss][k - 1] - a2s[ss][i];
        eee += P->internal_loop[u1];
      }

      if (c == energy + ggg[index[l] + k] + eee) {
        *p  = k;
        *q  = l;
        return 1;
      }
    }

  return 0;
}


/**
 *  backtrack an interior loop like enclosed g-quadruplex
 *  with closing pair (i,j) with underlying Lfold matrix
 *
 *  @param c      The total contribution the loop should resemble
 *  @param i      position i of enclosing pair
 *  @param j      position j of enclosing pair
 *  @param type   base pair type of enclosing pair (must be reverse type)
 *  @param S      integer encoded sequence
 *  @param ggg    triangular matrix containing g-quadruplex contributions
 *  @param p      here the 5' position of the gquad is stored
 *  @param q      here the 3' position of the gquad is stored
 *  @param P      the datastructure containing the precalculated contibutions
 *
 *  @return       1 on success, 0 if no gquad found
 */
INLINE PRIVATE int
backtrack_GQuad_IntLoop_L(int           c,
                          int           i,
                          int           j,
                          int           type,
                          short         *S,
                          int           **ggg,
                          int           maxdist,
                          int           *p,
                          int           *q,
                          vrna_param_t  *P)
{
  int   energy, dangles, k, l, maxl, minl, c0, l1;
  short si, sj;

  dangles = P->model_details.dangles;
  si      = S[i + 1];
  sj      = S[j - 1];
  energy  = 0;

  if (dangles == 2)
    energy += P->mismatchI[type][si][sj];

  if (type > 2)
    energy += P->TerminalAU;

  k = i + 1;
  if (S[k] == 3) {
    if (k < j - VRNA_GQUAD_MIN_BOX_SIZE) {
      minl  = j - i + k - MAXLOOP - 2;
      c0    = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      minl  = MAX2(c0, minl);
      c0    = j - 3;
      maxl  = k + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      maxl  = MIN2(c0, maxl);
      for (l = minl; l < maxl; l++) {
        if (S[l] != 3)
          continue;

        if (c == energy + ggg[k][l - k] + P->internal_loop[j - l - 1]) {
          *p  = k;
          *q  = l;
          return 1;
        }
      }
    }
  }

  for (k = i + 2;
       k < j - VRNA_GQUAD_MIN_BOX_SIZE;
       k++) {
    l1 = k - i - 1;
    if (l1 > MAXLOOP)
      break;

    if (S[k] != 3)
      continue;

    minl  = j - i + k - MAXLOOP - 2;
    c0    = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
    minl  = MAX2(c0, minl);
    c0    = j - 1;
    maxl  = k + VRNA_GQUAD_MAX_BOX_SIZE + 1;
    maxl  = MIN2(c0, maxl);
    for (l = minl; l < maxl; l++) {
      if (S[l] != 3)
        continue;

      if (c == energy + ggg[k][l - k] + P->internal_loop[l1 + j - l - 1]) {
        *p  = k;
        *q  = l;
        return 1;
      }
    }
  }

  l = j - 1;
  if (S[l] == 3)
    for (k = i + 4;
         k < j - VRNA_GQUAD_MIN_BOX_SIZE;
         k++) {
      l1 = k - i - 1;
      if (l1 > MAXLOOP)
        break;

      if (S[k] != 3)
        continue;

      if (c == energy + ggg[k][l - k] + P->internal_loop[l1]) {
        *p  = k;
        *q  = l;
        return 1;
      }
    }

  return 0;
}


INLINE PRIVATE int
backtrack_GQuad_IntLoop_L_comparative(int           c,
                                      int           i,
                                      int           j,
                                      unsigned int  *type,
                                      short         *S_cons,
                                      short         **S5,
                                      short         **S3,
                                      unsigned int  **a2s,
                                      int           **ggg,
                                      int           *p,
                                      int           *q,
                                      int           n_seq,
                                      vrna_param_t  *P)
{
  /*
   * The case that is handled here actually resembles something like
   * an interior loop where the enclosing base pair is of regular
   * kind and the enclosed pair is not a canonical one but a g-quadruplex
   * that should then be decomposed further...
   */
  int mm, dangle_model, k, l, maxl, minl, c0, l1, ss, tt, eee, u1, u2;

  dangle_model = P->model_details.dangles;

  mm = 0;
  for (ss = 0; ss < n_seq; ss++) {
    tt = type[ss];

    if (dangle_model == 2)
      mm += P->mismatchI[tt][S3[ss][i]][S5[ss][j]];

    if (tt > 2)
      mm += P->TerminalAU;
  }

  for (k = i + 2;
       k < j - VRNA_GQUAD_MIN_BOX_SIZE;
       k++) {
    if (S_cons[k] != 3)
      continue;

    l1 = k - i - 1;
    if (l1 > MAXLOOP)
      break;

    minl  = j - i + k - MAXLOOP - 2;
    c0    = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
    minl  = MAX2(c0, minl);
    c0    = j - 1;
    maxl  = k + VRNA_GQUAD_MAX_BOX_SIZE + 1;
    maxl  = MIN2(c0, maxl);
    for (l = minl; l < maxl; l++) {
      if (S_cons[l] != 3)
        continue;

      eee = 0;

      for (ss = 0; ss < n_seq; ss++) {
        u1  = a2s[ss][k - 1] - a2s[ss][i];
        u2  = a2s[ss][j - 1] - a2s[ss][l];
        eee += P->internal_loop[u1 + u2];
      }

      c0 = mm +
           ggg[k][l - k] +
           eee;

      if (c == c0) {
        *p  = k;
        *q  = l;
        return 1;
      }
    }
  }
  k = i + 1;
  if (S_cons[k] == 3) {
    if (k < j - VRNA_GQUAD_MIN_BOX_SIZE) {
      minl  = j - i + k - MAXLOOP - 2;
      c0    = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      minl  = MAX2(c0, minl);
      c0    = j - 3;
      maxl  = k + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      maxl  = MIN2(c0, maxl);
      for (l = minl; l < maxl; l++) {
        if (S_cons[l] != 3)
          continue;

        eee = 0;

        for (ss = 0; ss < n_seq; ss++) {
          u1  = a2s[ss][j - 1] - a2s[ss][l];
          eee += P->internal_loop[u1];
        }

        if (c == mm + ggg[k][l - k] + eee) {
          *p  = k;
          *q  = l;
          return 1;
        }
      }
    }
  }

  l = j - 1;
  if (S_cons[l] == 3) {
    for (k = i + 4; k < j - VRNA_GQUAD_MIN_BOX_SIZE; k++) {
      l1 = k - i - 1;
      if (l1 > MAXLOOP)
        break;

      if (S_cons[k] != 3)
        continue;

      eee = 0;

      for (ss = 0; ss < n_seq; ss++) {
        u1  = a2s[ss][k - 1] - a2s[ss][i];
        eee += P->internal_loop[u1];
      }

      if (c == mm + ggg[k][l - k] + eee) {
        *p  = k;
        *q  = l;
        return 1;
      }
    }
  }

  return 0;
}


PRIVATE INLINE
int
E_GQuad_IntLoop(int           i,
                int           j,
                int           type,
                short         *S,
                int           *ggg,
                int           *index,
                vrna_param_t  *P)
{
  int   energy, ge, dangles, p, q, l1, minq, maxq, c0;
  short si, sj;

  dangles = P->model_details.dangles;
  si      = S[i + 1];
  sj      = S[j - 1];
  energy  = 0;

  if (dangles == 2)
    energy += P->mismatchI[type][si][sj];

  if (type > 2)
    energy += P->TerminalAU;

  ge = INF;

  p = i + 1;
  if (S[p] == 3) {
    if (p < j - VRNA_GQUAD_MIN_BOX_SIZE) {
      minq  = j - i + p - MAXLOOP - 2;
      c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      minq  = MAX2(c0, minq);
      c0    = j - 3;
      maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      maxq  = MIN2(c0, maxq);
      for (q = minq; q < maxq; q++) {
        if (S[q] != 3)
          continue;

        c0  = energy + ggg[index[q] + p] + P->internal_loop[j - q - 1];
        ge  = MIN2(ge, c0);
      }
    }
  }

  for (p = i + 2;
       p < j - VRNA_GQUAD_MIN_BOX_SIZE;
       p++) {
    l1 = p - i - 1;
    if (l1 > MAXLOOP)
      break;

    if (S[p] != 3)
      continue;

    minq  = j - i + p - MAXLOOP - 2;
    c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
    minq  = MAX2(c0, minq);
    c0    = j - 1;
    maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
    maxq  = MIN2(c0, maxq);
    for (q = minq; q < maxq; q++) {
      if (S[q] != 3)
        continue;

      c0  = energy + ggg[index[q] + p] + P->internal_loop[l1 + j - q - 1];
      ge  = MIN2(ge, c0);
    }
  }

  q = j - 1;
  if (S[q] == 3)
    for (p = i + 4;
         p < j - VRNA_GQUAD_MIN_BOX_SIZE;
         p++) {
      l1 = p - i - 1;
      if (l1 > MAXLOOP)
        break;

      if (S[p] != 3)
        continue;

      c0  = energy + ggg[index[q] + p] + P->internal_loop[l1];
      ge  = MIN2(ge, c0);
    }

#if 0
  /* here comes the additional stuff for the odd dangle models */
  if (dangles % 1) {
    en1 = energy + P->dangle5[type][si];
    en2 = energy + P->dangle5[type][sj];
    en3 = energy + P->mismatchI[type][si][sj];

    /* first case with 5' dangle (i.e. j-1) onto enclosing pair */
    p = i + 1;
    if (S[p] == 3) {
      if (p < j - VRNA_GQUAD_MIN_BOX_SIZE) {
        minq  = j - i + p - MAXLOOP - 2;
        c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
        minq  = MAX2(c0, minq);
        c0    = j - 4;
        maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
        maxq  = MIN2(c0, maxq);
        for (q = minq; q < maxq; q++) {
          if (S[q] != 3)
            continue;

          c0  = en1 + ggg[index[q] + p] + P->internal_loop[j - q - 1];
          ge  = MIN2(ge, c0);
        }
      }
    }

    for (p = i + 2; p < j - VRNA_GQUAD_MIN_BOX_SIZE; p++) {
      l1 = p - i - 1;
      if (l1 > MAXLOOP)
        break;

      if (S[p] != 3)
        continue;

      minq  = j - i + p - MAXLOOP - 2;
      c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      minq  = MAX2(c0, minq);
      c0    = j - 2;
      maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      maxq  = MIN2(c0, maxq);
      for (q = minq; q < maxq; q++) {
        if (S[q] != 3)
          continue;

        c0  = en1 + ggg[index[q] + p] + P->internal_loop[l1 + j - q - 1];
        ge  = MIN2(ge, c0);
      }
    }

    q = j - 2;
    if (S[q] == 3)
      for (p = i + 4; p < j - VRNA_GQUAD_MIN_BOX_SIZE; p++) {
        l1 = p - i - 1;
        if (l1 > MAXLOOP)
          break;

        if (S[p] != 3)
          continue;

        c0  = en1 + ggg[index[q] + p] + P->internal_loop[l1 + 1];
        ge  = MIN2(ge, c0);
      }

    /* second case with 3' dangle (i.e. i+1) onto enclosing pair */
  }

#endif
  return ge;
}


PRIVATE INLINE
int
E_GQuad_IntLoop_comparative(int           i,
                            int           j,
                            unsigned int  *tt,
                            short         *S_cons,
                            short         **S5,
                            short         **S3,
                            unsigned int  **a2s,
                            int           *ggg,
                            int           *index,
                            int           n_seq,
                            vrna_param_t  *P)
{
  unsigned int  type;
  int           eee, energy, ge, p, q, l1, u1, u2, minq, maxq, c0, s;
  vrna_md_t     *md;

  md      = &(P->model_details);
  energy  = 0;

  for (s = 0; s < n_seq; s++) {
    type = tt[s];
    if (md->dangles == 2)
      energy += P->mismatchI[type][S3[s][i]][S5[s][j]];

    if (type > 2)
      energy += P->TerminalAU;
  }

  ge = INF;

  p = i + 1;
  if (S_cons[p] == 3) {
    if (p < j - VRNA_GQUAD_MIN_BOX_SIZE) {
      minq  = j - i + p - MAXLOOP - 2;
      c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      minq  = MAX2(c0, minq);
      c0    = j - 3;
      maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      maxq  = MIN2(c0, maxq);
      for (q = minq; q < maxq; q++) {
        if (S_cons[q] != 3)
          continue;

        eee = 0;

        for (s = 0; s < n_seq; s++) {
          u1  = a2s[s][j - 1] - a2s[s][q];
          eee += P->internal_loop[u1];
        }

        c0 = energy +
             ggg[index[q] + p] +
             eee;
        ge = MIN2(ge, c0);
      }
    }
  }

  for (p = i + 2;
       p < j - VRNA_GQUAD_MIN_BOX_SIZE;
       p++) {
    l1 = p - i - 1;
    if (l1 > MAXLOOP)
      break;

    if (S_cons[p] != 3)
      continue;

    minq  = j - i + p - MAXLOOP - 2;
    c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
    minq  = MAX2(c0, minq);
    c0    = j - 1;
    maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
    maxq  = MIN2(c0, maxq);
    for (q = minq; q < maxq; q++) {
      if (S_cons[q] != 3)
        continue;

      eee = 0;

      for (s = 0; s < n_seq; s++) {
        u1  = a2s[s][p - 1] - a2s[s][i];
        u2  = a2s[s][j - 1] - a2s[s][q];
        eee += P->internal_loop[u1 + u2];
      }

      c0 = energy +
           ggg[index[q] + p] +
           eee;
      ge = MIN2(ge, c0);
    }
  }

  q = j - 1;
  if (S_cons[q] == 3)
    for (p = i + 4;
         p < j - VRNA_GQUAD_MIN_BOX_SIZE;
         p++) {
      l1 = p - i - 1;
      if (l1 > MAXLOOP)
        break;

      if (S_cons[p] != 3)
        continue;

      eee = 0;

      for (s = 0; s < n_seq; s++) {
        u1  = a2s[s][p - 1] - a2s[s][i];
        eee += P->internal_loop[u1];
      }

      c0 = energy +
           ggg[index[q] + p] +
           eee;
      ge = MIN2(ge, c0);
    }

  return ge;
}


PRIVATE INLINE
int
E_GQuad_IntLoop_L_comparative(int           i,
                              int           j,
                              unsigned int  *tt,
                              short         *S_cons,
                              short         **S5,
                              short         **S3,
                              unsigned int  **a2s,
                              int           **ggg,
                              int           n_seq,
                              vrna_param_t  *P)
{
  unsigned int  type;
  int           eee, energy, ge, p, q, l1, u1, u2, minq, maxq, c0, s;
  vrna_md_t     *md;

  md      = &(P->model_details);
  energy  = 0;

  for (s = 0; s < n_seq; s++) {
    type = tt[s];
    if (md->dangles == 2)
      energy += P->mismatchI[type][S3[s][i]][S5[s][j]];

    if (type > 2)
      energy += P->TerminalAU;
  }

  ge = INF;

  p = i + 1;
  if (S_cons[p] == 3) {
    if (p < j - VRNA_GQUAD_MIN_BOX_SIZE) {
      minq  = j - i + p - MAXLOOP - 2;
      c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      minq  = MAX2(c0, minq);
      c0    = j - 3;
      maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      maxq  = MIN2(c0, maxq);
      for (q = minq; q < maxq; q++) {
        if (S_cons[q] != 3)
          continue;

        eee = 0;

        for (s = 0; s < n_seq; s++) {
          u1  = a2s[s][j - 1] - a2s[s][q];
          eee += P->internal_loop[u1];
        }

        c0 = energy +
             ggg[p][q - p] +
             eee;
        ge = MIN2(ge, c0);
      }
    }
  }

  for (p = i + 2;
       p < j - VRNA_GQUAD_MIN_BOX_SIZE;
       p++) {
    l1 = p - i - 1;
    if (l1 > MAXLOOP)
      break;

    if (S_cons[p] != 3)
      continue;

    minq  = j - i + p - MAXLOOP - 2;
    c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
    minq  = MAX2(c0, minq);
    c0    = j - 1;
    maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
    maxq  = MIN2(c0, maxq);
    for (q = minq; q < maxq; q++) {
      if (S_cons[q] != 3)
        continue;

      eee = 0;

      for (s = 0; s < n_seq; s++) {
        u1  = a2s[s][p - 1] - a2s[s][i];
        u2  = a2s[s][j - 1] - a2s[s][q];
        eee += P->internal_loop[u1 + u2];
      }

      c0 = energy +
           ggg[p][q - p] +
           eee;
      ge = MIN2(ge, c0);
    }
  }

  q = j - 1;
  if (S_cons[q] == 3)
    for (p = i + 4;
         p < j - VRNA_GQUAD_MIN_BOX_SIZE;
         p++) {
      l1 = p - i - 1;
      if (l1 > MAXLOOP)
        break;

      if (S_cons[p] != 3)
        continue;

      eee = 0;

      for (s = 0; s < n_seq; s++) {
        u1  = a2s[s][p - 1] - a2s[s][i];
        eee += P->internal_loop[u1];
      }

      c0 = energy +
           ggg[p][q - p] +
           eee;
      ge = MIN2(ge, c0);
    }

  return ge;
}


PRIVATE INLINE
int *
E_GQuad_IntLoop_exhaustive(int          i,
                           int          j,
                           int          **p_p,
                           int          **q_p,
                           int          type,
                           short        *S,
                           int          *ggg,
                           int          threshold,
                           int          *index,
                           vrna_param_t *P)
{
  int   energy, *ge, dangles, p, q, l1, minq, maxq, c0;
  short si, sj;
  int   cnt = 0;

  dangles = P->model_details.dangles;
  si      = S[i + 1];
  sj      = S[j - 1];
  energy  = 0;

  if (dangles == 2)
    energy += P->mismatchI[type][si][sj];

  if (type > 2)
    energy += P->TerminalAU;

  /* guess how many gquads are possible in interval [i+1,j-1] */
  *p_p  = (int *)vrna_alloc(sizeof(int) * 256);
  *q_p  = (int *)vrna_alloc(sizeof(int) * 256);
  ge    = (int *)vrna_alloc(sizeof(int) * 256);

  p = i + 1;
  if (S[p] == 3) {
    if (p < j - VRNA_GQUAD_MIN_BOX_SIZE) {
      minq  = j - i + p - MAXLOOP - 2;
      c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      minq  = MAX2(c0, minq);
      c0    = j - 3;
      maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      maxq  = MIN2(c0, maxq);
      for (q = minq; q < maxq; q++) {
        if (S[q] != 3)
          continue;

        c0 = energy + ggg[index[q] + p] + P->internal_loop[j - q - 1];
        if (c0 <= threshold) {
          ge[cnt]       = energy + P->internal_loop[j - q - 1];
          (*p_p)[cnt]   = p;
          (*q_p)[cnt++] = q;
        }
      }
    }
  }

  for (p = i + 2;
       p < j - VRNA_GQUAD_MIN_BOX_SIZE;
       p++) {
    l1 = p - i - 1;
    if (l1 > MAXLOOP)
      break;

    if (S[p] != 3)
      continue;

    minq  = j - i + p - MAXLOOP - 2;
    c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
    minq  = MAX2(c0, minq);
    c0    = j - 1;
    maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
    maxq  = MIN2(c0, maxq);
    for (q = minq; q < maxq; q++) {
      if (S[q] != 3)
        continue;

      c0 = energy + ggg[index[q] + p] + P->internal_loop[l1 + j - q - 1];
      if (c0 <= threshold) {
        ge[cnt]       = energy + P->internal_loop[l1 + j - q - 1];
        (*p_p)[cnt]   = p;
        (*q_p)[cnt++] = q;
      }
    }
  }

  q = j - 1;
  if (S[q] == 3)
    for (p = i + 4;
         p < j - VRNA_GQUAD_MIN_BOX_SIZE;
         p++) {
      l1 = p - i - 1;
      if (l1 > MAXLOOP)
        break;

      if (S[p] != 3)
        continue;

      c0 = energy + ggg[index[q] + p] + P->internal_loop[l1];
      if (c0 <= threshold) {
        ge[cnt]       = energy + P->internal_loop[l1];
        (*p_p)[cnt]   = p;
        (*q_p)[cnt++] = q;
      }
    }

  (*p_p)[cnt] = -1;

  return ge;
}


PRIVATE INLINE
int
E_GQuad_IntLoop_L(int           i,
                  int           j,
                  int           type,
                  short         *S,
                  int           **ggg,
                  int           maxdist,
                  vrna_param_t  *P)
{
  int   energy, ge, dangles, p, q, l1, minq, maxq, c0;
  short si, sj;

  dangles = P->model_details.dangles;
  si      = S[i + 1];
  sj      = S[j - 1];
  energy  = 0;

  if (dangles == 2)
    energy += P->mismatchI[type][si][sj];

  if (type > 2)
    energy += P->TerminalAU;

  ge = INF;

  p = i + 1;
  if (S[p] == 3) {
    if (p < j - VRNA_GQUAD_MIN_BOX_SIZE) {
      minq  = j - i + p - MAXLOOP - 2;
      c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      minq  = MAX2(c0, minq);
      c0    = j - 3;
      maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      maxq  = MIN2(c0, maxq);
      for (q = minq; q < maxq; q++) {
        if (S[q] != 3)
          continue;

        c0  = energy + ggg[p][q - p] + P->internal_loop[j - q - 1];
        ge  = MIN2(ge, c0);
      }
    }
  }

  for (p = i + 2;
       p < j - VRNA_GQUAD_MIN_BOX_SIZE;
       p++) {
    l1 = p - i - 1;
    if (l1 > MAXLOOP)
      break;

    if (S[p] != 3)
      continue;

    minq  = j - i + p - MAXLOOP - 2;
    c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
    minq  = MAX2(c0, minq);
    c0    = j - 1;
    maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
    maxq  = MIN2(c0, maxq);
    for (q = minq; q < maxq; q++) {
      if (S[q] != 3)
        continue;

      c0  = energy + ggg[p][q - p] + P->internal_loop[l1 + j - q - 1];
      ge  = MIN2(ge, c0);
    }
  }

  q = j - 1;
  if (S[q] == 3)
    for (p = i + 4;
         p < j - VRNA_GQUAD_MIN_BOX_SIZE;
         p++) {
      l1 = p - i - 1;
      if (l1 > MAXLOOP)
        break;

      if (S[p] != 3)
        continue;

      c0  = energy + ggg[p][q - p] + P->internal_loop[l1];
      ge  = MIN2(ge, c0);
    }

  return ge;
}


PRIVATE INLINE
FLT_OR_DBL
exp_E_GQuad_IntLoop(int               i,
                    int               j,
                    int               type,
                    short             *S,
                    FLT_OR_DBL        *G,
                    FLT_OR_DBL        *scale,
                    int               *index,
                    vrna_exp_param_t  *pf)
{
  int         k, l, minl, maxl, u, r;
  FLT_OR_DBL  q, qe;
  double      *expintern;
  short       si, sj;

  q         = 0;
  si        = S[i + 1];
  sj        = S[j - 1];
  qe        = (FLT_OR_DBL)pf->expmismatchI[type][si][sj];
  expintern = &(pf->expinternal[0]);

  if (type > 2)
    qe *= (FLT_OR_DBL)pf->expTermAU;

  k = i + 1;
  if (S[k] == 3) {
    if (k < j - VRNA_GQUAD_MIN_BOX_SIZE) {
      minl  = j - MAXLOOP - 1;
      u     = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      minl  = MAX2(u, minl);
      u     = j - 3;
      maxl  = k + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      maxl  = MIN2(u, maxl);
      for (l = minl; l < maxl; l++) {
        if (S[l] != 3)
          continue;

        if (G[index[k] - l] == 0.)
          continue;

        q += qe
             * G[index[k] - l]
             * (FLT_OR_DBL)expintern[j - l - 1]
             * scale[j - l + 1];
      }
    }
  }

  for (k = i + 2;
       k <= j - VRNA_GQUAD_MIN_BOX_SIZE;
       k++) {
    u = k - i - 1;
    if (u > MAXLOOP)
      break;

    if (S[k] != 3)
      continue;

    minl  = j - i + k - MAXLOOP - 2;
    r     = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
    minl  = MAX2(r, minl);
    maxl  = k + VRNA_GQUAD_MAX_BOX_SIZE + 1;
    r     = j - 1;
    maxl  = MIN2(r, maxl);
    for (l = minl; l < maxl; l++) {
      if (S[l] != 3)
        continue;

      if (G[index[k] - l] == 0.)
        continue;

      q += qe
           * G[index[k] - l]
           * (FLT_OR_DBL)expintern[u + j - l - 1]
           * scale[u + j - l + 1];
    }
  }

  l = j - 1;
  if (S[l] == 3)
    for (k = i + 4; k <= j - VRNA_GQUAD_MIN_BOX_SIZE; k++) {
      u = k - i - 1;
      if (u > MAXLOOP)
        break;

      if (S[k] != 3)
        continue;

      if (G[index[k] - l] == 0.)
        continue;

      q += qe
           * G[index[k] - l]
           * (FLT_OR_DBL)expintern[u]
           * scale[u + 2];
    }

  return q;
}


PRIVATE INLINE
FLT_OR_DBL
exp_E_GQuad_IntLoop_comparative(int               i,
                                int               j,
                                unsigned int      *tt,
                                short             *S_cons,
                                short             **S5,
                                short             **S3,
                                unsigned int      **a2s,
                                FLT_OR_DBL        *G,
                                FLT_OR_DBL        *scale,
                                int               *index,
                                int               n_seq,
                                vrna_exp_param_t  *pf)
{
  unsigned int  type;
  int           k, l, minl, maxl, u, u1, u2, r, s;
  FLT_OR_DBL    q, qe, qqq;
  double        *expintern;
  vrna_md_t     *md;

  q         = 0;
  qe        = 1.;
  md        = &(pf->model_details);
  expintern = &(pf->expinternal[0]);

  for (s = 0; s < n_seq; s++) {
    type = tt[s];
    if (md->dangles == 2)
      qe *= (FLT_OR_DBL)pf->expmismatchI[type][S3[s][i]][S5[s][j]];

    if (type > 2)
      qe *= (FLT_OR_DBL)pf->expTermAU;
  }

  k = i + 1;
  if (S_cons[k] == 3) {
    if (k < j - VRNA_GQUAD_MIN_BOX_SIZE) {
      minl  = j - MAXLOOP - 1;
      u     = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      minl  = MAX2(u, minl);
      u     = j - 3;
      maxl  = k + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      maxl  = MIN2(u, maxl);
      for (l = minl; l < maxl; l++) {
        if (S_cons[l] != 3)
          continue;

        if (G[index[k] - l] == 0.)
          continue;

        qqq = 1.;

        for (s = 0; s < n_seq; s++) {
          u1  = a2s[s][j - 1] - a2s[s][l];
          qqq *= expintern[u1];
        }

        q += qe *
             G[index[k] - l] *
             qqq *
             scale[j - l + 1];
      }
    }
  }

  for (k = i + 2;
       k <= j - VRNA_GQUAD_MIN_BOX_SIZE;
       k++) {
    u = k - i - 1;
    if (u > MAXLOOP)
      break;

    if (S_cons[k] != 3)
      continue;

    minl  = j - i + k - MAXLOOP - 2;
    r     = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
    minl  = MAX2(r, minl);
    maxl  = k + VRNA_GQUAD_MAX_BOX_SIZE + 1;
    r     = j - 1;
    maxl  = MIN2(r, maxl);
    for (l = minl; l < maxl; l++) {
      if (S_cons[l] != 3)
        continue;

      if (G[index[k] - l] == 0.)
        continue;

      qqq = 1.;

      for (s = 0; s < n_seq; s++) {
        u1  = a2s[s][k - 1] - a2s[s][i];
        u2  = a2s[s][j - 1] - a2s[s][l];
        qqq *= expintern[u1 + u2];
      }

      q += qe *
           G[index[k] - l] *
           qqq *
           scale[u + j - l + 1];
    }
  }

  l = j - 1;
  if (S_cons[l] == 3)
    for (k = i + 4; k <= j - VRNA_GQUAD_MIN_BOX_SIZE; k++) {
      u = k - i - 1;
      if (u > MAXLOOP)
        break;

      if (S_cons[k] != 3)
        continue;

      if (G[index[k] - l] == 0.)
        continue;

      qqq = 1.;

      for (s = 0; s < n_seq; s++) {
        u1  = a2s[s][k - 1] - a2s[s][i];
        qqq *= expintern[u1];
      }

      q += qe *
           G[index[k] - l] *
           qqq *
           scale[u + 2];
    }

  return q;
}


/**
 * @}
 */


#endif
