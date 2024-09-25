/*
 *  ViennaRNA/utils/structures.c
 *
 *  Various functions to convert, parse, encode secondary structures
 *
 *  c  Ivo L Hofacker, Walter Fontana, Ronny Lorenz
 *              Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/structures/dotbracket.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE INLINE void
flatten_brackets(char       *string,
                 const char pair[3],
                 const char target[3]);


PRIVATE void
assign_elements_pair(short  *pt,
                     int    i,
                     int    j,
                     char   *elements);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC char *
vrna_db_pack(const char *struc)
{
  /* 5:1 compression using base 3 encoding */
  int           i, j, l, pi;
  unsigned char *packed;

  l       = (int)strlen(struc);
  packed  = (unsigned char *)vrna_alloc(((l + 4) / 5 + 1) * sizeof(unsigned char));

  j = i = pi = 0;
  while (i < l) {
    register int p;
    for (p = pi = 0; pi < 5; pi++) {
      p *= 3;
      switch (struc[i]) {
        case '(':
        case '\0':
          break;
        case ')':
          p++;
          break;
        case '.':
          p += 2;
          break;
        default:
          vrna_log_warning("vrna_db_pack: "
                               "illegal character %c at position %d in structure\n%s",
                               struc[i],
                               i + 1,
                               struc);
          return NULL;
      }
      if (i < l)
        i++;
    }
    packed[j++] = (unsigned char)(p + 1); /* never use 0, so we can use
                                           * strcmp()  etc. */
  }
  packed[j] = '\0';                       /* for str*() functions */
  return (char *)packed;
}


PUBLIC char *
vrna_db_unpack(const char *packed)
{
  /* 5:1 compression using base 3 encoding */
  int                 i, j, l;
  char                *struc;
  unsigned const char *pp;
  char                code[3] = {
    '(', ')', '.'
  };

  l     = (int)strlen(packed);
  pp    = (const unsigned char *)packed;
  struc = (char *)vrna_alloc((l * 5 + 1) * sizeof(char));   /* up to 4 byte extra */

  for (i = j = 0; i < l; i++) {
    register int p, c, k;

    p = (int)pp[i] - 1;
    for (k = 4; k >= 0; k--) {
      c             = p % 3;
      p             /= 3;
      struc[j + k]  = code[c];
    }
    j += 5;
  }
  struc[j--] = '\0';
  /* strip trailing ( */
  while ((j >= 0) &&
         (struc[j] == '('))
    struc[j--] = '\0';

  return struc;
}


PUBLIC char *
vrna_db_pk_remove(const char    *structure,
                  unsigned int  options)
{
  char  *s;
  short *pt_pk, *pt;

  s = NULL;

  if (structure) {
    pt_pk = vrna_ptable_from_string(structure, options & VRNA_BRACKETS_ANY);
    pt    = vrna_pt_pk_remove(pt_pk, options);
    s     = vrna_db_from_ptable(pt);

    free(pt_pk);
    free(pt);
  }

  return s;
}


PUBLIC char *
vrna_db_from_ptable(const short *ptable)
{
  char          *dotbracket           = NULL;
  const char    *bracket_open_avail   = "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  const char    *bracket_close_avail  = ")]}>abcdefghijklmnopqrstuvwxyz";
  short         *pt;
  unsigned int  n, i, stack_cnt, *stack, bracket_count, recheck;

  if (ptable) {
    n = (unsigned int)ptable[0];
    if (n > 0) {
      pt = (short *)vrna_alloc(sizeof(short) * (n + 1));
      pt = memcpy(pt, ptable, sizeof(short) * (n + 1));

      /* check for validity, i.e. both pairing positions must specify pairing partner */
      for (i = 1; i <= n; i++) {
        if (((unsigned int)pt[i] > i) &&
            ((unsigned int)pt[pt[i]] != i))
          return NULL;
      }

      /* prepare everything */
      dotbracket = (char *)vrna_alloc((n + 1) * sizeof(char));
      memset(dotbracket, '.', n);

      bracket_count = 0;
      recheck       = 1;
      stack         = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n + 1));
      stack_cnt     = 0;

      while (recheck) {
        recheck = 0;

        for (i = 1; i <= n; i++) {
          if (pt[i] > i) {
            /* check for clash */
            if ((stack_cnt > 0) &&
                (pt[i] > stack[stack_cnt - 1])) {
              recheck = 1;
              continue;
            } else {
              stack[stack_cnt++] = pt[i];
              dotbracket[i - 1] = bracket_open_avail[bracket_count];
              dotbracket[pt[i] - 1] = bracket_close_avail[bracket_count];
            }
          } else if (pt[i] != 0) {
            if (stack_cnt == 0)
              continue;

            if (i == stack[stack_cnt - 1]) {
              /* remove pair from pair table */
              pt[i] = pt[pt[i]] = 0;
              stack_cnt--;
            }
          }
        }

        stack_cnt = 0;
        bracket_count++;

        if (bracket_count >= 30) {
          vrna_log_warning("Not enough bracket types available in vrna_db_from_ptable()! Skipping remaining base pairs!");
          break;
        }
      }

      dotbracket[i - 1] = '\0';
      free(stack);
      free(pt);
    }
  }

  return dotbracket;
}


PUBLIC void
vrna_db_flatten(char          *string,
                unsigned int  options)
{
  vrna_db_flatten_to(string, "()", options);
}


PUBLIC void
vrna_db_flatten_to(char         *string,
                   const char   target[3],
                   unsigned int options)
{
  if (string) {
    if (options & VRNA_BRACKETS_RND)
      flatten_brackets(string, "()", target);

    if (options & VRNA_BRACKETS_ANG)
      flatten_brackets(string, "<>", target);

    if (options & VRNA_BRACKETS_CLY)
      flatten_brackets(string, "{}", target);

    if (options & VRNA_BRACKETS_SQR)
      flatten_brackets(string, "<>", target);

    if (options & VRNA_BRACKETS_ALPHA) {
      char pairs[3];

      for (int i = 65; i < 91; i++) {
        pairs[0]  = (char)i;
        pairs[1]  = (char)(i + 32);
        pairs[2]  = '\0';
        flatten_brackets(string, pairs, target);
      }
    }
  }
}


PUBLIC char *
vrna_db_from_WUSS(const char *wuss)
{
  char          *db, *tmp;
  short         *pt;
  unsigned int  i, p, q, pos, L, l[3], num_gq, n;

  db = NULL;

  if (wuss) {
    n       = strlen(wuss);
    num_gq  = 0;

    /*
     *  Note, in WUSS notation, matching pairs of (), <>, {}, [] are allowed but must not
     *  cross! Thus, we can simply flatten all brackets to ().
     */
    tmp = (char *)vrna_alloc(sizeof(char) * (n + 1));
    tmp = (char *)memcpy(tmp, wuss, n + 1);

    vrna_db_flatten(tmp, VRNA_BRACKETS_DEFAULT);

    /* now convert flattened structure string to pair-table (removes pseudo-knots) */
    pt = vrna_ptable_from_string(tmp, VRNA_BRACKETS_RND);

    /* convert back to dot-bracket (replaces all special characters for unpaired positions) */
    db = vrna_db_from_ptable(pt);

    /* check for G-Quadruplexes, annotated as G quartets */
    q = 0;
    while ((pos = vrna_gq_parse(wuss + q, &L, l)) > 0) {
      num_gq++;
      q += pos;

      if ((num_gq == 1) &&
          (4 * L + l[0] + l[1] + l[2] > pos)) {
        /* G-Quadruplex wraps around n,1 junction */
        p = n + pos + 1 - 4 * L - l[0] - l[1] - l[2] + 1;

        for (i = 0; i < L; i++) {
          unsigned int p1, p2, p3, p4;
          p1  = p + i;
          p2  = p1 + L + l[0];
          p3  = p2 + L + l[1];
          p4  = p3 + L + l[2];

          if (p4 > n)
            p4 = ((p4 - 1) % n) + 1;

          if (p3 > n)
            p3 = ((p3 - 1) % n) + 1;

          if (p2 > n)
            p2 = ((p2 - 1) % n) + 1;

          if (p1 > n)
            p1 = ((p1 - 1) % n) + 1;

          db[p1 - 1] = db[p2 - 1] = db[p3 - 1] = db[p4 - 1] = VRNA_GQUAD_DB_SYMBOL;

          if (i + 1 == L) /* last position of G-Quadruplex */
            db[p4 - 1] = VRNA_GQUAD_DB_SYMBOL_END;
        }
      } else {
        p = q - 4 * L - l[0] - l[1] - l[2] + 1;

        if (q > n)
          break;

        /* re-insert G-Quadruplex */
        for (i = 0; i < L; i++) {
          db[p + i - 1]                                 = VRNA_GQUAD_DB_SYMBOL;
          db[p + L + l[0] + i - 1]                      = VRNA_GQUAD_DB_SYMBOL;
          db[p + (2 * L) + l[0] + l[1] + i - 1]         = VRNA_GQUAD_DB_SYMBOL;
          db[p + (3 * L) + l[0] + l[1] + l[2] + i - 1]  = VRNA_GQUAD_DB_SYMBOL;
        }
      }

//      q++;
    }

    free(pt);
    free(tmp);
  }

  return db;
}


PUBLIC char
vrna_bpp_symbol(const float *x)
{
  /*  if( ((x[1]-x[2])*(x[1]-x[2]))<0.1&&x[0]<=0.677) return '|'; */
  if (x[0] > 0.667)
    return '.';

  if (x[1] > 0.667)
    return '(';

  if (x[2] > 0.667)
    return ')';

  if ((x[1] + x[2]) > x[0]) {
    if ((x[1] / (x[1] + x[2])) > 0.667)
      return '{';

    if ((x[2] / (x[1] + x[2])) > 0.667)
      return '}';
    else
      return '|';
  }

  if (x[0] > (x[1] + x[2]))
    return ',';

  return ':';
}


PUBLIC char *
vrna_db_from_probs(const FLT_OR_DBL *p,
                   unsigned int     length)
{
  int   i, j, *index;
  float P[3];    /* P[][0] unpaired, P[][1] upstream p, P[][2] downstream p */
  char  *s;

  s = NULL;

  if (p) {
    index = vrna_idx_row_wise(length);
    s     = (char *)vrna_alloc(sizeof(char) * (length + 1));

    for (j = 1; j <= length; j++) {
      P[0]  = 1.0;
      P[1]  = P[2] = 0.0;
      for (i = 1; i < j; i++) {
        P[2]  += (float)p[index[i] - j];  /* j is paired downstream */
        P[0]  -= (float)p[index[i] - j];  /* j is unpaired */
      }
      for (i = j + 1; i <= length; i++) {
        P[1]  += (float)p[index[j] - i];  /* j is paired upstream */
        P[0]  -= (float)p[index[j] - i];  /* j is unpaired */
      }
      s[j - 1] = vrna_bpp_symbol(P);
    }
    s[length] = '\0';
    free(index);
  }

  return s;
}


PUBLIC char *
vrna_pairing_tendency(vrna_fold_compound_t *fc)
{
  unsigned int  i, j, n;
  int           *idx;
  FLT_OR_DBL    *p, pp;
  float         P[3];    /* P[][0] unpaired, P[][1] upstream p, P[][2] downstream p */
  vrna_md_t     *md;
  char                      *s;
  vrna_smx_csr(FLT_OR_DBL)  *p_gq;

  s = NULL;

  if ((fc) &&
      (fc->exp_matrices) &&
      (fc->exp_matrices->probs)) {
    n   = fc->length;
    idx = fc->iindx;
    p   = fc->exp_matrices->probs;
    md  = &(fc->exp_params->model_details);
    s   = (char *)vrna_alloc(sizeof(char) * (n + 1));

    if ((md->circ) &&
        (md->gquad) &&
        (fc->exp_matrices->p_gq))
      p_gq = fc->exp_matrices->p_gq;
    else
      p_gq = NULL;

    for (j = 1; j <= n; j++) {
      P[0]  = 1.0;
      P[1]  = P[2] = 0.0;

      for (i = 1; i < j; i++) {
        P[2]  += (float)p[idx[i] - j];  /* j is paired downstream */
        P[0]  -= (float)p[idx[i] - j];  /* j is unpaired */
      }

      for (i = j + 1; i <= n; i++) {
        P[1]  += (float)p[idx[j] - i];  /* j is paired upstream */
        P[0]  -= (float)p[idx[j] - i];  /* j is unpaired */
      }

      if (p_gq) {
        /* do something about the gquads wrapping around the n,1 junction */
        for (i = 1; i < j; i++) {
#ifndef VRNA_DISABLE_C11_FEATURES
          pp    = vrna_smx_csr_get(p_gq, j, i, 0.);
#else
          pp    = vrna_smx_csr_FLT_OR_DBL_get(p_gq, j, i, 0.);
#endif
          P[1]  += (float)pp;  /* j is paired downstream */
          P[0]  -= (float)pp;  /* j is unpaired */
        }

        for (i = j + 1; i <= n; i++) {
#ifndef VRNA_DISABLE_C11_FEATURES
          pp    = vrna_smx_csr_get(p_gq, i, j, 0.);
#else
          pp    = vrna_smx_csr_FLT_OR_DBL_get(p_gq, i, j, 0.);
#endif
          P[2]  += (float)pp;  /* j is paired downstream */
          P[0]  -= (float)pp;  /* j is unpaired */
        }
      }

      s[j - 1] = vrna_bpp_symbol(P);
    }

    s[n] = '\0';
  }

  return s;
}


PUBLIC void
vrna_letter_structure(char            *structure,
                      vrna_bp_stack_t *bp,
                      unsigned int    length)
{
  int   n, k, x, y;
  char  alpha[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

  if (length > 0) {
    memset(structure, '.', length);
    structure[length] = '\0';

    for (n = 0, k = 1; k <= bp[0].i; k++) {
      y = bp[k].j;
      x = bp[k].i;
      if ((x - 1 > 0) && (y + 1 <= length)) {
        if (structure[x - 2] != ' ' && structure[y] == structure[x - 2]) {
          structure[x - 1]  = structure[x - 2];
          structure[y - 1]  = structure[x - 1];
          continue;
        }
      }

      if ((structure[x] != ' ') && (structure[y - 2] == structure[x])) {
        structure[x - 1]  = structure[x];
        structure[y - 1]  = structure[x - 1];
        continue;
      }

      n++;
      structure[x - 1]  = alpha[n - 1];
      structure[y - 1]  = alpha[n - 1];
    }
  }
}


PUBLIC char *
vrna_db_from_bps(vrna_bps_t   bp_stack,
                 unsigned int length)
{
  size_t        k;
  unsigned int  i, j, tmp;
  int           temp;
  char          *structure;
  vrna_bp_t     bp;

  structure = NULL;

  if (bp_stack) {
    structure = vrna_alloc(sizeof(char) * (length + 1));

    if (length > 0)
      memset(structure, '.', length);

    structure[length] = '\0';

    for (k = 0; k < vrna_bps_size(bp_stack); k++) {
      bp = vrna_bps_at(bp_stack, k);
      i = bp.i;
      j = bp.j;
      if (i > length)
        i -= length;

      if (j > length)
        j -= length;

      if (i > j) {
        temp  = i;
        i     = j;
        j     = temp;
      }

      if (i == j) {
        /* Gquad bonds are marked as bp[i].i == bp[i].j */
        if (bp.L > 0) {
          vrna_db_insert_gq(structure, i, bp.L, bp.l, length);
        } else {
          structure[i - 1] = '+';
        }
      } else {
        /* the following ones are regular base pairs */
        structure[i - 1]  = '(';
        structure[j - 1]  = ')';
      }
    }
  }

  return structure;
}


PUBLIC char *
vrna_db_from_bp_stack(vrna_bp_stack_t *bp,
                      unsigned int    length)
{
  int   k, i, j, temp;
  char  *structure;

  structure = NULL;

  if (bp) {
    structure = vrna_alloc(sizeof(char) * (length + 1));

    if (length > 0)
      memset(structure, '.', length);

    structure[length] = '\0';

    for (k = 1; k <= bp[0].i; k++) {
      i = bp[k].i;
      j = bp[k].j;
      if (i > length)
        i -= length;

      if (j > length)
        j -= length;

      if (i > j) {
        temp  = i;
        i     = j;
        j     = temp;
      }

      if (i == j) {
        /* Gquad bonds are marked as bp[i].i == bp[i].j */
        structure[i - 1] = '+';
      } else {
        /* the following ones are regular base pairs */
        structure[i - 1]  = '(';
        structure[j - 1]  = ')';
      }
    }
  }

  return structure;
}


PUBLIC char *
vrna_db_from_plist(vrna_ep_t    *pairs,
                   unsigned int n)
{
  vrna_ep_t *ptr;
  char      *structure = NULL;

  if ((n > 0) &&
      (pairs)) {
    structure = (char *)vrna_alloc(sizeof(char) * (n + 1));
    memset(structure, '.', n);
    structure[n] = '\0';

    for (ptr = pairs; (*ptr).i; ptr++) {
      if (((*ptr).i < n) && ((*ptr).j <= n)) {
        structure[(*ptr).i - 1] = '(';
        structure[(*ptr).j - 1] = ')';
      }
    }
  }

  return structure;
}


PUBLIC char *
vrna_db_to_element_string(const char *structure)
{
  char  *elements;
  int   n, i;
  short *pt;

  elements = NULL;

  if (structure) {
    n         = (int)strlen(structure);
    pt        = vrna_ptable(structure);
    elements  = (char *)vrna_alloc(sizeof(char) * (n + 1));

    for (i = 1; i <= n; i++) {
      if (!pt[i]) {
        /* mark nucleotides in exterior loop */
        elements[i - 1] = 'e';
      } else {
        assign_elements_pair(pt, i, pt[i], elements);
        i = pt[i];
      }
    }

    elements[n] = '\0';
    free(pt);
  }

  return elements;
}


PUBLIC unsigned int
vrna_gq_parse(const char    *db_string,
              unsigned int  *L,
              unsigned int  l[3])
{
  unsigned int i, n, end, start, stop, pre, tetrad, stacks, LL[4];

  end = 0;

  if ((db_string == NULL) ||
      (L == NULL))
    return 0;

  /*
   *  scan along the dot-bracket string to identify the first character that
   *  indicates a part of a G-Quadruplex
   */

  n       = strlen(db_string);
  LL[0]   = LL[1] = LL[2] = LL[3] = 0;
  stop    = 0;
  *L      = 0;
  l[0]    = l[1] = l[2] = 0;
  tetrad  = 0;

  for (i = 0;
       (db_string[i]) &&
       (db_string[i] != VRNA_GQUAD_DB_SYMBOL) &&
       (db_string[i] != VRNA_GQUAD_DB_SYMBOL_END);
       i++);

  if (db_string[i]) {
    pre = i;

    /* we've encountered some piece that suspicially looks like a G-Quadruplex */
    if (db_string[i] == VRNA_GQUAD_DB_SYMBOL_END) {
      stop  = 1;
      LL[0] = 1;
    }

    while (stop == 0) {
      start = i;
      while (db_string[++i] == VRNA_GQUAD_DB_SYMBOL)
        /* stop consuming gquad symbols if we already know the number of stacks */
        if ((tetrad > 1) &&
            (i - start == LL[tetrad - 1]))
          break;

      end         = i;
      LL[tetrad]  = end - start;

      if (db_string[i] == VRNA_GQUAD_DB_SYMBOL_END)
        LL[tetrad]++;

      if ((tetrad > 1) &&
          (LL[tetrad] != LL[tetrad - 1])) {
        vrna_log_debug("unequal stack sizes (%u vs. %u) in G-quadruplex",
                       LL[tetrad - 1],
                       LL[tetrad]);
        return 0;
      }

      /* check whether we stopped due to linker or end symbol */
      if (db_string[i] == VRNA_GQUAD_DB_SYMBOL_END) {
        stop = i + 1;
        break;
      }

      if (tetrad == 3)
        break; /* no more linkers */

      /*
       * next character must be linker!
       * count number of linker nucleotides
       */
      while (db_string[i++] == '.');
      l[tetrad] = i - end - 1;

      if (l[tetrad] == 0) {
        vrna_log_debug("zero-length linker in G-Quadruplex");
        return 0;
      } else {
        i--;
      }

      if (db_string[i] == '\0') {
        return 0;
      }

      tetrad++;
    }

    if (stop > 0) {
      if (tetrad < 3) {
        /* move currently parsed lengths to correct positions */
        unsigned int s;
        for (s = 0; s <= tetrad; s++)
          LL[3 - s] = LL[tetrad - s];

        for (; s <= 3; s++)
          LL[3 - s] = 0;

        for (s = 0; s < tetrad; s++)
          l[2 - s] = l[tetrad - s - 1];

        l[2 - s++] = pre;

        for (; s < 3; s++)
          l[2 - s] = 0;
      }

      /* only continue parsing if we did not find a complete G-Quadruplex yet */
      if (LL[0] != LL[3]) {
        i = n - 1;
        /* check whether the end of the structure is a continuation of a g-run or a linker */
        if (db_string[i] == '.') {
          while (db_string[--i] == '.');
          l[2 - tetrad] += n - i - 1;
          tetrad++;
        } else if (pre > 0) {
          tetrad++;
        }

        while (tetrad < 4) {
          start = i;

          while (db_string[--i] == VRNA_GQUAD_DB_SYMBOL) {
            /* stop consuming gquad symbols if we already know the number of stacks */
            if ((tetrad == 3) &&
                (start - i + LL[3 - tetrad] >= LL[3]))
              break;
          }

          end             = i;
          LL[3 - tetrad]  += start - end;

          if ((tetrad > 0) &&
              (LL[3 - tetrad] != LL[4 - tetrad])) {
            vrna_log_debug("unequal stack sizes (%u vs. %u) in G-quadruplex",
                           LL[3 - tetrad],
                           LL[4 - tetrad]);
            return 0;
          }

          if (tetrad == 3)
            break; /* no more linkers */

          /*
           * next character must be linker!
           * count number of linker nucleotides
           */
          while (db_string[i--] == '.');
          l[2 - tetrad] = end - i - 1;

          if (l[2 - tetrad] == 0) {
            vrna_log_debug("zero-length linker in G-Quadruplex");
            return 0;
          } else {
            i++;
          }

          if (i <= stop + 1) {
            return 0;
          }

          tetrad++;
        }
      }

      end = stop;
    }
  }

  /* last sanity check */
  if ((LL[0] == LL[1]) &&
      (LL[1] == LL[2]) &&
      (LL[2] == LL[3])) {
    *L = LL[0];
  } else {
    vrna_log_debug("Unequal G-quadruplex stack tetrad lengths (%u, %u, %u, %d)",
                   LL[0],
                   LL[1],
                   LL[2],
                   LL[3]);
    return 0;
  }

  return end;
}


PUBLIC void
vrna_db_insert_gq(char          *db,
                  unsigned int  i,
                  unsigned int  L,
                  unsigned int l[3],
                  unsigned int  n)
{
  if (db) {
    if (n == 0)
      n = strlen(db);

    if (4 * L + l[0] + l[1] + l[2] <= n) {
      unsigned int j, ll;
      for (ll = 0; ll < L; ll++) {
        j     = (i + ll - 1) % (n);
        db[j] = VRNA_GQUAD_DB_SYMBOL;
      }

      i += L + l[0];
      for (ll = 0; ll < L; ll++) {
        j     = (i + ll - 1) % (n);
        db[j] = VRNA_GQUAD_DB_SYMBOL;
      }

      i += L + l[1];
      for (ll = 0; ll < L; ll++) {
        j     = (i + ll - 1) % (n);
        db[j] = VRNA_GQUAD_DB_SYMBOL;
      }

      i += L + l[2];
      for (ll = 0; ll < L - 1; ll++) {
        j     = (i + ll - 1) % (n);
        db[j] = VRNA_GQUAD_DB_SYMBOL;
      }
      j     = (i + ll - 1) % (n);
      db[j] = VRNA_GQUAD_DB_SYMBOL_END;
    }
  }
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE INLINE void
flatten_brackets(char       *string,
                 const char pair[3],
                 const char target[3])
{
  unsigned int i;

  for (i = 0; string[i] != '\0'; i++) {
    if (string[i] == pair[0])
      string[i] = target[0];
    else if (string[i] == pair[1])
      string[i] = target[1];
  }
}


PRIVATE void
assign_elements_pair(short  *pt,
                     int    i,
                     int    j,
                     char   *elements)
{
  int p, k, num_pairs;

  num_pairs = 0;
  /* first, determine the number of pairs (i,j) is enclosing */
  for (k = i + 1; k < j; k++) {
    if (k < pt[k]) {
      num_pairs++;
      k = pt[k];
    }
  }

  switch (num_pairs) {
    case 0:   /* hairpin loop */
      elements[i - 1] = elements[j - 1] = 'H';
      for (k = i + 1; k < j; k++)
        elements[k - 1] = 'h';
      break;

    case 1:   /* internal loop */
      elements[i - 1] = elements[j - 1] = 'I';
      p               = 0;
      for (k = i + 1; k < j; k++) {
        if (!pt[k]) {
          elements[k - 1] = 'i';
        } else {
          p = k;
          k = pt[k];
        }
      }

      if (p)
        assign_elements_pair(pt, p, pt[p], elements);

      break;

    default:  /* multibranch loop */
      elements[i - 1] = elements[j - 1] = 'M';
      for (k = i + 1; k < j; k++) {
        if (!pt[k]) {
          elements[k - 1] = 'm';
        } else {
          assign_elements_pair(pt, k, pt[k], elements);
          k = pt[k];
        }
      }
      break;
  }
}


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
 * ###########################################
 * # deprecated functions below              #
 *###########################################
 */


PUBLIC char *
pack_structure(const char *struc)
{
  return vrna_db_pack(struc);
}


PUBLIC char *
unpack_structure(const char *packed)
{
  return vrna_db_unpack(packed);
}


PUBLIC void
parenthesis_structure(char            *structure,
                      vrna_bp_stack_t *bp,
                      int             length)
{
  char *s = vrna_db_from_bp_stack(bp, length);

  strncpy(structure, s, length + 1);
  free(s);
}


PUBLIC void
letter_structure(char             *structure,
                 vrna_bp_stack_t  *bp,
                 int              length)
{
  vrna_letter_structure(structure, bp, length);
}


PUBLIC void
parenthesis_zuker(char            *structure,
                  vrna_bp_stack_t *bp,
                  int             length)
{
  char *s = vrna_db_from_bp_stack(bp, length);

  strncpy(structure, s, length + 1);
  free(s);
}


PUBLIC char
bppm_symbol(const float *x)
{
  return vrna_bpp_symbol(x);
}


PUBLIC void
bppm_to_structure(char          *structure,
                  FLT_OR_DBL    *p,
                  unsigned int  length)
{
  char *s = vrna_db_from_probs((const FLT_OR_DBL *)p, length);

  memcpy(structure, s, length);
  structure[length] = '\0';
  free(s);
}


PUBLIC int
parse_gquad(const char  *struc,
            int         *L,
            int         l[3])
{
  unsigned int ret, LL, ll[3];

  ret = vrna_gq_parse(struc, &LL, ll);

  if (ret) {
    *L    = (int)LL;
    l[0]  = (int)ll[0];
    l[1]  = (int)ll[1];
    l[2]  = (int)ll[2];
  }

  return ret;
}


#endif
