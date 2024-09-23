/*
 * gquad.c
 *
 * Ronny Lorenz 2012
 *
 * ViennaRNA Package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/params/constants.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/alignments.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/loops/gquad.h"
#include "ViennaRNA/eval/gquad.h"

#include "ViennaRNA/loops/gquad_intern.h"


#ifndef INLINE
#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif
#endif


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

/**
 *  MFE callback for process_gquad_enumeration()
 */
PRIVATE
void
gquad_mfe(unsigned int   i,
          unsigned int   L,
          unsigned int   *l,
          void  *data,
          void  *P,
          void  *NA,
          void  *NA2);


/**
 * Partition function callback for process_gquad_enumeration()
 */
PRIVATE
void
gquad_pf(unsigned int  i,
         unsigned int  L,
         unsigned int  *l,
         void *data,
         void *P,
         void *NA,
         void *NA2);


PRIVATE void
gquad_pf_ali(unsigned int  i,
             unsigned int  L,
             unsigned int  *l,
             void *data,
             void *helper,
             void *NA,
             void *NA2);


/**
 * MFE (alifold) callback for process_gquad_enumeration()
 */
PRIVATE
void
gquad_mfe_ali(unsigned int   i,
              unsigned int   L,
              unsigned int   *l,
              void  *data,
              void  *helper,
              void  *NA,
              void  *NA2);


/**
 * MFE (alifold) callback for process_gquad_enumeration()
 * with seperation of free energy and penalty contribution
 */
PRIVATE
void
gquad_mfe_ali_en(unsigned int  i,
                 unsigned int  L,
                 unsigned int  *l,
                 void *data,
                 void *helper,
                 void *NA,
                 void *NA2);


PRIVATE
void
gquad_count(unsigned int   i,
            unsigned int   L,
            unsigned int   *l,
            void  *data,
            void  *NA,
            void  *NA2,
            void  *NA3);


PRIVATE
void
gquad_count_layers(unsigned int  i,
                   unsigned int  L,
                   unsigned int  *l,
                   void *data,
                   void *NA,
                   void *NA2,
                   void *NA3);


/* other useful static functions */

PRIVATE void
gquad_mfe_ali_pos(unsigned int   i,
                  unsigned int   L,
                  unsigned int   *l,
                  void  *data,
                  void  *P,
                  void  *Lmfe,
                  void  *lmfe);


/*
 #########################################
 # BEGIN OF PUBLIC FUNCTION DEFINITIONS  #
 #      (all available in RNAlib)        #
 #########################################
 */

/* 4. Parsing */

PUBLIC plist *
get_plist_gquad_from_db(const char  *structure,
                        float       pr)
{
  unsigned int  x, L, l[3], n, size, ge, ee, actual_size, gb;
  plist         *pl;

  actual_size = 0;
  ge          = 0;
  n           = 2;
  size        = strlen(structure);
  pl          = (plist *)vrna_alloc(n * size * sizeof(plist));

  while ((ee = vrna_gq_parse(structure + ge, &L, l)) > 0) {
    ge  += ee;

    if (4 * L + l[0] + l[1] + l[2] > ee) {
      gb = size + ge - L * 4 - l[0] - l[1] - l[2] + 1;
    } else {
      gb = ge - L * 4 - l[0] - l[1] - l[2] + 1;
    }

    /* add pseudo-base pair enclosing gquad */
    if (actual_size >= n * size - 5) {
      n   *= 2;
      pl  = (plist *)vrna_realloc(pl, n * size * sizeof(plist));
    }

    pl[actual_size].i       = gb;
    pl[actual_size].j       = ge;
    pl[actual_size].p       = pr;
    pl[actual_size++].type  = VRNA_PLIST_TYPE_GQUAD;

    for (x = 0; x < L; x++) {
      if (actual_size >= n * size - 5) {
        n   *= 2;
        pl  = (plist *)vrna_realloc(pl, n * size * sizeof(plist));
      }

      pl[actual_size].i       = (gb + x - 1) % (size) + 1;
      pl[actual_size].j       = (ge + x - L + 1 - 1) % (size) + 1;
      pl[actual_size].p       = pr;
      pl[actual_size++].type  = VRNA_PLIST_TYPE_TRIPLE;

      pl[actual_size].i       = (gb + x - 1) % (size) + 1;
      pl[actual_size].j       = (gb + x + l[0] + L - 1) % (size) + 1;
      pl[actual_size].p       = pr;
      pl[actual_size++].type  = VRNA_PLIST_TYPE_TRIPLE;

      pl[actual_size].i       = (gb + x + l[0] + L - 1) % (size) + 1;
      pl[actual_size].j       = (ge + x - 2 * L - l[2] + 1 - 1) % (size) + 1;
      pl[actual_size].p       = pr;
      pl[actual_size++].type  = VRNA_PLIST_TYPE_TRIPLE;

      pl[actual_size].i       = (ge + x - 2 * L - l[2] + 1 - 1) % (size) + 1;
      pl[actual_size].j       = (ge + x - L + 1 - 1) % (size) + 1;
      pl[actual_size].p       = pr;
      pl[actual_size++].type  = VRNA_PLIST_TYPE_TRIPLE;
    }
  }

  pl[actual_size].i   = pl[actual_size].j = 0;
  pl[actual_size++].p = 0;
  pl                  = (plist *)vrna_realloc(pl, actual_size * sizeof(plist));
  return pl;
}


PUBLIC int
get_gquad_count(short *S,
                int   i,
                int   j)
{
  unsigned int p, q, *gg;
  int           counter;

  gg = get_g_islands_sub(S, (unsigned int)i, (unsigned int)j);
  counter = 0;

  FOR_EACH_GQUAD(p, q, (unsigned int)i, (unsigned int)j)
  process_gquad_enumeration(gg, p, q,
                            &gquad_count,
                            (void *)(&counter),
                            NULL,
                            NULL,
                            NULL);

  gg += (unsigned int)i - 1;
  free(gg);
  return counter;
}


PUBLIC int
get_gquad_layer_count(short *S,
                      int   i,
                      int   j)
{
  unsigned int  p, q, *gg;
  int           counter;

  gg = get_g_islands_sub(S, (unsigned int)i, (unsigned int)j);
  counter = 0;

  FOR_EACH_GQUAD(p, q, (unsigned int)i, (unsigned int)j)
  process_gquad_enumeration(gg, p, q,
                            &gquad_count_layers,
                            (void *)(&counter),
                            NULL,
                            NULL,
                            NULL);

  gg += (unsigned int)i - 1;
  free(gg);
  return counter;
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
 #########################################
 # BEGIN OF PRIVATE FUNCTION DEFINITIONS #
 #          (internal use only)          #
 #########################################
 */
PRIVATE void
gquad_mfe(unsigned int   i,
          unsigned int   L,
          unsigned int   *l,
          void  *data,
          void  *P,
          void  *NA,
          void  *NA2)
{
  int cc = ((vrna_param_t *)P)->gquad[L][l[0] + l[1] + l[2]];

  if (cc < *((int *)data))
    *((int *)data) = cc;
}


PRIVATE
void
gquad_count(unsigned int   i,
            unsigned int   L,
            unsigned int   *l,
            void  *data,
            void  *NA,
            void  *NA2,
            void  *NA3)
{
  *((int *)data) += 1;
}


PRIVATE
void
gquad_count_layers(unsigned int  i,
                   unsigned int  L,
                   unsigned int  *l,
                   void *data,
                   void *NA,
                   void *NA2,
                   void *NA3)
{
  *((unsigned int *)data) += L;
}


PRIVATE void
gquad_pf(unsigned int  i,
         unsigned int  L,
         unsigned int  *l,
         void *data,
         void *pf,
         void *NA,
         void *NA2)
{
  *((FLT_OR_DBL *)data) += ((vrna_exp_param_t *)pf)->expgquad[L][l[0] + l[1] + l[2]];
}


PRIVATE void
gquad_pf_ali(unsigned int  i,
             unsigned int  L,
             unsigned int  *l,
             void *data,
             void *helper,
             void *NA,
             void *NA2)
{
  const short             **S;
  const unsigned int      **a2s;
  unsigned int            n, s, n_seq;
  FLT_OR_DBL              penalty;
  vrna_exp_param_t        *pf;
  struct gquad_ali_helper *gq_help;

  gq_help = (struct gquad_ali_helper *)helper;
  S       = gq_help->S;
  a2s     = gq_help->a2s;
  n       = gq_help->length;
  n_seq   = gq_help->n_seq;
  pf      = gq_help->pf;

  *((FLT_OR_DBL *)data) = vrna_exp_E_consensus_gquad(L,
                                                     l,
                                                     pf,
                                                     i,
                                                     n,
                                                     n_seq,
                                                     S,
                                                     a2s);
}


PRIVATE void
gquad_mfe_ali(unsigned int   i,
              unsigned int   L,
              unsigned int   *l,
              void  *data,
              void  *helper,
              void  *NA,
              void  *NA2)
{
  int en[2], cc;

  en[0] = en[1] = INF;

  CHECK_GQUAD(L, l, return );

  gquad_mfe_ali_en(i, L, l, (void *)(&(en[0])), helper, NULL, NULL);
  if (en[1] != INF) {
    cc = en[0] + en[1];
    if (cc < *((int *)data))
      *((int *)data) = cc;
  }
}


PRIVATE void
gquad_mfe_ali_en(unsigned int  i,
                 unsigned int  L,
                 unsigned int  *l,
                 void *data,
                 void *helper,
                 void *NA,
                 void *NA2)
{
  const short             **S;
  const unsigned int      **a2s;
  unsigned int            s, n_seq, n;
  int                     en[2], cc, dd, u1, u2, u3;
  vrna_param_t            *P;
  struct gquad_ali_helper *gq_help;

  gq_help = (struct gquad_ali_helper *)helper;
  S       = gq_help->S;
  a2s     = gq_help->a2s;
  n       = gq_help->length;
  n_seq   = gq_help->n_seq;
  P       = gq_help->P;

  vrna_E_consensus_gquad(L,
                         l,
                         i,
                         n,
                         n_seq,
                         S,
                         a2s,
                         P,
                         en);

  if (en[1] != INF) {
    cc  = en[0] + en[1];
    dd  = ((int *)data)[0] + ((int *)data)[1];
    if (cc < dd) {
      ((int *)data)[0]  = en[0];
      ((int *)data)[1]  = en[1];
    }
  }
}


/*
 *###########################################
 *# deprecated functions below              #
 *###########################################
 */


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

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
