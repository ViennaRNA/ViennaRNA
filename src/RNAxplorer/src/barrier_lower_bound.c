/*
 *
 * Compute a lower bound for a refolding energy barrier using 2D representation
 * of the landscape, and application of the Dijkstra algorithm on the resulting
 * graph
 *
 * (c) Ronny Lorenz
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/utils/structures.h>
#include <ViennaRNA/2Dfold.h>

#include "barrier_lower_bound.h"

typedef struct nb_t {
  float B;
  int   n_lu;
  int   n_ld;
  int   n_ru;
  int   n_rd;
  int   k;
  int   l;
  int   prev_idx;

  float en;
  char  *s;
} nb_t;

void printBarrier(float B,
                  float E,
                  char  *s);


void
printBarrier(float  B,
             float  E,
             char   *s)
{
  if (s)
    printf("Estimated energy barrier is %6.2f\nSaddle: %s (%6.2f)\n", B, s, E);
  else
    printf("Estimated energy barrier is %6.2f\n", B);
}


/**
 * Compute a shortest path via Dijkstra algorithm trough the 2D-Landscape.
 * @param seq - the RNA sequence
 * @param md - model details
 * @param s1 - dot-bracket structure
 * @param s2 - dot-bracket structure
 * @param maximum_distance1 - initial value for space allocation. Will be computed within this method.
 * @param maximum_distance2 - initial value for space allocation. Will be computed within this method.
 * @return the energybarrier of the path.
 */
float
barrier_estimate_2D(char      *seq,
                    vrna_md_t *md,
                    char      *s1,
                    char      *s2,
                    int       maximum_distance1,
                    int       maximum_distance2)
{
  short *pt1, *pt2;

  pt1 = vrna_ptable(s1);
  pt2 = vrna_ptable(s2);
  int   i, n = pt1[0];
  /* compute symmetrical difference between both structures */
  int   a = 0;
  int   b = 0;
  int   c = 0;
  for (i = 1; i <= n; i++) {
    if (pt1[i] == pt2[i]) {
      if (pt1[i] != 0)
        c++;
    } else {
      if (i < pt1[i])
        a++;

      if (i < pt2[i])
        b++;
    }
  }

  free(pt1);
  free(pt2);

  printf("%d:%d ... %d:%d\n", a + c, a + maximum_distance1, b + c, b + maximum_distance2);
  int                   maxD1 = (a + c < a + maximum_distance1) ? a + maximum_distance1 : a + c;
  int                   maxD2 = (b + c < b + maximum_distance2) ? b + maximum_distance2 : a + c;

  vrna_fold_compound_t  *vc     = vrna_fold_compound_TwoD(seq, s1, s2, md, VRNA_OPTION_MFE);
  vrna_sol_TwoD_t       *mfe_s  = vrna_mfe_TwoD(vc, maxD1, maxD2);


  /* make a lucky guess for the real max distances */

  float                 mfe_s1, mfe_s2;
  int                   map_s1, map_s2;

  int                   max_k   = 0;
  int                   min_k   = INF;
  int                   *max_l  = (int *)vrna_alloc(sizeof(int) * (maxD1 + 1));
  int                   *min_l  = (int *)vrna_alloc(sizeof(int) * (maxD1 + 1));
  for (i = 0; i < (maxD1 + 1); i++) {
    max_l[i]  = 0;
    min_l[i]  = INF;
  }
  int                   **mapping = (int **)vrna_alloc(sizeof(int *) * (maxD1 + 1));
  for (i = 0; i < (maxD1 + 1); i++)
    mapping[i] = (int *)vrna_alloc(sizeof(int) * (maxD2 + 1));

  /* get some statistics of the 2D fold output */
  int                   number_of_states = 0;
  for (i = number_of_states = 0; mfe_s[i].k != INF; i++) {
    int k = mfe_s[i].k;
    int l = mfe_s[i].l;

    if (k == -1)
      continue;

    max_k     = MAX2(max_k, k);
    min_k     = MIN2(min_k, k);
    max_l[k]  = MAX2(max_l[k], l);
    min_l[k]  = MIN2(min_l[k], l);

    if (k == 0) {
      map_s1  = i;
      mfe_s1  = mfe_s[i].en;
    }

    if (l == 0) {
      map_s2  = i;
      mfe_s2  = mfe_s[i].en;
    }

    mapping[k][l] = i;
    number_of_states++;
  }

  /* begin actual graph node construction */
  nb_t  *nodes      = (nb_t *)vrna_alloc(sizeof(nb_t) * (number_of_states + 1));
  int   map_source  = (mfe_s1 < mfe_s2) ? map_s1 : map_s2;
  int   map_target  = (mfe_s2 <= mfe_s1) ? map_s1 : map_s2;
  float min12       = MIN2(mfe_s1, mfe_s2);

  for (i = 0; mfe_s[i].k != INF; i++) {
    int k = mfe_s[i].k;
    int l = mfe_s[i].l;

    if (k == -1)
      continue;

    nodes[i].en       = mfe_s[i].en;
    nodes[i].s        = mfe_s[i].s;
    nodes[i].k        = k;
    nodes[i].l        = l;
    nodes[i].n_lu     = nodes[i].n_ld = nodes[i].n_ru = nodes[i].n_rd = -1;
    nodes[i].B        = (float)INF / 100.;
    nodes[i].prev_idx = -1;

    /* set the four neighbor indices, if they are within range */
    if (k - 1 >= min_k) {
      if ((l - 1 >= min_l[k - 1]) && (l - 1 <= max_l[k - 1]))
        nodes[i].n_ld = mapping[k - 1][l - 1];

      if ((l + 1 <= max_l[k - 1]) && (l + 1 >= min_l[k - 1]))
        nodes[i].n_lu = mapping[k - 1][l + 1];
    }

    if (k + 1 <= max_k) {
      if ((l - 1 >= min_l[k + 1]) && (l - 1 <= max_l[k + 1]))
        nodes[i].n_rd = mapping[k + 1][l - 1];

      if ((l + 1 <= max_l[k + 1]) && (l + 1 >= min_l[k + 1]))
        nodes[i].n_ru = mapping[k + 1][l + 1];
    }
  }

  free(mfe_s);
  free(max_l);
  free(min_l);

  nb_t  *vertices     = nodes;
  int   vertex_count  = number_of_states;
  int   shift         = 0;

  vertices[map_source].B = 0.;

  /* search for optimal path using Dijkstra algorithm */
  do {
    /* put node with smallest barrier to front */
    int min_idx = 0;
    int min_b   = vertices[0].B;
    for (i = 1; i < vertex_count; i++) {
      if (vertices[i].B < min_b) {
        min_idx = i;
        min_b   = vertices[i].B;
      }
    }
    if (min_idx != 0) {
      /* swap entries if vertex with smallest barrier was not in front */
      nb_t tmp;
      memcpy(&tmp, vertices, sizeof(nb_t));
      memcpy(vertices, vertices + min_idx, sizeof(nb_t));
      memcpy(vertices + min_idx, &tmp, sizeof(nb_t));

      /* update all vertices that had either of the swapped nodes as neighbors */
      for (i = 0; i < vertex_count; i++) {
        if (vertices[i].n_ld == 0)
          vertices[i].n_ld = min_idx;
        else if (vertices[i].n_ld == min_idx)
          vertices[i].n_ld = 0;

        if (vertices[i].n_lu == 0)
          vertices[i].n_lu = min_idx;
        else if (vertices[i].n_lu == min_idx)
          vertices[i].n_lu = 0;

        if (vertices[i].n_rd == 0)
          vertices[i].n_rd = min_idx;
        else if (vertices[i].n_rd == min_idx)
          vertices[i].n_rd = 0;

        if (vertices[i].n_ru == 0)
          vertices[i].n_ru = min_idx;
        else if (vertices[i].n_ru == min_idx)
          vertices[i].n_ru = 0;
      }

      /* did we move the source or target vertex? */
      if (min_idx == map_target) {
        map_target = 0;
        break;
      }

      if (min_idx == map_source)
        map_source = 0;

      if (map_target == 0)
        map_target = min_idx;
    }

    /* now lets update the neighboring vertices */
    if (vertices[0].n_ld > 0) {
      float tmp_b = MAX2(vertices[0].B, vertices[vertices[0].n_ld].en - min12);
      if (tmp_b < vertices[vertices[0].n_ld].B) {
        vertices[vertices[0].n_ld].B        = tmp_b;
        vertices[vertices[0].n_ld].prev_idx = shift;
      }
    }

    if (vertices[0].n_lu > 0) {
      float tmp_b = MAX2(vertices[0].B, vertices[vertices[0].n_lu].en - min12);
      if (tmp_b < vertices[vertices[0].n_lu].B) {
        vertices[vertices[0].n_lu].B        = tmp_b;
        vertices[vertices[0].n_lu].prev_idx = shift;
      }
    }

    if (vertices[0].n_rd > 0) {
      float tmp_b = MAX2(vertices[0].B, vertices[vertices[0].n_rd].en - min12);
      if (tmp_b < vertices[vertices[0].n_rd].B) {
        vertices[vertices[0].n_rd].B        = tmp_b;
        vertices[vertices[0].n_rd].prev_idx = shift;
      }
    }

    if (vertices[0].n_ru > 0) {
      float tmp_b = MAX2(vertices[0].B, vertices[vertices[0].n_ru].en - min12);
      if (tmp_b < vertices[vertices[0].n_ru].B) {
        vertices[vertices[0].n_ru].B        = tmp_b;
        vertices[vertices[0].n_ru].prev_idx = shift;
      }
    }

    vertices++;
    vertex_count--;
    map_target--;
    shift++;
    /* decrease vertex index of neighboring nodes for remaining set of vertices */
    for (i = 0; i < vertex_count; i++) {
      vertices[i].n_ld--;
      vertices[i].n_lu--;
      vertices[i].n_rd--;
      vertices[i].n_ru--;
    }
  } while (vertex_count > 0);

  /* now backtrack the optimal path */
  int   curr      = map_target + shift;
  float B         = nodes[curr].B;
  float E_saddle  = nodes[curr].en;
  char  *s_saddle = nodes[curr].s;

  while (curr != 0) {
    if (nodes[curr].en > E_saddle) {
      E_saddle  = nodes[curr].en;
      s_saddle  = nodes[curr].s;
    }

    printf("%d %d -> %6.2f\n", nodes[curr].k, nodes[curr].l, nodes[curr].B);
    curr = nodes[curr].prev_idx;
  }
  printf("%d %d -> %6.2f\n", nodes[0].k, nodes[0].l, nodes[0].B);

  printBarrier(B, E_saddle, s_saddle);

  vrna_fold_compound_free(vc);
  for (i = 0; i < number_of_states; i++)
    free(nodes[i].s);
  free(nodes);
  for (i = 0; i < (maxD1 + 1); i++)
    free(mapping[i]);
  free(mapping);

  return B;
}
