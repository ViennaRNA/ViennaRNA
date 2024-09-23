/**
 * This file is a container for all plotting layout algorithms
 *
 *  c Ronny Lorenz
 *    The ViennaRNA Package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/structures/dotbracket.h"

#include "ViennaRNA/plotting/layouts.h"

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

PUBLIC int rna_plot_type = VRNA_PLOT_TYPE_DEFAULT;

#endif


#ifndef PI
#define  PI       3.141592654
#endif
#define  PIHALF       PI / 2.


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE struct vrna_plot_layout_s *
rna_layout(const char   *structure,
           unsigned int plot_type,
           void         *options);


PRIVATE int
coords_simple(const short *pt,
              float       **x,
              float       **y);


PRIVATE void
loop(const short  *pair_table,
     int          i,
     int          j,
     float        *angle,
     int          *stack_size,
     int          *loop_size,
     int          *stk,
     int          *lp);


PRIVATE int
coords_circular(const short *pt,
                float       **x,
                float       **y);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC struct vrna_plot_layout_s *
vrna_plot_layout(const char   *structure,
                 unsigned int plot_type)
{
  if (structure)
    return rna_layout(structure, plot_type, NULL);

  return NULL;
}


PUBLIC struct vrna_plot_layout_s *
vrna_plot_layout_simple(const char *structure)
{
  return vrna_plot_layout(structure, VRNA_PLOT_TYPE_SIMPLE);
}


#ifdef VRNA_WITH_NAVIEW_LAYOUT
PUBLIC struct vrna_plot_layout_s *
vrna_plot_layout_naview(const char *structure)
{
  return vrna_plot_layout(structure, VRNA_PLOT_TYPE_NAVIEW);
}
#endif

PUBLIC struct vrna_plot_layout_s *
vrna_plot_layout_circular(const char *structure)
{
  return vrna_plot_layout(structure, VRNA_PLOT_TYPE_CIRCULAR);
}


PUBLIC struct vrna_plot_layout_s *
vrna_plot_layout_turtle(const char *structure)
{
  return vrna_plot_layout(structure, VRNA_PLOT_TYPE_TURTLE);
}


PUBLIC struct vrna_plot_layout_s *
vrna_plot_layout_puzzler(const char                   *structure,
                         vrna_plot_options_puzzler_t  *options)
{
  if (structure)
    return rna_layout(structure, VRNA_PLOT_TYPE_PUZZLER, (void *)options);

  return NULL;
}


PUBLIC void
vrna_plot_layout_free(struct vrna_plot_layout_s *layout)
{
  if (layout) {
    free(layout->x);
    free(layout->y);
    free(layout->arcs);
    free(layout);
  }
}


PUBLIC int
vrna_plot_coords(const char *structure,
                 float      **x,
                 float      **y,
                 int        plot_type)
{
  if (structure) {
    int   ret = 0;
    short *pt = vrna_ptable(structure);

    ret = vrna_plot_coords_pt(pt, x, y, plot_type);

    free(pt);

    return ret;
  }

  if (x)
    *x = NULL;

  if (y)
    *y = NULL;

  return 0;
}


PUBLIC int
vrna_plot_coords_pt(const short *pt,
                    float       **x,
                    float       **y,
                    int         plot_type)
{
  if ((pt) && (x) && (y)) {
    switch (plot_type) {
      case VRNA_PLOT_TYPE_SIMPLE:
        return coords_simple(pt, x, y);

#ifdef VRNA_WITH_NAVIEW_LAYOUT
      case VRNA_PLOT_TYPE_NAVIEW:
        return vrna_plot_coords_naview_pt(pt, x, y);
#endif

      case VRNA_PLOT_TYPE_CIRCULAR:
        return coords_circular(pt, x, y);

      case VRNA_PLOT_TYPE_TURTLE:
        return vrna_plot_coords_turtle_pt(pt, x, y, NULL);

      case VRNA_PLOT_TYPE_PUZZLER:
        return vrna_plot_coords_puzzler_pt(pt, x, y, NULL, NULL);

#ifdef VRNA_WITH_NAVIEW_LAYOUT
      default:
        return vrna_plot_coords_naview_pt(pt, x, y);
#else
      case VRNA_PLOT_TYPE_NAVIEW:
        vrna_log_warning("Naview layout algorithm not available in this version of ViennaRNA Package");
        break;

      default:
        return vrna_plot_coords_puzzler_pt(pt, x, y, NULL, NULL);
#endif

    }
  }

  if (x)
    *x = NULL;

  if (y)
    *y = NULL;

  return 0;
}


PUBLIC int
vrna_plot_coords_simple(const char  *structure,
                        float       **x,
                        float       **y)
{
  return vrna_plot_coords(structure, x, y, VRNA_PLOT_TYPE_SIMPLE);
}


PUBLIC int
vrna_plot_coords_simple_pt(const short  *pt,
                           float        **x,
                           float        **y)
{
  return vrna_plot_coords_pt(pt, x, y, VRNA_PLOT_TYPE_SIMPLE);
}


PUBLIC int
vrna_plot_coords_circular(const char  *structure,
                          float       **x,
                          float       **y)
{
  return vrna_plot_coords(structure, x, y, VRNA_PLOT_TYPE_CIRCULAR);
}


PUBLIC int
vrna_plot_coords_circular_pt(const short  *pt,
                             float        **x,
                             float        **y)
{
  return vrna_plot_coords_pt(pt, x, y, VRNA_PLOT_TYPE_CIRCULAR);
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE struct vrna_plot_layout_s *
rna_layout(const char   *structure,
           unsigned int plot_type,
           void         *options)
{
  struct vrna_plot_layout_s *layout;
  short                     *pt, *pt_g;
  unsigned int              i, n, ee, Lg, l[3];
  int                       xmin, xmax, ymin, ymax, gb, ge, r;


  n               = strlen(structure);
  layout          = (struct vrna_plot_layout_s *)vrna_alloc(sizeof(struct vrna_plot_layout_s));
  layout->type    = plot_type;
  layout->length  = n;
  layout->x       = NULL;
  layout->y       = NULL;
  layout->arcs    = NULL;

  /* convert dot-bracket string to pair table */
  pt    = vrna_ptable(structure);
  pt_g  = vrna_ptable_copy(pt);

  /* account for possible G-Quadruplexes in dot-bracket string */
  ge = 0;
  while ((ee = vrna_gq_parse(structure + ge, &Lg, l)) > 0) {
    ge  += ee;
    if (4 * Lg + l[0] + l[1] + l[2] > ee) {
      gb = n + ge - 4 * Lg - l[0] - l[1] - l[2] + 1;
    } else {
      gb  = ge - Lg * 4 - l[0] - l[1] - l[2] + 1;
    }

    /* add pseudo-base pair encloding gquad */
    for (i = 0; i < Lg; i++) {
      unsigned int ii, jj;
      ii = (gb + i - 1) % n + 1;
      jj = (n + ge - i - 1) % n + 1;
      pt_g[ii]  = jj;
      pt_g[jj]  = ii;
    }
  }

  switch (plot_type) {
    case VRNA_PLOT_TYPE_SIMPLE:
      i = coords_simple(pt_g,
                        &(layout->x),
                        &(layout->y));
      break;

    case VRNA_PLOT_TYPE_CIRCULAR:
      r = 3 * n;
      i = coords_circular(pt_g,
                          &(layout->x),
                          &(layout->y));

      for (i = 0; i < n; i++) {
        layout->x[i]  *= r;
        layout->x[i]  += r;
        layout->y[i]  *= r;
        layout->y[i]  += r;
      }
      break;

    case VRNA_PLOT_TYPE_TURTLE:
      i = vrna_plot_coords_turtle_pt(pt,
                                     &(layout->x),
                                     &(layout->y),
                                     &(layout->arcs));
      break;

    case VRNA_PLOT_TYPE_PUZZLER:
      i = vrna_plot_coords_puzzler_pt(pt,
                                      &(layout->x),
                                      &(layout->y),
                                      &(layout->arcs),
                                      (vrna_plot_options_puzzler_t *)options);
      break;

#ifdef VRNA_WITH_NAVIEW_LAYOUT
    default:
      i = vrna_plot_coords_naview_pt(pt_g,
                                     &(layout->x),
                                     &(layout->y));
      break;
#else
    case VRNA_PLOT_TYPE_NAVIEW:
      i = 0;
      vrna_log_warning("Naview layout algorithm not available in this version of ViennaRNA Package");
      break;

    default:
      i = vrna_plot_coords_puzzler_pt(pt,
                                      &(layout->x),
                                      &(layout->y),
                                      &(layout->arcs),
                                      (vrna_plot_options_puzzler_t *)options);
      break;

#endif

  }

  if (i != n) {
    vrna_log_warning("strange things happening in vrna_plot_layout*()...");
    layout->bbox[0] = layout->bbox[1] = layout->bbox[2] = layout->bbox[3] = 0;
  } else {
    /* adjust bounding box coordinates */
    xmin  = xmax = layout->x[0];
    ymin  = ymax = layout->y[0];

    for (i = 1; i < n; i++) {
      xmin  = layout->x[i] < xmin ? layout->x[i] : xmin;
      xmax  = layout->x[i] > xmax ? layout->x[i] : xmax;
      ymin  = layout->y[i] < ymin ? layout->y[i] : ymin;
      ymax  = layout->y[i] > ymax ? layout->y[i] : ymax;
    }

    layout->bbox[0] = xmin;
    layout->bbox[1] = ymin;
    layout->bbox[2] = xmax;
    layout->bbox[3] = ymax;
  }

  free(pt);
  free(pt_g);

  return layout;
}


PRIVATE int
coords_simple(const short *pt,
              float       **x,
              float       **y)
{
  float INIT_ANGLE  = 0.;     /* initial bending angle */
  float INIT_X      = 100.;   /* coordinate of first digit */
  float INIT_Y      = 100.;   /* see above */
  float RADIUS      = 15.;

  int   i, length;
  float alpha;

  float *angle;
  int   *loop_size, *stack_size;
  int   lp, stk;

  length      = pt[0];
  angle       = (float *)vrna_alloc((length + 5) * sizeof(float));
  loop_size   = (int *)vrna_alloc(16 + (length / 5) * sizeof(int));
  stack_size  = (int *)vrna_alloc(16 + (length / 5) * sizeof(int));
  lp          = stk = 0;

  *x  = (float *)vrna_alloc(sizeof(float) * (length + 1));
  *y  = (float *)vrna_alloc(sizeof(float) * (length + 1));

  loop(pt, 0, length, angle, stack_size, loop_size, &stk, &lp);

  loop_size[lp] -= 2;     /* correct for cheating with function loop */

  alpha   = INIT_ANGLE;
  (*x)[0] = INIT_X;
  (*y)[0] = INIT_Y;

  for (i = 1; i <= length; i++) {
    (*x)[i] = (*x)[i - 1] + RADIUS * cos(alpha);
    (*y)[i] = (*y)[i - 1] + RADIUS * sin(alpha);
    alpha   += PI - angle[i + 1];
  }

  free(angle);
  free(loop_size);
  free(stack_size);

  return length;
}


/*
 *  i, j are the positions AFTER the last pair of a stack; i.e
 *  i-1 and j+1 are paired.
 */
PRIVATE void
loop(const short  *pair_table,
     int          i,
     int          j,
     float        *angle,
     int          *stack_size,
     int          *loop_size,
     int          *stk,
     int          *lp)
{
  int count = 2;            /*
                             * counts the VERTICES of a loop polygon; that's
                             *   NOT necessarily the number of unpaired bases!
                             *   Upon entry the loop has already 2 vertices, namely
                             *   the pair i-1/j+1.
                             */

  int   r = 0, bubble = 0;  /* bubble counts the unpaired digits in loops */

  int   i_old, partner, k, l, start_k, start_l, fill, ladder;
  int   begin, v, diff;
  float polygon;

  short *remember;

  remember = (short *)vrna_alloc((3 + (j - i) / 5) * 2 * sizeof(short));

  i_old = i - 1, j++;         /* j has now been set to the partner of the
                               * previous pair for correct while-loop
                               * termination.  */
  while (i != j) {
    partner = pair_table[i];
    if ((!partner) || (i == 0)) {
      i++, count++, bubble++;
    } else {
      count         += 2;
      k             = i, l = partner; /* beginning of stack */
      remember[++r] = k;
      remember[++r] = l;
      i             = partner + 1; /* next i for the current loop */

      start_k = k, start_l = l;
      ladder  = 0;
      do
        k++, l--, ladder++;        /* go along the stack region */
      while ((pair_table[k] == l) && (pair_table[k] > k));

      fill = ladder - 2;
      if (ladder >= 2) {
        angle[start_k + 1 + fill] += PIHALF;  /*  Loop entries and    */
        angle[start_l - 1 - fill] += PIHALF;  /*  exits get an        */
        angle[start_k]            += PIHALF;  /*  additional PI/2.    */
        angle[start_l]            += PIHALF;  /*  Why ? (exercise)    */
        if (ladder > 2) {
          for (; fill >= 1; fill--) {
            angle[start_k + fill] = PI;     /*  fill in the angles  */
            angle[start_l - fill] = PI;     /*  for the backbone    */
          }
        }
      }

      stack_size[++(*stk)] = ladder;
      if (k <= l)
        loop(pair_table, k, l, angle, stack_size, loop_size, stk, lp);
    }
  }
  polygon       = PI * (count - 2) / (float)count; /* bending angle in loop polygon */
  remember[++r] = j;
  begin         = i_old < 0 ? 0 : i_old;
  for (v = 1; v <= r; v++) {
    diff = remember[v] - begin;
    for (fill = 0; fill <= diff; fill++)
      angle[begin + fill] += polygon;
    if (v > r)
      break;

    begin = remember[++v];
  }
  loop_size[++(*lp)] = bubble;
  free(remember);
}


PRIVATE int
coords_circular(const short *pt,
                float       **x,
                float       **y)
{
  unsigned int  length = (unsigned int)pt[0];
  unsigned int  i;
  float         d = 2 * PI / length;

  *x  = (float *)vrna_alloc(sizeof(float) * (length + 1));
  *y  = (float *)vrna_alloc(sizeof(float) * (length + 1));

  for (i = 0; i < length; i++) {
    (*x)[i] = cos(i * d - PIHALF);
    (*y)[i] = sin(i * d - PIHALF);
  }

  return length;
}


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
 *###########################################
 *# deprecated functions below              #
 *###########################################
 */
PUBLIC int
simple_xy_coordinates(short *pair_table,
                      float *x,
                      float *y)
{
  if ((pair_table) && (x) && (y)) {
    int   ret, n;
    float *xx, *yy;

    n   = pair_table[0];
    ret = coords_simple(pair_table, &xx, &yy);

    memcpy(x, xx, sizeof(float) * (n + 1));
    memcpy(y, yy, sizeof(float) * (n + 1));

    free(xx);
    free(yy);
  }

  return 0;
}


PUBLIC int
simple_circplot_coordinates(short *pair_table,
                            float *x,
                            float *y)
{
  if ((pair_table) && (x) && (y)) {
    int   ret, n;
    float *xx, *yy;

    n   = pair_table[0];
    ret = coords_circular(pair_table, &xx, &yy);

    memcpy(x, xx, sizeof(float) * (n + 1));
    memcpy(y, yy, sizeof(float) * (n + 1));

    free(xx);
    free(yy);
  }

  return 0;
}


#endif
