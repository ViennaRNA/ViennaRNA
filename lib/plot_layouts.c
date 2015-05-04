/**
 * This file is a container for all plotting layout algorithms
 *
 *  c Ronny Lorenz
 *    The ViennaRNA Package
 */

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"
#include "plot_layouts.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define PUBLIC
#define PRIVATE static

PUBLIC  int     rna_plot_type = 1;  /* 0 = simple, 1 = naview, 2 = circular plot */

PRIVATE float   *angle;
PRIVATE int     *loop_size, *stack_size;
PRIVATE int     lp, stk;

PRIVATE void  loop(int i, int j, short *pair_table);

#ifdef _OPENMP
/* NOTE: all threadprivate variables are uninitialized when entering a thread! */
#pragma omp threadprivate(angle, loop_size, stack_size, lp, stk)
#endif

/*---------------------------------------------------------------------------*/

PUBLIC int simple_xy_coordinates(short *pair_table, float *x, float *y)
{
  float INIT_ANGLE=0.;     /* initial bending angle */
  float INIT_X = 100.;     /* coordinate of first digit */
  float INIT_Y = 100.;     /* see above */
  float RADIUS =  15.;

  int i, length;
  float  alpha;

  length = pair_table[0];
  angle =      (float*) space( (length+5)*sizeof(float) );
  loop_size  =   (int*) space( 16+(length/5)*sizeof(int) );
  stack_size =   (int*) space( 16+(length/5)*sizeof(int) );
  lp = stk = 0;
  loop(0, length+1, pair_table);
  loop_size[lp] -= 2;     /* correct for cheating with function loop */

  alpha = INIT_ANGLE;
  x[0]  = INIT_X;
  y[0]  = INIT_Y;

  for (i = 1; i <= length; i++) {
    x[i] = x[i-1]+RADIUS*cos(alpha);
    y[i] = y[i-1]+RADIUS*sin(alpha);
    alpha += PI-angle[i+1];
  }
  free(angle);
  free(loop_size);
  free(stack_size);

  return length;

}

/*---------------------------------------------------------------------------*/

PRIVATE void loop(int i, int j, short *pair_table)
             /* i, j are the positions AFTER the last pair of a stack; i.e
                i-1 and j+1 are paired. */
{
  int    count = 2;   /* counts the VERTICES of a loop polygon; that's
                           NOT necessarily the number of unpaired bases!
                           Upon entry the loop has already 2 vertices, namely
                           the pair i-1/j+1.  */

  int    r = 0, bubble = 0; /* bubble counts the unpaired digits in loops */

  int    i_old, partner, k, l, start_k, start_l, fill, ladder;
  int    begin, v, diff;
  float  polygon;

  short *remember;

  remember = (short *) space((1+(j-i)/5)*2*sizeof(short));

  i_old = i-1, j++;         /* j has now been set to the partner of the
                               previous pair for correct while-loop
                               termination.  */
  while (i != j) {
    partner = pair_table[i];
    if ((!partner) || (i==0))
      i++, count++, bubble++;
    else {
      count += 2;
      k = i, l = partner;    /* beginning of stack */
      remember[++r] = k;
      remember[++r] = l;
      i = partner+1;         /* next i for the current loop */

      start_k = k, start_l = l;
      ladder = 0;
      do {
        k++, l--, ladder++;        /* go along the stack region */
      }
      while (pair_table[k] == l);

      fill = ladder-2;
      if (ladder >= 2) {
        angle[start_k+1+fill] += PIHALF;   /*  Loop entries and    */
        angle[start_l-1-fill] += PIHALF;   /*  exits get an        */
        angle[start_k]        += PIHALF;   /*  additional PI/2.    */
        angle[start_l]        += PIHALF;   /*  Why ? (exercise)    */
        if (ladder > 2) {
          for (; fill >= 1; fill--) {
            angle[start_k+fill] = PI;    /*  fill in the angles  */
            angle[start_l-fill] = PI;    /*  for the backbone    */
          }
        }
      }
      stack_size[++stk] = ladder;
      loop(k, l, pair_table);
    }
  }
  polygon = PI*(count-2)/(float)count; /* bending angle in loop polygon */
  remember[++r] = j;
  begin = i_old < 0 ? 0 : i_old;
  for (v = 1; v <= r; v++) {
    diff  = remember[v]-begin;
    for (fill = 0; fill <= diff; fill++)
      angle[begin+fill] += polygon;
    if (v > r)
      break;
    begin = remember[++v];
  }
  loop_size[++lp] = bubble;
  free(remember);
}

/*---------------------------------------------------------------------------*/

PUBLIC int simple_circplot_coordinates(short *pair_table, float *x, float *y){
  unsigned int  length = (unsigned int) pair_table[0];
  unsigned int  i;
  float         d = 2*PI/length;
  for(i=0; i < length; i++){
    x[i] = cos(i * d - PI/2);
    y[i] = sin(i * d - PI/2);
  }
  return length;
}
