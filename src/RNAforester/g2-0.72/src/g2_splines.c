/*****************************************************************************
**  Copyright (C) 1998-2001  Ljubomir Milanovic & Horst Wagner
**  Copyright (C) 1999-2006  Tijs Michels
**  This file is part of the g2 library
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Lesser General Public
**  License as published by the Free Software Foundation; either
**  version 2.1 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Lesser General Public License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License along with this library; if not, write to the Free Software
**  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
******************************************************************************/
/*
 * g2_splines.c
 * Tijs Michels
 * tijs@users.sourceforge.net
 * 06/16/99 : initial release
 * 19/02/06 : eliminated duplicates by using pointers to functions
 */

#include <math.h>
#include <stdio.h>
#include "g2.h"
#include "g2_util.h"

/**
 * \ingroup graphic
 * \defgroup splines splines
 */

typedef void calc_f( int, const double *, int,    double *);
typedef void calc_d( int, const double *, double, double *);

static calc_f g2_c_spline;
static calc_f g2_c_b_spline;
static calc_d g2_c_raspln;
static calc_f g2_c_para_3;
static calc_f g2_c_para_5;

static void g2_p_spline(int id, int n, const double *points, int o, calc_f *f)
{
   const int m = (n-1)*o+1;
   double * const sxy = (double *) g2_malloc(m*2*sizeof(double));

   (*f)(n, points, m, sxy);
   g2_poly_line(id, m, sxy);

   g2_free(sxy);
}

static void g2_p_filled_spline(int id, int n, const double *points, int o, calc_f *f)
{
   const int m = (n-1)*o+2;
   double * const sxy = (double *) g2_malloc(m*2*sizeof(double));

   (*f)(n, points, m, sxy);
   sxy[m+m-2] = points[n+n-2];
   sxy[m+m-1] = points[1];
   g2_filled_polygon(id, m, sxy);
   g2_free(sxy);
}

static void g2_split(int n, const double *points, double *x, double *y)
{
   int i;
   for (i = 0; i < n; i++) {
      x[i] = points[i+i];
      y[i] = points[i+i+1];
   }
}

#define eps 1.e-12

void g2_c_spline(int n, const double *points, int m, double *sxy)

/*
 *	FUNCTIONAL DESCRIPTION:
 *
 *	Compute a curve of m points (sx[j],sy[j])
 *	-- j being a positive integer < m --
 *	passing through the n data points (x[i],y[i])
 *	-- i being a positive integer < n --
 *	supplied by the user.
 *	The procedure to determine sy[j] involves
 *	Young's method of successive over-relaxation.
 *
 *	FORMAL ARGUMENTS:
 *
 *	n			number of data points
 *	points			data points (x[i],y[i])
 *	m			number of interpolated points; m = (n-1)*o+1
 *				for o curve points for every data point
 *	sxy			interpolated points (sx[j],sy[j])
 *
 *	IMPLICIT INPUTS:	NONE
 *	IMPLICIT OUTPUTS:	NONE
 *	SIDE EFFECTS:		NONE
 *
 *	REFERENCES:
 *
 *	1. Ralston and Wilf, Mathematical Methods for Digital Computers,
 *	   Vol. II, John Wiley and Sons, New York 1967, pp. 156-158.
 *	2. Greville, T.N.E., Ed., Proceedings of An Advanced Seminar
 *	   Conducted by the Mathematics Research Center, U.S. Army,
 *	   University of Wisconsin, Madison. October 7-9, 1968. Theory
 *	   and Applications of Spline Functions, Academic Press,
 *	   New York / London 1969, pp. 156-167.
 *
 *	AUTHORS:
 *
 *	Josef Heinen	04/06/88	<J.Heinen@KFA-Juelich.de>
 *	Tijs Michels	06/16/99	<t.michels@vimec.nl>
 */

{
   int i, j;
   double *x, *y, *g, *h;
   double k, u, delta_g;

   if (n < 3) {
      fputs("\nERROR calling function \"g2_c_spline\":\n"
	    "number of data points input should be at least three\n", stderr);
      return;
   }
   if ((m-1)%(n-1)) {
      fputs("\nWARNING from function \"g2_c_spline\":\n"
	    "number of curve points output for every data point input "
	    "is not an integer\n", stderr);
   }

   x = (double *) g2_malloc(n*4*sizeof(double));
   y = x + n;
   g = y + n;
   h = g + n; /* for the constant copy of g */
   g2_split(n, points, x, y);

   n--; /* last value index */
   k =  x[0]; /* look up once */
   u = (x[n] - k) / (m - 1); /* calculate step outside loop */
   for (j = 0; j < m; j++)	sxy[j+j] = j * u + k; /* x-coordinates */

   for (i = 1; i < n; i++) {
      g[i] = 2. * ((y[i+1] - y[i]) / (x[i+1] - x[i]) -
		   (y[i] - y[i-1]) / (x[i] - x[i-1]))
	/ (x[i+1] - x[i-1]); /* whereas g[i] will later be changed repeatedly */
      h[i] = 1.5 * g[i];     /* copy h[i] of g[i] will remain constant */
   }

   k = 0.;

   do {
      for (u = 0., i = 1; i < n; i++) {
	 delta_g = .5 * (x[i] - x[i-1]) / (x[i+1] - x[i-1]);
	 delta_g = (h[i] -
		    g[i] -
		    g[i-1] * delta_g -      /* 8. - 4 * sqrt(3.) */
		    g[i+1] * (.5 - delta_g)) * 1.0717967697244907832;
	 g[i] += delta_g;

	 if (fabs(delta_g) > u) u = fabs(delta_g);
      }	/* On loop termination u holds the largest delta_g. */

      if (k == 0.)	k = u * eps;
	/* Only executed once, at the end of pass one. So k preserves
	 * the largest delta_g of pass one, multiplied by eps.
	 */
   } while (u > k);

   m += m, i = 1, j = 0;
   do {
      u = sxy[j++]; /* x-coordinate */

      while (x[i] < u)	i++;

      if (--i > n)	i = n;

      k = (u - x[i]) / (x[i+1] - x[i]); /* calculate outside loop */
      sxy[j++] = y[i] +
	(y[i+1] - y[i]) * k +
	(u - x[i]) * (u - x[i+1]) *
	((2. - k) * g[i] +
	 (1. + k) * g[i+1]) / 6.; /* y-coordinate */
   } while (j < m);
   g2_free(x);
}

/**
 *
 * Using Young's method of successive over-relaxation,
 * plot a spline curve with \a o interpolated points per data point.
 * So the larger \a o, the more fluent the curve.
 *
 * \param dev device id
 * \param n number of data points (not the size of buffer \a points)
 * \param points buffer of \a n data points x1, y1, ... x\a n, y\a n
 * \param o number of interpolated points per data point
 *
 * \ingroup splines
 */
void g2_spline(int dev, int n, double *points, int o)
{
   g2_p_spline(dev, n, points, o, g2_c_spline);
}

/**
 *
 * Using Young's method of successive over-relaxation,
 * plot a filled spline curve with \a o interpolated points per data point.
 * So the larger \a o, the more fluent the curve.
 *
 * \param dev device id
 * \param n number of data points (not the size of buffer \a points)
 * \param points buffer of \a n data points x1, y1, ... x\a n, y\a n
 * \param o number of interpolated points per data point
 *
 * \ingroup splines
 */
void g2_filled_spline(int dev, int n, double *points, int o)
{
   g2_p_filled_spline(dev, n, points, o, g2_c_spline);
}

void g2_c_b_spline(int n, const double *points, int m, double *sxy)

/*
 * g2_c_b_spline takes n input points. It uses parameter t
 * to compute sx(t) and sy(t) respectively
 */

{
   int i, j;
   double *x, *y;
   double t, bl1, bl2, bl3, bl4;
   double interval, xi_3, yi_3, xi, yi;

   if (n < 3) {
      fputs("\nERROR calling function \"g2_c_b_spline\":\n"
	    "number of data points input should be at least three\n", stderr);
      return;
   }
   x = (double *) g2_malloc(n*2*sizeof(double));
   y = x + n;
   g2_split(n, points, x, y);

   m--; /* last value index */
   n--; /* last value index */
   interval = (double)n / m;

   for (m += m, i = 2, j = 0; i <= n+1; i++) {
      if (i == 2) {
	 xi_3 = 2 * x[0] - x[1];
	 yi_3 = 2 * y[0] - y[1];
      } else {
	 xi_3 = x[i-3];
	 yi_3 = y[i-3];
      }
      if (i == n+1) {
	 xi = 2 * x[n] - x[n-1];
	 yi = 2 * y[n] - y[n-1];
      } else {
	 xi = x[i];
	 yi = y[i];
      }

      t = fmod(j * interval, 1.);

      while (t < 1. && j < m) {
	 bl1 = (1. - t);
	 bl2 = t * t;	/* t^2 */
	 bl4 = t * bl2;	/* t^3 */
	 bl3 = bl4 - bl2;

	 bl1 = bl1 * bl1 * bl1;
	 bl2 = 3. * (bl3 - bl2) + 4.;
	 bl3 = 3. * (  t - bl3) + 1.;

	 sxy[j++] = (bl1 * xi_3 + bl2 * x[i-2] + bl3 * x[i-1] + bl4 * xi) / 6.; /* x-coordinate */
	 sxy[j++] = (bl1 * yi_3 + bl2 * y[i-2] + bl3 * y[i-1] + bl4 * yi) / 6.; /* y-coordinate */

	 t += interval;
      }
   }
   sxy[m]   = x[n];
   sxy[m+1] = y[n];
   g2_free(x);
}

/**
 *
 * Plot a b-spline curve with \a o interpolated points per data point.
 * So the larger \a o, the more fluent the curve.
 * For most averaging purposes, this is the right spline.
 *
 * \param dev device id
 * \param n number of data points (not the size of buffer \a points)
 * \param points buffer of \a n data points x1, y1, ... x\a n, y\a n
 * \param o number of interpolated points per data point
 *
 * \ingroup splines
 */
void g2_b_spline(int dev, int n, double *points, int o)
{
   g2_p_spline(dev, n, points, o, g2_c_b_spline);
}

/**
 *
 * Plot a filled b-spline curve with \a o interpolated points per data point.
 * So the larger \a o, the more fluent the curve.
 * For most averaging purposes, this is the right spline.
 *
 * \param dev device id
 * \param n number of data points (not the size of buffer \a points)
 * \param points buffer of \a n data points x1, y1, ... x\a n, y\a n
 * \param o number of interpolated points per data point
 *
 * \ingroup splines
 */
void g2_filled_b_spline(int dev, int n, double *points, int o)
{
   g2_p_filled_spline(dev, n, points, o, g2_c_b_spline);
}

/*
 *	FUNCTION g2_c_raspln
 *
 *	FUNCTIONAL DESCRIPTION:
 *
 *	This function draws a piecewise cubic polynomial through
 *	the specified data points. The (n-1) cubic polynomials are
 *	basically parametric cubic Hermite polynomials through the
 *	n specified data points with tangent values at the data
 *	points determined by a weighted average of the slopes of
 *	the secant lines. A tension parameter "tn" is provided to
 *	adjust the length of the tangent vector at the data points.
 *	This allows the "roundness" of the curve to be adjusted.
 *	For further information and references on this technique see:
 *
 *	D. Kochanek and R. Bartels, Interpolating Splines With Local
 *	Tension, Continuity and Bias Control, Computer Graphics,
 *	18(1984)3.
 *
 *	AUTHORS:
 *
 *	Dennis Mikkelson	distributed in GPLOT	Jan 7, 1988	F77
 *	Tijs Michels		t.michels@vimec.nl	Jun 7, 1999	C
 *
 *	FORMAL ARGUMENTS:
 *
 *	n	number of data points, n > 2
 *	points	double array holding the x and y-coords of the data points
 *	tn	double parameter in [0.0, 2.0]. When tn = 0.0,
 *		the curve through the data points is very rounded.
 *		As tn increases the curve is gradually pulled tighter.
 *		When tn = 2.0, the curve is essentially a polyline
 *		through the given data points.
 *	sxy	double array holding the coords of the spline curve
 *
 *	IMPLICIT INPUTS:	NONE
 *	IMPLICIT OUTPUTS:	NONE
 *	SIDE EFFECTS:		NONE
 */

#define nb 40
/*
 * Number of straight connecting lines of which each polynomial consists.
 * So between one data point and the next, (nb-1) points are placed.
 */

void g2_c_raspln(int n, const double *points, double tn, double *sxy)
{
   int i, j;
   double bias, tnFactor, tangentL1, tangentL2;
   double D1x, D1y, D2x, D2y, t1x, t1y, t2x, t2y;
   double h1[nb+1];	/* Values of the Hermite basis functions */
   double h2[nb+1];	/* at nb+1 evenly spaced points in [0,1] */
   double h3[nb+1];
   double h4[nb+1];

   double * const x = (double *) g2_malloc(n*2*sizeof(double));
   double * const y = x + n;
   g2_split(n, points, x, y);

/*
 * First, store the values of the Hermite basis functions in a table h[ ]
 * so no time is wasted recalculating them
 */
   for (i = 0; i < nb+1; i++) {
      double t, tt, ttt;
      t = (double) i / nb;
      tt  = t * t;
      ttt = t * tt;
      h1[i] =  2. * ttt - 3. * tt + 1.;
      h2[i] = -2. * ttt + 3. * tt;
      h3[i] =       ttt - 2. * tt + t;
      h4[i] =       ttt -      tt;
   }

/*
 * Set local tnFactor based on input parameter tn
 */
   if (tn <= 0.) {
      tnFactor = 2.;
      fputs("g2_c_raspln: Using Tension Factor 0.0: very rounded", stderr);
   }
   else if (tn >= 2.) {
      tnFactor = 0.;
      fputs("g2_c_raspln: Using Tension Factor 2.0: not rounded at all", stderr);
   }
   else			tnFactor = 2. - tn;

   D1x = D1y = 0.; /* first point has no preceding point */
   for (j = 0; j < n - 2; j++) {
      t1x = x[j+1] - x[j];
      t1y = y[j+1] - y[j];
      t2x = x[j+2] - x[j+1];
      t2y = y[j+2] - y[j+1];
      tangentL1 = t1x * t1x + t1y * t1y;
      tangentL2 = t2x * t2x + t2y * t2y;
      if (tangentL1 + tangentL2 == 0) bias = .5;
      else bias = tangentL2 / (tangentL1 + tangentL2);
      D2x = tnFactor * (bias  * t1x + (1 - bias) * t2x);
      D2y = tnFactor * (bias  * t1y + (1 - bias) * t2y);
      for (i = 0; i < nb; i++) {
	sxy[2 * nb * j + i + i] =
	   h1[i] * x[j] + h2[i] * x[j+1] + h3[i] * D1x + h4[i] * D2x;
	sxy[2 * nb * j + i + i + 1] =
	   h1[i] * y[j] + h2[i] * y[j+1] + h3[i] * D1y + h4[i] * D2y;
      }
      D1x = D2x; /* store as preceding point in */
      D1y = D2y; /* the next pass */
   }

/*
 * Do the last subinterval as a special case since no point follows the
 * last point
 */
   for (i = 0; i < nb+1; i++) {
      sxy[2 * nb * (n-2) + i + i] =
	h1[i] * x[n-2] + h2[i] * x[n-1] + h3[i] * D1x;
      sxy[2 * nb * (n-2) + i + i + 1] =
	h1[i] * y[n-2] + h2[i] * y[n-1] + h3[i] * D1y;
   }
   g2_free(x);
}

/**
 *
 * Plot a piecewise cubic polynomial with adjustable roundness
 * through the given data points.
 * Each Hermite polynomial between two data points is made up of 40 lines.
 * Tension factor \a tn must be between 0.0 (very rounded)
 * and 2.0 (not rounded at all, i.e. essentially a \ref g2_poly_line "polyline").
 *
 * \param dev device id
 * \param n number of data points (not the size of buffer \a points)
 * \param points buffer of \a n data points x1, y1, ... x\a n, y\a n
 * \param tn tension factor in the range [0.0, 2.0]
 *
 * \ingroup splines
 */
void g2_raspln(int dev, int n, double *points, double tn)
{
   const int m = (n-1)*nb+1;
   double * const sxy = (double *) g2_malloc(m*2*sizeof(double)); /* coords of the entire spline curve */

   g2_c_raspln(n, points, tn, sxy);
   g2_poly_line(dev, m, sxy);

   g2_free(sxy);
}

/**
 *
 * Plot a filled piecewise cubic polynomial with adjustable roundness
 * through the given data points.
 * Each Hermite polynomial between two data points is made up of 40 lines.
 * Tension factor \a tn must be between 0.0 (very rounded)
 * and 2.0 (not rounded at all, i.e. essentially a \ref g2_poly_line "polyline").
 *
 * \param dev device id
 * \param n number of data points (not the size of buffer \a points)
 * \param points buffer of \a n data points x1, y1, ... x\a n, y\a n
 * \param tn tension factor in the range [0.0, 2.0]
 *
 * \ingroup splines
 */
void g2_filled_raspln(int dev, int n, double *points, double tn)
{
   const int m = (n-1)*nb+2;
   double * const sxy = (double *) g2_malloc(m*2*sizeof(double)); /* coords of the entire spline curve */

   g2_c_raspln(n, points, tn, sxy);
   sxy[m+m-2] = points[n+n-2];
   sxy[m+m-1] = points[1];
   g2_filled_polygon(dev, m, sxy);

   g2_free(sxy);
}

/* ---- And now for a rather different approach ---- */

/*
 *	FUNCTION g2_c_newton
 *
 *	FUNCTIONAL DESCRIPTION:
 *
 *	Use Newton's Divided Differences to calculate an interpolation
 *	polynomial through the specified data points.
 *	This function is called by
 *		g2_c_para_3 and
 *		g2_c_para_5.
 *
 *	Dennis Mikkelson	distributed in GPLOT	Jan  5, 1988	F77
 *	Tijs Michels		t.michels@vimec.nl	Jun 16, 1999	C
 *
 *	FORMAL ARGUMENTS:
 *
 *	n	number of entries in c1 and c2, 4 <= n <= MaxPts
 *		for para_3	(degree 3)	n = 4
 *		for para_5	(degree 5)	n = 6
 *		for para_i	(degree i)	n = (i + 1)
 *	c1	double array holding at most MaxPts values giving the
 *		first  coords of the points to be interpolated
 *	c2	double array holding at most MaxPts values giving the
 *		second coords of the points to be interpolated
 *	o	number of points at which the interpolation
 *		polynomial is to be evaluated
 *	xv	double array holding o points at which to
 *		evaluate the interpolation polynomial
 *	yv	double array holding upon return the values of the
 *		interpolation polynomial at the corresponding points in xv
 *
 *		yv is the OUTPUT
 *
 *	IMPLICIT INPUTS:	NONE
 *	IMPLICIT OUTPUTS:	NONE
 *	SIDE EFFECTS:		NONE
 */

#define MaxPts 21
#define xstr(s) __str(s)
#define __str(s) #s

/*
 * Maximum number of data points allowed
 * 21 would correspond to a polynomial of degree 20
 */

static void g2_c_newton(int n, const double *c1, const double *c2,
			int o, const double *xv, double *yv)
{
   int i, j;
   double p, s;
   double ddt[MaxPts][MaxPts];		/* Divided Difference Table */

   if (n < 4) {
      fputs("g2_c_newton: Error! Less than 4 points passed "
	    "to function g2_c_newton\n", stderr);
      return;
   }

   if (n > MaxPts) {
      fputs("g2_c_newton: Error! More than " xstr(MaxPts) " points passed "
	    "to function g2_c_newton\n", stderr);
      return;
   }

/* First, build the divided difference table */

   for (i = 0; i < n; i++)	ddt[i][0] = c2[i];
   for (j = 1; j < n; j++) {
      for (i = 0; i < n - j; i++)
	ddt[i][j] = (ddt[i+1][j-1] - ddt[i][j-1]) / (c1[i+j] - c1[i]);
   }

/* Next, evaluate the polynomial at the specified points */

   for (i = 0; i < o; i++) {
      for (p = 1., s = ddt[0][0], j = 1; j < n; j++) {
	 p *= xv[i] - c1[j-1];
	 s += p * ddt[0][j];
      }
      yv[i] = s;
   }
}

/*
 *	FUNCTION: g2_c_para_3
 *
 *	FUNCTIONAL DESCRIPTION:
 *
 *	This function draws a piecewise parametric interpolation
 *	polynomial of degree 3 through the specified data points.
 *	The effect is similar to that obtained using DISSPLA to
 *	draw a curve after a call to the DISSPLA routine PARA3.
 *	The curve is parameterized using an approximation to the
 *	curve's arc length. The basic interpolation is done
 *	using function g2_c_newton.
 *
 *	Dennis Mikkelson	distributed in GPLOT	Jan  7, 1988	F77
 *	Tijs Michels		t.michels@vimec.nl	Jun 17, 1999	C
 *
 *	FORMAL ARGUMENTS:
 *
 *	n	number of data points through which to draw the curve
 *	points	double array containing the x and y-coords of the data points
 *
 *	IMPLICIT INPUTS:	NONE
 *	IMPLICIT OUTPUTS:	NONE
 *	SIDE EFFECTS:		NONE
 */

/*
 * #undef  nb
 * #define nb 40
 * Number of straight connecting lines of which each polynomial consists.
 * So between one data point and the next, (nb-1) points are placed.
 */

void g2_c_para_3(int n, const double *points, int m, double *sxy)
{
#define dgr	(3+1)
#define nb2	(nb*2)
   int i, j;
   double x1t, y1t;
   double o, step;
   double X[nb2];		/* x-coords of the current curve piece */
   double Y[nb2];		/* y-coords of the current curve piece */
   double t[dgr];		/* data point parameter values */
   double Xpts[dgr];		/* x-coords data point subsection */
   double Ypts[dgr];		/* y-coords data point subsection */
   double s[nb2];		/* parameter values at which to interpolate */
   m = m;			/* dummy argument; stop compiler complaints */

   /* Do first TWO subintervals first */

   g2_split(dgr, points, Xpts, Ypts);

   t[0] = 0.;
   for (i = 1; i < dgr; i++) {
      x1t = Xpts[i] - Xpts[i-1];
      y1t = Ypts[i] - Ypts[i-1];
      t[i] = t[i-1] + sqrt(x1t * x1t + y1t * y1t);
   }

   step = t[2] / nb2;
   for (i = 0; i < nb2; i++)	s[i] = i * step;

   g2_c_newton(dgr, t, Xpts, nb2, s, X);
   g2_c_newton(dgr, t, Ypts, nb2, s, Y);
   for (i = 0; i < nb2; i++) {
      sxy[i+i]   = X[i];
      sxy[i+i+1] = Y[i];
   }

   /* Next, do later central subintervals */

   for (j = 1; j < n - dgr + 1; j++) {
      g2_split(dgr, points + j + j, Xpts, Ypts);

      for (i = 1; i < dgr; i++) {
	 x1t = Xpts[i] - Xpts[i-1];
	 y1t = Ypts[i] - Ypts[i-1];
	 t[i] = t[i-1] + sqrt(x1t * x1t + y1t * y1t);
      }

      o = t[1]; /* look up once */
      step = (t[2] - o) / nb;
      for (i = 0; i < nb; i++)	s[i] = i * step + o;

      g2_c_newton(dgr, t, Xpts, nb, s, X);
      g2_c_newton(dgr, t, Ypts, nb, s, Y);

      for (i = 0; i < nb; i++) {
	 sxy[(j + 1) * nb2 + i + i]     = X[i];
	 sxy[(j + 1) * nb2 + i + i + 1] = Y[i];
      }
   }

   /* Now do last subinterval */

   o = t[2];
   step = (t[3] - o) / nb;
   for (i = 0; i < nb; i++)	s[i] = i * step + o;

   g2_c_newton(dgr, t, Xpts, nb, s, X);
   g2_c_newton(dgr, t, Ypts, nb, s, Y);

   for (i = 0; i < nb; i++) {
      sxy[(n - dgr + 2) * nb2 + i + i]     = X[i];
      sxy[(n - dgr + 2) * nb2 + i + i + 1] = Y[i];
   }
   sxy[(n - 1) * nb2]     = points[n+n-2];
   sxy[(n - 1) * nb2 + 1] = points[n+n-1];
}

/**
 *
 * Using Newton's Divided Differences method, plot a piecewise
 * parametric interpolation polynomial of degree 3
 * through the given data points.
 *
 * \param dev device id
 * \param n number of data points (not the size of buffer \a points)
 * \param points buffer of \a n data points x1, y1, ... x\a n, y\a n
 *
 * \ingroup splines
 */
void g2_para_3(int dev, int n, double *points)
{
   g2_p_spline(dev, n, points, nb, g2_c_para_3);
}

/**
 *
 * Using Newton's Divided Differences method, plot a filled piecewise
 * parametric interpolation polynomial of degree 3
 * through the given data points.
 *
 * \param dev device id
 * \param n number of data points (not the size of buffer \a points)
 * \param points buffer of \a n data points x1, y1, ... x\a n, y\a n
 *
 * \ingroup splines
 */
void g2_filled_para_3(int dev, int n, double *points)
{
   g2_p_filled_spline(dev, n, points, nb, g2_c_para_3);
}

/*
 *	FUNCTION: g2_c_para_5
 *
 *	As g2_c_para_3, but now plot a polynomial of degree 5
 */

/*
 * Number of straight connecting lines of which each polynomial consists.
 * So between one data point and the next, (nb-1) points are placed.
 */

void g2_c_para_5(int n, const double *points, int m, double *sxy)
{
#undef	dgr
#define dgr	(5+1)
#define nb3	(nb*3)
   int i, j;
   double x1t, y1t;
   double o, step;
   double X[nb3];		/* x-coords of the current curve piece */
   double Y[nb3];		/* y-coords of the current curve piece */
   double t[dgr];		/* data point parameter values */
   double Xpts[dgr];		/* x-coords data point subsection */
   double Ypts[dgr];		/* y-coords data point subsection */
   double s[nb3];		/* parameter values at which to interpolate */
   m = m;			/* dummy argument; stop compiler complaints */

   /* Do first THREE subintervals first */

   g2_split(dgr, points, Xpts, Ypts);

   t[0] = 0.;
   for (i = 1; i < dgr; i++) {
      x1t = Xpts[i] - Xpts[i-1];
      y1t = Ypts[i] - Ypts[i-1];
      t[i] = t[i-1] + sqrt(x1t * x1t + y1t * y1t);
   }

   step = t[3] / nb3;
   for (i = 0; i < nb3; i++)	s[i] = i * step;

   g2_c_newton(dgr, t, Xpts, nb3, s, X);
   g2_c_newton(dgr, t, Ypts, nb3, s, Y);
   for (i = 0; i < nb3; i++) {
      sxy[i+i]   = X[i];
      sxy[i+i+1] = Y[i];
   }

   /* Next, do later central subintervals */

   for (j = 1; j < n - dgr + 1; j++) {
      g2_split(dgr, points + j + j, Xpts, Ypts);

      for (i = 1; i < dgr; i++) {
	 x1t = Xpts[i] - Xpts[i-1];
	 y1t = Ypts[i] - Ypts[i-1];
	 t[i] = t[i-1] + sqrt(x1t * x1t + y1t * y1t);
      }

      o = t[2]; /* look up once */
      step = (t[3] - o) / nb;
      for (i = 0; i < nb; i++)	s[i] = i * step + o;

      g2_c_newton(dgr, t, Xpts, nb, s, X);
      g2_c_newton(dgr, t, Ypts, nb, s, Y);

      for (i = 0; i < nb; i++) {
	 sxy[(j + 2) * nb2 + i + i]     = X[i];
	 sxy[(j + 2) * nb2 + i + i + 1] = Y[i];
      }
   }

   /* Now do last TWO subinterval */

   o = t[3];
   step = (t[5] - o) / nb2;
   for (i = 0; i < nb2; i++)	s[i] = i * step + o;

   g2_c_newton(dgr, t, Xpts, nb2, s, X);
   g2_c_newton(dgr, t, Ypts, nb2, s, Y);

   for (i = 0; i < nb2; i++) {
      sxy[(n - dgr + 3) * nb2 + i + i]     = X[i];
      sxy[(n - dgr + 3) * nb2 + i + i + 1] = Y[i];
   }
   sxy[(n - 1) * nb2]     = points[n+n-2];
   sxy[(n - 1) * nb2 + 1] = points[n+n-1];
}

/**
 *
 * Using Newton's Divided Differences method, plot a piecewise
 * parametric interpolation polynomial of degree 5
 * through the given data points.
 *
 * \param dev device id
 * \param n number of data points (not the size of buffer \a points)
 * \param points buffer of \a n data points x1, y1, ... x\a n, y\a n
 *
 * \ingroup splines
 */
void g2_para_5(int dev, int n, double *points)
{
   g2_p_spline(dev, n, points, nb, g2_c_para_5);
}

/**
 *
 * Using Newton's Divided Differences method, plot a filled piecewise
 * parametric interpolation polynomial of degree 5
 * through the given data points.
 *
 * \param dev device id
 * \param n number of data points (not the size of buffer \a points)
 * \param points buffer of \a n data points x1, y1, ... x\a n, y\a n
 *
 * \ingroup splines
 */
void g2_filled_para_5(int dev, int n, double *points)
{
   g2_p_filled_spline(dev, n, points, nb, g2_c_para_5);
}
