/**
 * \file plot_layouts.h
 *
 * \brief Secondary structure plot layout algorithms
 *
 *  c Ronny Lorenz
 *    The ViennaRNA Package
 */
#ifndef __VIENNA_RNA_PACKAGE_PLOT_LAYOUTS_H__
#define __VIENNA_RNA_PACKAGE_PLOT_LAYOUTS_H__

#include "data_structures.h"
#include "naview.h"

#ifndef PI
#define  PI       3.141592654
#endif
#define  PIHALF       PI/2.


/**
 *  \brief Definition of Plot type <i>simple</i>
 *
 *  This is the plot type definition for several RNA structure plotting functions telling
 *  them to use <b>Simple</b> plotting algorithm
 *
 *  \see rna_plot_type, PS_rna_plot_a(), PS_rna_plot(), svg_rna_plot(), gmlRNA(), ssv_rna_plot(), xrna_plot()
 */
#define VRNA_PLOT_TYPE_SIMPLE     0

/**
 *  \brief Definition of Plot type <i>Naview</i>
 *
 *  This is the plot type definition for several RNA structure plotting functions telling
 *  them to use <b>Naview</b> plotting algorithm
 *
 *  \see rna_plot_type, PS_rna_plot_a(), PS_rna_plot(), svg_rna_plot(), gmlRNA(), ssv_rna_plot(), xrna_plot()
 */
#define VRNA_PLOT_TYPE_NAVIEW     1

/**
 *  \brief Definition of Plot type <i>Circular</i>
 *
 *  This is the plot type definition for several RNA structure plotting functions telling
 *  them to produce a <b>Circular plot</b>
 *
 *  \see rna_plot_type, PS_rna_plot_a(), PS_rna_plot(), svg_rna_plot(), gmlRNA(), ssv_rna_plot(), xrna_plot()
 */
#define VRNA_PLOT_TYPE_CIRCULAR   2


/**
 *  \brief Switch for changing the secondary structure layout algorithm
 *
 *  Current possibility are 0 for a simple radial drawing or 1 for the modified
 *  radial drawing taken from the \e naview program of \ref bruccoleri_88 "Bruccoleri & Heinrich (1988)".
 *
 *  \see #VRNA_PLOT_TYPE_SIMPLE, #VRNA_PLOT_TYPE_NAVIEW, #VRNA_PLOT_TYPE_CIRCULAR
 *
 */
extern int rna_plot_type;

/**
 *  \brief Calculate nucleotide coordinates for secondary structure plot the <i>Simple way</i>
 *
 *  \see make_pair_table(), rna_plot_type, simple_circplot_coordinates(), naview_xy_coordinates(), PS_rna_plot_a(),
 *  PS_rna_plot, svg_rna_plot()
 *
 *  \param  pair_table  The pair table of the secondary structure
 *  \param  x           a pointer to an array with enough allocated space to hold the x coordinates
 *  \param  y           a pointer to an array with enough allocated space to hold the y coordinates
 *  \return             length of sequence on success, 0 otherwise
 */
int simple_xy_coordinates(short *pair_table, float *X, float *Y);

/**
 *  \brief Calculate nucleotide coordinates for <i>Circular Plot</i>
 *
 *  This function calculates the coordinates of nucleotides mapped in equal distancies onto a unit circle.
 *  Additional to the actual R<sup>2</sup> coordinates the function provides a radius scaling factor
 *  for all nucleotides that are involved in a base pairing. This scaling factor may be used for
 *  computing coordinates for a second tangential point \f$P^t\f$ needed to draw quadratic bezier curves.
 *  E.g. \f$ P^{t}_x[i] = X[i] * R[i]\f$ and \f$P^{t}_y[i] = Y[i] * R[i]\f$
 *
 *  \see make_pair_table(), rna_plot_type, simple_xy_coordinates(), naview_xy_coordinates(), PS_rna_plot_a(),
 *  PS_rna_plot, svg_rna_plot()
 *
 *  \param  pair_table  The pair table of the secondary structure
 *  \param  x           a pointer to an array with enough allocated space to hold the x coordinates
 *  \param  y           a pointer to an array with enough allocated space to hold the y coordinates
 *  \param  r           a pointer to an array with enough allocated space to hold the scaling factors
 *  \return             length of sequence on success, 0 otherwise
 */
int simple_circplot_coordinates(short *pair_table, float *x, float *y, float *r);


#endif
