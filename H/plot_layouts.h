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
 *  \note To provide thread safety please do not rely on this global variable in future implementations
 *  but pass a plot type flag directly to the function that decides which layout algorithm it may use!
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
 *  \param  X           a pointer to an array with enough allocated space to hold the x coordinates
 *  \param  Y           a pointer to an array with enough allocated space to hold the y coordinates
 *  \return             length of sequence on success, 0 otherwise
 */
int simple_xy_coordinates(short *pair_table,
                          float *X,
                          float *Y);

/**
 *  \brief Calculate nucleotide coordinates for <i>Circular Plot</i>
 *
 *  This function calculates the coordinates of nucleotides mapped in equal distancies onto a unit circle.
 *
 *  \note In order to draw nice arcs using quadratic bezier curves that connect base pairs one may calculate
 *  a second tangential point \f$P^t\f$ in addition to the actual R<sup>2</sup> coordinates.
 *  the simplest way to do so may be to compute a radius scaling factor \f$rs\f$ in the interval \f$[0,1]\f$ that
 *  weights the proportion of base pair span to the actual length of the sequence. This scaling factor
 *  can then be used to calculate the coordinates for \f$P^t\f$, i.e. \f$ P^{t}_x[i] = X[i] * rs\f$
 *  and \f$P^{t}_y[i] = Y[i] * rs\f$.
 *
 *  \see make_pair_table(), rna_plot_type, simple_xy_coordinates(), naview_xy_coordinates(), PS_rna_plot_a(),
 *  PS_rna_plot, svg_rna_plot()
 *
 *  \param  pair_table  The pair table of the secondary structure
 *  \param  x           a pointer to an array with enough allocated space to hold the x coordinates
 *  \param  y           a pointer to an array with enough allocated space to hold the y coordinates
 *  \return             length of sequence on success, 0 otherwise
 */
int simple_circplot_coordinates(short *pair_table,
                                float *x,
                                float *y);


#endif
