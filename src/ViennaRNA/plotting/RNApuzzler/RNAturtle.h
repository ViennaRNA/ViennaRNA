#ifndef RNATURTLE_H
#define RNATURTLE_H

/**
 * Compute layout using RNAturtle algorithm
 */
int
vrna_plot_coords_turtle(const char  *structure,
                        float       **x,
                        float       **y,
                        double      **arc_coords);


int
vrna_plot_coords_turtle_pt(short const *const pair_table,
                           float              **x,
                           float              **y,
                           double             **arc_coords);


#endif
