#ifndef __VIENNA_RNA_PACKAGE_PS_DOT_H__
#define __VIENNA_RNA_PACKAGE_PS_DOT_H__

#include "data_structures.h"

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif

/** \file PS_dot.h **/

/* routines from PS_dot.c */
int PS_rna_plot(char *string, char *structure, char *file);
/* write PostScript drawing of structure to file */
int PS_rna_plot_a(char *string, char *structure, char *file, char *pre, char *post);
/* write PostScript drawing of structure to file with annotation */
int gmlRNA(char *string, char *structure, char *ssfile, char option);
/* structure drawing in gml */
int ssv_rna_plot(char *string, char *structure, char *ssfile);
/*write coord file for SStructView */
int svg_rna_plot(char *string, char *structure, char *ssfile);
/*write RNAplot in SVG */
int xrna_plot(char *string, char *structure, char *ssfile);
/*write .ss file for further editing in XRNA */

/**
*** Wrapper to PS_dot_plot_list
***
*** \note DO NOT USE THIS FUNCTION ANYMORE AS IT IS NOT THREADSAFE
***
*** \deprecated This function is deprecated and will be removed soon! Use \ref PS_dot_plot_list() instead!
**/
int DEPRECATED(PS_dot_plot(char *string, char *file));

extern int rna_plot_type;   /* 0= simple coordinates, 1= naview */

int PS_color_dot_plot(char *string, cpair *pi, char *filename);
int PS_color_dot_plot_turn(char *seq, cpair *pi, char *filename, int winSize);

int PS_dot_plot_list(char *seq, char *filename, struct plist *pl,
         struct plist *mf, char *comment);
int PS_dot_plot_turn(char *seq, struct plist *pl, char *filename,
         int winSize);
int PS_color_aln(const char *structure, const char *filename,
			const char *seqs[], const char *names[]);

int simple_xy_coordinates(short *pair_table, float *X, float *Y);

#endif
