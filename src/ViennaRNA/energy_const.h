#ifndef __VIENNA_RNA_PACKAGE_ENERGY_CONST_H__
#define __VIENNA_RNA_PACKAGE_ENERGY_CONST_H__

#include <limits.h>

/**
 *  \file energy_const.h
 *  energy constants
 */

/** The gas constant */
#define GASCONST 1.98717  /* in [cal/K] */
/** 0 deg Celsius in Kelvin */
#define K0  273.15
/** Infinity as used in minimization routines */
#define INF 10000000 /* (INT_MAX/10) */

#define EMAX (INF/10)
/** forbidden */
#define FORBIDDEN 9999
/** bonus contribution */
#define BONUS 10000
/** The number of distinguishable base pairs */
#define NBPAIRS 7
/** The minimum loop length */
#define TURN 3
/** The maximum loop length */
#define MAXLOOP 30

#define   VRNA_GQUAD_MISMATCH_PENALTY   300   /* penalty for incompatible nucleotides in an alignment that destruct a gquad layer */
#define   VRNA_GQUAD_MISMATCH_NUM_ALI   1   /* maximum number of mismatching sequences in the alignment when gquad should be formed */

#endif
