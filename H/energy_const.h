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
#define INF (INT_MAX/10)

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

#endif
