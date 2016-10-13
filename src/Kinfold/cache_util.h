/*
  Last changed Time-stamp: <2006-10-03 09:37:34 xtof>
  c  Christoph Flamm and Ivo L Hofacker
  {xtof,ivo}@tbi.univie.ac.at
  Kinfold: $Name:  $
  $Id: cache_util.h,v 1.2 2006/10/04 12:45:12 xtof Exp $
*/

#ifndef CACHE_UTIL_H
#define CACHE_UTIL_H

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

typedef struct {
  char *structure;
  int top;           /* number of neighbors */
  int lmin;          /* is a local minimum ? */
  double flux;       /* sum of rates */
  double energy;     /* energy of this structure */
  short *neighbors;  
  float *rates;
  double *energies;
} cache_entry;

extern cache_entry *lookup_cache (char *x);
extern int write_cache (cache_entry *x);
void kill_cache(void);

#endif
