/*
  Last changed Time-stamp: <2005-02-16 10:10:56 ivo>
  c  Christoph Flamm and Ivo L Hofacker
  {xtof,ivo}@tbi.univie.ac.at
  Kinfold: $Name:  $
  $Id: nachbar.h,v 1.2 2005/02/16 17:00:48 ivo Exp $
*/

#ifndef NACHBAR_H
#define NACHBAR_H

/* used in baum.c */
extern void ini_nbList(int chords);
extern void update_nbList(int i,int j, int iE);

/* used in main.c */
extern int sel_nb(void);
extern void clean_up_nbList(void);

extern void grow_chain(void);
#endif
