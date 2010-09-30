/*
      minimum free energy
      RNA secondary structure with
      basepair distance d to reference structure prediction
      
*/
#ifndef __VIENNA_RNA_PACKAGE_TWO_D_FOLD_H__
#define __VIENNA_RNA_PACKAGE_TWO_D_FOLD_H__

/**
*** \file 2Dfold.h
*** \brief Compute the minimum free energy (MFE) for secondary structures with a 
*** basepair distance to two fixed reference structures
**/

#include "data_structures.h"

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif

/**
*** \brief Get a structure of type TwoDfold_vars prefilled with current global settings
***
*** This function returns a datastructure of type TwoDfold_vars.
*** The data fields inside the TwoDfold_vars are prefilled by global settings and all memory
*** allocations necessary to start a computation are already done for the convenience of the user
***
*** \note Make sure that the reference structures are compatible with the sequence according to Watson-Crick- and Wobble-base pairing
***
*** \see destroy_TwoDfold_variables(), TwoDfold(), TwoDfold_circ
***
*** \param seq          The RNA sequence
*** \param structure1   The first reference structure in dot-bracket notation
*** \param structure2   The second reference structure in dot-bracket notation
*** \param circ         A switch to indicate the assumption to fold a circular instead of linear RNA (0=OFF, 1=ON)
*** \returns            A datastructure prefilled with folding options and allocated memory
**/
TwoDfold_vars     *get_TwoDfold_variables(const char *seq, const char *structure1, const char *structure2, int circ);

/**
*** \brief Destroy a TwoDfold_vars datastructure without memory loss
***
*** This function free's all allocated memory that depends on the datastructure given.
***
*** \see get_TwoDfold_variables()
***
*** \param our_variables  A pointer to the datastructure to be destroyed 
**/
void              destroy_TwoDfold_variables(TwoDfold_vars *our_variables);

/**
***
**/
TwoDfold_solution **TwoDfold(const char *string, char *structure1, char *structure2);
/**
***
**/
TwoDfold_solution **TwoDfold_bound(TwoDfold_vars *our_variables, int distance1, int distance2);

/**
***
**/
TwoDfold_solution **TwoDfold_circ(const char *string, char *structure1, char *structure2);
/**
***
**/
TwoDfold_solution **TwoDfold_circ_bound(TwoDfold_vars *our_variables, int maxDistance1, int maxDistance2);

/* function for Pathfinder.c */
/**
***
**/
char *TwoDfold_backtrack_f5(unsigned int j, unsigned int k, unsigned int l, TwoDfold_vars *vars);

/* deprecated */
DEPRECATED(TwoDfold_solution **TwoDfold_bounded(const char *string, char *structure1, char *structure2, int maxDistance1, int maxDistance2));
/* deprecated */
DEPRECATED(TwoDfold_solution **TwoDfold_circ_bounded(const char *string, char *structure1, char *structure2, int maxDistance1, int maxDistance2));

#endif
