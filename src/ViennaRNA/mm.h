#ifndef VIENNA_RNA_PACKAGE_MM_H
#define VIENNA_RNA_PACKAGE_MM_H

/**
 *
 *  @file mm.h
 *  @ingroup subopt_and_representatives
 *  @brief Several Maximum Matching implementations
 *
 *  This file contains the declarations for several maximum matching implementations
 */

#include <ViennaRNA/fold_compound.h>

int
vrna_maximum_matching(vrna_fold_compound_t *fc);


int
vrna_maximum_matching_simple(const char *sequence);


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

unsigned int
maximumMatching(const char *string);


unsigned int *
maximumMatchingConstraint(const char  *string,
                          short       *ptable);


unsigned int *
maximumMatching2Constraint(const char *string,
                           short      *ptable,
                           short      *ptable2);


#endif

#endif
