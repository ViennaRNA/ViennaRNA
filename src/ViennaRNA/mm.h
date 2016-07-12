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


unsigned int  maximumMatching(const char *string);

unsigned int *maximumMatchingConstraint(const char *string,
                                        short *ptable);

unsigned int *maximumMatching2Constraint( const char *string,
                                          short *ptable,
                                          short *ptable2);

#endif
