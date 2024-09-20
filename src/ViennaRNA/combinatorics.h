#ifndef VIENNA_RNA_PACKAGE_COMBINATORICS_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_COMBINATORICS_DEPRECATED_H

/**
 *  @file       combinatorics.h
 *  @brief      Deprecated include file for combinatorics utils API
 *  @deprecated Use ViennaRNA/combinatorics/basic.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/combinatorics.h>! Use <ViennaRNA/combinatorics/basic.h> instead!"
# endif
#include <ViennaRNA/combinatorics/basic.h>
#endif

#endif
