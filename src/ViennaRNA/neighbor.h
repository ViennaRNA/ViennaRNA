#ifndef VIENNA_RNA_PACKAGE_NEIGHBOR_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_NEIGHBOR_DEPRECATED_H

/**
 *  @file ViennaRNA/neighbor.h
 *  @brief      Use ViennaRNA/landscape/neighbor.h instead
 *  @deprecated Use ViennaRNA/landscape/neighbor.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/neighbor.h>! Use <ViennaRNA/landscape/neighbor.h> instead!"
# endif
#include <ViennaRNA/landscape/neighbor.h>
#include <ViennaRNA/landscape/move.h>
#endif

#endif
