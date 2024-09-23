#ifndef VIENNA_RNA_PACKAGE_UTILS_ALIGNMENTS_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_UTILS_ALIGNMENTS_DEPRECATED_H

/**
 *  @file       ViennaRNA/utils/alignments.h
 *  @brief      Deprecated include file for multiple sequence alignment utils API
 *  @deprecated Use ViennaRNA/sequences/alignments.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/utils/alignments.h>! Use <ViennaRNA/sequences/alignments.h> instead!"
# endif
#include <ViennaRNA/sequences/alignments.h>
#endif

#endif
