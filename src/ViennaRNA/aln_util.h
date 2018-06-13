#ifndef VIENNA_RNA_PACKAGE_ALN_UTILS_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_ALN_UTILS_DEPRECATED_H

/**
 *  @file ViennaRNA/aln_util.h
 *  @brief      Use ViennaRNA/utils/alignments.h instead
 *  @deprecated Use ViennaRNA/utils/alignments.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/aln_util.h>! Use <ViennaRNA/utils/alignments.h> instead!"
# endif
#include <ViennaRNA/utils/alignments.h>
#endif

#endif
