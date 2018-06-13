#ifndef VIENNA_RNA_PACKAGE_STRUCTURE_UTILS_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_STRUCTURE_UTILS_DEPRECATED_H

/**
 *  @file ViennaRNA/structure_utils.h
 *  @brief      Use ViennaRNA/utils/structures.h instead
 *  @deprecated Use ViennaRNA/utils/structures.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/structure_utils.h>! Use <ViennaRNA/utils/structures.h> instead!"
# endif
#include <ViennaRNA/utils/structures.h>
#endif

#endif
