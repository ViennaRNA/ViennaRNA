#ifndef VIENNA_RNA_PACKAGE_DATA_STRUCTURES_LIST_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_DATA_STRUCTURES_LIST_DEPRECATED_H

/**
 *  @file ViennaRNA/list.h
 *  @brief      Use ViennaRNA/datastructures/lists.h instead
 *  @deprecated Use ViennaRNA/datastructures/lists.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/list.h>! Use <ViennaRNA/datastructures/lists.h> instead!"
# endif
#include <ViennaRNA/datastructures/lists.h>
#endif

#endif
