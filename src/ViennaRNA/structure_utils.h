#ifndef VIENNA_RNA_PACKAGE_STRUCTURE_UTILS_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_STRUCTURE_UTILS_DEPRECATED_H

/**
 *  @file ViennaRNA/structure_utils.h
 *  @brief      Use ViennaRNA/utils/structures.h instead
 *  @deprecated Use ViennaRNA/utils/structures.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/structure_utils.h>! Use one of the headers in ViennaRNA/structures/ instead!"
# endif
#include <ViennaRNA/structures/benchmark.h>
#include <ViennaRNA/structures/dotbracket.h>
#include <ViennaRNA/structures/helix.h>
#include <ViennaRNA/structures/metrics.h>
#include <ViennaRNA/structures/pairtable.h>
#include <ViennaRNA/structures/problist.h>
#include <ViennaRNA/structures/shapes.h>
#include <ViennaRNA/structures/tree.h>
#include <ViennaRNA/structures/utils.h>
#endif

#endif
