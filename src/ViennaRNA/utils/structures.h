#ifndef VIENNA_RNA_PACKAGE_UTILS_STRUCTURES_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_UTILS_STRUCTURES_DEPRECATED_H

/**
 *  @file       ViennaRNA/utils/structures.h
 *  @brief      Deprecated include file for structure utilities
 *  @deprecated Use an appropriate header from ViennaRNA/structures/ instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/utils/structures.h>! Use one of the headers in ViennaRNA/structures/ instead!"
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
