#ifndef VIENNA_RNA_PACKAGE_WALK_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_WALK_DEPRECATED_H

/**
 *  @file ViennaRNA/walk.h
 *  @brief      Use ViennaRNA/landscape/walk.h instead
 *  @deprecated Use ViennaRNA/landscape/walk.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/walk.h>! Use <ViennaRNA/landscape/walk.h> instead!"
# endif
#include <ViennaRNA/landscape/walk.h>
#include <ViennaRNA/landscape/neighbor.h>
#endif

#endif
