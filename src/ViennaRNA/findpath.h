#ifndef VIENNA_RNA_PACKAGE_FINDPATH_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_FINDPATH_DEPRECATED_H

/**
 *  @file ViennaRNA/findpath.h
 *  @brief      Use ViennaRNA/landscape/findpath.h instead
 *  @deprecated Use ViennaRNA/landscape/findpath.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/findpath.h>! Use <ViennaRNA/landscape/findpath.h> instead!"
# endif
#include <ViennaRNA/landscape/findpath.h>
#endif

#endif
