#ifndef VIENNA_RNA_PACKAGE_INVERSE_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_INVERSE_DEPRECATED_H

/**
 *  @file       inverse.h
 *  @brief      Deprecated include file for inverse folding API
 *  @deprecated Use ViennaRNA/inverse/basic.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/inverse.h>! Use <ViennaRNA/inverse/basic.h> instead!"
# endif
#include <ViennaRNA/inverse/basic.h>
#endif

#endif
