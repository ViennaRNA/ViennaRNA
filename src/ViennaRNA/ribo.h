#ifndef VIENNA_RNA_PACKAGE_RIBOSUM_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_RIBOSUM_DEPRECATED_H

/**
 *  @file       ribo.h
 *  @brief      Deprecated include file for ribosum scores
 *  @deprecated Use ViennaRNA/params/ribosum.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/ribo.h>! Use <ViennaRNA/params/ribosum.h> instead!"
# endif
#include <ViennaRNA/params/ribosum.h>
#endif

#endif
