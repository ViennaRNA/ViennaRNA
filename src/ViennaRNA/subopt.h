#ifndef VIENNA_RNA_PACKAGE_SUBOPT_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_SUBOPT_DEPRECATED_H

/**
 *  @file       subopt.h
 *  @brief      Deprecated include file for subopt API
 *  @deprecated Use ViennaRNA/subopt/basic.h and ViennaRNA/subopt/wuchty.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/subopt.h>! Use <ViennaRNA/subopt/basic.h> and <ViennaRNA/subopt/wuchty.h> instead!"
# endif
#include <ViennaRNA/subopt/wuchty.h>
#endif

#endif
