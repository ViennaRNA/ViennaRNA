#ifndef VIENNA_RNA_PACKAGE_LOOPS_HAIRPIN_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_LOOPS_HAIRPIN_DEPRECATED_H

/**
 *  @file       ViennaRNA/loops/hairpin.h
 *  @brief      Deprecated include file for hairpin energy evaluation function declarations
 *  @deprecated Use ViennaRNA/eval/hairpin.h and ViennaRNA/backtrack/hairpin.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/loops/hairpin.h>! Use <ViennaRNA/eval/hairpin.h> and <ViennaRNA/backtrack/hairpin.h> instead!"
# endif
#include <ViennaRNA/eval/hairpin.h>
#include <ViennaRNA/backtrack/hairpin.h>
#endif

#endif
