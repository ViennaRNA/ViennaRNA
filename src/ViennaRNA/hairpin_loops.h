#ifndef VIENNA_RNA_PACKAGE_LOOPS_HAIRPIN_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_LOOPS_HAIRPIN_DEPRECATED_H

/**
 *  @file ViennaRNA/hairpin_loops.h
 *  @brief      Use ViennaRNA/loops/hairpin.h instead
 *  @deprecated Use ViennaRNA/loops/hairpin.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/hairpin_loops.h>! Use <ViennaRNA/loops/hairpin.h> instead!"
# endif
#include <ViennaRNA/loops/hairpin.h>
#endif

#endif
