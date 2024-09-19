#ifndef VIENNA_RNA_PACKAGE_LOOPS_EXTERNAL_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_LOOPS_EXTERNAL_DEPRECATED_H

/**
 *  @file ViennaRNA/exterior_loops.h
 *  @brief      Use ViennaRNA/loops/external.h instead
 *  @deprecated Use ViennaRNA/eval/external.h or ViennaRNA/eval/external.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/exterior_loops.h>! Use <ViennaRNA/eval/exterior.h> or <ViennaRNA/mfe/exterior.h> instead!"
# endif
#include <ViennaRNA/loops/external.h>
#endif

#endif
