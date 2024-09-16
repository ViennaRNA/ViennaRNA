#ifndef VIENNA_RNA_PACKAGE_EVAL_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_EVAL_DEPRECATED_H

/**
 *  @file eval.h
 *  @brief      Deprecated include file for structure energy evaluation function declarations
 *  @deprecated Use ViennaRNA/eval/structures.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/eval.h>! Use <ViennaRNA/eval/structures.h> instead!"
# endif
#include <ViennaRNA/eval/structures.h>
#endif

#endif
