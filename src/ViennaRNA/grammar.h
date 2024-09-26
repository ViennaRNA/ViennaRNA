#ifndef VIENNA_RNA_PACKAGE_GRAMMAR_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_GRAMMAR_DEPRECATED_H

/**
 *  @file       grammar.h
 *  @brief      Deprecated include file for grammar extension API
 *  @deprecated Use ViennaRNA/grammar/basic.h, ViennaRNA/grammar/mfe.h or ViennaRNA/grammar/partfunc.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/grammar.h>! Use <ViennaRNA/grammar/basic.h>, <ViennaRNA/grammar/mfe.h> or <ViennaRNA/grammar/partfunc.h> instead!"
# endif
#include <ViennaRNA/grammar/basic.h>
#include <ViennaRNA/grammar/mfe.h>
#include <ViennaRNA/grammar/partfunc.h>
#endif

#endif
