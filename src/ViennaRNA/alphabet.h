#ifndef VIENNA_RNA_PACKAGE_ALPHABET_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_ALPHABET_DEPRECATED_H

/**
 *  @file       alphabet.h
 *  @brief      Deprecated include file for sequence alphabet API
 *  @deprecated Use ViennaRNA/sequences/alphabet.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/alphabet.h>! Use <ViennaRNA/sequences/alphabet.h> instead!"
# endif
#include <ViennaRNA/sequences/alphabet.h>
#endif

#endif
