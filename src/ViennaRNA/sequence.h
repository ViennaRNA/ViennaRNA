#ifndef VIENNA_RNA_PACKAGE_SEQUENCE_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_SEQUENCE_DEPRECATED_H

/**
 *  @file       sequence.h
 *  @brief      Deprecated include file for sequence API
 *  @deprecated Use ViennaRNA/sequences/sequence.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/sequence.h>! Use <ViennaRNA/sequences/sequence.h> instead!"
# endif
#include <ViennaRNA/sequences/sequence.h>
#endif

#endif
