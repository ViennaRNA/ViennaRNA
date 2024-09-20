#ifndef VIENNA_RNA_PACKAGE_EQUILIBRIUM_PROBS_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_EQUILIBRIUM_PROBS_DEPRECATED_H

/**
 *  @file       equilibrium_probs.h
 *  @brief      Deprecated include file for equilibrium probabilities API
 *  @deprecated Use ViennaRNA/probabilities/basepairs.h and ViennaRNA/probabilities/structures.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/equilibrium_probs.h>! Use <ViennaRNA/probabilities/basepairs.h> and <ViennaRNA/probabilities/structures.h> instead!"
# endif
#include <ViennaRNA/probabilities/basepairs.h>
#include <ViennaRNA/probabilities/structures.h>
#endif

#endif
