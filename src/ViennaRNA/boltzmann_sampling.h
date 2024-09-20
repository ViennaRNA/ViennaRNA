#ifndef VIENNA_RNA_PACKAGE_SAMPLING_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_SAMPLING_DEPRECATED_H

/**
 *  @file       boltzmann_sampling.h
 *  @brief      Deprecated include file for Boltzmann Sampling API
 *  @deprecated Use ViennaRNA/sampling/basic.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/boltzmann_sampling.h>! Use <ViennaRNA/sampling/basic.h> instead!"
# endif
#include <ViennaRNA/sampling/basic.h>
#endif

#endif
