#ifndef VIENNA_RNA_PACKAGE_MEA_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_MEA_DEPRECATED_H

/**
 *  @file       MEA.h
 *  @brief      Deprecated include file for MEA API
 *  @deprecated Use ViennaRNA/structures/mea.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/MEA.h>! Use <ViennaRNA/structures/mea.h> instead!"
# endif
#include <ViennaRNA/structures/mea.h>
#endif

#endif
