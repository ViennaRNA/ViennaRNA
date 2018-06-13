#ifndef VIENNA_RNA_PACKAGE_PARAMS_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_PARAMS_DEPRECATED_H

/**
 *  @file ViennaRNA/params.h
 *  @brief      Use ViennaRNA/params/basic.h instead
 *  @deprecated Use ViennaRNA/params/basic.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/params.h>! Use <ViennaRNA/params/basic.h> instead!"
# endif
#include <ViennaRNA/params/basic.h>
#endif

#endif
