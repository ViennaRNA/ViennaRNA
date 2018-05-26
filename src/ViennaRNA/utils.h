#ifndef VIENNA_RNA_PACKAGE_UTILS_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_UTILS_DEPRECATED_H

/**
 *  @file ViennaRNA/utils.h
 *  @brief      Use ViennaRNA/utils/basic.h instead
 *  @deprecated Use ViennaRNA/utils/basic.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/utils.h>! Use <ViennaRNA/utils/basic.h> instead!"
# endif
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/utils/strings.h>
#include <ViennaRNA/utils/structures.h>
#include <ViennaRNA/io/utils.h>
#include <ViennaRNA/alphabet.h>
#endif

#endif
