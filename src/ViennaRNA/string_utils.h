#ifndef VIENNA_RNA_PACKAGE_STRING_UTILS_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_STRING_UTILS_DEPRECATED_H

/**
 *  @file ViennaRNA/string_utils.h
 *  @brief      Use ViennaRNA/utils/strings.h instead
 *  @deprecated Use ViennaRNA/utils/strings.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/string_utils.h>! Use <ViennaRNA/utils/strings.h> instead!"
# endif
#include <ViennaRNA/utils/strings.h>
#endif

#endif
