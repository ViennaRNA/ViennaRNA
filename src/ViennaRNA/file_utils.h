#ifndef VIENNA_RNA_PACKAGE_FILE_UTILS_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_FILE_UTILS_DEPRECATED_H

/**
 *  @file ViennaRNA/file_utils.h
 *  @brief      Use ViennaRNA/io/utils.h instead
 *  @deprecated Use ViennaRNA/io/utils.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/file_utils.h>! Use <ViennaRNA/io/utils.h> instead!"
# endif
#include <ViennaRNA/io/utils.h>
#endif

#endif
