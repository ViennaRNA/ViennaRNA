#ifndef VIENNA_RNA_PACKAGE_FILE_FORMATS_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_FILE_FORMATS_DEPRECATED_H

/**
 *  @file ViennaRNA/file_formats.h
 *  @brief      Use ViennaRNA/io/file_formats.h instead
 *  @deprecated Use ViennaRNA/io/file_formats.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/file_formats.h>! Use <ViennaRNA/io/file_formats.h> instead!"
# endif
#include <ViennaRNA/io/file_formats.h>
#include <ViennaRNA/io/file_formats_msa.h>
#endif

#endif
