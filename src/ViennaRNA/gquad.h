#ifndef VIENNA_RNA_PACKAGE_GQUAD_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_GQUAD_DEPRECATED_H

/**
 *  @file gquad.h
 *  @brief      Deprecated include file for G-quadruplex declarations
 *  @deprecated Use ViennaRNA/loops/gquad.h and ViennaRNA/loops/gquad.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/gquad.h>! Use <ViennaRNA/loops/gquad.h> and <ViennaRNA/backtrack/gquad.h> instead!"
# endif
#include <ViennaRNA/loops/gquad.h>
#include <ViennaRNA/backtrack/gquad.h>
#endif

#endif
