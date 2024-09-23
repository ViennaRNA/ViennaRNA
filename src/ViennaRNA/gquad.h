#ifndef VIENNA_RNA_PACKAGE_GQUAD_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_GQUAD_DEPRECATED_H

/**
 *  @file gquad.h
 *  @brief      Deprecated include file for G-quadruplex declarations
 *  @deprecated Use ViennaRNA/eval/gquad.h and other header files instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/gquad.h>! Use <ViennaRNA/eval/gquad.h> and other header files instead!"
# endif
#include <ViennaRNA/eval/gquad.h>
#include <ViennaRNA/mfe/gquad.h>
#include <ViennaRNA/backtrack/gquad.h>
#include <ViennaRNA/partfunc/gquad.h>
#endif

#endif
