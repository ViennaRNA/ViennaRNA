#ifndef VIENNA_RNA_PACKAGE_LOOPS_EXTLOOP_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_LOOPS_EXTLOOP_DEPRECATED_H

/**
 *  @file       ViennaRNA/loops/exterior.h
 *  @brief      Deprecated include file for exterior loop energy evaluation function declarations
 *  @deprecated Use ViennaRNA/eval/exterior.h, ViennaRNA/mfe/exterior.h, ViennaRNA/backtrack/exterior.h, and ViennaRNA/partfunc/exterior.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/loops/exterior.h>! Use <ViennaRNA/eval/exterior.h>, <ViennaRNA/mfe/exterior.h>, <ViennaRNA/backtrack/exterior.h>, and <ViennaRNA/partfunc/exterior.h> instead!"
# endif
#include <ViennaRNA/eval/exterior.h>
#include <ViennaRNA/mfe/exterior.h>
#include <ViennaRNA/backtrack/exterior.h>
#include <ViennaRNA/partfunc/exterior.h>
#endif

#endif
