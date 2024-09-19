#ifndef VIENNA_RNA_PACKAGE_LOOPS_INTLOOP_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_LOOPS_INTLOOP_DEPRECATED_H

/**
 *  @file       ViennaRNA/loops/internal.h
 *  @brief      Deprecated include file for internal loop energy evaluation function declarations
 *  @deprecated Use ViennaRNA/eval/internal.h, ViennaRNA/mfe/internal.h, ViennaRNA/backtrack/internal.h, and ViennaRNA/partfunc/internal.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/loops/internal.h>! Use <ViennaRNA/eval/internal.h>, <ViennaRNA/mfe/internal.h>, <ViennaRNA/backtrack/internal.h>, and <ViennaRNA/partfunc/internal.h> instead!"
# endif
#include <ViennaRNA/eval/internal.h>
#include <ViennaRNA/mfe/internal.h>
#include <ViennaRNA/backtrack/internal.h>
#include <ViennaRNA/partfunc/internal.h>
#endif

#endif
