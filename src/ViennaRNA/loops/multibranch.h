#ifndef VIENNA_RNA_PACKAGE_LOOPS_MULTIBRANCH_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_LOOPS_MULTIBRANCH_DEPRECATED_H

/**
 *  @file       ViennaRNA/loops/multibranch.h
 *  @brief      Deprecated include file for multibranch loop energy evaluation function declarations
 *  @deprecated Use ViennaRNA/eval/multibranch.h, ViennaRNA/mfe/multibranch.h, ViennaRNA/backtrack/multibranch.h, and ViennaRNA/partfunc/multibranch.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/loops/multibranch.h>! Use <ViennaRNA/eval/multibranch.h>, <ViennaRNA/mfe/multibranch.h>, <ViennaRNA/backtrack/multibranch.h>, and <ViennaRNA/partfunc/multibranch.h> instead!"
# endif
#include <ViennaRNA/eval/multibranch.h>
#include <ViennaRNA/mfe/multibranch.h>
#include <ViennaRNA/backtrack/multibranch.h>
#include <ViennaRNA/partfunc/multibranch.h>
#endif

#endif
