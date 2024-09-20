#ifndef VIENNA_RNA_PACKAGE_PF_MULTIFOLD_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_PF_MULTIFOLD_DEPRECATED_H

/**
 *  @file       pf_multifold.h
 *  @brief      Deprecated include file for multi-strand partition function API
 *  @deprecated Use ViennaRNA/partfunc/multifold.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/pf_multifold.h>! Use <ViennaRNA/partfunc/multifold.h> and <ViennaRNA/backtrack/global.h> instead!"
# endif
#include <ViennaRNA/partfunc/multifold.h>
#endif

#endif
