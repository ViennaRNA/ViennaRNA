#ifndef VIENNA_RNA_PACKAGE_PARTFUNC_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_PARTFUNC_DEPRECATED_H

/**
 *  @file       part_func.h
 *  @brief      Deprecated include file for partition function API
 *  @deprecated Use ViennaRNA/partfunc/global.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/part_func.h>! Use <ViennaRNA/partfunc/global.h> instead!"
# endif
#include <ViennaRNA/partfunc/global.h>
#include <ViennaRNA/probabilities/basepairs.h>
#include <ViennaRNA/sampling/basic.h>
#endif

#endif
