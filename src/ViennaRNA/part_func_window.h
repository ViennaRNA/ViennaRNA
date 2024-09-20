#ifndef VIENNA_RNA_PACKAGE_PF_WINDOW_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_PF_WINDOW_DEPRECATED_H

/**
 *  @file       part_func_window.h
 *  @brief      Deprecated include file for sliding-window partition function API
 *  @deprecated Use ViennaRNA/partfunc/local.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/part_func_window.h>! Use <ViennaRNA/partfunc/local.h> and <ViennaRNA/backtrack/global.h> instead!"
# endif
#include <ViennaRNA/partfunc/local.h>
#endif

#endif
