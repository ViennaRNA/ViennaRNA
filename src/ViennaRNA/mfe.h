#ifndef VIENNA_RNA_PACKAGE_MFE_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_MFE_DEPRECATED_H

/**
 *  @file       mfe.h
 *  @brief      Deprecated include file for MFE API
 *  @deprecated Use ViennaRNA/mfe/global.h and ViennaRNA/backtrack/global.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/mfe.h>! Use <ViennaRNA/mfe/global.h> and <ViennaRNA/backtrack/global.h> instead!"
# endif
#include <ViennaRNA/mfe/global.h>
#include <ViennaRNA/backtrack/global.h>
#endif

#endif
