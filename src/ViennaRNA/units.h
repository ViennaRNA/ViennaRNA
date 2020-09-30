#ifndef VIENNA_RNA_PACKAGE_UNITS_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_UNITS_DEPRECATED_H

/**
 *  @file ViennaRNA/units.h
 *  @brief      Use ViennaRNA/utils/units.h instead
 *  @deprecated Use ViennaRNA/utils/units.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/units.h>! Use <ViennaRNA/utils/units.h> instead!"
# endif
#include <ViennaRNA/utils/units.h>
#endif

#endif
