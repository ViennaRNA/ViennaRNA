#ifndef VIENNA_RNA_PACKAGE_PLOT_NAVIEW_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_PLOT_NAVIEW_DEPRECATED_H

/**
 *  @file ViennaRNA/naview.h
 *  @brief      Use ViennaRNA/plotting/naview/naview.h instead
 *  @deprecated Use ViennaRNA/plotting/naview/naview.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#   warning "Including deprecated header file <ViennaRNA/naview.h>! Use <ViennaRNA/plotting/naview.h> instead!"
# endif
# ifdef VRNA_WITH_NAVIEW_LAYOUT
#   include <ViennaRNA/plotting/naview/naview.h>
# else
#   warning "Naview Layout algorithm is not available in this version!"
# endif
#endif

#endif
