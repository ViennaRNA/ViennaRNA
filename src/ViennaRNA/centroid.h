#ifndef VIENNA_RNA_PACKAGE_CENTROID_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_CENTROID_DEPRECATED_H

/**
 *  @file       centroid.h
 *  @brief      Deprecated include file for centroid API
 *  @deprecated Use ViennaRNA/structures/centroid.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/centroid.h>! Use <ViennaRNA/structures/centroid.h> instead!"
# endif
#include <ViennaRNA/structures/centroid.h>
#endif

#endif
