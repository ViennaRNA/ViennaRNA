#ifndef VIENNA_RNA_PACKAGE_PKPLEX_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_PKPLEX_DEPRECATED_H

/**
 *  @file ViennaRNA/utils.h
 *  @brief      Use ViennaRNA/utils/basic.h instead
 *  @deprecated Use ViennaRNA/utils/basic.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/PKplex.h>! Use <ViennaRNA/pk_plex.h> instead!"
# endif

#ifdef VRNA_WARN_DEPRECATED
# if defined(__clang__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated("", msg)))
# elif defined(__GNUC__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated(msg)))
# else
#  define DEPRECATED(func, msg) func
# endif
#else
# define DEPRECATED(func, msg) func
#endif

#include <ViennaRNA/datastructures/basic.h>


DEPRECATED(dupVar *
PKLduplexfold_XS(const char *s1,
                 const int  **access_s1,
                 int  penalty,
                 int  max_interaction_length,
                 int  delta),
          "Use vrna_pk_plex() instead!");

#include <ViennaRNA/pk_plex.h>

#endif

#endif
