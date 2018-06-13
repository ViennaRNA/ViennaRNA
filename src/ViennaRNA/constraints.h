#ifndef VIENNA_RNA_PACKAGE_CONSTRAINTS_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_CONSTRAINTS_DEPRECATED_H

/**
 *  @file constraints.h
 *  @brief      Use ViennaRNA/constraints/basic.h instead
 *  @deprecated Use ViennaRNA/constraints/basic.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/constraints.h>! Use <ViennaRNA/constraints/basic.h> instead!"
# endif
#include <ViennaRNA/constraints/basic.h>
#include <ViennaRNA/constraints/hard.h>
#include <ViennaRNA/constraints/soft.h>
#include <ViennaRNA/constraints/SHAPE.h>
#include <ViennaRNA/constraints/ligand.h>
#endif

#endif
