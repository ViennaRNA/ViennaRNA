#ifndef VIENNA_RNA_PACKAGE_COMMANDS_DEPRECATED_H
#define VIENNA_RNA_PACKAGE_COMMANDS_DEPRECATED_H

/**
 *  @file       commands.h
 *  @brief      Deprecated include file for commands API
 *  @deprecated Use ViennaRNA/io/commands.h instead
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
# ifdef VRNA_WARN_DEPRECATED
#warning "Including deprecated header file <ViennaRNA/commands.h>! Use <ViennaRNA/io/commands.h> instead!"
# endif
#include <ViennaRNA/io/commands.h>
#endif

#endif
