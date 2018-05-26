/** \file data_structures.c **/

/*
 *                Data structure creation/destruction
 *
 *                This file contains everything which is necessary to
 *                obtain and destroy datastructures used in the folding
 *                recurrences throughout the ViennaRNA package
 *
 *                c Ronny Lorenz
 *
 *                ViennaRNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "ViennaRNA/utils/basic.h"

/*
 #################################
 # PRIVATE MACROS                #
 #################################
 */

/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */

#ifndef VRNA_DISABLE_C11_FEATURES
PUBLIC void
vrna_C11_features(void)
{
  __asm("nop");
}


#endif
