/*
 * gquad.c
 *
 * Ronny Lorenz 2012
 *
 * ViennaRNA Package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/params/constants.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/alignments.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/loops/gquad.h"

#include "ViennaRNA/loops/gquad_intern.h"


#ifndef INLINE
#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif
#endif


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */


/*
 #########################################
 # BEGIN OF PUBLIC FUNCTION DEFINITIONS  #
 #      (all available in RNAlib)        #
 #########################################
 */

/********************************
 * Here are the G-quadruplex
 * dynamic programming matrices
 ********************************/

/*
 #########################################
 # BEGIN OF PRIVATE FUNCTION DEFINITIONS #
 #          (internal use only)          #
 #########################################
 */

/*
 *###########################################
 *# deprecated functions below              #
 *###########################################
 */


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

#endif
