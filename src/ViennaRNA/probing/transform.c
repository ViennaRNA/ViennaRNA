/* SHAPE reactivity transform */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/params/default.h"
#include "ViennaRNA/params/constants.h" /* defines MINPSCORE */
#include "ViennaRNA/datastructures/array.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/io/utils.h"
#include "ViennaRNA/io/file_formats.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/probing/basic.h"

#include "ViennaRNA/probing/transform.h"


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
 * Transform reactivity values
 */
PRIVATE FLT_OR_DBL *
reactivity_transform(unsigned int n,
                     const double *reactivity,
                     vrna_probing_transform_f trans,
                     void *options);


/*
 * Keep the given reactivity value
 */
PRIVATE double
reactivity_trans_identity(double r, void *options);


/*
 * Set negative reactivity as missing
 */
PRIVATE double
reactivity_trans_neg_ignore(double r, void *options);


/*
 * Set negative reactivity as zero
 */
PRIVATE double
reactivity_trans_neg_to_zero(double r, void *options);


/*
 * Log(1+x) for x > -1, ignore otherwise
 */
PRIVATE double
reactivity_trans_log1p(double r, void *options);

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */

PUBLIC vrna_probing_transform_f
vrna_reactivity_trans_default(unsigned int flag)
{
  switch (flag) {
    case VRNA_PROBING_METHOD_EDDY2014_2:
      return reactivity_trans_log1p;
  }
  return reactivity_trans_neg_ignore;
}


PUBLIC vrna_probing_transform_f
vrna_reactivity_trans_method(unsigned int flag)
{
  switch (flag) {
    case VRNA_REACTIVITY_TRANS_DEFAULT:
      return NULL;
    case VRNA_REACTIVITY_TRANS_IDEN:
      return reactivity_trans_identity;
    case VRNA_REACTIVITY_TRANS_NEG_IGNORE:
      return reactivity_trans_neg_ignore;
    case VRNA_REACTIVITY_TRANS_NEG_ZERO:
      return reactivity_trans_neg_to_zero;
    case VRNA_REACTIVITY_TRANS_LOG1P:
      return reactivity_trans_log1p;
  }
  return NULL;
}


PUBLIC double *
vrna_reactivity_transform(unsigned int n,
                          const double *reactivity,
                          vrna_probing_transform_f trans,
                          void *options)
{
  /* init the transformed reactivity array */
  double *a = NULL;

  if ((reactivity) &&
      (trans)) {
    a     = (double *)vrna_alloc(sizeof(double) * (n + 1));
    a[0]  = VRNA_REACTIVITY_MISSING;

    for (size_t i = 1; i <= n; ++i)
      a[i] = (reactivity[i] == VRNA_REACTIVITY_MISSING) ? VRNA_REACTIVITY_MISSING : trans(reactivity[i], options);
  }

  return a;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE double
reactivity_trans_identity(double r, void *options)
{
  return r;
}


PRIVATE double
reactivity_trans_neg_ignore(double r, void *options)
{
  return r < 0 ? VRNA_REACTIVITY_MISSING : r;
}


PRIVATE double
reactivity_trans_neg_to_zero(double r, void *options)
{
  return  r < 0 ? 0 : r;
}

PRIVATE double
reactivity_trans_log1p(double r, void *options){
  return r <= -1 ? VRNA_REACTIVITY_MISSING : log(1+r);
}
