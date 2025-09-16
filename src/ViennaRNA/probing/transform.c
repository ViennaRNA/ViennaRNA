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


typedef struct {
  double slope;
  double intersect;
  double oob;
} linear_transform_param_t;

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

PRIVATE double
transform_linear(double value,
                 void   *options);


PRIVATE double
transform_linear_log(double value,
                     void   *options);


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
    case VRNA_REACTIVITY_TRANS_LINEAR_MODEL:
      return transform_linear;
    case VRNA_REACTIVITY_TRANS_LINEAR_LOG_MODEL:
      return transform_linear_log;
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


PRIVATE double
transform_linear(double value,
                 void   *options)
{
  linear_transform_param_t *o = (linear_transform_param_t *)options;

  return o->intersect + o->slope * value;
}


PRIVATE double
transform_linear_log(double value,
                     void   *options)
{
  linear_transform_param_t *o = (linear_transform_param_t *)options;

  return (value > 0.) ? o->intersect + o->slope * log(value) : o->oob;
}

typedef struct {
  double        lower_bound;
  double        upper_bound;
  unsigned char map_to_ub;
  unsigned char map_to_lb;
} domain_transform_param_t;


PRIVATE double
transform_domain(double v,
                 void   *options)
{
  double vv = v;
  domain_tranform_param_t *o = (domain_transform_param_t *)options;

  if (v < o->lower_bound)
    return (o->map_to_lb) ? o->lower_bound : VRNA_REACTIVITY_MISSING;
  else if (v > o->upper_bound)
    return (o->map_to_ub) ? o->upper_bound : VRNA_READTIVITY_MISSING;

  return vv;
}


typedef struct {
  vrna_array(double)  threshold_source;
  vrna_array(double)  threshold_target;
  unsigned char       project;
} bin_transform_param_t;


PRIVATE double
transform_bin(double  v,
              void    *options)
{
  bin_transform_param_t *o = (bin_transform_param_t *)options;

  if (v < o->threshold_source[0])
    return o->threshold_target[0];

  size_t  start, end;
  start = 0;
  end   = vrna_array_size(o->threshold_source) - 1;

  /* binary search for bin the value will be mapped to */
  while (start <= end) {
    size_t mid = start + (end - start) / 2;

    if ((v >= o->threshold_source[mid]) &&
        (v < o->threshold_source[mid + 1])) {
      if (o->project) {
        double diff_source = o->threshold_source[mid + 1] - o->threshold_source[mid];
        double diff_target = o->threshold_target[mid + 1] - o->threshold_target[mid];
        return  (v - o->threshold_source[mid]) /
                diff_source *
                diff_target +
                o->threshold_target[mid];
      } else {
        return o->threshold_target[mid];
      }
    }

    if (v > o->threshold_source[mid])
      start = mid + 1;
    else
      end = mid - 1;
  }

  return o->threshold_target[vrna_array_size(o->threshold_target) - 1];
}
