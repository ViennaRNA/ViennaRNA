/* SHAPE reactivity transform */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/data/transform.h"

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
PUBLIC double *
vrna_data_lin_transform(const double        *data,
                        size_t              data_size,
                        vrna_math_fun_dbl_f     trans,
                        vrna_math_fun_dbl_opt_t trans_options,
                        double              domain[4],
                        double              oob_value,
                        unsigned int        options)
{
  /* init the transformed reactivity array */
  double  *a  = NULL;
  size_t  i, cnt, *vs = NULL;

  if ((data) &&
      (data_size > 0)) {
    a = (double *)vrna_alloc(sizeof(double) * data_size);
    a = memcpy(a, data, sizeof(double) * data_size);

    if (domain == NULL)
      options &= ~VRNA_TRANSFORM_ENFORCE_DOMAINS;

    if (options & VRNA_TRANSFORM_ENFORCE_DOMAIN_SOURCE) {
      vs = (size_t *)vrna_alloc(sizeof(size_t) * data_size);

      for (i = cnt = 0; i < data_size; ++i) {
        if (a[i] < domain[0]) {
          if (options & VRNA_TRANSFORM_MAP_SOURCE_LOW) {
            a[i] = domain[0];
            vs[cnt++] = i;
          } else {
            a[i] = oob_value;
          }
        } else if (a[i] > domain[1]) {
          if (options & VRNA_TRANSFORM_MAP_SOURCE_HIGH) {
            a[i] = domain[1];
            vs[cnt++] = i;
          } else {
            a[i] = oob_value;
          }
        } else {
          vs[cnt++] = i;
        }
      }

      if (trans)
        for (i = 0; i < cnt; ++i)
          a[vs[i]] = trans(a[vs[i]], trans_options);

      free(vs);
    } else if (trans) {
      for (size_t i = 0; i < data_size; ++i)
        a[i] = trans(a[i], trans_options);
    }

    if (options & VRNA_TRANSFORM_ENFORCE_DOMAIN_TARGET) {
      for (i = 0; i < data_size; ++i) {
        if (a[i] != oob_value) {
          if (a[i] < domain[2]) {
            if (options & VRNA_TRANSFORM_MAP_TARGET_LOW)
              a[i] = domain[2];
            else
              a[i] = oob_value;
          } else if (a[i] > domain[3]) {
            if (options & VRNA_TRANSFORM_MAP_TARGET_HIGH)
              a[i] = domain[3];
            else
              a[i] = oob_value;
          }
        }
      }
    }
  }

  return a;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
