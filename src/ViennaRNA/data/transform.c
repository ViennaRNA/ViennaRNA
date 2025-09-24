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
vrna_data_lin_transform(const double              *data,
                        size_t                    data_size,
                        vrna_data_lin_trans_f     trans,
                        vrna_data_lin_trans_opt_t options)
{
  /* init the transformed reactivity array */
  double *a = NULL;

  if ((data) &&
      (data_size > 0)) {
    a = (double *)vrna_alloc(sizeof(double) * data_size);
    a = memcpy(a, data, sizeof(double) * data_size);

    if (trans)
      for (size_t i = 0; i < data_size; ++i)
        a[i] = trans(a[i], options);
  }

  return a;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
