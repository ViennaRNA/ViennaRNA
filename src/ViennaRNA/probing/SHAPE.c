/* SHAPE reactivity data handling */

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

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/io/utils.h"
#include "ViennaRNA/io/file_formats.h"
#include "ViennaRNA/probing/basic.h"
#include "ViennaRNA/probing/strategy_deigan.h"
#include "ViennaRNA/probing/strategy_zarringhalam.h"
#include "ViennaRNA/probing/strategy_eddy.h"

#include "ViennaRNA/probing/SHAPE.h"



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

PRIVATE double **
load_n_distribute(unsigned int  n_seq,
                  unsigned int  *ns,
                  const char    **sequences,
                  const char    **file_names,
                  const int     *file_name_association,
                  unsigned int  options);

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */


PUBLIC int
vrna_sc_add_SHAPE_zarringhalam(vrna_fold_compound_t *vc,
                               const double         *reactivities,
                               double               b,
                               double               default_value,
                               const char           *shape_conversion,
                               unsigned int         options)
{
  int ret;

  ret = 0; /* error */

  if ((vc) &&
      (reactivities && (vc->type == VRNA_FC_TYPE_SINGLE))) {
    vrna_probing_data_t d = vrna_probing_data_zarringhalam(reactivities,
                                                           vc->length,
                                                           b,
                                                           shape_conversion,
                                                           default_value);
    ret = vrna_sc_probing(vc, d);
    vrna_probing_data_free(d);
  }

  return ret;
}


PUBLIC int
vrna_sc_add_SHAPE_deigan(vrna_fold_compound_t *vc,
                         const double         *reactivities,
                         double               m,
                         double               b,
                         unsigned int         options)
{
  int                 ret = 0;
  vrna_probing_data_t d;

  if ((vc) &&
      (reactivities)) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        d = vrna_probing_data_deigan(reactivities, vc->length, m, b);
        ret = vrna_sc_probing(vc, d);
        vrna_probing_data_free(d);
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        vrna_log_warning("vrna_sc_add_SHAPE_deigan() not implemented for comparative prediction! "
                             "Use vrna_sc_add_SHAPE_deigan_ali() instead!");
        break;
    }
  }

  return ret;
}


PUBLIC int
vrna_sc_add_SHAPE_deigan_ali(vrna_fold_compound_t *vc,
                             const char           **shape_files,
                             const int            *shape_file_association,
                             double               m,
                             double               b,
                             unsigned int         options)
{
  int           ret;
  unsigned int  s;
  double        **r;

  ret = 0;

  if ((vc) &&
      (vc->type == VRNA_FC_TYPE_COMPARATIVE)) {
    r = load_n_distribute(vc->n_seq,
                          vc->alignment->gapfree_size,
                          (const char **)vc->alignment->gapfree_seq,
                          shape_files,
                          shape_file_association,
                          VRNA_PROBING_DATA_CHECK_SEQUENCE);

    vrna_probing_data_t d = vrna_probing_data_deigan_comparative((const double **)r,
                                                                  vc->alignment->gapfree_size,
                                                                  vc->n_seq,
                                                                  &m,
                                                                  &b,
                                                                  VRNA_PROBING_METHOD_MULTI_PARAMS_0);
    ret = vrna_sc_probing(vc, d);
    vrna_probing_data_free(d);

    for (s = 0; s < vc->n_seq; s++)
      free(r[s]);

    free(r);
  }

  return ret;
}



PUBLIC int
vrna_sc_add_SHAPE_eddy_2(vrna_fold_compound_t *fc,
                         const double         *reactivities,
                         int                  unpaired_nb,
                         const double         *unpaired_data,
                         int                  paired_nb,
                         const double         *paired_data)
{
  int ret;

  ret = 0; /* error */

  if ((fc) &&
      (reactivities) &&
      (fc->type == VRNA_FC_TYPE_SINGLE)) {
    vrna_probing_data_t d = vrna_probing_data_eddy(reactivities,
                                                   fc->length,
                                                   fc->params->temperature,
                                                   VRNA_PROBING_STRATEGY_EDDY_OPTIONS_DEFAULT,
                                                   unpaired_data,
                                                   unpaired_nb,
                                                   paired_data,
                                                   paired_nb);
    ret = vrna_sc_probing(fc, d);
    vrna_probing_data_free(d);
  }

  return ret;
}


PRIVATE double **
load_n_distribute(unsigned int  n_seq,
                  unsigned int  *ns,
                  const char    **sequences,
                  const char    **file_names,
                  const int     *file_name_association,
                  unsigned int  options)
{
  char          *sequence;
  unsigned int  s, ss;
  double        *values, **r;

  r = NULL;

  if ((ns) &&
      (file_names) &&
      (file_name_association)) {
    r = (double **)vrna_alloc(sizeof(double *) * n_seq);

    for (s = 0; file_name_association[s] >= 0; s++) {
      ss = file_name_association[s]; /* actual sequence number in alignment */

      if (ss >= n_seq) {
        vrna_log_warning("Failed to associate probing data file \"%s\" with sequence %d in alignment! "
                         "Omitting data since alignment has only %d sequences!",
                         file_names[s],
                         ss,
                         n_seq);
        continue;
      }

      sequence  = vrna_alloc(sizeof(char) * (ns[ss] + 1));
      values    = vrna_alloc(sizeof(double) * (ns[ss] + 1));

      if (vrna_file_SHAPE_read(file_names[s], ns[ss], -1, sequence, values)) {
        r[ss] = values;

        if ((sequence) &&
            (sequences) &&
            (options & VRNA_PROBING_DATA_CHECK_SEQUENCE)) {
          /* double check information by comparing the sequence read from */
          if (strcmp(sequence, sequences[ss]))
            vrna_log_warning("Input sequence %d differs from sequence provided via probing data file!\n%s\n%s",
                             file_name_association[s] + 1,
                             sequences[ss],
                             sequence);

        }
      } else {
        vrna_log_warning("Failed to open probing data file \"%d\"! "
                         "No data will be used for sequence %d.",
                         s,
                         ss + 1);
      }

      free(sequence);
    }
  }

  return r;
}
