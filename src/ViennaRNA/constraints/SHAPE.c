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

#include "ViennaRNA/params/default.h"
#include "ViennaRNA/params/constants.h" /* defines MINPSCORE */
#include "ViennaRNA/datastructures/array.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/utils/alignments.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/io/utils.h"
#include "ViennaRNA/io/file_formats.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/constraints/probing.h"
#include "ViennaRNA/constraints/SHAPE.h"



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
PRIVATE void
sc_parse_parameters(const char  *string,
                    char        c1,
                    char        c2,
                    float       *v1,
                    float       *v2);


PRIVATE FLT_OR_DBL
conversion_deigan(double  reactivity,
                  double  m,
                  double  b);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC void
vrna_constraints_add_SHAPE(vrna_fold_compound_t *vc,
                           const char           *shape_file,
                           const char           *shape_method,
                           const char           *shape_conversion,
                           int                  verbose,
                           unsigned int         constraint_type)
{
  float             p1, p2;
  char              method;
  char              *sequence;
  double            *values;
  int               i, length = vc->length;
  vrna_probing_data_t d = NULL;

  if (!vrna_sc_SHAPE_parse_method(shape_method, &method, &p1, &p2)) {
    vrna_log_warning("Method for SHAPE reactivity data conversion not recognized!");
    return;
  }

  if (verbose) {
    if (method != 'W') {
      if (method == 'Z') {
        vrna_log_info("Using SHAPE method '%c' with parameter p1=%f", method, p1);
      } else {
        vrna_log_info("Using SHAPE method '%c' with parameters p1=%f and p2=%f",
                          method,
                          p1,
                          p2);
      }
    }
  }

  sequence  = vrna_alloc(sizeof(char) * (length + 1));
  values    = vrna_alloc(sizeof(double) * (length + 1));
  vrna_file_SHAPE_read(shape_file, length, method == 'W' ? 0 : -1, sequence, values);

  switch (method) {
    case 'D':
      d = vrna_probing_data_Deigan2009(values,
                                     length,
                                     p1,
                                     p2);
      break;

    case 'Z':
      d = vrna_probing_data_Zarringhalam2012(values,
                                           length,
                                           p1,
                                           shape_conversion,
                                           VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_probability);
      break;

    case 'W':
      FLT_OR_DBL *v = vrna_alloc(sizeof(FLT_OR_DBL) * (length + 1));
      for (i = 0; i < length; i++)
        v[i] = values[i];

      vrna_sc_set_up(vc, v, constraint_type);
      free(v);
      free(values);
      free(sequence);

      return;
  }

  (void)vrna_sc_probing(vc, d);
  vrna_probing_data_free(d);

  free(values);
  free(sequence);
}


PUBLIC void
vrna_constraints_add_SHAPE_ali(vrna_fold_compound_t *vc,
                               const char           *shape_method,
                               const char           **shape_files,
                               const int            *shape_file_association,
                               int                  verbose,
                               unsigned int         constraint_type)
{
  float p1, p2;
  char  method;

  if (!vrna_sc_SHAPE_parse_method(shape_method, &method, &p1, &p2)) {
    vrna_log_warning("Method for SHAPE reactivity data conversion not recognized!");
    return;
  }

  if (method != 'D') {
    vrna_log_warning("SHAPE method %c not implemented for comparative prediction!",
                         method);
    vrna_log_warning("Ignoring SHAPE reactivity data!");
    return;
  } else {
    if (verbose) {
      vrna_log_info("Using SHAPE method '%c' with parameters p1=%f and p2=%f",
                        method,
                        p1,
                        p2);
    }

    vrna_sc_add_SHAPE_deigan_ali(vc, shape_files, shape_file_association, p1, p2, constraint_type);
    return;
  }
}


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
    vrna_probing_data_t d = vrna_probing_data_Zarringhalam2012(reactivities,
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
  int         ret = 0;

  if ((vc) &&
      (reactivities)) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        vrna_probing_data_t d = vrna_probing_data_Deigan2009(reactivities, vc->length, m, b);
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
  FILE          *fp;
  float         reactivity, *reactivities, weight;
  char          *line, nucleotide, *sequence;
  int           s, i, r, n_data, position, n_seq, ret;
  FLT_OR_DBL    **contributions, energy;
  unsigned int  **a2s;

  ret = 0;

  if (vc && (vc->type == VRNA_FC_TYPE_COMPARATIVE)) {
    n_seq = vc->n_seq;
    a2s   = vc->a2s;

    vrna_sc_init(vc);

    /* count number of SHAPE data available for this alignment */
    for (n_data = s = 0; shape_file_association[s] != -1; s++) {
      if (shape_file_association[s] >= n_seq)
        continue;

      /* try opening the shape data file */
      if ((fp = fopen(shape_files[s], "r"))) {
        fclose(fp);
        n_data++;
      }
    }

    weight = (n_data > 0) ? ((float)n_seq / (float)n_data) : 0.;

    /* collect contributions for the sequences in the alignment */
    contributions = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (n_seq));

    for (s = 0; shape_file_association[s] != -1; s++) {
      int ss = shape_file_association[s]; /* actual sequence number in alignment */

      if (ss >= n_seq) {
        vrna_log_warning("Failed to associate SHAPE file \"%s\" with sequence %d in alignment! "
                             "Alignment has only %d sequences!",
                             shape_files[s],
                             ss,
                             n_seq);
        continue;
      }

      /* read the shape file */
      if (!(fp = fopen(shape_files[s], "r"))) {
        vrna_log_warning("Failed to open SHAPE data file \"%d\"! "
                             "No shape data will be used for sequence %d.",
                             s,
                             ss + 1);
      } else {
        reactivities  = (float *)vrna_alloc(sizeof(float) * (vc->length + 1));
        sequence      = (char *)vrna_alloc(sizeof(char) * (vc->length + 1));

        /* initialize reactivities with missing data for entire alignment length */
        for (i = 1; i <= vc->length; i++)
          reactivities[i] = -1.;

        while ((line = vrna_read_line(fp))) {
          r = sscanf(line, "%d %c %f", &position, &nucleotide, &reactivity);
          if (r) {
            if (position <= 0) {
              vrna_log_warning("SHAPE data for position %d outside alignment!", position);
            } else if (position > vc->length) {
              vrna_log_warning("SHAPE data for position %d outside alignment!", position);
            } else {
              switch (r) {
                case 1:
                  nucleotide = 'N';
                /* fall through */
                case 2:
                  reactivity = -1.;
                /* fall through */
                default:
                  sequence[position - 1]  = nucleotide;
                  reactivities[position]  = reactivity;
                  break;
              }
            }
          }

          free(line);
        }
        fclose(fp);

        sequence[vc->length] = '\0';

        /* double check information by comparing the sequence read from */
        char *tmp_seq = vrna_seq_ungapped(vc->sequences[shape_file_association[s]]);
        if (strcmp(tmp_seq, sequence))
          vrna_log_warning("Input sequence %d differs from sequence provided via SHAPE file!",
                               shape_file_association[s] + 1);

        free(tmp_seq);

        /*  begin preparation of the pseudo energies */
        /*  beware of the fact that energy_stack will be accessed through a2s[s] array,
         *  hence pseudo_energy might be gap-free (default)
         */
        int gaps, is_gap;
        contributions[ss] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (vc->length + 1));
        for (gaps = 0, i = 1; i <= vc->length; i++) {
          is_gap  = (vc->sequences[ss][i - 1] == '-') ? 1 : 0;
          energy  =
            ((i - gaps > 0) && !(is_gap)) ? conversion_deigan(reactivities[i - gaps], m,
                                                              b) * weight : 0.;

          if (vc->params->model_details.oldAliEn)
            contributions[ss][i] = energy;
          else if (!is_gap)
            contributions[ss][a2s[ss][i]] = energy;

          gaps += is_gap;
        }

        free(reactivities);
      }
    }

    ret = vrna_sc_set_stack_comparative(vc, (const FLT_OR_DBL **)contributions, options);

    for (s = 0; s < n_seq; s++)
      free(contributions[s]);

    free(contributions);
  }

  return ret;
}


PUBLIC int
vrna_sc_SHAPE_parse_method(const char *method_string,
                           char       *method,
                           float      *param_1,
                           float      *param_2)
{
  const char *params = method_string + 1;

  *param_1  = 0;
  *param_2  = 0;

  if (!method_string || !method_string[0])
    return 0;

  *method = method_string[0];

  switch (method_string[0]) {
    case 'Z':
      *param_1 = 0.89;
      sc_parse_parameters(params, 'b', '\0', param_1, NULL);
      break;

    case 'D':
      *param_1  = 1.8;
      *param_2  = -0.6;
      sc_parse_parameters(params, 'm', 'b', param_1, param_2);
      break;

    case 'W':
      break;

    default:
      *method = 0;
      return 0;
  }

  return 1;
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
      (unpaired_data) &&
      (paired_data) &&
      (fc->type == VRNA_FC_TYPE_SINGLE)) {
    vrna_probing_data_t d = vrna_probing_data_Eddy2014_2(reactivities,
                                                         fc->length,
                                                         unpaired_data,
                                                         unpaired_nb,
                                                         paired_data,
                                                         paired_nb);
    ret = vrna_sc_probing(fc, d);
    vrna_probing_data_free(d);
  }

  return ret;
}


PRIVATE void
sc_parse_parameters(const char  *string,
                    char        c1,
                    char        c2,
                    float       *v1,
                    float       *v2)
{
  char        *fmt;
  const char  warning[] = "SHAPE method parameters not recognized! Using default parameters!";
  int         r;

  assert(c1);
  assert(v1);

  if (!string || !(*string))
    return;

  if (c2 == 0 || v2 == NULL) {
    fmt = vrna_strdup_printf("%c%%f", c1);
    r   = sscanf(string, fmt, v1);

    if (!r)
      vrna_log_warning(warning);

    free(fmt);

    return;
  }

  fmt = vrna_strdup_printf("%c%%f%c%%f", c1, c2);
  r   = sscanf(string, fmt, v1, v2);

  if (r != 2) {
    free(fmt);
    fmt = vrna_strdup_printf("%c%%f", c1);
    r   = sscanf(string, fmt, v1);

    if (!r) {
      free(fmt);
      fmt = vrna_strdup_printf("%c%%f", c2);
      r   = sscanf(string, fmt, v2);

      if (!r)
        vrna_log_warning(warning);
    }
  }

  free(fmt);
}


PRIVATE FLT_OR_DBL
conversion_deigan(double  reactivity,
                  double  m,
                  double  b)
{
  return reactivity < 0 ? 0. : (FLT_OR_DBL)(m * log(reactivity + 1) + b);
}
