#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/utils/log.h"
#include <ViennaRNA/utils/units.h>
#include "ViennaRNA/io/file_formats.h"
#include "ViennaRNA/probing/basic.h"
#include "ViennaRNA/probing/strategy_deigan.h"
#include "ViennaRNA/probing/strategy_eddy.h"
#include "ViennaRNA/probing/strategy_zarringhalam.h"
#include "ViennaRNA/probing/SHAPE.h"


#include "probing_data_helpers.h"


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


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC void
probing_data_free(probing_data_t *dat) {
  if (dat) {
    for (size_t i = 0; i < vrna_array_size(dat->files); ++i)
      free(dat->files[i]);

    vrna_array_free(dat->files);

    for (size_t i = 0; i < vrna_array_size(dat->strategies); ++i)
      free(dat->strategies[i]);

    vrna_array_free(dat->strategies);

    for (size_t i = 0; i < vrna_array_size(dat->preprocessing); ++i)
      free(dat->preprocessing[i]);

    vrna_array_free(dat->preprocessing);

    for (size_t i = 0; i < vrna_array_size(dat->prior_unpaired); ++i)
      free(dat->prior_unpaired[i]);

    vrna_array_free(dat->prior_unpaired);

    for (size_t i = 0; i < vrna_array_size(dat->prior_paired); ++i)
      free(dat->prior_paired[i]);

    vrna_array_free(dat->prior_paired);

    free(dat);
  }
}


PUBLIC int
apply_probing_data(vrna_fold_compound_t *fc,
                   probing_data_t       *d)
{
  char                method, *sequence;
  unsigned int        length;
  size_t              prior_unpaired_size, prior_paired_size;
  float               p1, p2;
  double              *values, *prior_unpaired, *prior_paired;
  vrna_probing_data_t pd;

  length = fc->length;

  for (size_t i = 0; i < vrna_array_size(d->strategies); i++) {
    pd = NULL;

    /* parse the strategy and corresponding arguments */
    if (!vrna_sc_SHAPE_parse_method(d->strategies[i] ? d->strategies[i] : "D", &method, &p1, &p2)) {
      vrna_log_warning("Method for SHAPE reactivity data conversion not recognized!");
      return 1;
    }

    sequence  = vrna_alloc(sizeof(char) * (length + 1));
    values    = vrna_alloc(sizeof(double) * (length + 1));

    if (vrna_file_SHAPE_read(d->files[i],
                             length,
                             method == 'W' ? 0 : -1,
                             sequence,
                             values)) {
      switch (method) {
        case 'D':
          pd = vrna_probing_data_deigan(values,
                                       length,
                                       p1,
                                       p2);
          if (pd)
            vrna_log_info("Probing data set %ld: file = %s, strategy = Deigan et al. 2009, m = %g, b = %g, pre-processing = %s",
                          i,
                          d->files[i],
                          p1,
                          p2,
                          d->preprocessing[i] ? d->preprocessing[i] : "None");

          break;

        case 'E':
          prior_unpaired_size = prior_paired_size = 0;
          prior_unpaired      = prior_paired      = NULL;

          if (d->prior_unpaired[i])
            prior_unpaired = vrna_file_doubles_read(d->prior_unpaired[i],
                                                    0,
                                                    &prior_unpaired_size);

          if (d->prior_paired[i])
            prior_paired = vrna_file_doubles_read(d->prior_paired[i],
                                                    0,
                                                    &prior_paired_size);

          pd = vrna_probing_data_eddy(values,
                                     length,
                                     p1,
                                     (p1 > -300.) ? VRNA_PROBING_STRATEGY_EDDY_NO_TEMPERATURE_RESCALING : VRNA_PROBING_STRATEGY_EDDY_OPTIONS_DEFAULT,
                                     prior_unpaired,
                                     prior_unpaired_size,
                                     prior_paired,
                                     prior_paired_size);

          if (pd)
            vrna_log_info("Probing data set %ld: file = %s, strategy = Eddy 2014, t = %g, pre-processing = %s",
                          i,
                          d->files[i],
                          (p1 > -300.) ? p1 : fc->params->temperature,
                          d->preprocessing[i] ? d->preprocessing[i] : "None");

          free(prior_unpaired);
          free(prior_paired);

          break;

        case 'W':
          /* this case won't trigger vrna_sc_probing() since we add soft constraints here already */
          for (i = 1; i < length; i++)
            (void)vrna_sc_add_up(fc, i, (FLT_OR_DBL)values[i], VRNA_OPTION_DEFAULT);

          vrna_log_info("Probing data set %ld: file = %s, strategy = Wuchty et al. 2012, pre-processing = %s",
                        i,
                        d->files[i],
                        d->preprocessing[i] ? d->preprocessing[i] : "None");

          break;

        case 'Z':
          pd = vrna_probing_data_zarringhalam(values,
                                             length,
                                             p1,
                                             d->preprocessing[i],
                                             VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_probability);

          if (pd)
            vrna_log_info("Probing data set %ld: file = %s, strategy = Zarringhalam 2012, b = %g, pre-processing = %s",
                          i,
                          d->files[i],
                          p1,
                          d->preprocessing[i] ? d->preprocessing[i] : "None");

          break;

        default:
          vrna_log_warning("Probing data set %ld not recognized", i);
          break;
      }
    }

    (void)vrna_sc_probing(fc, pd);

    vrna_probing_data_free(pd);
    free(values);
    free(sequence);
  }

  return 0;
}


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
  FLT_OR_DBL        *v;
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
      d = vrna_probing_data_deigan(values,
                                   length,
                                   p1,
                                   p2);
      break;

    case 'Z':
      d = vrna_probing_data_zarringhalam(values,
                                         length,
                                         p1,
                                         shape_conversion,
                                         VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_probability);
      break;

    case 'W':
      v = vrna_alloc(sizeof(FLT_OR_DBL) * (length + 1));
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
      *param_1 = VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_beta;
      sc_parse_parameters(params, 'b', '\0', param_1, NULL);
      break;

    case 'E':
      *param_1 = -300.;

      sc_parse_parameters(params, 't', '\0', param_1, NULL);

      if (*param_1 == -300.)
        sc_parse_parameters(params, 'k', '\0', param_1, NULL);

      if (*param_1 == -300.)
        sc_parse_parameters(params, 'c', '\0', param_1, NULL);
      else
        *param_1 = (float)vrna_convert_temperature((float)(*param_1),
                                                   VRNA_UNIT_K,
                                                   VRNA_UNIT_DEG_C);
      break;

    case 'D':
      *param_1  = VRNA_PROBING_METHOD_DEIGAN2009_DEFAULT_m;
      *param_2  = VRNA_PROBING_METHOD_DEIGAN2009_DEFAULT_b;
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


PUBLIC double **
vrna_probing_data_load_n_distribute(unsigned int  n_seq,
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


PUBLIC probing_data_t *
extract_probing_options(int     argc,
                        char    *argv[],
                        size_t  expected_num_data)
{
  char  *file_r, *file_tmp, *file_pu, *file_pp, *strategy, *preprocess;

  probing_data_t *probing_data_p = (probing_data_t *)vrna_alloc(sizeof(probing_data_t));

  probing_data_p->count = 0;
  vrna_array_init_size(probing_data_p->files, expected_num_data);
  vrna_array_init_size(probing_data_p->strategies, expected_num_data);
  vrna_array_init_size(probing_data_p->preprocessing, expected_num_data);
  vrna_array_init_size(probing_data_p->prior_unpaired, expected_num_data);
  vrna_array_init_size(probing_data_p->prior_paired, expected_num_data);

  file_r = file_pu = file_pp = file_tmp = NULL;
  strategy = preprocess = NULL;

  /* Go through argument options and collect all probing data */
  for (size_t i = 1; i < argc; i++) {
    if (strncmp(argv[i], "--sp-", 5) == 0) {
      /* extract argument for this option */
      char *cmd, *arg;
      cmd = argv[i];

      if ((arg = strchr(cmd, '=')) != NULL) {
        arg++; /* remove leading '=' */
      } else if (i + 1 < argc) {
        arg = argv[++i];
      } else {
        vrna_log_error("Command line option %s is missing its argument!");
      }

      /* now, do something */
      switch (cmd[5]) {
        case 'd':
          /* data file */
          if (file_tmp) {
            /*
             *  we read a new data file after already having a data file from the previous round
             *  without any further specification of the strategy, so we assume default strategy
             *  and add the previous data file with default strategy to our list
             */
            vrna_array_append(probing_data_p->files, file_tmp);
            vrna_array_append(probing_data_p->strategies, NULL);
            vrna_array_append(probing_data_p->preprocessing, preprocess);
            vrna_array_append(probing_data_p->prior_unpaired, NULL);
            vrna_array_append(probing_data_p->prior_paired, NULL);
            /* reset everything */
            file_r = file_pu = file_pp = strategy = preprocess = NULL;
          } else if ((file_r) &&
                     (strategy) &&
                     (strategy[0] != 'E')) {
            /* we already have a strategy other than Eddy 2014, so let us store everything */
            vrna_array_append(probing_data_p->files, file_r);
            vrna_array_append(probing_data_p->strategies, strategy);
            vrna_array_append(probing_data_p->preprocessing, preprocess);
            vrna_array_append(probing_data_p->prior_unpaired, NULL);
            vrna_array_append(probing_data_p->prior_paired, NULL);
            file_r = file_pu = file_pp = strategy = preprocess = NULL;
          }

          /* store current file name */
          file_tmp = strdup(arg);

          break;

        case 's':
          /* strategy */

          /* first, check whether we are expecting to use Eddy 2014 approach and may use user-defined prior distributions */
          if ((strategy) &&
              (strategy[0] == 'E')) {
            if (arg[0] == 'P') { /* last data file was prior data set */
              if (arg[1] == 'u') {
                free(file_pu);
                file_pu   = file_tmp;
                file_tmp  = NULL;
              } else if (arg[1] == 'p') {
                free(file_pp);
                file_pp   = file_tmp;
                file_tmp  = NULL;
              } else {
                vrna_log_error("Unrecognzed prior data set!");
                free(file_tmp);
                file_tmp  = NULL;
              }
              /* break out of case */
              break;
            } else {
              /* last data file was a new probing data set, so store data for previously read Eddy strategy */
              vrna_array_append(probing_data_p->files, file_r);
              vrna_array_append(probing_data_p->strategies, strategy);
              vrna_array_append(probing_data_p->preprocessing, preprocess);
              vrna_array_append(probing_data_p->prior_unpaired, file_pu);
              vrna_array_append(probing_data_p->prior_paired, file_pp);
              file_r = file_pu = file_pp = strategy = preprocess = NULL;
            }
          }

          /* overwrite/store strategy */
          free(strategy);
          strategy  = strdup(arg);

          /* last data file contains actual probing data */
          file_r    = file_tmp;
          file_tmp  = NULL;

          break;

        case 'p':
          /* preprocessing */

          /* overwrite/store pre-processing option */
          free(preprocess);
          preprocess = strdup(arg);
          break;

        default:
          break;
      }
    }
  }

  if (file_r) {
    vrna_array_append(probing_data_p->files, file_r);
    vrna_array_append(probing_data_p->strategies, strategy);
    vrna_array_append(probing_data_p->preprocessing, preprocess);
    vrna_array_append(probing_data_p->prior_unpaired, file_pu);
    vrna_array_append(probing_data_p->prior_paired, file_pp);
  }

  if (file_tmp) {
    vrna_array_append(probing_data_p->files, file_tmp);
    vrna_array_append(probing_data_p->strategies, NULL);
    vrna_array_append(probing_data_p->preprocessing, preprocess);
    vrna_array_append(probing_data_p->prior_unpaired, NULL);
    vrna_array_append(probing_data_p->prior_paired, NULL);
  }

  return probing_data_p;
}

