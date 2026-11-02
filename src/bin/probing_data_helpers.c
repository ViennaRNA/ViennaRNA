#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/utils/log.h"
#include <ViennaRNA/utils/units.h>
#include "ViennaRNA/datastructures/array.h"
#include "ViennaRNA/io/file_formats.h"
#include "ViennaRNA/probing/basic.h"
#include "ViennaRNA/probing/strategy_deigan.h"
#include "ViennaRNA/probing/strategy_eddy.h"
#include "ViennaRNA/probing/strategy_zarringhalam.h"
#include "ViennaRNA/probing/SHAPE.h"
#include "ViennaRNA/probing/crosslink.h"


#include "probing_data_helpers.h"

typedef struct {
  char    key;
  float   value;
} kv_pair;


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE vrna_array(kv_pair)
sc_parse_parameters(const char        *string,
                    vrna_array(char)  *aux_flags);

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
                   probing_data_t       *d,
                   const char           *allowed_methods)
{
  const char          *default_strategies = "DEZW";
  char                method, *sequence;
  unsigned int        length, options;
  size_t              prior_unpaired_size, prior_paired_size;
  float               p1, p2, p3;
  double              *values, *prior_unpaired, *prior_paired;
  vrna_probing_data_t pd;

  if (allowed_methods == NULL)
    allowed_methods = default_strategies;

  length = fc->length;

  for (size_t i = 0; i < vrna_array_size(d->strategies); i++) {
    pd = NULL;
    p1 = 0.;
    p2 = 0.;

    /* parse the strategy and corresponding arguments */
    method = d->strategies[i] ? d->strategies[i][0] : 'D';

    if (strchr(allowed_methods, method) == NULL) {
      vrna_log_error("Probing strategy \"%c\" not allowed in this context, skipping data! Available strategies are: \"%s\"",
                     method,
                     allowed_methods);
      continue;
    }

    /* parse parameters, if any */
    vrna_array(char) aux_flags;
    vrna_array(kv_pair) parameters;
    parameters = sc_parse_parameters(d->strategies[i] ? d->strategies[i] + 1: "",
                                     &aux_flags);

    /* process parameters */
    switch (method) {
      case 'Z':
        p1 = VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_beta;
        if (parameters) {
          for (size_t i = 0; i < vrna_array_size(parameters); i++)
            if (parameters[i].key == 'b') {
              p1 = parameters[i].value;
              break;
            }
        }
        break;

      case 'E':
        p1 = -300.;
        if (parameters) {
          for (size_t i = 0; i < vrna_array_size(parameters); i++)
            if ((parameters[i].key == 't') ||
                (parameters[i].key == 'k')) {
              p1 = (float)vrna_convert_temperature(parameters[i].value,
                                                         VRNA_UNIT_K,
                                                         VRNA_UNIT_DEG_C);
              break;
            } else if (parameters[i].key == 'c') {
              p1 = parameters[i].value;
              break;
            }
        }
        break;

      case 'D':
        p1  = VRNA_PROBING_METHOD_DEIGAN2009_DEFAULT_m;
        p2  = VRNA_PROBING_METHOD_DEIGAN2009_DEFAULT_b;
        if (parameters) {
          for (size_t i = 0; i < vrna_array_size(parameters); i++)
            if (parameters[i].key == 'b') {
              p2 = parameters[i].value;
            } else if (parameters[i].key == 'm') {
              p1 = parameters[i].value;
              break;
            }
        }
        break;

      case 'W':
        break;

      case 'X':
        /* Cross-linking data */
        p1      = fc->params->model_details.temperature; /* default temperature */
        p2      = VRNA_PROBING_METHOD_CROSSLINK_DEFAULT_SCALING; /* default scaling factor */
        p3      = VRNA_PROBING_METHOD_CROSSLINK_DEFAULT_LIMIT; /* limit in energy */
        options = VRNA_PROBING_METHOD_CROSSLINK_DEFAULT_OPTIONS;

        if (aux_flags) {
          for (size_t i = 0; i < vrna_array_size(aux_flags); i++) {
            if (aux_flags[i] == 'E')
              options |= VRNA_PROBING_METHOD_CROSSLINK_INPUT_ENERGY;  /* input file states energies instead of probabilities */
            else if (aux_flags[i] == 'P')
              options &= ~VRNA_PROBING_METHOD_CROSSLINK_INPUT_ENERGY; /* input file states probabilities */
          }
        }

        if (parameters) {
          for (size_t i = 0; i < vrna_array_size(parameters); i++)
            if ((parameters[i].key == 't') ||
                (parameters[i].key == 'k')) {
              p1 = (float)vrna_convert_temperature(parameters[i].value,
                                                         VRNA_UNIT_K,
                                                         VRNA_UNIT_DEG_C);
            } else if (parameters[i].key == 'c') {
              p1 = parameters[i].value;
            } else if (parameters[i].key == 's') {
              p2 = parameters[i].value;
            } else if (parameters[i].key == 'e') {
              p3 = parameters[i].value;
              options &= ~VRNA_PROBING_METHOD_CROSSLINK_LIMIT_PROBABILITY; /* limit is probability */
            } else if (parameters[i].key == 'p') {
              p3 = parameters[i].value;
              options |= VRNA_PROBING_METHOD_CROSSLINK_LIMIT_PROBABILITY; /* limit is probability */
            }
        }

        vrna_log_info("Add crosslinking data from \"%s\" with T=%3.2f, scale=%g, limit=%g, options = %u",
                      d->files[i],
                      p1,
                      p2,
                      p3,
                      options);

        if (!vrna_probing_data_xlink_file(fc,
                                          d->files[i],
                                          p1,
                                          p2,
                                          p3,
                                          options))
          vrna_log_error("failed to load crosslink data from \"%s\"\n",
                         d->files[i]);

        break;
    }

    if (method == 'X') {
      vrna_array_free(parameters);
      vrna_array_free(aux_flags);
      /* do not proceed with the usual probing data integration below */
      continue;
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

          vrna_log_info("Probing data set %ld: file = %s, strategy = Washietl et al. 2012, pre-processing = %s",
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
          vrna_log_warning("Probing data set/strategy %ld not recognized", i);
          break;
      }
    }

    (void)vrna_sc_probing(fc, pd);

    vrna_probing_data_free(pd);
    vrna_array_free(parameters);
    vrna_array_free(aux_flags);
    free(values);
    free(sequence);
  }

  return 0;
}


PUBLIC void
constraints_add_SHAPE(vrna_fold_compound_t *vc,
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

  if ((!vrna_sc_SHAPE_parse_method(shape_method, &method, &p1, &p2)) ||
      (method == 'X')) {
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
constraints_add_SHAPE_ali(vrna_fold_compound_t *vc,
                          const char           *shape_method,
                          const char           **shape_files,
                          const int            *shape_file_association,
                          int                  verbose,
                          unsigned int         constraint_type)
{
  float p1, p2;
  char  method;

  if ((!vrna_sc_SHAPE_parse_method(shape_method, &method, &p1, &p2)) ||
      (method == 'X')) {
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

/* Parse key-value pairs with key = char and value = float from a string */
PRIVATE vrna_array(kv_pair)
sc_parse_parameters(const char        *string,
                    vrna_array(char)  *aux_flags)
{
  vrna_array(kv_pair) parameters;

  parameters = NULL;
  *aux_flags = NULL;

  if (!string || !(*string))
    return parameters;

  size_t l    = strlen(string);
  size_t pos  = 0;

  vrna_array_init(parameters);
  vrna_array_init(*aux_flags);

  while (pos < l) {
    int   consumed = 0;
    char  p;
    float v;
    /* loop through the string and try extracting as many parameters as possible */
    if (sscanf(string + pos, " %c %f%n", &p, &v, &consumed) == 2) {
      vrna_array_append(parameters,
                        ((kv_pair){
                          .key = p,
                          .value = v
                        }));
      pos += (size_t)consumed;
    } else if (sscanf(string + pos, " %c = %f%n", &p, &v, &consumed) == 2) {
      vrna_array_append(parameters,
                        ((kv_pair){
                          .key = p,
                          .value = v
                        }));
      pos += (size_t)consumed;
    } else if (sscanf(string + pos, " %c : %f%n", &p, &v, &consumed) == 2) {
      vrna_array_append(parameters,
                        ((kv_pair){
                          .key = p,
                          .value = v
                        }));
      pos += (size_t)consumed;
    } else if ((sscanf(string + pos, " %[A-Z]%n", &p, &consumed) == 1)) {
      vrna_array_append(*aux_flags, p);
      pos += (size_t)consumed;
    } else {
      vrna_log_warning("can't parse any more parameters from \"%s\"", string + pos);
      break;
    }
  }

  return parameters;
}


PUBLIC double **
probing_data_load_n_distribute(unsigned int  n_seq,
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
  for (size_t i = 1; i < (size_t)argc; i++) {
    if (strncmp(argv[i], "--sp-", 5) == 0) {
      /* extract argument for this option */
      char *cmd, *arg;
      cmd = argv[i];

      if ((arg = strchr(cmd, '=')) != NULL) {
        arg++; /* remove leading '=' */
      } else if (i + 1 < (size_t)argc) {
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

