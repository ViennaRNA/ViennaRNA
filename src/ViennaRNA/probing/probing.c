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
#include "ViennaRNA/sequences/alignments.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/io/utils.h"
#include "ViennaRNA/io/file_formats.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/probing/basic.h"


#define gaussian(u) (1 / (sqrt(2 * PI)) * exp(-u * u / 2))


struct vrna_probing_data_s {
  unsigned int method;
  vrna_array(double)    params1;
  vrna_array(double)    params2;
  vrna_array(double *)  reactivities;
  vrna_array(double *)  datas1;
  vrna_array(double *)  datas2;
};


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
PRIVATE int
apply_Deigan2009_method(vrna_fold_compound_t        *fc,
                        struct vrna_probing_data_s  *data);


PRIVATE int
apply_Zarringhalam2012_method(vrna_fold_compound_t        *fc,
                              struct vrna_probing_data_s  *data);


PRIVATE int
apply_Washietl2012_method(vrna_fold_compound_t        *fc,
                          struct vrna_probing_data_s  *data);


PRIVATE int
apply_Eddy2014_method(vrna_fold_compound_t        *fc,
                      struct vrna_probing_data_s  *data);


PRIVATE FLT_OR_DBL
get_msa_weight(unsigned int n_seq,
               double       **datasets);


PRIVATE FLT_OR_DBL
conversion_deigan(double  reactivity,
                  double  m,
                  double  b);


PRIVATE FLT_OR_DBL
conversion_zarringhalam_up(double       beta,
                           double       *pr,
                           unsigned int i);


PRIVATE FLT_OR_DBL
conversion_zarringhalam_bp(double       beta,
                           double       *pr,
                           unsigned int i,
                           unsigned int j);


PRIVATE FLT_OR_DBL
conversion_eddy_up(double       kT,
                   double       *data,
                   unsigned int i);


PRIVATE FLT_OR_DBL
conversion_eddy_bp(double       kT,
                   double       *data,
                   unsigned int i,
                   unsigned int j);


/* PDF of x using Gaussian KDE
 * n is data number
 * h is bandwidth
 */
PRIVATE FLT_OR_DBL
gaussian_kde_pdf(double       x,
                 unsigned int n,
                 float        h,
                 const double *data);


PRIVATE FLT_OR_DBL
exp_pdf(double  x,
        double  lambda) VRNA_UNUSED;


/* We use same format as scitpy */
PRIVATE FLT_OR_DBL
gev_pdf(double  x,
        double  c,
        double  loc,
        double  scale) VRNA_UNUSED;


/*
 * Bandwidth for univariate KDE with Scott factor as in scipy
 * bandwidth = Scott facter * std with ddof = 1
 */
PRIVATE FLT_OR_DBL
bandwidth(unsigned int  n,
          const double  *data);


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
PUBLIC int
vrna_sc_probing(vrna_fold_compound_t  *fc,
                vrna_probing_data_t   data)
{
  int ret = 0;

  if ((fc) && (data)) {
    switch (data->method) {
      case VRNA_PROBING_METHOD_DEIGAN2009:
        ret = apply_Deigan2009_method(fc, data);
        break;

      case VRNA_PROBING_METHOD_ZARRINGHALAM2012:
        ret = apply_Zarringhalam2012_method(fc, data);
        break;

      case VRNA_PROBING_METHOD_WASHIETL2012:
        ret = apply_Washietl2012_method(fc, data);
        break;

      case VRNA_PROBING_METHOD_EDDY2014_2:
        ret = apply_Eddy2014_method(fc, data);
        break;

      default:
        break;
    }
  }

  return ret;
}


PUBLIC struct vrna_probing_data_s *
vrna_probing_data_Deigan2009(const double *reactivities,
                             unsigned int n,
                             double       m,
                             double       b)
{
  struct vrna_probing_data_s *d = NULL;

  if (reactivities)
    d = vrna_probing_data_Deigan2009_comparative(&reactivities,
                                                 &n,
                                                 1,
                                                 &m,
                                                 &b,
                                                 VRNA_PROBING_METHOD_MULTI_PARAMS_0);

  return d;
}


PUBLIC struct vrna_probing_data_s *
vrna_probing_data_Deigan2009_comparative(const double       **reactivities,
                                         const unsigned int *n,
                                         unsigned int       n_seq,
                                         double             *ms,
                                         double             *bs,
                                         unsigned int       multi_params)
{
  struct vrna_probing_data_s  *d = NULL;
  double                      m, b;

  if ((reactivities) && (n)) {
    m = (ms) ? *ms : VRNA_PROBING_METHOD_DEIGAN2009_DEFAULT_m;
    b = (bs) ? *bs : VRNA_PROBING_METHOD_DEIGAN2009_DEFAULT_b;

    if (((ms == NULL) && (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_1)) ||
        ((bs == NULL) && (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_2))) {
      /* if multi_params != 0, either ms or bs or both must be provided! */
      return d;
    }

    d = (struct vrna_probing_data_s *)vrna_alloc(sizeof(struct vrna_probing_data_s));

    d->method = VRNA_PROBING_METHOD_DEIGAN2009;
    vrna_array_init_size(d->params1, n_seq);
    vrna_array_init_size(d->params2, n_seq);
    vrna_array_init_size(d->reactivities, n_seq);

    for (unsigned int i = 0; i < n_seq; i++) {
      if (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_1)
        m = ms[i];

      if (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_2)
        b = bs[i];

      vrna_array_append(d->params1, m);
      vrna_array_append(d->params2, b);

      if (reactivities[i]) {
        /* init and store reactivity data */
        vrna_array(FLT_OR_DBL)  a;
        vrna_array_init_size(a, n[i] + 1);
        for (unsigned int j = 0; j <= n[i]; j++)
          vrna_array_append(a, (FLT_OR_DBL)reactivities[i][j]);

        vrna_array_append(d->reactivities, a);
      } else {
        vrna_array_append(d->reactivities, NULL);
      }
    }

    vrna_array_init(d->datas1);
    vrna_array_init(d->datas2);
  }

  return d;
}


PUBLIC struct vrna_probing_data_s *
vrna_probing_data_Zarringhalam2012(const double *reactivities,
                                   unsigned int n,
                                   double       beta,
                                   const char   *pr_conversion,
                                   double       pr_default)
{
  struct vrna_probing_data_s *d = NULL;

  if (reactivities)
    d = vrna_probing_data_Zarringhalam2012_comparative(&reactivities,
                                                       &n,
                                                       1,
                                                       &beta,
                                                       &pr_conversion,
                                                       &pr_default,
                                                       VRNA_PROBING_METHOD_MULTI_PARAMS_0);

  return d;
}


PUBLIC struct vrna_probing_data_s *
vrna_probing_data_Zarringhalam2012_comparative(const double **reactivities,
                                               unsigned int *n,
                                               unsigned int n_seq,
                                               double       *betas,
                                               const char   **pr_conversions,
                                               double       *pr_defaults,
                                               unsigned int multi_params)
{
  struct vrna_probing_data_s  *d = NULL;
  double                      beta;
  const char                  *pr_conversion;
  double                      pr_default;

  if (reactivities) {
    beta          = (betas) ? *betas : VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_beta;
    pr_conversion = (pr_conversions) ? *pr_conversions : VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_conversion;
    pr_default    = (pr_defaults) ? *pr_defaults : VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_probability;

    if (((betas == NULL) && (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_1)) ||
        ((pr_conversions == NULL) && (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_2)) ||
        ((pr_defaults == NULL) && (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_3))) {
      /* if multi_params != 0, betas must be provided! */
      return d;
    }

    d = (struct vrna_probing_data_s *)vrna_alloc(sizeof(struct vrna_probing_data_s));

    d->method = VRNA_PROBING_METHOD_ZARRINGHALAM2012;
    vrna_array_init_size(d->params1, n_seq);
    vrna_array_init(d->params2);
    vrna_array_init_size(d->reactivities, n_seq);
    vrna_array_init_size(d->datas1, n_seq);

    for (unsigned int i = 0; i < n_seq; i++) {
      if (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_1)
        beta = betas[i];

      vrna_array_append(d->params1, beta);

      if (reactivities[i]) {
        /* init and store reactivity data */
        vrna_array(FLT_OR_DBL)  a;
        vrna_array_init_size(a, n[i] + 1);
        for (unsigned int j = 0; j <= n[i]; j++)
          vrna_array_append(a, (FLT_OR_DBL)reactivities[i][j]);

        vrna_array_append(d->reactivities, a);

        /* prepare probability data according to pr_conversion strategy */
        vrna_array(FLT_OR_DBL)  pr;
        vrna_array_init_size(pr, n[i] + 1);
        for (unsigned int j = 0; j <= n[i]; j++)
          vrna_array_append(pr, (FLT_OR_DBL)reactivities[i][j]);

        if (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_2)
          pr_conversion = pr_conversions[i];

        if (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_3)
          pr_default = pr_defaults[i];

        vrna_sc_SHAPE_to_pr(pr_conversion, pr, n[i], pr_default);
        vrna_array_append(d->datas1, pr);
      } else {
        vrna_array_append(d->reactivities, NULL);
        vrna_array_append(d->datas1, NULL);
      }
    }

    vrna_array_init(d->datas2);
  }

  return d;
}


PUBLIC struct vrna_probing_data_s *
vrna_probing_data_Eddy2014_2(const double *reactivities,
                             unsigned int n,
                             const double *unpaired_data,
                             unsigned int unpaired_len,
                             const double *paired_data,
                             unsigned int paired_len)
{
  struct vrna_probing_data_s *d = NULL;

  if (reactivities)
    d = vrna_probing_data_Eddy2014_2_comparative(&reactivities,
                                                 &n,
                                                 1,
                                                 &unpaired_data,
                                                 &unpaired_len,
                                                 &paired_data,
                                                 &paired_len,
                                                 VRNA_PROBING_METHOD_MULTI_PARAMS_0);

  return d;
}


PUBLIC struct vrna_probing_data_s *
vrna_probing_data_Eddy2014_2_comparative(const double **reactivities,
                                         unsigned int *n,
                                         unsigned int n_seq,
                                         const double **unpaired_datas,
                                         unsigned int *unpaired_lens,
                                         const double **paired_datas,
                                         unsigned int *paired_lens,
                                         unsigned int multi_params)
{
  struct vrna_probing_data_s  *d = NULL;
  double                      unpaired_h, paired_h;

  if ((reactivities) &&
      (unpaired_datas) &&
      (paired_datas) &&
      (unpaired_datas[0]) &&
      (paired_datas[0])) {
    d = (struct vrna_probing_data_s *)vrna_alloc(sizeof(struct vrna_probing_data_s));

    d->method = VRNA_PROBING_METHOD_EDDY2014_2;
    vrna_array_init(d->params1);
    vrna_array_init(d->params2);
    vrna_array_init_size(d->reactivities, n_seq);
    vrna_array_init_size(d->datas1, n_seq);
    vrna_array_init_size(d->datas2, n_seq);

    /* prepare first probabilities */
    unpaired_h  = bandwidth(unpaired_lens[0], unpaired_datas[0]);
    paired_h    = bandwidth(paired_lens[0], paired_datas[0]);

    for (unsigned int i = 0; i < n_seq; i++) {
      if (reactivities[i]) {
        /* init and store reactivity data */
        vrna_array(FLT_OR_DBL)  a;
        vrna_array_init_size(a, n[i] + 1);
        for (unsigned int j = 0; j <= n[i]; j++)
          vrna_array_append(a, (FLT_OR_DBL)reactivities[i][j]);

        vrna_array_append(d->reactivities, a);

        /* kernel-density probability computations */
        vrna_array(FLT_OR_DBL) unpaired;
        vrna_array(FLT_OR_DBL) paired;

        vrna_array_init_size(unpaired, n[i] + 1);
        vrna_array_init_size(paired, n[i] + 1);

        /* Compute bandwidth? */
        if (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_1)
          unpaired_h = bandwidth(unpaired_lens[i], unpaired_datas[i]);

        if (multi_params & VRNA_PROBING_METHOD_MULTI_PARAMS_2)
          paired_h = bandwidth(paired_lens[i], paired_datas[i]);

        /* convert and add */
        vrna_array_append(unpaired, 0);
        vrna_array_append(paired, 0.);

        for (unsigned int j = 1; j <= n[i]; ++j) {
          vrna_array_append(unpaired,
                            log(gaussian_kde_pdf(reactivities[i][j], unpaired_lens[i], unpaired_h,
                                                 unpaired_datas[i])));
          vrna_array_append(paired,
                            log(gaussian_kde_pdf(reactivities[i][j], paired_lens[i], paired_h,
                                                 paired_datas[i])));
        }

        vrna_array_append(d->datas1, unpaired);
        vrna_array_append(d->datas2, paired);
      } else {
        vrna_array_append(d->reactivities, NULL);
        vrna_array_append(d->datas1, NULL);
        vrna_array_append(d->datas2, NULL);
      }
    }
  }

  return d;
}


PUBLIC void
vrna_probing_data_free(struct vrna_probing_data_s *d)
{
  if (d) {
    /* free all reactivity data */
    for (unsigned int i = 0; i < vrna_array_size(d->reactivities); i++)
      vrna_array_free(d->reactivities[i]);
    vrna_array_free(d->reactivities);

    /* free parameters */
    vrna_array_free(d->params1);
    vrna_array_free(d->params2);

    /* free auxiliary data */
    for (unsigned int i = 0; i < vrna_array_size(d->datas1); i++)
      vrna_array_free(d->datas1[i]);

    vrna_array_free(d->datas1);

    for (unsigned int i = 0; i < vrna_array_size(d->datas2); i++)
      vrna_array_free(d->datas2[i]);

    vrna_array_free(d->datas2);

    free(d);
  }
}


PUBLIC int
vrna_sc_SHAPE_to_pr(const char  *shape_conversion,
                    double      *values,
                    int         length,
                    double      default_value)
{
  int *indices;
  int i, j;
  int index;
  int ret = 1;

  if (!shape_conversion || !(*shape_conversion) || length <= 0)
    return 0;

  if (*shape_conversion == 'S')
    return 1;

  indices = vrna_alloc(sizeof(int) * (length + 1));
  for (i = 1, j = 0; i <= length; ++i) {
    if (values[i] < 0)
      values[i] = default_value;
    else
      indices[j++] = i;
  }

  if (*shape_conversion == 'M') {
    double  max;
    double  map_info[4][2] = { {
                                 0.25, 0.35
                               },
                               {
                                 0.30, 0.55
                               },
                               {
                                 0.70, 0.85
                               },
                               {
                                 0, 1
                               } };

    max = values[1];
    for (i = 2; i <= length; ++i)
      max = MAX2(max, values[i]);
    map_info[3][0] = max;

    for (i = 0; indices[i]; ++i) {
      double  lower_source  = 0;
      double  lower_target  = 0;

      index = indices[i];

      if (values[index] == 0)
        continue;

      for (j = 0; j < 4; ++j) {
        if (values[index] > lower_source && values[index] <= map_info[j][0]) {
          double  diff_source = map_info[j][0] - lower_source;
          double  diff_target = map_info[j][1] - lower_target;
          values[index] = (values[index] - lower_source) / diff_source * diff_target + lower_target;
          break;
        }

        lower_source  = map_info[j][0];
        lower_target  = map_info[j][1];
      }
    }
  } else if (*shape_conversion == 'C') {
    float cutoff = 0.25;
    int   i;

    sscanf(shape_conversion + 1, "%f", &cutoff);

    for (i = 0; indices[i]; ++i) {
      index         = indices[i];
      values[index] = values[index] < cutoff ? 0 : 1;
    }
  } else if (*shape_conversion == 'L' || *shape_conversion == 'O') {
    int   i;
    float slope     = (*shape_conversion == 'L') ? 0.68 : 1.6;
    float intercept = (*shape_conversion == 'L') ? 0.2 : -2.29;

    sc_parse_parameters(shape_conversion + 1, 's', 'i', &slope, &intercept);

    for (i = 0; indices[i]; ++i) {
      double v;
      index = indices[i];

      v             = (*shape_conversion == 'L') ? values[index] : log(values[index]);
      values[index] = MAX2(MIN2((v - intercept) / slope, 1), 0);
    }
  } else {
    ret = 0;
  }

  free(indices);

  return ret;
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


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE int
apply_Deigan2009_method(vrna_fold_compound_t        *fc,
                        struct vrna_probing_data_s  *data)
{
  unsigned int  i, s, n, **a2s;
  int           ret;
  FLT_OR_DBL    *vs, **cvs, weight;

  ret = 0;
  n   = fc->length;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      if ((vrna_array_size(data->reactivities) > 0) &&
          (n <= vrna_array_size(data->reactivities[0]))) {
        vs = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));

        /* first convert the values according to provided slope and intercept values */
        for (i = 1; i <= n; ++i)
          vs[i] = conversion_deigan(data->reactivities[0][i], data->params1[0], data->params2[0]);

        /* always store soft constraints in plain format */
        ret = vrna_sc_set_stack(fc, (const FLT_OR_DBL *)vs, VRNA_OPTION_DEFAULT);

        free(vs);
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      if (vrna_array_size(data->reactivities) >= fc->n_seq) {
        a2s = fc->a2s;

        /* collect contributions for the sequences in the alignment */
        cvs = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (fc->n_seq));

        weight = get_msa_weight(fc->n_seq, data->reactivities);

        for (s = 0; s < fc->n_seq; s++) {
          if (data->reactivities[s] != NULL) {
            n = a2s[s][fc->length];

            /*  begin preparation of the pseudo energies */
            cvs[s] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));

            for (i = 1; i <= n; i++) {
              cvs[s][i] = conversion_deigan(data->reactivities[s][i],
                                            data->params1[s],
                                            data->params2[s]) *
                          weight;
            }
          }
        }

        ret = vrna_sc_set_stack_comparative(fc, (const FLT_OR_DBL **)cvs, VRNA_OPTION_DEFAULT);
      }

      break;
  }

  return ret;
}


PRIVATE int
apply_Zarringhalam2012_method(vrna_fold_compound_t        *fc,
                              struct vrna_probing_data_s  *data)
{
  unsigned int  i, j, n, s, **a2s;
  int           ret;
  FLT_OR_DBL    *up, **bp, **ups, ***bps, weight;
  double        *pr, b;

  n   = fc->length;
  ret = 0;
  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      if ((vrna_array_size(data->datas1) > 0) &&
          (n <= vrna_array_size(data->datas1[0]))) {
        /*  now, convert probabilities into pseudo free energies for unpaired, and
         *  paired nucleotides
         */
        pr  = data->datas1[0];
        b   = data->params1[0];

        up  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
        bp  = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (n + 1));

        for (i = 1; i <= n; ++i) {
          bp[i] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
          up[i] = conversion_zarringhalam_up(b, pr, i);

          for (j = i + 1; j <= n; ++j)
            bp[i][j] = conversion_zarringhalam_bp(b, pr, i, j);
        }

        /* add the pseudo energies as soft constraints */
        vrna_sc_set_up(fc, (const FLT_OR_DBL *)up, VRNA_OPTION_DEFAULT);
        vrna_sc_set_bp(fc, (const FLT_OR_DBL **)bp, VRNA_OPTION_DEFAULT);

        /* clean up memory */
        for (i = 1; i <= n; ++i)
          free(bp[i]);
        free(bp);
        free(up);

        ret = 1; /* success */
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      if (vrna_array_size(data->reactivities) >= fc->n_seq) {
        a2s = fc->a2s;

        /* compute weight for each set of probing data */
        weight = get_msa_weight(fc->n_seq, data->reactivities);

        ups = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * fc->n_seq);
        bps = (FLT_OR_DBL ***)vrna_alloc(sizeof(FLT_OR_DBL * *) * fc->n_seq);
        for (s = 0; s < fc->n_seq; s++) {
          if (data->reactivities[s]) {
            /* prepare memory */
            ups[s]  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (fc->length + 1));
            bps[s]  = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (fc->length + 1));

            for (i = 1; i <= fc->length; ++i)
              bps[s][i] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (fc->length + 1));


            /*  now, convert probabilities into pseudo free energies for unpaired, and
             *  paired nucleotides
             */
            n   = a2s[s][fc->length];
            pr  = data->datas1[s];
            b   = data->params1[s] * weight;

            /* currently, unpaired probabilities have to be stored in sequence coordinates */
            for (i = 1; i <= n; ++i)
              ups[s][i] = conversion_zarringhalam_up(b, pr, i);

            /* currently, paired probabilities have to be stored in alignment coordinates */
            for (i = 1; i <= fc->length; ++i)
              for (j = i + 1; j <= fc->length; ++j)
                bps[s][i][j] = conversion_zarringhalam_bp(b, pr, a2s[s][i], a2s[s][j]);
          }
        }

        /* add the pseudo energies as soft constraints */
        (void)vrna_sc_set_up_comparative(fc, (const FLT_OR_DBL **)ups, VRNA_OPTION_DEFAULT);
        ret = vrna_sc_set_bp_comparative(fc, (const FLT_OR_DBL ***)bps, VRNA_OPTION_DEFAULT);

        /* clean up memory */
        for (s = 0; s < fc->n_seq; s++) {
          if (data->reactivities[s]) {
            for (i = 1; i <= fc->length; ++i)
              free(bps[s][i]);
            free(bps[s]);
            free(ups[s]);
          }
        }
        free(bps);
        free(ups);
      }

      break;
  }

  return ret;
}


PRIVATE int
apply_Washietl2012_method(vrna_fold_compound_t        *fc,
                          struct vrna_probing_data_s  *data VRNA_UNUSED)
{
  int ret;

  ret = 0;
  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      break;
  }

  return ret;
}


PRIVATE int
apply_Eddy2014_method(vrna_fold_compound_t        *fc,
                      struct vrna_probing_data_s  *data)
{
  unsigned int  i, j, n, s, **a2s;
  int           ret;
  FLT_OR_DBL    weight;
  double        kT;

  n   = fc->length;
  ret = 0;

  kT = GASCONST * ((fc->params)->temperature + K0) / 1000; /* in kcal/mol */

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      if ((vrna_array_size(data->reactivities) > 0) &&
          (n <= vrna_array_size(data->reactivities[0])) &&
          (vrna_array_size(data->datas1) > 0) &&
          (n <= vrna_array_size(data->datas1[0])) &&
          (vrna_array_size(data->datas2) > 0) &&
          (n <= vrna_array_size(data->datas2[0]))) {
        /* add for unpaired position */
        for (i = 1; i <= n; ++i)
          vrna_sc_add_up(fc, i, conversion_eddy_up(kT, data->datas1[0], i), VRNA_OPTION_DEFAULT);

        /* add for paired position */
        for (i = 1; i <= n; ++i)
          for (j = i + 1; j <= n; ++j)
            vrna_sc_add_bp(fc,
                           i,
                           j,
                           conversion_eddy_bp(kT, data->datas2[0], i, j),
                           VRNA_OPTION_DEFAULT);

        ret = 1; /* success */
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      if (vrna_array_size(data->reactivities) >= fc->n_seq) {
        a2s = fc->a2s;

        /* compute weight for each set of probing data */
        weight = get_msa_weight(fc->n_seq, data->reactivities);

        kT *= weight;

        for (s = 0; s < fc->n_seq; s++) {
          if (data->reactivities[s]) {
            /*  now, convert probabilities into pseudo free energies for unpaired, and
             *  paired nucleotides
             */
            n = a2s[s][fc->length];

            /* add for unpaired position */
            for (i = 1; i <= n; ++i)
              vrna_sc_add_up_comparative_seq(fc,
                                             s,
                                             i,
                                             conversion_eddy_up(kT, data->datas1[s], i),
                                             VRNA_OPTION_DEFAULT);

            /*
             * add for paired position
             * currently, paired probabilities have to be stored in alignment coordinates
             */
            for (i = 1; i <= fc->length; ++i)
              for (j = i + 1; j <= fc->length; ++j)
                vrna_sc_add_bp_comparative_seq(fc,
                                               s,
                                               i,
                                               j,
                                               conversion_eddy_bp(kT, data->datas2[s], a2s[s][i],
                                                                  a2s[s][j]),
                                               VRNA_OPTION_DEFAULT);
          }
        }
      }

      break;
  }

  return ret;
}


/*
 *  computes the individual weight for each set of probing data
 *  by simply spreading-out the probing data contributions to
 *  all sequences where no probing data is provided for.
 */
PRIVATE FLT_OR_DBL
get_msa_weight(unsigned int n_seq,
               double       **datasets)
{
  unsigned int  s, n_data;
  FLT_OR_DBL    weight = 0;

  for (n_data = s = 0; s < n_seq; s++)
    if (datasets[s] != NULL)
      n_data++;

  if (n_data != 0)
    weight = ((FLT_OR_DBL)n_seq / (FLT_OR_DBL)n_data);

  return weight;
}


PRIVATE FLT_OR_DBL
conversion_deigan(double  reactivity,
                  double  m,
                  double  b)
{
  return reactivity < 0 ? 0. : (FLT_OR_DBL)(m * log(reactivity + 1) + b);
}


PRIVATE FLT_OR_DBL
conversion_zarringhalam_up(double       beta,
                           double       *pr,
                           unsigned int i)
{
  return beta * fabs(pr[i] - 1.);
}


PRIVATE FLT_OR_DBL
conversion_zarringhalam_bp(double       beta,
                           double       *pr,
                           unsigned int i,
                           unsigned int j)
{
  return beta * (pr[i] + pr[j]);
}


PRIVATE FLT_OR_DBL
conversion_eddy_up(double       kT,
                   double       *data,
                   unsigned int i)
{
  return -kT * data[i];
}


PRIVATE FLT_OR_DBL
conversion_eddy_bp(double       kT,
                   double       *data,
                   unsigned int i,
                   unsigned int j)
{
  return -kT * (data[i] + data[j]);
}


/*
 * ****************
 *      Eddy
 *****************
 */
PRIVATE FLT_OR_DBL
gaussian_kde_pdf(double       x,
                 unsigned int n,
                 float        h,
                 const double *data)
{
  FLT_OR_DBL total;

  total = 0.;
  for (unsigned int i = 0; i < n; i++)
    total += gaussian((x - data[i]) / h);
  return total / (n * h);
}


PRIVATE FLT_OR_DBL
exp_pdf(double  x,
        double  lambda)
{
  return x < 0 ? 0. : (FLT_OR_DBL)(lambda * exp(-lambda * x));
}


PRIVATE FLT_OR_DBL
gev_pdf(double  x,
        double  c,
        double  loc,
        double  scale)
{
  FLT_OR_DBL s, t;

  s = (FLT_OR_DBL)((x - loc) / scale);
  t = c * s;

  if (c == 0) {
    return exp(-s) * exp(-exp(-s));
  } else if (t < 1) {
    return pow(1 - t, 1 / c - 1) * exp(-pow(1 - t, 1 / c));
  } else {
    return 0;
  }
}


PRIVATE FLT_OR_DBL
bandwidth(unsigned int  n,
          const double  *data)
{
  double factor, mu, std;

  mu  = 0.;
  std = 0.;

  factor = (pow(n, -1. / 5));

  for (unsigned int i = 0; i < n; i++)
    mu += data[i];
  mu /= n;

  for (unsigned int i = 0; i < n; i++)
    std += (data[i] - mu) * (data[i] - mu);
  std = sqrt(std / (n - 1));
  return (FLT_OR_DBL)(factor * std);
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
