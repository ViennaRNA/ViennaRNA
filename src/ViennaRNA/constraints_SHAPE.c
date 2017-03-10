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

#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/energy_const.h" /* defines MINPSCORE */
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/aln_util.h"
#include "ViennaRNA/file_formats.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/constraints_SHAPE.h"

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


PRIVATE void
sc_add_stack_en_mfe(vrna_fold_compound_t  *vc,
                    const FLT_OR_DBL      *constraints,
                    unsigned int          options);


PRIVATE void
prepare_Boltzmann_weights_stack(vrna_fold_compound_t *vc);


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
  float   p1, p2;
  char    method;
  char    *sequence;
  double  *values;
  int     i, length = vc->length;

  if (!vrna_sc_SHAPE_parse_method(shape_method, &method, &p1, &p2)) {
    vrna_message_warning("Method for SHAPE reactivity data conversion not recognized!");
    return;
  }

  if (verbose) {
    if (method != 'W') {
      if (method == 'Z')
        vrna_message_info(stderr, "Using SHAPE method '%c' with parameter p1=%f", method, p1);
      else
        vrna_message_info(stderr, "Using SHAPE method '%c' with parameters p1=%f and p2=%f", method, p1, p2);
    }
  }

  sequence  = vrna_alloc(sizeof(char) * (length + 1));
  values    = vrna_alloc(sizeof(double) * (length + 1));
  vrna_file_SHAPE_read(shape_file, length, method == 'W' ? 0 : -1, sequence, values);

  if (method == 'D') {
    (void)vrna_sc_add_SHAPE_deigan(vc, (const double *)values, p1, p2, constraint_type);
  } else if (method == 'Z') {
    (void)vrna_sc_add_SHAPE_zarringhalam(vc, (const double *)values, p1, 0.5, shape_conversion, constraint_type);
  } else {
    assert(method == 'W');
    FLT_OR_DBL *v = vrna_alloc(sizeof(FLT_OR_DBL) * (length + 1));
    for (i = 0; i < length; i++)
      v[i] = values[i];

    vrna_sc_set_up(vc, v, constraint_type);

    free(v);
  }

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
    vrna_message_warning("Method for SHAPE reactivity data conversion not recognized!");
    return;
  }

  if (verbose) {
    if (method != 'W') {
      if (method == 'Z')
        vrna_message_info(stderr,
                          "Using SHAPE method '%c' with parameter p1=%f",
                          method, p1);
      else
        vrna_message_info(stderr,
                          "Using SHAPE method '%c' with parameters p1=%f and p2=%f",
                          method, p1, p2);
    }
  }

  if (method == 'D') {
    vrna_sc_add_SHAPE_deigan_ali(vc, shape_files, shape_file_association, p1, p2, constraint_type);
    return;
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
    double  map_info[4][2] = { { 0.25, 0.35 },
                               { 0.30, 0.55 },
                               { 0.70, 0.85 },
                               { 0,    1    } };

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


PUBLIC int
vrna_sc_add_SHAPE_zarringhalam(vrna_fold_compound_t *vc,
                               const double         *reactivities,
                               double               b,
                               double               default_value,
                               const char           *shape_conversion,
                               unsigned int         options)
{
  int         i, j, n, ret;
  double      *pr;
  FLT_OR_DBL  *up, **bp;
  vrna_md_t   *md;

  ret = 0; /* error */

  if (vc && reactivities && (vc->type == VRNA_FC_TYPE_SINGLE)) {
    n   = vc->length;
    md  = &(vc->params->model_details);

    /* first we copy over the reactivities to convert them into probabilities later on */
    pr = (double *)vrna_alloc(sizeof(double) * (n + 1));
    for (i = 0; i <= n; i++)
      pr[i] = reactivities[i];

    if (vrna_sc_SHAPE_to_pr(shape_conversion, pr, n, default_value)) {
      /*  now, convert them into pseudo free energies for unpaired, and
       *  paired nucleotides
       */
      up  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
      bp  = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (n + 1));
      for (i = 1; i <= n; ++i) {
        up[i] = b * fabs(pr[i] - 1);
        bp[i] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
        for (j = i + md->min_loop_size + 1; j <= n; ++j)
          bp[i][j] = b * (pr[i] + pr[j]);
      }

      /* add the pseudo energies as soft constraints */
      vrna_sc_set_up(vc, (const FLT_OR_DBL *)up, options);
      vrna_sc_set_bp(vc, (const FLT_OR_DBL **)bp, options);

      /* clean up memory */
      for (i = 1; i <= n; ++i)
        free(bp[i]);
      free(bp);
      free(up);

      ret = 1; /* success */
    }

    free(pr);
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
  int         i;
  FLT_OR_DBL  *values;

  if (vc && (vc->type == VRNA_FC_TYPE_SINGLE)) {
    if (reactivities) {
      values = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (vc->length + 1));

      /* first convert the values according to provided slope and intercept values */
      for (i = 1; i <= vc->length; ++i)
        values[i] = reactivities[i] < 0 ? 0. : (FLT_OR_DBL)(m * log(reactivities[i] + 1) + b);

      /* always store soft constraints in plain format */
      sc_add_stack_en_mfe(vc, (const FLT_OR_DBL *)values, options);
      free(values);
    }

    if (options & VRNA_OPTION_PF)
      prepare_Boltzmann_weights_stack(vc);

    return 1; /* success */
  }

  return 0; /* error */
}


PUBLIC int
vrna_sc_add_SHAPE_deigan_ali(vrna_fold_compound_t *vc,
                             const char           **shape_files,
                             const int            *shape_file_association,
                             double               m,
                             double               b,
                             unsigned int         options)
{
  FILE            *fp;
  float           reactivity, *reactivities, e1, weight;
  char            *line, nucleotide, *sequence;
  int             s, i, p, r, n_data, position, *pseudo_energies, n_seq;
  unsigned short  **a2s;

  if (vc && (vc->type == VRNA_FC_TYPE_COMPARATIVE)) {
    n_seq = vc->n_seq;
    a2s   = vc->a2s;

    vrna_sc_init(vc);

    /* count number of SHAPE data available for this alignment */
    for (n_data = s = 0; shape_file_association[s] != -1; s++) {
      if (shape_file_association[s] >= n_seq)
        continue;

      /* try opening the shape data file */
      if (fp = fopen(shape_files[s], "r")) {
        fclose(fp);
        n_data++;
      }
    }

    weight = (n_data > 0) ? ((float)n_seq / (float)n_data) : 0.;

    for (s = 0; shape_file_association[s] != -1; s++) {
      int ss = shape_file_association[s]; /* actual sequence number in alignment */

      if (ss >= n_seq) {
        vrna_message_warning("SHAPE file association exceeds sequence number in alignment");
        continue;
      }

      /* read the shape file */
      if (!(fp = fopen(shape_files[s], "r"))) {
        vrna_message_warning("SHAPE data file %d could not be opened. No shape data will be used.", s);
      } else {
        reactivities  = (float *)vrna_alloc(sizeof(float) * (vc->length + 1));
        sequence      = (char *)vrna_alloc(sizeof(char) * (vc->length + 1));

        /* initialize reactivities with missing data for entire alignment length */
        for (i = 1; i <= vc->length; i++)
          reactivities[i] = -1.;

        while ((line = vrna_read_line(fp))) {
          r = sscanf(line, "%d %c %f", &position, &nucleotide, &reactivity);
          if (r) {
            if ((position <= 0) || (position > vc->length))
              vrna_message_error("provided shape data outside of sequence scope");

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

          free(line);
        }
        fclose(fp);

        sequence[vc->length] = '\0';

        /* double check information by comparing the sequence read from */
        char *tmp_seq = get_ungapped_sequence(vc->sequences[shape_file_association[s]]);
        if (strcmp(tmp_seq, sequence))
          vrna_message_warning("Input sequence %d differs from sequence provided via SHAPE file!\n", shape_file_association[s]);

        free(tmp_seq);

        /* convert reactivities to pseudo energies */
        for (i = 1; i <= vc->length; i++) {
          if (reactivities[i] < 0)
            reactivities[i] = 0.;
          else
            reactivities[i] = m * log(reactivities[i] + 1.) + b; /* this should be a value in kcal/mol */

          /* weight SHAPE data derived pseudo energies for this alignment */
          reactivities[i] *= weight;
        }

        /*  begin actual storage of the pseudo energies */
        /*  beware of the fact that energy_stack will be accessed through a2s[s] array,
         *  hence pseudo_energy might be gap-free (default)
         */
        /* ALWAYS store soft constraints in plain format */
        int energy, cnt, gaps, is_gap;
        pseudo_energies = (int *)vrna_alloc(sizeof(int) * (vc->length + 1));
        for (gaps = cnt = 0, i = 1; i <= vc->length; i++) {
          is_gap  = (vc->sequences[ss][i - 1] == '-') ? 1 : 0;
          energy  = ((i - gaps > 0) && !(is_gap)) ? (int)roundf(reactivities[i - gaps] * 100.) : 0;

          if (vc->params->model_details.oldAliEn) {
            pseudo_energies[i] = energy;
            cnt++;
          } else if (!is_gap) {
            /* store gap-free */
            pseudo_energies[a2s[ss][i]] = energy;
            cnt++;
          }

          gaps += is_gap;
        }

        /* resize to actual number of entries */
        pseudo_energies           = vrna_realloc(pseudo_energies, sizeof(int) * (cnt + 2));
        vc->scs[ss]->energy_stack = pseudo_energies;

        if (options & VRNA_OPTION_PF) {
          FLT_OR_DBL *exp_pe = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (vc->length + 1));
          for (i = 0; i <= vc->length; i++)
            exp_pe[i] = 1.;

          for (p = 0, i = 1; i <= vc->length; i++) {
            e1 = (i - p > 0) ? reactivities[i - p] : 0.;
            if (vc->sequences[ss][i - 1] == '-') {
              p++;
              e1 = 0.;
            }

            exp_pe[i] = (FLT_OR_DBL)exp(-(e1 * 1000.) / vc->exp_params->kT);
          }
          vc->scs[ss]->exp_energy_stack = exp_pe;
        }

        free(reactivities);
      }
    }

    return 1; /* success */
  } else {
    return 0; /* error */
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
      vrna_message_warning(warning);

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
        vrna_message_warning(warning);
    }
  }

  free(fmt);
}


PRIVATE void
sc_add_stack_en_mfe(vrna_fold_compound_t  *vc,
                    const FLT_OR_DBL      *constraints,
                    unsigned int          options)
{
  int i;

  if (!vc->sc)
    vrna_sc_init(vc);

  if (!vc->sc->energy_stack)
    vc->sc->energy_stack = (int *)vrna_alloc(sizeof(int) * (vc->length + 1));

  for (i = 1; i <= vc->length; ++i)
    vc->sc->energy_stack[i] = (int)roundf(constraints[i] * 100.);
}


PRIVATE void
prepare_Boltzmann_weights_stack(vrna_fold_compound_t *vc)
{
  int       i;
  vrna_sc_t *sc = vc->sc;

  if (sc->energy_stack) {
    if (!sc->exp_energy_stack) {
      sc->exp_energy_stack = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (vc->length + 1));
      for (i = 0; i <= vc->length; ++i)
        sc->exp_energy_stack[i] = 1.;
    }

    for (i = 1; i <= vc->length; ++i)
      sc->exp_energy_stack[i] = (FLT_OR_DBL)exp(-(sc->energy_stack[i] * 10.) / vc->exp_params->kT);
  }
}
