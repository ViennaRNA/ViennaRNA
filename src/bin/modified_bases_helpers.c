#include <stdlib.h>
#include <string.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"


#include "ViennaRNA/constraints/sc_cb_intern.h"
#include "ViennaRNA/static/energy_parameter_sets.h"

#include "modified_bases_helpers.h"


size_t **
mod_positions_seq_prepare(char                *sequence,
                          vrna_sc_mod_param_t *params,
                          int                 verbose,
                          size_t              *param_set_num)
{
  size_t **mod_positions = NULL;

  *param_set_num = 0;

  if (sequence) {
    if (params)
      for (size_t i = 0; params[i] != NULL; i++)
        (*param_set_num)++;

    if (*param_set_num > 0)
      mod_positions = vrna_alloc(sizeof(size_t *) * *param_set_num);

    if (params) {
      /* split input sequences at default delimiter '&' */
      char **sequences = vrna_strsplit(sequence, NULL);

      size_t i = 0;
      for (vrna_sc_mod_param_t *ptr = params; *ptr != NULL; ptr++, i++) {
        size_t strand_num = 0;
        size_t total_length = 0;

        for (char **s = sequences; *s; s++) {
          size_t *mpos = vrna_strchr(*s, (int)(*ptr)->one_letter_code, 0);

          if (mpos) {
            /* increase size of mod_positions */
            size_t old_size = 0;
            size_t new_size = mpos[0];

            if (mod_positions[i]) {
              old_size = mod_positions[i][0];
              new_size += mod_positions[i][0];
            }

            mod_positions[i] = vrna_realloc(mod_positions[i], sizeof(size_t) * (new_size + 1));

            /* update size of mod_positions */
            mod_positions[i][0] = new_size;

            /* copy-over modified positions in this strand */
            for (size_t j = 1; j <= mpos[0]; j++) {
              mod_positions[i][old_size + j] = mpos[j] + total_length;

              /* change modified base to fallback */
              sequence[mod_positions[i][old_size + j] - 1 + strand_num] = (*ptr)->fallback;

              if (verbose)
                printf("Found modified base %c at position %d\n",
                       (*ptr)->one_letter_code,
                       mod_positions[i][old_size + j]);

            }
          }

          total_length += strlen(*s);
          strand_num++;

          free(mpos);
        }
      }

      /* release memory of individual strands */
      for (char **s = sequences; *s; s++)
        free(*s);
      free(sequences);
    }

  }

  return mod_positions;
}


void
mod_bases_apply(vrna_fold_compound_t  *fc,
                size_t                param_set_num,
                size_t                **mod_positions,
                vrna_sc_mod_param_t   *params)
{
  /* apply modified base support if requested */
  if (param_set_num > 0) {
    size_t        i, j;
    unsigned int  *modification_sites;

    i                   = 0;
    modification_sites  = vrna_alloc(sizeof(unsigned int) * (fc->length + 1));

    if (params) {
      for (vrna_sc_mod_param_t *ptr = params; *ptr != NULL; ptr++, i++) {
        if (mod_positions[i][0] > 0) {
          for (j = 1; j <= mod_positions[i][0]; j++)
            modification_sites[j - 1] = mod_positions[i][j];

          modification_sites[j - 1] = 0;
          vrna_sc_mod(fc, *ptr, modification_sites, VRNA_SC_MOD_DEFAULT);
        }

        free(mod_positions[i]);
        mod_positions[i] = NULL;
      }
    }

    free(modification_sites);
    free(mod_positions);
  }
}


vrna_sc_mod_param_t *
mod_params_collect_from_string(const char           *string,
                               size_t               *num_params,
                               vrna_sc_mod_param_t  *mod_params,
                               vrna_md_t            *md)
{
  if (string) {
    mod_params =
      vrna_realloc(mod_params, sizeof(vrna_sc_mod_param_t) * (*num_params + strlen(string) + 1));

    for (const char *ptr = string; *ptr != '\0'; ptr++) {
      switch (*ptr) {
        case '7':
          mod_params[(*num_params)++] = vrna_sc_mod_read_from_json(
            (const char *)parameter_set_rna_mod_7DA_parameters,
            md);
          break;
        case 'I':
          mod_params[(*num_params)++] = vrna_sc_mod_read_from_json(
            (const char *)parameter_set_rna_mod_inosine_parameters,
            md);
          break;
        case '6':
          mod_params[(*num_params)++] = vrna_sc_mod_read_from_json(
            (const char *)parameter_set_rna_mod_m6A_parameters,
            md);
          break;
        case 'P':
          mod_params[(*num_params)++] = vrna_sc_mod_read_from_json(
            (const char *)parameter_set_rna_mod_pseudouridine_parameters,
            md);
          break;
        case '1':
          mod_params[(*num_params)++] = vrna_sc_mod_read_from_json(
            (const char *)parameter_set_rna_mod_n1methylpseudouridine_parameters,
            md);
          break;
        case '9':
          mod_params[(*num_params)++] = vrna_sc_mod_read_from_json(
            (const char *)parameter_set_rna_mod_purine_parameters,
            md);
          break;
        case 'D':
          mod_params[(*num_params)++] = vrna_sc_mod_read_from_json(
            (const char *)parameter_set_rna_mod_dihydrouridine_parameters,
            md);
          break;
        default:
          break;
      }
    }
    mod_params[*num_params] = NULL;
  }

  return mod_params;
}


vrna_sc_mod_param_t *
mod_params_collect_from_files(const char          **filenames,
                              unsigned int        file_num,
                              size_t              *num_params,
                              vrna_sc_mod_param_t *mod_params,
                              vrna_md_t           *md)
{
  if ((file_num > 0) &&
      (filenames != NULL)) {
    size_t shift = *num_params;

    *num_params += file_num;

    mod_params = vrna_realloc(mod_params, sizeof(vrna_sc_mod_param_t) * (*num_params + 1));

    for (size_t i = 0; i < file_num; i++)
      mod_params[shift + i] = vrna_sc_mod_read_from_jsonfile(filenames[i], md);

    mod_params[*num_params] = NULL;
  }

  return mod_params;
}
