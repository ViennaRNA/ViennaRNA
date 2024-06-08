#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "json/json.h"

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/fold_compound.h"
#include "ViennaRNA/datastructures/string.h"
#include "ViennaRNA/io/utils.h"
#include "ViennaRNA/constraints/soft_special.h"

#include "ViennaRNA/constraints/sc_cb_intern.h"


#ifndef INLINE
# ifdef __GNUC__
#   define INLINE inline
# else
#   define INLINE
# endif
#endif

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

PRIVATE unsigned int
parse_stacks(JsonNode   *dom,
             const char *identifier,
             const char *bases,
             size_t     (*ptypes)[MAX_ALPHABET][MAX_ALPHABET],
             int        (*storage)[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET]);


PRIVATE unsigned int
parse_mismatch(JsonNode   *dom,
               const char *identifier,
               const char *bases,
               size_t     (*ptypes)[MAX_ALPHABET][MAX_ALPHABET],
               vrna_md_t  *md,
               int        (*storage)[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET]);


PRIVATE unsigned int
parse_dangles(JsonNode    *dom,
              const char  *identifier,
              const char  *bases,
              size_t      (*ptypes)[MAX_ALPHABET][MAX_ALPHABET],
              vrna_md_t   *md,
              int         (*storage)[MAX_PAIRS][MAX_ALPHABET]);


PRIVATE unsigned int
parse_terminal(JsonNode   *dom,
               const char *identifier,
               const char *bases,
               size_t     (*ptypes)[MAX_ALPHABET][MAX_ALPHABET],
               int        (*storage)[MAX_PAIRS]);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_sc_mod_param_t
vrna_sc_mod_read_from_jsonfile(const char *filename,
                               vrna_md_t  *md)
{
  char                *ptr;
  FILE                *param_file;
  vrna_sc_mod_param_t params;

  params      = NULL;
  param_file  = fopen(filename, "r");

  if (param_file) {
    vrna_string_t param_content = vrna_string_make("");

    while ((ptr = vrna_read_line(param_file))) {
      param_content = vrna_string_append_cstring(param_content, ptr);
      free(ptr);
    }

    fclose(param_file);

    params = vrna_sc_mod_read_from_json(param_content, md);

    if (!params)
      vrna_log_warning("JSON content could not be read from file \"%s\"",
                           filename);

    vrna_string_free(param_content);
  }

  return params;
}


PUBLIC vrna_sc_mod_param_t
vrna_sc_mod_read_from_json(const char *json,
                           vrna_md_t  *md_p)
{
  char                *ptr, bases[8] = "_ACGUTM";
  vrna_md_t           *md, md_default;
  vrna_sc_mod_param_t parameters = NULL;

  if (json) {
    if (!json_validate(json)) {
      vrna_log_warning("JSON content is not valid\n");
      return NULL;
    }

    JsonNode *dom = json_decode(json);
    if (md_p) {
      md = md_p;
    } else {
      vrna_md_set_default(&md_default);
      md = &md_default;
    }

    if (dom) {
      JsonNode *e, *entry, *mod_data;

      parameters = (struct vrna_sc_mod_param_s *)vrna_alloc(sizeof(struct vrna_sc_mod_param_s));

      parameters->name                = NULL;
      parameters->available           = 0;
      parameters->num_ptypes          = 0;
      parameters->one_letter_code     = '\0';
      parameters->pairing_partners[0] = '\0';
      parameters->unmodified          = '\0';

      mod_data = json_find_member(dom, "modified_base");

      if ((mod_data) &&
          (e = json_find_member(mod_data, "name")) &&
          (e->tag == JSON_STRING))
        parameters->name = strdup(e->string_);

      /* use one-letter code as specified in json file */
      if ((mod_data) &&
          (e = json_find_member(mod_data, "one_letter_code")) &&
          (e->tag == JSON_STRING) &&
          (strlen(e->string_) == 1))
        parameters->one_letter_code = bases[6] = toupper(e->string_[0]);

      if ((mod_data) &&
          (e = json_find_member(mod_data, "unmodified")) &&
          (e->tag == JSON_STRING) &&
          (strlen(e->string_) == 1) &&
          (ptr = strchr(bases, e->string_[0]))) {
        parameters->unmodified = toupper(e->string_[0]);
        size_t enc = ptr - &(bases[0]);
        if (enc > 4)
          enc--;

        parameters->unmodified_encoding = enc;
      }

      if ((mod_data) &&
          (e = json_find_member(mod_data, "fallback")) &&
          (e->tag == JSON_STRING) &&
          (strlen(e->string_) == 1) &&
          (ptr = strchr(bases, e->string_[0]))) {
        parameters->fallback = toupper(e->string_[0]);
        size_t enc = ptr - &(bases[0]);
        if (enc > 4)
          enc--;

        parameters->fallback_encoding = enc;
      }

      size_t cnt = 0;

      if ((mod_data) &&
          (e = json_find_member(mod_data, "pairing_partners")) &&
          (e->tag == JSON_ARRAY)) {
        json_foreach(entry, e) {
          if ((entry->tag == JSON_STRING) &&
              (strlen(entry->string_) == 1) &&
              (ptr = strchr(bases, entry->string_[0]))) {
            size_t enc = ptr - &(bases[0]);
            if (enc > 4)
              enc--;

            parameters->ptypes[5][enc]                    = ++(parameters->num_ptypes);
            parameters->ptypes[enc][5]                    = ++(parameters->num_ptypes);
            parameters->pairing_partners[cnt]             = entry->string_[0];
            parameters->pairing_partners_encoding[cnt++]  = enc;
          }
        }
      }

      parameters->pairing_partners[cnt] = '\0';

      if (parse_stacks(dom, "stacking_energies", &bases[0], &(parameters->ptypes),
                       &(parameters->stack_dG)) > 0)
        parameters->available |= MOD_PARAMS_STACK_dG;

      if (parse_stacks(dom, "stacking_enthalpies", &bases[0], &(parameters->ptypes),
                       &(parameters->stack_dH)) > 0)
        parameters->available |= MOD_PARAMS_STACK_dH;

      if (parse_mismatch(dom, "mismatch_energies", &bases[0], &(parameters->ptypes), md,
                         &(parameters->mismatch_dG)) > 0)
        parameters->available |= MOD_PARAMS_MISMATCH_dG;

      if (parse_mismatch(dom, "mismatch_enthalpies", &bases[0], &(parameters->ptypes), md,
                         &(parameters->mismatch_dH)) > 0)
        parameters->available |= MOD_PARAMS_MISMATCH_dH;

      if (parse_terminal(dom, "terminal_energies", &bases[0], &(parameters->ptypes),
                         &(parameters->terminal_dG)) > 0)
        parameters->available |= MOD_PARAMS_TERMINAL_dG;

      if (parse_terminal(dom, "terminal_enthalpies", &bases[0], &(parameters->ptypes),
                         &(parameters->terminal_dH)) > 0)
        parameters->available |= MOD_PARAMS_TERMINAL_dH;

      if (parse_dangles(dom, "dangle5_energies", &bases[0], &(parameters->ptypes), md,
                        &(parameters->dangle5_dG)) > 0)
        parameters->available |= MOD_PARAMS_DANGLES_dG;

      if (parse_dangles(dom, "dangle5_enthalpies", &bases[0], &(parameters->ptypes), md,
                        &(parameters->dangle5_dH)) > 0)
        parameters->available |= MOD_PARAMS_DANGLES_dH;

      if (parse_dangles(dom, "dangle3_energies", &bases[0], &(parameters->ptypes), md,
                        &(parameters->dangle3_dG)) > 0)
        parameters->available |= MOD_PARAMS_DANGLES_dG;

      if (parse_dangles(dom, "dangle3_enthalpies", &bases[0], &(parameters->ptypes), md,
                        &(parameters->dangle3_dH)) > 0)
        parameters->available |= MOD_PARAMS_DANGLES_dH;

      json_delete(dom);
    }
  }

  return parameters;
}


PUBLIC void
vrna_sc_mod_parameters_free(vrna_sc_mod_param_t params)
{
  if (params) {
    free(params->name);
    free(params);
  }
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE unsigned int
parse_stacks(JsonNode   *dom,
             const char *identifier,
             const char *bases,
             size_t     (*ptypes)[MAX_ALPHABET][MAX_ALPHABET],
             int        (*storage)[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET])
{
  unsigned char num_params = 0;
  char          *enc_ptr;
  size_t        i;
  unsigned int  enc[5] = {
    0, 0, 0, 0, 0
  };
  JsonNode      *entry, *e, *mod_data;

  /* go through storage and initialize */
  for (size_t i = 0; i < MAX_PAIRS; i++)
    for (size_t k = 0; k < MAX_ALPHABET; k++)
      for (size_t l = 0; l < MAX_ALPHABET; l++)
        (*storage)[i][k][l] = INF;

  if (!(mod_data = json_find_member(dom, "modified_base")))
    mod_data = dom;

  if ((e = json_find_member(mod_data, identifier)) &&
      (e->tag == JSON_OBJECT)) {
    json_foreach(entry, e) {
      if ((entry->key) &&
          (entry->tag == JSON_NUMBER) &&
          (strlen(entry->key) == 4)) {
        /* encode sequence */
        for (i = 0; i < 4; i++) {
          if (!(enc_ptr = strchr(&bases[0], entry->key[i]))) {
            vrna_log_warning("Unrecognized character in \"%s\" base: %s\n",
                                 identifier,
                                 entry->key);
            break;
          }

          enc[i] = enc_ptr - &bases[0];
          if (enc[i] > 4)
            enc[i]--;
        }

        if (i == 4) {
          num_params++;

          if ((enc[0] == 5) || (enc[2] == 5))
            (*storage)[(*ptypes)[enc[0]][enc[2]]][enc[3]][enc[1]] = (int)(entry->number_ * 100.);
          else if ((enc[1] == 5) || (enc[3] == 5))
            (*storage)[(*ptypes)[enc[3]][enc[1]]][enc[0]][enc[2]] = (int)(entry->number_ * 100.);
          else
            num_params--;
        }
      }
    }
  }

  return num_params;
}


PRIVATE unsigned int
parse_mismatch(JsonNode   *dom,
               const char *identifier,
               const char *bases,
               size_t     (*ptypes)[MAX_ALPHABET][MAX_ALPHABET],
               vrna_md_t  *md,
               int        (*storage)[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET])
{
  unsigned char num_params = 0;
  char          *enc_ptr;
  size_t        i;
  unsigned int  enc[5] = {
    0, 0, 0, 0, 0
  };
  JsonNode      *entry, *e, *mod_data;

  /* go through storage and initialize */
  for (size_t i = 0; i < MAX_PAIRS; i++)
    for (size_t k = 0; k < MAX_ALPHABET; k++)
      for (size_t l = 0; l < MAX_ALPHABET; l++)
        (*storage)[i][k][l] = INF;

  if (!(mod_data = json_find_member(dom, "modified_base")))
    mod_data = dom;

  if ((e = json_find_member(mod_data, identifier)) &&
      (e->tag == JSON_OBJECT)) {
    json_foreach(entry, e) {
      if ((entry->key) &&
          (entry->tag == JSON_NUMBER) &&
          (strlen(entry->key) == 4)) {
        /* encode sequence */
        for (i = 0; i < 4; i++) {
          if (!(enc_ptr = strchr(&bases[0], entry->key[i]))) {
            vrna_log_warning("Unrecognized character in \"%s\" base: %s\n",
                                 identifier,
                                 entry->key);
            break;
          }

          enc[i] = enc_ptr - &bases[0];
          if (enc[i] > 4)
            enc[i]--;
        }

        if (i == 4) {
          num_params++;

          if ((enc[0] == 5) || (enc[2] == 5))
            (*storage)[NBPAIRS +
                       (*ptypes)[enc[0]][enc[2]]][enc[1]][enc[3]] = (int)(entry->number_ * 100.);
          else if ((enc[1] == 5) || (enc[3] == 5))
            (*storage)[md->pair[enc[0]][enc[2]]][enc[1]][enc[3]] = (int)(entry->number_ * 100.);
          else
            num_params--;
        }
      }
    }
  }

  return num_params;
}


PRIVATE unsigned int
parse_dangles(JsonNode    *dom,
              const char  *identifier,
              const char  *bases,
              size_t      (*ptypes)[MAX_ALPHABET][MAX_ALPHABET],
              vrna_md_t   *md,
              int         (*storage)[MAX_PAIRS][MAX_ALPHABET])
{
  unsigned char num_params = 0;
  char          *enc_ptr;
  size_t        i;
  unsigned int  enc[5] = {
    0, 0, 0, 0, 0
  };
  JsonNode      *entry, *e, *mod_data;

  /* go through storage and initialize */
  for (size_t i = 0; i < MAX_PAIRS; i++)
    for (size_t k = 0; k < MAX_ALPHABET; k++)
      (*storage)[i][k] = INF;

  if (!(mod_data = json_find_member(dom, "modified_base")))
    mod_data = dom;

  if ((e = json_find_member(mod_data, identifier)) &&
      (e->tag == JSON_OBJECT)) {
    json_foreach(entry, e) {
      if ((entry->key) &&
          (entry->tag == JSON_NUMBER) &&
          (strlen(entry->key) == 3)) {
        /* encode sequence */
        for (i = 0; i < 3; i++) {
          if (!(enc_ptr = strchr(&bases[0], entry->key[i]))) {
            vrna_log_warning("Unrecognized character in \"%s\" base: %s\n",
                                 identifier,
                                 entry->key);
            break;
          }

          enc[i] = enc_ptr - &bases[0];
          if (enc[i] > 4)
            enc[i]--;
        }

        if (i == 3) {
          num_params++;
          if ((enc[0] == 5) || (enc[1] == 5))
            (*storage)[NBPAIRS + (*ptypes)[enc[0]][enc[1]]][enc[2]] = (int)(entry->number_ * 100.);
          else if (enc[2] == 5)
            (*storage)[md->pair[enc[0]][enc[1]]][enc[2]] = (int)(entry->number_ * 100.);
          else
            num_params--;
        }
      }
    }
  }

  return num_params;
}


PRIVATE unsigned int
parse_terminal(JsonNode   *dom,
               const char *identifier,
               const char *bases,
               size_t     (*ptypes)[MAX_ALPHABET][MAX_ALPHABET],
               int        (*storage)[MAX_PAIRS])
{
  unsigned char num_params = 0;
  char          *enc_ptr;
  size_t        i;
  unsigned int  enc[5] = {
    0, 0, 0, 0, 0
  };
  JsonNode      *entry, *e, *mod_data;

  /* go through storage and initialize */
  for (size_t i = 0; i < MAX_PAIRS; i++)
    (*storage)[i] = INF;

  if (!(mod_data = json_find_member(dom, "modified_base")))
    mod_data = dom;

  if ((e = json_find_member(mod_data, identifier)) &&
      (e->tag == JSON_OBJECT)) {
    json_foreach(entry, e) {
      if ((entry->key) &&
          (entry->tag == JSON_NUMBER) &&
          (strlen(entry->key) == 2)) {
        /* encode sequence */
        for (i = 0; i < 2; i++) {
          if (!(enc_ptr = strchr(&bases[0], entry->key[i]))) {
            vrna_log_warning("Unrecognized character in \"%s\" base: %s\n",
                                 identifier,
                                 entry->key);
            break;
          }

          enc[i] = enc_ptr - &bases[0];
          if (enc[i] > 4)
            enc[i]--;
        }

        if (i == 2) {
          num_params++;

          if ((enc[0] == 5) || (enc[1] == 5))
            (*storage)[(*ptypes)[enc[0]][enc[1]]] = (int)(entry->number_ * 100.);
          else
            num_params--;
        }
      }
    }
  }

  return num_params;
}
