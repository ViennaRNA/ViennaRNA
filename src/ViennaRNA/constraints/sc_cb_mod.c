#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/fold_compound.h"
#include "ViennaRNA/datastructures/string.h"
#include "ViennaRNA/datastructures/array.h"
#include "ViennaRNA/io/utils.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/params/default.h"
#include "json/json.h"

#ifndef INLINE
# ifdef __GNUC__
#   define INLINE inline
# else
#   define INLINE
# endif
#endif

#define MOD_PARAMS_STACK_dG     (1 << 0)
#define MOD_PARAMS_STACK_dH     (1 << 1)
#define MOD_PARAMS_MISMATCH_dG  (1 << 2)
#define MOD_PARAMS_MISMATCH_dH  (1 << 3)
#define MOD_PARAMS_TERMINAL_dG  (1 << 4)
#define MOD_PARAMS_TERMINAL_dH  (1 << 5)
#define MOD_PARAMS_DANGLES_dG   (1 << 6)
#define MOD_PARAMS_DANGLES_dH   (1 << 7)

#define DEBUG
#define MAX_ALPHABET  (6)
#define MAX_PAIRS     (NBPAIRS + 1 + 25)


/* a container to store the data read from a json parameter file */
typedef struct {
  unsigned int available;

  char    *name;
  char    one_letter_code;
  char    unmodified;
  char    pairing_partners[7];
  unsigned int  pairing_partners_encoding[7];
  unsigned int  unmodified_encoding;

  size_t  num_ptypes;
  size_t  ptypes[MAX_ALPHABET][MAX_ALPHABET];

  int stack_dG[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET];
  int stack_dH[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET];

  int dangle5_dG[MAX_PAIRS][MAX_ALPHABET];
  int dangle5_dH[MAX_PAIRS][MAX_ALPHABET];
  int dangle3_dG[MAX_PAIRS][MAX_ALPHABET];
  int dangle3_dH[MAX_PAIRS][MAX_ALPHABET];

  int mismatch_dG[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET];
  int mismatch_dH[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET];

  int terminal_dG[MAX_PAIRS];
  int terminal_dH[MAX_PAIRS];
} energy_parameters;

/* the actual data structure passed around while evaluating */
typedef struct {
  short   *enc;
  size_t  ptypes[MAX_ALPHABET][MAX_ALPHABET];

  int stack_diff[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET];

  int dangle5_diff[MAX_PAIRS][MAX_ALPHABET];
  int dangle3_diff[MAX_PAIRS][MAX_ALPHABET];

  int mismatch_diff[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET];

  int terminal_diff[MAX_PAIRS];
} energy_corrections;


PRIVATE unsigned int
parse_stacks(JsonNode         *dom,
           const char       *identifier,
           const char       *bases,
           size_t           (*ptypes)[MAX_ALPHABET][MAX_ALPHABET],
           int              (*storage)[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET])
{
  unsigned char num_params = 0;
  char          *enc_ptr;
  size_t        i;
  unsigned int  enc[5] = { 0, 0, 0, 0, 0 };
  JsonNode      *entry, *e;

  /* go through storage and initialize */
  for (size_t i = 0; i < MAX_PAIRS; i++)
    for (size_t k = 0; k < MAX_ALPHABET; k++)
      for (size_t l = 0; l < MAX_ALPHABET; l++)
        (*storage)[i][k][l] = INF;

  if ((e = json_find_member(dom, identifier)) &&
      (e->tag == JSON_OBJECT)) {
    json_foreach(entry, e) {
      if ((entry->key) &&
          (entry->tag == JSON_NUMBER) &&
          (strlen(entry->key) == 4)) {
        /* encode sequence */
        for (i = 0; i < 4; i++) {
          if (!(enc_ptr = strchr(&bases[0], entry->key[i]))) {
            vrna_message_warning("Unrecognized character in \"%s\" base: %s\n", identifier, entry->key);
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
parse_mismatch(JsonNode         *dom,
               const char       *identifier,
               const char       *bases,
               size_t           (*ptypes)[MAX_ALPHABET][MAX_ALPHABET],
               vrna_md_t        *md,
               int              (*storage)[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET])
{
  unsigned char num_params = 0;
  char          *enc_ptr;
  size_t        i;
  unsigned int  enc[5] = { 0, 0, 0, 0, 0 };
  JsonNode      *entry, *e;

  /* go through storage and initialize */
  for (size_t i = 0; i < MAX_PAIRS; i++)
    for (size_t k = 0; k < MAX_ALPHABET; k++)
      for (size_t l = 0; l < MAX_ALPHABET; l++)
        (*storage)[i][k][l] = INF;

  if ((e = json_find_member(dom, identifier)) &&
      (e->tag == JSON_OBJECT)) {
    json_foreach(entry, e) {
      if ((entry->key) &&
          (entry->tag == JSON_NUMBER) &&
          (strlen(entry->key) == 4)) {
        /* encode sequence */
        for (i = 0; i < 4; i++) {
          if (!(enc_ptr = strchr(&bases[0], entry->key[i]))) {
            vrna_message_warning("Unrecognized character in \"%s\" base: %s\n", identifier, entry->key);
            break;
          }
          enc[i] = enc_ptr - &bases[0];
          if (enc[i] > 4)
            enc[i]--;
        }

        if (i == 4) {
          num_params++;

          if ((enc[0] == 5) || (enc[2] == 5))
            (*storage)[NBPAIRS + (*ptypes)[enc[0]][enc[2]]][enc[1]][enc[3]] = (int)(entry->number_ * 100.);
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
parse_dangles(JsonNode         *dom,
              const char       *identifier,
              const char       *bases,
              size_t           (*ptypes)[MAX_ALPHABET][MAX_ALPHABET],
              vrna_md_t        *md,
              int              (*storage)[MAX_PAIRS][MAX_ALPHABET])
{
  unsigned char num_params = 0;
  char          *enc_ptr;
  size_t        i;
  unsigned int  enc[5] = { 0, 0, 0, 0, 0 };
  JsonNode      *entry, *e;

  /* go through storage and initialize */
  for (size_t i = 0; i < MAX_PAIRS; i++)
    for (size_t k = 0; k < MAX_ALPHABET; k++)
      (*storage)[i][k] = INF;

  if ((e = json_find_member(dom, identifier)) &&
      (e->tag == JSON_OBJECT)) {
    json_foreach(entry, e) {
      if ((entry->key) &&
          (entry->tag == JSON_NUMBER) &&
          (strlen(entry->key) == 3)) {
        /* encode sequence */
        for (i = 0; i < 3; i++) {
          if (!(enc_ptr = strchr(&bases[0], entry->key[i]))) {
            vrna_message_warning("Unrecognized character in \"%s\" base: %s\n", identifier, entry->key);
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
parse_terminal(JsonNode         *dom,
               const char       *identifier,
               const char       *bases,
               size_t           (*ptypes)[MAX_ALPHABET][MAX_ALPHABET],
               int              (*storage)[MAX_PAIRS])
{
  unsigned char num_params = 0;
  char          *enc_ptr;
  size_t        i;
  unsigned int  enc[5] = { 0, 0, 0, 0, 0 };
  JsonNode      *entry, *e;

  /* go through storage and initialize */
  for (size_t i = 0; i < MAX_PAIRS; i++)
    (*storage)[i] = INF;

  if ((e = json_find_member(dom, identifier)) &&
      (e->tag == JSON_OBJECT)) {
    json_foreach(entry, e) {
      if ((entry->key) &&
          (entry->tag == JSON_NUMBER) &&
          (strlen(entry->key) == 2)) {
        /* encode sequence */
        for (i = 0; i < 2; i++) {
          if (!(enc_ptr = strchr(&bases[0], entry->key[i]))) {
            vrna_message_warning("Unrecognized character in \"%s\" base: %s\n", identifier, entry->key);
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

PRIVATE energy_parameters *
read_params_from_json_file(const char *filename, vrna_md_t *md) {
  energy_parameters *parameters = NULL;

  FILE *param_file = fopen(filename, "r");
  if (param_file) {
    char *ptr;
    vrna_string_t param_content = vrna_string_make("");

    while ((ptr = vrna_read_line(param_file))) {
      param_content = vrna_string_append_cstring(param_content, ptr);
      free(ptr);
    };
    fclose(param_file);

    if (!json_validate(param_content)) {
      vrna_message_warning("JSON file \"%s\" not valid\n", filename);
      return NULL;
    }

    JsonNode *dom = json_decode(param_content);
    char bases[8] = "_ACGUTM";
    vrna_string_t pairing_partners = vrna_string_make("");

    if (dom) {
      JsonNode *e, *entry, *elem;

      parameters = (energy_parameters *)vrna_alloc(sizeof(energy_parameters));

      parameters->name  = NULL;
      parameters->available = 0;
      parameters->num_ptypes = 0;
      parameters->one_letter_code = '\0';
      parameters->pairing_partners[0] = '\0';
      parameters->unmodified = '\0';

      if ((e = json_find_member(dom, "modified_base")) &&
          (e = json_find_member(e, "name")) &&
          (e->tag == JSON_STRING)) {
        parameters->name = strdup(e->string_);
      }

      /* use one-letter code as specified in json file */
      if ((e = json_find_member(dom, "modified_base")) &&
          (e = json_find_member(e, "one_letter_code")) &&
          (e->tag == JSON_STRING) &&
          (strlen(e->string_) == 1)) {
        parameters->one_letter_code = bases[6] = toupper(e->string_[0]);
      }

      if ((e = json_find_member(dom, "modified_base")) &&
          (e = json_find_member(e, "unmodified")) &&
          (e->tag == JSON_STRING) &&
          (strlen(e->string_) == 1) &&
          (ptr = strchr(bases, e->string_[0]))) {
        parameters->unmodified = toupper(e->string_[0]);
        size_t enc = ptr - &(bases[0]);
        if (enc > 4)
          enc--;
        parameters->unmodified_encoding = enc;
      }

      size_t  cnt = 0;

      if ((e = json_find_member(dom, "modified_base")) &&
          (e = json_find_member(e, "pairing_partners")) &&
          (e->tag == JSON_ARRAY)) {
        json_foreach(entry, e) {
          if ((entry->tag == JSON_STRING) &&
              (strlen(entry->string_) == 1) &&
              (ptr = strchr(bases, entry->string_[0]))) {
            pairing_partners = vrna_string_append_cstring(pairing_partners, entry->string_);

            size_t enc = ptr - &(bases[0]);
            if (enc > 4)
              enc--;

            parameters->ptypes[5][enc] = ++(parameters->num_ptypes);
            parameters->ptypes[enc][5] = ++(parameters->num_ptypes);
            parameters->pairing_partners[cnt]             = entry->string_[0];
            parameters->pairing_partners_encoding[cnt++]  = enc;
          }
        }
      }

      parameters->pairing_partners[cnt] = '\0';

      if (parse_stacks(dom, "stacking_energies", &bases[0], &(parameters->ptypes), &(parameters->stack_dG)) > 0)
        parameters->available |= MOD_PARAMS_STACK_dG;

      if (parse_stacks(dom, "stacking_enthalpies", &bases[0], &(parameters->ptypes), &(parameters->stack_dH)) > 0)
        parameters->available |= MOD_PARAMS_STACK_dH;

      if (parse_mismatch(dom, "mismatch_energies", &bases[0], &(parameters->ptypes), md, &(parameters->mismatch_dG)) > 0)
        parameters->available |= MOD_PARAMS_MISMATCH_dG;

      if (parse_mismatch(dom, "mismatch_enthalpies", &bases[0], &(parameters->ptypes), md, &(parameters->mismatch_dH)) > 0)
        parameters->available |= MOD_PARAMS_MISMATCH_dH;

      if (parse_terminal(dom, "terminal_energies", &bases[0], &(parameters->ptypes), &(parameters->terminal_dG)) > 0)
        parameters->available |= MOD_PARAMS_TERMINAL_dG;

      if (parse_terminal(dom, "terminal_enthalpies", &bases[0], &(parameters->ptypes), &(parameters->terminal_dH)) > 0)
        parameters->available |= MOD_PARAMS_TERMINAL_dH;

      if (parse_dangles(dom, "dangle5_energies", &bases[0], &(parameters->ptypes), md, &(parameters->dangle5_dG)) > 0)
        parameters->available |= MOD_PARAMS_DANGLES_dG;

      if (parse_dangles(dom, "dangle5_enthalpies", &bases[0], &(parameters->ptypes), md, &(parameters->dangle5_dH)) > 0)
        parameters->available |= MOD_PARAMS_DANGLES_dH;

      if (parse_dangles(dom, "dangle3_energies", &bases[0], &(parameters->ptypes), md, &(parameters->dangle3_dG)) > 0)
        parameters->available |= MOD_PARAMS_DANGLES_dG;

      if (parse_dangles(dom, "dangle3_enthalpies", &bases[0], &(parameters->ptypes), md, &(parameters->dangle3_dH)) > 0)
        parameters->available |= MOD_PARAMS_DANGLES_dH;

      json_delete(dom);
    }
    vrna_string_free(param_content);
  }

  return parameters;
}


PRIVATE INLINE void
init_stacks(energy_parameters   *params,
            energy_corrections  *diffs,
            vrna_param_t        *P)
{
  unsigned int  i, si, sj, enc_unmod, enc_pp, tt, pair_MP, pair_PM;
  vrna_md_t     *md     = &(P->model_details);
  double        tempf   = (md->temperature + K0) / (37. + K0);
  char          nt[MAX_ALPHABET] = {'\0', 'A', 'C','G','U', 'M'};
  int           e, (*dG)[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET], (*dH)[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET];

  enc_unmod = params->unmodified_encoding;
  dG        = &(params->stack_dG);
  dH        = &(params->stack_dH);

  nt[5] = params->one_letter_code;

  if (params->available & MOD_PARAMS_STACK_dG) {
    for (i = 1; i <= params->num_ptypes; i += 2) {
      enc_pp = params->pairing_partners_encoding[(i - 1)/ 2];
      /* pair type of unmodified version as encoded in RNAlib */
      pair_MP = md->pair[enc_unmod][enc_pp];
      pair_PM = md->pair[enc_pp][enc_unmod];

      if (pair_MP == 0)
        pair_MP = 7;
      if (pair_PM == 0)
        pair_PM = 7;

      for (si = 1; si < MAX_ALPHABET; si++) {
        for (sj = 1; sj < MAX_ALPHABET; sj++) {
          if (si == 5) {
            if (sj == 5) {
              tt = md->pair[enc_unmod][enc_unmod];
            } else {
              tt = md->pair[sj][enc_unmod];
            }
          } else if (sj == 5) {
            tt = md->pair[enc_unmod][si];
          } else {
            tt = md->pair[sj][si];
          }

          if (tt == 0)
            tt = 7;

          if ((*dG)[i][sj][si] != INF) {
            if (params->available & MOD_PARAMS_STACK_dH)
              e = (*dH)[i][sj][si] - ((*dH)[i][sj][si] - (*dG)[i][sj][si]) * tempf;
            else
              e = (*dG)[i][sj][si];

            diffs->stack_diff[i][sj][si]      = e - P->stack[pair_MP][tt];
#ifdef DEBUG
            printf("d_stack(%c%c, %c%c) = %d = %d - %d\n", nt[5], nt[enc_pp], nt[sj], nt[si],  diffs->stack_diff[i][sj][si], e, P->stack[pair_MP][tt]);
#endif
          }

          if ((*dG)[i + 1][sj][si] != INF) {
            if (params->available & MOD_PARAMS_STACK_dH)
              e = (*dH)[i + 1][sj][si] - ((*dH)[i + 1][sj][si] - (*dG)[i + 1][sj][si]) * tempf;
            else
              e = (*dG)[i + 1][sj][si];

            diffs->stack_diff[i + 1][sj][si]  = e - P->stack[pair_PM][tt];
#ifdef DEBUG
            printf("d_stack(%c%c, %c%c) = %d = %d - %d\n", nt[enc_pp], nt[5], nt[sj], nt[si], diffs->stack_diff[i + 1][sj][si], e, P->stack[pair_PM][tt]);
#endif
          }
        }
      }
    }
  }
}

PRIVATE INLINE void
init_mismatches(energy_parameters   *params,
                energy_corrections  *diffs,
                vrna_param_t        *P)
{
  unsigned int  i, si, sj, enc_unmod, enc_pp, siu, sju, pair_MP, pair_PM;
  vrna_md_t     *md     = &(P->model_details);
  double        tempf   = (md->temperature + K0) / (37. + K0);
  char          nt[MAX_ALPHABET] = {'\0', 'A', 'C','G','U', 'M'};
  int           e, (*dG)[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET], (*dH)[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET];

#ifdef DEBUG
  char          bp[3] = {0}, *bpairs[7] = {"NN", "CG", "GC", "GU", "UG", "AU", "UA"};
#endif

  enc_unmod = params->unmodified_encoding;
  dG        = &(params->mismatch_dG);
  dH        = &(params->mismatch_dH);

  nt[5] = params->one_letter_code;

  if (params->available & MOD_PARAMS_MISMATCH_dG) {
    /* process all closing pairs without modified bases */
    for (i = 1; i <= params->num_ptypes + NBPAIRS; i += 2) {
      if (i > NBPAIRS) {
        /* closing pair consists of at least one modified base */
        enc_pp = params->pairing_partners_encoding[(i - NBPAIRS - 1)/ 2];
        /* pair type of unmodified version as encoded in RNAlib */
        pair_MP = md->pair[enc_unmod][enc_pp];
        pair_PM = md->pair[enc_pp][enc_unmod];

#ifdef DEBUG
        bp[0] = nt[5];
        bp[1] = nt[enc_pp];
#endif

      } else {
        /* closing pair is regular base pair without modified bases */
        pair_MP = i;
        pair_PM = i + 1;
#ifdef DEBUG
        bp[0] = bpairs[i][0];
        bp[1] = bpairs[i][1];
#endif
      }

      if (pair_MP == 0)
        pair_MP = 7;
      if (pair_PM == 0)
        pair_PM = 7;

      for (si = 1; si < MAX_ALPHABET; si++) {
        for (sj = 1; sj < MAX_ALPHABET; sj++) {
          if (si == 5)
            siu = enc_unmod;
          else
            siu = si;
          if (sj == 5)
            sju = enc_unmod;
          else
            sju = sj;

          if ((*dG)[i][si][sj] != INF) {
            if (params->available & MOD_PARAMS_MISMATCH_dH)
              e = (*dH)[i][si][sj] - ((*dH)[i][si][sj] - (*dG)[i][si][sj]) * tempf;
            else
              e = (*dG)[i][si][sj];

            /*
             * take mismatch energies from multiloops as they do not contain anything
             * but mismatch contributions. Also note, that multiloop mismatches are
             * stored such that unpaired bases are outside of pair, in contrast to
             * considering the unpaired bases enclosed. Thus, the pair and 5' and 3'
             * unpaired bases must be rotatated
             */
            diffs->mismatch_diff[i][si][sj] = e - P->mismatchM[pair_PM][sju][siu];
#ifdef DEBUG
            printf("d_mm(%c%c, %c, %c) = %d = %d - %d\n", bp[0], bp[1], nt[si], nt[sj],  diffs->mismatch_diff[i][si][sj], e, P->mismatchM[pair_PM][sju][siu]);
#endif
          }

          if ((*dG)[i + 1][si][sj] != INF) {
            if (params->available & MOD_PARAMS_MISMATCH_dH)
              e = (*dH)[i + 1][si][sj] - ((*dH)[i + 1][si][sj] - (*dG)[i + 1][si][sj]) * tempf;
            else
              e = (*dG)[i + 1][si][sj];

            /*
             * take mismatch energies from multiloops as they do not contain anything
             * but mismatch contributions. Also note, that multiloop mismatches are
             * stored such that unpaired bases are outside of pair, in contrast to
             * considering the unpaired bases enclosed. Thus, the pair and 5' and 3'
             * unpaired bases must be rotatated
             */
            diffs->mismatch_diff[i + 1][si][sj]  = e - P->mismatchM[pair_MP][sju][siu];
#ifdef DEBUG
            printf("d_mm(%c%c, %c, %c) = %d = %d - %d\n", bp[1], bp[0], nt[si], nt[sj], diffs->mismatch_diff[i + 1][si][sj], e, P->mismatchM[pair_MP][sju][siu]);
#endif
          }
        }
      }
    }
  }
}

PRIVATE INLINE void
init_dangles(energy_parameters   *params,
             energy_corrections  *diffs,
             vrna_param_t        *P)
{
  unsigned int  i, si, sj, enc_unmod, enc_pp, siu, sju, pair_MP, pair_PM;
  vrna_md_t     *md     = &(P->model_details);
  double        tempf   = (md->temperature + K0) / (37. + K0);
  char          nt[MAX_ALPHABET] = {'\0', 'A', 'C','G','U', 'M'};
  int           e, (*dG5)[MAX_PAIRS][MAX_ALPHABET], (*dH5)[MAX_PAIRS][MAX_ALPHABET], (*dG3)[MAX_PAIRS][MAX_ALPHABET], (*dH3)[MAX_PAIRS][MAX_ALPHABET];

#ifdef DEBUG
  char          bp[3] = {0}, *bpairs[7] = {"NN", "CG", "GC", "GU", "UG", "AU", "UA"};
#endif

  enc_unmod = params->unmodified_encoding;
  dG5       = &(params->dangle5_dG);
  dH5       = &(params->dangle5_dH);
  dG3       = &(params->dangle3_dG);
  dH3       = &(params->dangle3_dH);

  nt[5] = params->one_letter_code;

  if (params->available & MOD_PARAMS_DANGLES_dG) {
    /* process all closing pairs without modified bases */
    for (i = 1; i <= params->num_ptypes + NBPAIRS; i += 2) {
      if (i > NBPAIRS) {
        /* closing pair consists of at least one modified base */
        enc_pp = params->pairing_partners_encoding[(i - NBPAIRS - 1)/ 2];
        /* pair type of unmodified version as encoded in RNAlib */
        pair_MP = md->pair[enc_unmod][enc_pp];
        pair_PM = md->pair[enc_pp][enc_unmod];

#ifdef DEBUG
        bp[0] = nt[5];
        bp[1] = nt[enc_pp];
#endif

      } else {
        /* closing pair is regular base pair without modified bases */
        pair_MP = i;
        pair_PM = i + 1;
#ifdef DEBUG
        bp[0] = bpairs[i][0];
        bp[1] = bpairs[i][1];
#endif
      }

      if (pair_MP == 0)
        pair_MP = 7;
      if (pair_PM == 0)
        pair_PM = 7;

      for (si = 1; si < MAX_ALPHABET; si++) {
        if (si == 5)
          siu = enc_unmod;
        else
          siu = si;

        if ((*dG5)[i][si] != INF) {
          if (params->available & MOD_PARAMS_DANGLES_dH)
            e = (*dH5)[i][si] - ((*dH5)[i][si] - (*dG5)[i][si]) * tempf;
          else
            e = (*dG5)[i][si];

          /*
           * Note, that dangling ends arestored such that unpaired bases are
           * outside of pair, in contrast to considering the unpaired bases
           * enclosed. Thus, the pair must be rotatated
           */
          diffs->dangle5_diff[i][si] = e - P->dangle5[pair_PM][siu];
#ifdef DEBUG
          printf("d_d5(%c%c, %c) = %d = %d - %d\n", bp[0], bp[1], nt[si],  diffs->dangle5_diff[i][si], e, P->dangle5[pair_PM][siu]);
#endif
        }

        if ((*dG3)[i][si] != INF) {
          if (params->available & MOD_PARAMS_DANGLES_dH)
            e = (*dH3)[i][si] - ((*dH3)[i][si] - (*dG3)[i][si]) * tempf;
          else
            e = (*dG3)[i][si];

          /*
           * Note, that dangling ends arestored such that unpaired bases are
           * outside of pair, in contrast to considering the unpaired bases
           * enclosed. Thus, the pair must be rotatated
           */
          diffs->dangle3_diff[i][si] = e - P->dangle3[pair_PM][siu];
#ifdef DEBUG
          printf("d_d3(%c%c, %c) = %d = %d - %d\n", bp[0], bp[1], nt[si],  diffs->dangle3_diff[i][si], e, P->dangle3[pair_PM][siu]);
#endif
        }

        if ((*dG5)[i + 1][si] != INF) {
          if (params->available & MOD_PARAMS_DANGLES_dH)
            e = (*dH5)[i + 1][si] - ((*dH5)[i + 1][si] - (*dG5)[i + 1][si]) * tempf;
          else
            e = (*dG5)[i + 1][si];

          /*
           * Note, that dangling ends arestored such that unpaired bases are
           * outside of pair, in contrast to considering the unpaired bases
           * enclosed. Thus, the pair must be rotatated
           */
          diffs->dangle5_diff[i + 1][si] = e - P->dangle5[pair_MP][siu];
#ifdef DEBUG
          printf("d_d5(%c%c, %c) = %d = %d - %d\n", bp[1], bp[0], nt[si],  diffs->dangle5_diff[i + 1][si], e, P->dangle5[pair_MP][siu]);
#endif
        }

        if ((*dG3)[i + 1][si] != INF) {
          if (params->available & MOD_PARAMS_DANGLES_dH)
            e = (*dH3)[i + 1][si] - ((*dH3)[i + 1][si] - (*dG3)[i + 1][si]) * tempf;
          else
            e = (*dG3)[i + 1][si];

          /*
           * Note, that dangling ends arestored such that unpaired bases are
           * outside of pair, in contrast to considering the unpaired bases
           * enclosed. Thus, the pair must be rotatated
           */
          diffs->dangle3_diff[i + 1][si] = e - P->dangle3[pair_MP][siu];
#ifdef DEBUG
          printf("d_d3(%c%c, %c) = %d = %d - %d\n", bp[1], bp[0], nt[si],  diffs->dangle3_diff[i + 1][si], e, P->dangle3[pair_MP][siu]);
#endif
        }
      }
    }
  }
}

PRIVATE INLINE void
init_terminal(energy_parameters   *params,
              energy_corrections  *diffs,
              vrna_param_t        *P)
{
  unsigned int  i, si, sj, enc_unmod, enc_pp, tt, pair_MP, pair_PM;
  vrna_md_t     *md     = &(P->model_details);
  double        tempf   = (md->temperature + K0) / (37. + K0);
  char          nt[MAX_ALPHABET] = {'\0', 'A', 'C','G','U', 'M'};
  int           e, (*dG)[MAX_PAIRS], (*dH)[MAX_PAIRS], Terminal_unmod;

  enc_unmod = params->unmodified_encoding;
  dG        = &(params->terminal_dG);
  dH        = &(params->terminal_dH);

  nt[5] = params->one_letter_code;

  if (params->available & MOD_PARAMS_TERMINAL_dG) {
    for (i = 1; i <= params->num_ptypes; i += 2) {
      enc_pp = params->pairing_partners_encoding[(i - 1)/ 2];
      tt     = md->pair[enc_unmod][enc_pp];
      Terminal_unmod = (tt > 2) ? P->TerminalAU : 0; /* terminal AU only for AU and GU pairs */

      if ((*dG)[i] != INF) {
        if (params->available & MOD_PARAMS_TERMINAL_dH)
          e = (*dH)[i] - ((*dH)[i] - (*dG)[i]) * tempf;
        else
          e = (*dG)[i];

        diffs->terminal_diff[i] = e - Terminal_unmod;
#ifdef DEBUG
        printf("d_term(%c%c) = %d = %d - %d\n", nt[5], nt[enc_pp], diffs->terminal_diff[i], e, Terminal_unmod);
#endif
      }

      if ((*dG)[i + 1] != INF) {
        if (params->available & MOD_PARAMS_TERMINAL_dH)
          e = (*dH)[i + 1] - ((*dH)[i + 1] - (*dG)[i + 1]) * tempf;
        else
          e = (*dG)[i + 1];

        diffs->terminal_diff[i + 1] = e - Terminal_unmod;
#ifdef DEBUG
        printf("d_term(%c%c) = %d = %d - %d\n", nt[enc_pp], nt[5], diffs->terminal_diff[i + 1], e, Terminal_unmod);
#endif
      }
    }
  }
}


PRIVATE INLINE int
mismatch(vrna_fold_compound_t *fc,
         unsigned int i,
         unsigned int j,
         energy_corrections *data)
{
  short *enc      = data->enc;
  unsigned int tt = data->ptypes[enc[i]][enc[j]];
  vrna_md_t *md   = &(fc->params->model_details);

  if (tt == 0)
    /* if we don't know the base pair, it must be canonical */
    tt = md->pair[enc[i]][enc[j]];
  else
    tt += NBPAIRS;

  if (j > 1) {
    if (i < fc->length)
      return data->mismatch_diff[tt][enc[i + 1]][enc[j - 1]];
    else
      return data->dangle5_diff[tt][enc[j - 1]];
  } else if (i < fc->length) {
    return data->dangle3_diff[tt][enc[i + 1]];
  }

  return 0;
}


PRIVATE INLINE int
terminal(unsigned int i,
         unsigned int j,
         energy_corrections *data)
{
  short *enc      = data->enc;
  unsigned int tt = data->ptypes[enc[i]][enc[j]];

  return data->terminal_diff[tt];
}


/* hairpin loop correction including terminalAU terms */
PRIVATE INLINE int
sc_PAIR_HP_terminal(vrna_fold_compound_t  *fc,
                int                   i,
      int                   j,
      int                   k,
      int                   l,
      void                  *d)
{
  return terminal(i, j,
                  (energy_corrections *)d);
}

/* hairpin loop correction including mismatch terms */
PRIVATE INLINE int
sc_PAIR_HP_mismatch(vrna_fold_compound_t  *fc,
            int                   i,
      int                   j,
      int                   k,
      int                   l,
      void                  *d)
{
  return mismatch(fc,
                  i, j,
                  (energy_corrections *)d);
}

/* hairpin loop correction including mismatch and terminalAU terms */
PRIVATE int
sc_PAIR_HP(vrna_fold_compound_t  *fc,
      int                   i,
      int                   j,
      int                   k,
      int                   l,
      void                  *d)
{
  return sc_PAIR_HP_mismatch(fc, i, j, k, l, d) +
         sc_PAIR_HP_terminal(fc, i, j, k, l, d);
}


PRIVATE INLINE int
sc_PAIR_IL_stack(vrna_fold_compound_t *fc,
                 int  i,
                 int  j,
                 int  k,
                 int  l,
                 void *d)
{
  if ((i + 1 == k) &&
      (l == j - 1)) {
    energy_corrections *data = (energy_corrections *)d;
    short *enc = data->enc;

    unsigned int enc_i = enc[i];
    unsigned int enc_j = enc[j];
    unsigned int enc_k = enc[k];
    unsigned int enc_l = enc[l];
    unsigned int t1 = data->ptypes[enc_i][enc_j];
    unsigned int t2 = data->ptypes[enc_l][enc_k];
    if (t1 != 0)
      return data->stack_diff[t1][enc_l][enc_k];
    else if (t2 != 0)
      return data->stack_diff[t2][enc_i][enc_j];
  }

  return 0;
}

PRIVATE INLINE int
sc_PAIR_IL_terminal(vrna_fold_compound_t *fc,
                    int  i,
                    int  j,
                    int  k,
                    int  l,
                    void *d)
{
  if ((i + 1 < k) ||
      (l + 1 < j))
    return terminal(i, j, (energy_corrections *)d) +
           terminal(l, k, (energy_corrections *)d);

  return 0;
}

PRIVATE INLINE int
sc_PAIR_IL_mismatch(vrna_fold_compound_t *fc,
                    int  i,
                    int  j,
                    int  k,
                    int  l,
                    void *d)
{
  if (((k - i - 1) > 2) &&
      ((j - l - 1) > 2))
    return mismatch(fc, i, j, (energy_corrections *)d) +
           mismatch(fc, l, k, (energy_corrections *)d);

  return 0;
}

PRIVATE INLINE int
sc_PAIR_IL_mismatch_terminal(vrna_fold_compound_t *fc,
                             int  i,
                             int  j,
                             int  k,
                             int  l,
                             void *d)
{
  return sc_PAIR_IL_mismatch(fc, i, j, k, l, d) +
         sc_PAIR_IL_terminal(fc, i, j, k, l, d);
}

PRIVATE INLINE int
sc_PAIR_IL_stack_terminal(vrna_fold_compound_t *fc,
                          int  i,
                          int  j,
                          int  k,
                          int  l,
                          void *d)
{
  return sc_PAIR_IL_stack(fc, i, j, k, l, d) +
         sc_PAIR_IL_terminal(fc, i, j, k, l, d);
}

PRIVATE INLINE int
sc_PAIR_IL_stack_mismatch(vrna_fold_compound_t *fc,
                          int  i,
                          int  j,
                          int  k,
                          int  l,
                          void *d)
{
  return sc_PAIR_IL_stack(fc, i, j, k, l, d) +
         sc_PAIR_IL_mismatch(fc, i, j, k, l, d);
}

PRIVATE INLINE int
sc_PAIR_IL(vrna_fold_compound_t *fc,
               int  i,
               int  j,
               int  k,
               int  l,
               void *d)
{
  return sc_PAIR_IL_stack(fc, i, j, k, l, d) +
         sc_PAIR_IL_mismatch(fc, i, j, k, l, d) +
         sc_PAIR_IL_terminal(fc, i, j, k, l, d);
}

PRIVATE INLINE int
sc_PAIR_ML_terminal(vrna_fold_compound_t  *fc,
      int                   i,
      int                   j,
      int                   k,
      int                   l,
      void                  *d)
{
  return terminal(i, j,
                  (energy_corrections *)d);
}


PRIVATE INLINE int
sc_PAIR_ML_mismatch(vrna_fold_compound_t  *fc,
      int                   i,
      int                   j,
      int                   k,
      int                   l,
      void                  *d)
{
  return mismatch(fc,
                  i, j,
                  (energy_corrections *)d);
}


PRIVATE INLINE int
sc_PAIR_ML(vrna_fold_compound_t  *fc,
      int                   i,
      int                   j,
      int                   k,
      int                   l,
      void                  *d)
{
  return sc_PAIR_ML_mismatch(fc, i, j, k, l, d) +
         sc_PAIR_ML_terminal(fc, i, j, k, l, d);
}


PRIVATE INLINE int
sc_STEM_terminal(vrna_fold_compound_t  *fc,
            int                   i,
            int                   j,
            int                   k,
            int                   l,
            void                  *d)
{
  return terminal(l, k,
                  (energy_corrections *)d);
}


PRIVATE INLINE int
sc_STEM_mismatch(vrna_fold_compound_t  *fc,
            int                   i,
            int                   j,
            int                   k,
            int                   l,
            void                  *d)
{
  return mismatch(fc,
                  l, k,
                  (energy_corrections *)d);
}


PRIVATE INLINE int
sc_STEM(vrna_fold_compound_t  *fc,
            int                   i,
            int                   j,
            int                   k,
            int                   l,
            void                  *d)
{
  return sc_STEM_mismatch(fc, i, j, k, l, d) +
         sc_STEM_terminal(fc, i, j, k, l, d);
}


PRIVATE INLINE int
sc_EXT_STEM_EXT_terminal(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j,
                int                   k,
                int                   l,
                void                  *d)
{
  return terminal(k, i,
                  (energy_corrections *)d);
}


PRIVATE INLINE int
sc_EXT_STEM_EXT_mismatch(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j,
                int                   k,
                int                   l,
                void                  *d)
{
  return mismatch(fc,
                  k, i,
                  (energy_corrections *)d);
}


PRIVATE INLINE int
sc_EXT_STEM_EXT(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j,
                int                   k,
                int                   l,
                void                  *d)
{
  return sc_EXT_STEM_EXT_mismatch(fc, i, j, k, l, d) +
         sc_EXT_STEM_EXT_terminal(fc, i, j, k, l, d);
}


PRIVATE INLINE int
sc_EXT_EXT_STEM_terminal(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j,
                int                   k,
                int                   l,
                void                  *d)
{
  return terminal(j, l,
                  (energy_corrections *)d);
}


PRIVATE INLINE int
sc_EXT_EXT_STEM_mismatch(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j,
                int                   k,
                int                   l,
                void                  *d)
{
  return mismatch(fc,
                  j, l,
                  (energy_corrections *)d);
}


PRIVATE INLINE int
sc_EXT_EXT_STEM(vrna_fold_compound_t  *fc,
                int                   i,
                int                   j,
                int                   k,
                int                   l,
                void                  *d)
{
  return sc_EXT_EXT_STEM_mismatch(fc, i, j, k, l, d) +
         sc_EXT_EXT_STEM_terminal(fc, i, j, k, l, d);
}


PRIVATE INLINE int
sc_EXT_STEM_OUTSIDE_terminal(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   k,
                    int                   l,
                    void                  *d)
{
  return terminal(l, k,
                  (energy_corrections *)d);
}


PRIVATE INLINE int
sc_EXT_STEM_OUTSIDE_mismatch(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   k,
                    int                   l,
                    void                  *d)
{
  return mismatch(fc,
                  l, k,
                  (energy_corrections *)d);
}


PRIVATE INLINE int
sc_EXT_STEM_OUTSIDE(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   k,
                    int                   l,
                    void                  *d)
{
  return sc_EXT_STEM_OUTSIDE_mismatch(fc, i, j, k, l, d) +
         sc_EXT_STEM_OUTSIDE_terminal(fc, i, j, k, l, d);
}


PRIVATE INLINE int
sc_ML_ML_STEM_terminal(vrna_fold_compound_t  *fc,
              int                   i,
              int                   j,
              int                   k,
              int                   l,
              void                  *d)
{
  return terminal(j, l,
                  (energy_corrections *)d);
}


PRIVATE INLINE int
sc_ML_ML_STEM_mismatch(vrna_fold_compound_t  *fc,
              int                   i,
              int                   j,
              int                   k,
              int                   l,
              void                  *d)
{
  return mismatch(fc,
                  j, l,
                  (energy_corrections *)d);
}


PRIVATE INLINE int
sc_ML_ML_STEM(vrna_fold_compound_t  *fc,
              int                   i,
              int                   j,
              int                   k,
              int                   l,
              void                  *d)
{
  return sc_ML_ML_STEM_mismatch(fc, i, j, k, l, d) +
         sc_ML_ML_STEM_terminal(fc, i, j, k, l, d);
}


PRIVATE void
free_energy_corrections(void *d)
{
  energy_corrections *diff = (energy_corrections *)d;
  free(diff->enc);
  free(diff);
}


PUBLIC void
vrna_sc_mod_json(vrna_fold_compound_t *fc,
                 const char           *json_file,
                 unsigned int         *modification_sites) {

  if ((fc) &&
      (json_file) &&
      (modification_sites)) {
    energy_parameters *params = read_params_from_json_file(json_file, &(fc->params->model_details));

    if (params) {
      vrna_md_t *md = &(fc->params->model_details);
      char bases[8] = "_ACGUTM";
      bases[6] = params->one_letter_code;

      energy_corrections *diffs = (energy_corrections *)vrna_alloc(sizeof(energy_corrections));

      /* copy ptypes */
      memcpy(&(diffs->ptypes[0][0]), &(params->ptypes[0][0]), sizeof(params->ptypes));

      diffs->enc = (short *)vrna_alloc(sizeof(short) * (fc->length + 2));
      memcpy(diffs->enc, fc->sequence_encoding, sizeof(short) * (fc->length + 1));

      for (size_t i = 0; modification_sites[i]; i++) {
        unsigned int msite = modification_sites[i];
        if (msite > fc->length) {
          vrna_message_warning("modification site %u after sequence length (%u)",
                               msite,
                               fc->length);
          continue;
        }

        if (fc->sequence_encoding[msite] != params->unmodified_encoding) {
          vrna_message_warning("modification site %u lists wrong unmodified base %c (should be %c)",
                               msite,
                               bases[fc->sequence_encoding[msite]],
                               params->unmodified);
          continue;
        }

        diffs->enc[msite] = 5;

        unsigned int *sn = fc->strand_number;

        /* allow for all pairing partners specified in the input */
        for (unsigned int j = 1; j < msite; j++) {
          if ((sn[msite] != sn[j]) ||
              ((msite - j - 1) >= md->min_loop_size))
            for (unsigned int cnt = 0; cnt < params->num_ptypes / 2; cnt++) {
              unsigned int pp_enc = params->pairing_partners_encoding[cnt];
              if (fc->sequence_encoding[j] == pp_enc) {
                vrna_hc_add_bp(fc,
                               j,
                               msite,
                               VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS | VRNA_CONSTRAINT_CONTEXT_NO_REMOVE);
              }
            }
        }
        for (unsigned int j = msite + 1; j <= fc->length; j++) {
          if ((sn[msite] != sn[j]) ||
              ((j - msite - 1) >= md->min_loop_size))
            for (unsigned int cnt = 0; cnt < params->num_ptypes / 2; cnt++) {
              unsigned int pp_enc = params->pairing_partners_encoding[cnt];
              if (fc->sequence_encoding[j] == pp_enc) {
                vrna_hc_add_bp(fc,
                               msite,
                               j,
                               VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS | VRNA_CONSTRAINT_CONTEXT_NO_REMOVE);
              }
            }
        }
      }

      init_stacks(params, diffs, fc->params);
      init_terminal(params, diffs, fc->params);
      init_mismatches(params, diffs, fc->params);
      init_dangles(params, diffs, fc->params);

      unsigned int available = params->available;

      /* bind callbacks depending on what data is provided for this modified base */
      if (available & MOD_PARAMS_TERMINAL_dG) {
        if (available & MOD_PARAMS_MISMATCH_dG) {
          vrna_sc_multi_cb_add(fc,
                               &sc_PAIR_HP,
                               NULL,
                               (void *)diffs,
                               &free_energy_corrections,
                               VRNA_DECOMP_PAIR_HP);

          vrna_sc_multi_cb_add(fc,
                               (available & MOD_PARAMS_STACK_dG) ? &sc_PAIR_IL : sc_PAIR_IL_mismatch_terminal,
                               NULL,
                               (void *)diffs,
                               NULL,
                               VRNA_DECOMP_PAIR_IL);

          vrna_sc_multi_cb_add(fc,
                               &sc_PAIR_ML,
                               NULL,
                               (void *)diffs,
                               NULL,
                               VRNA_DECOMP_PAIR_ML);

          vrna_sc_multi_cb_add(fc,
                               &sc_STEM,
                               NULL,
                               (void *)diffs,
                               NULL,
                               VRNA_DECOMP_EXT_STEM);

          vrna_sc_multi_cb_add(fc,
                               &sc_EXT_STEM_EXT,
                               NULL,
                               (void *)diffs,
                               NULL,
                               VRNA_DECOMP_EXT_STEM_EXT);

          vrna_sc_multi_cb_add(fc,
                               &sc_EXT_EXT_STEM,
                               NULL,
                               (void *)diffs,
                               NULL,
                               VRNA_DECOMP_EXT_EXT_STEM);

          vrna_sc_multi_cb_add(fc,
                               &sc_EXT_STEM_OUTSIDE,
                               NULL,
                               (void *)diffs,
                               NULL,
                               VRNA_DECOMP_EXT_STEM_OUTSIDE);

          vrna_sc_multi_cb_add(fc,
                               &sc_STEM,
                               NULL,
                               (void *)diffs,
                               NULL,
                               VRNA_DECOMP_ML_STEM);

          vrna_sc_multi_cb_add(fc,
                               &sc_ML_ML_STEM,
                               NULL,
                               (void *)diffs,
                               NULL,
                               VRNA_DECOMP_ML_ML_STEM);
        } else {
          /* no mismatch energies available */
          vrna_sc_multi_cb_add(fc,
                               &sc_PAIR_HP_terminal,
                               NULL,
                               (void *)diffs,
                               &free_energy_corrections,
                               VRNA_DECOMP_PAIR_HP);

          vrna_sc_multi_cb_add(fc,
                               (available & MOD_PARAMS_STACK_dG) ? &sc_PAIR_IL_stack_terminal : &sc_PAIR_IL_terminal,
                               NULL,
                               (void *)diffs,
                               NULL,
                               VRNA_DECOMP_PAIR_IL);

          vrna_sc_multi_cb_add(fc,
                               &sc_PAIR_ML_terminal,
                               NULL,
                               (void *)diffs,
                               NULL,
                               VRNA_DECOMP_PAIR_ML);

          vrna_sc_multi_cb_add(fc,
                               &sc_STEM_terminal,
                               NULL,
                               (void *)diffs,
                               NULL,
                               VRNA_DECOMP_EXT_STEM);

          vrna_sc_multi_cb_add(fc,
                               &sc_EXT_STEM_EXT_terminal,
                               NULL,
                               (void *)diffs,
                               NULL,
                               VRNA_DECOMP_EXT_STEM_EXT);

          vrna_sc_multi_cb_add(fc,
                               &sc_EXT_EXT_STEM_terminal,
                               NULL,
                               (void *)diffs,
                               NULL,
                               VRNA_DECOMP_EXT_EXT_STEM);

          vrna_sc_multi_cb_add(fc,
                               &sc_EXT_STEM_OUTSIDE_terminal,
                               NULL,
                               (void *)diffs,
                               NULL,
                               VRNA_DECOMP_EXT_STEM_OUTSIDE);

          vrna_sc_multi_cb_add(fc,
                               &sc_STEM_terminal,
                               NULL,
                               (void *)diffs,
                               NULL,
                               VRNA_DECOMP_ML_STEM);

          vrna_sc_multi_cb_add(fc,
                               &sc_ML_ML_STEM_terminal,
                               NULL,
                               (void *)diffs,
                               NULL,
                               VRNA_DECOMP_ML_ML_STEM);
        }
      } else if (available & MOD_PARAMS_MISMATCH_dG) {
        /* no terminalAU-like parameters available */
        vrna_sc_multi_cb_add(fc,
                             &sc_PAIR_HP_mismatch,
                             NULL,
                             (void *)diffs, &free_energy_corrections,
                             VRNA_DECOMP_PAIR_HP);

        vrna_sc_multi_cb_add(fc,
                             (available & MOD_PARAMS_STACK_dG) ? &sc_PAIR_IL_stack_mismatch : &sc_PAIR_IL_mismatch,
                             NULL,
                             (void *)diffs,
                             NULL,
                             VRNA_DECOMP_PAIR_IL);

        vrna_sc_multi_cb_add(fc,
                             &sc_PAIR_ML_mismatch,
                             NULL,
                             (void *)diffs,
                             NULL,
                             VRNA_DECOMP_PAIR_ML);

        vrna_sc_multi_cb_add(fc,
                             &sc_STEM_mismatch,
                             NULL,
                             (void *)diffs,
                             NULL,
                             VRNA_DECOMP_EXT_STEM);

        vrna_sc_multi_cb_add(fc,
                             &sc_EXT_STEM_EXT_mismatch,
                             NULL,
                             (void *)diffs,
                             NULL,
                             VRNA_DECOMP_EXT_STEM_EXT);

        vrna_sc_multi_cb_add(fc,
                             &sc_EXT_EXT_STEM_mismatch,
                             NULL,
                             (void *)diffs,
                             NULL,
                             VRNA_DECOMP_EXT_EXT_STEM);

        vrna_sc_multi_cb_add(fc,
                             &sc_EXT_STEM_OUTSIDE_mismatch,
                             NULL,
                             (void *)diffs,
                             NULL,
                             VRNA_DECOMP_EXT_STEM_OUTSIDE);

        vrna_sc_multi_cb_add(fc,
                             &sc_STEM_mismatch,
                             NULL,
                             (void *)diffs,
                             NULL,
                             VRNA_DECOMP_ML_STEM);

        vrna_sc_multi_cb_add(fc,
                             &sc_ML_ML_STEM_mismatch,
                             NULL,
                             (void *)diffs,
                             NULL,
                             VRNA_DECOMP_ML_ML_STEM);
      } else if (available & MOD_PARAMS_STACK_dG) {
        /* just stacking parameters available */
        vrna_sc_multi_cb_add(fc,
                             &sc_PAIR_IL_stack,
                             NULL,
                             (void *)diffs,
                             &free_energy_corrections,
                             VRNA_DECOMP_PAIR_IL);
      }

    }
  }
}
