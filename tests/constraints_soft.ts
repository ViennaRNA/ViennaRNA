#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/string_utils.h>
#include <ViennaRNA/constraints_soft.h>

#suite Constraints

#tcase  SoftConstraints

#test test_vrna_sc_add_up_simple

  int i, j;

  char *seq = vrna_random_string(10, "ACGU");

  vrna_fold_compound_t *fc = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);

  for(i = 1; i <= fc->length; i++)
    vrna_sc_add_up(fc, i, -1. * i, VRNA_OPTION_DEFAULT);

  vrna_sc_t *sc = fc->sc;

  ck_assert(sc != NULL);
  ck_assert(sc->energy_up != NULL);

  for(i = 1; i <= fc->length; i++){
    int counter = 0;
    for(j = 1; i + j - 1 <= fc->length; j++){
      counter += (j+i-1) * -100;
      ck_assert_int_eq(sc->energy_up[i][j], counter);
    }
  }

  vrna_sc_remove(fc);
  ck_assert(fc->sc == NULL);

  /* clean up */
  vrna_fold_compound_free(fc);
  free(seq);


#test test_vrna_sc_add_up_addition

  int i, e;

  char *seq = vrna_random_string(10, "ACGU");

  vrna_fold_compound_t *fc = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);

  for(e = 0, i = 1; i <= fc->length; i++){
    vrna_sc_add_up(fc, 1, -1. * i, VRNA_OPTION_DEFAULT);
    e += -100 * i;
  }
  vrna_sc_t *sc = fc->sc;

  ck_assert(sc != NULL);
  ck_assert(sc->energy_up != NULL);

  for(i = 1; i <= fc->length; i++){
    ck_assert_int_eq(sc->energy_up[1][i], e);
  }

  vrna_sc_remove(fc);
  ck_assert(fc->sc == NULL);

  /* clean up */
  vrna_fold_compound_free(fc);
  free(seq);


#test test_vrna_sc_add_up_addition_extended

  int i, j, e;

  char *seq = vrna_random_string(10, "ACGU");

  vrna_fold_compound_t *fc = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);

  for(i = 1, e = 0; i <= fc->length; i++){
    e += -100 * i;
    for(j = 1; j <= fc->length; j++){
      vrna_sc_add_up(fc, j, -1. * i, VRNA_OPTION_DEFAULT);
    }
  }
  vrna_sc_t *sc = fc->sc;

  ck_assert(sc != NULL);
  ck_assert(sc->energy_up != NULL);

  for(i = 1; i <= fc->length; i++){
    for(j = 1; i + j - 1 <= fc->length; j++){
      ck_assert_int_eq(sc->energy_up[i][j], j * e);
    }
  }

  vrna_sc_remove(fc);
  ck_assert(fc->sc == NULL);

  /* clean up */
  vrna_fold_compound_free(fc);
  free(seq);


#test test_vrna_sc_add_bp

  int i, j;

  char *seq = vrna_random_string(10, "ACGU");

  vrna_fold_compound_t *fc = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);

  for(i = 1; i < fc->length; i++)
    for(j = i + 1; j <= fc->length; j++)
      vrna_sc_add_bp(fc, i, j, -1. * (i + j), VRNA_OPTION_DEFAULT);

  vrna_sc_t *sc = fc->sc;

  ck_assert(sc != NULL);
  ck_assert(sc->energy_bp != NULL);

  for(i = 1; i < fc->length; i++){
    for(j = i + 1; j <= fc->length; j++){
      ck_assert_int_eq(sc->energy_bp[fc->jindx[j] + i], -100 * (i + j));
    }
  }

  vrna_sc_remove(fc);
  ck_assert(fc->sc == NULL);

  /* clean up */
  vrna_fold_compound_free(fc);
  free(seq);


#test test_vrna_sc_add_bp_addition

  int i, j, c, num;

  num = 10;

  char *seq = vrna_random_string(10, "ACGU");

  vrna_fold_compound_t *fc = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);

  for(c = 0; c < num; c++)
    for(i = 1; i < fc->length; i++)
      for(j = i + 1; j <= fc->length; j++)
        vrna_sc_add_bp(fc, i, j, -1. * (i + j), VRNA_OPTION_DEFAULT);

  vrna_sc_t *sc = fc->sc;

  ck_assert(sc != NULL);
  ck_assert(sc->energy_bp != NULL);

  for(i = 1; i < fc->length; i++){
    for(j = i + 1; j <= fc->length; j++){
      ck_assert_int_eq(sc->energy_bp[fc->jindx[j] + i], -100 * num * (i + j));
    }
  }

  vrna_sc_remove(fc);
  ck_assert(fc->sc == NULL);

  /* clean up */
  vrna_fold_compound_free(fc);
  free(seq);


#test test_vrna_sc_add_bp_removal

  int i, j, c, num, num_r;

  num   = 10;
  num_r = 3;

  char *seq = vrna_random_string(10, "ACGU");

  vrna_fold_compound_t *fc = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);

  for(c = 0; c < num; c++)
    for(i = 1; i < fc->length; i++)
      for(j = i + 1; j <= fc->length; j++)
        vrna_sc_add_bp(fc, i, j, -1. * (i + j), VRNA_OPTION_DEFAULT);

  for(c = 0; c < num_r; c++)
    for(i = 1; i < fc->length; i++)
      for(j = i + 1; j <= fc->length; j++)
        vrna_sc_add_bp(fc, i, j, 1. * (i + j), VRNA_OPTION_DEFAULT);

  vrna_sc_t *sc = fc->sc;

  ck_assert(sc != NULL);
  ck_assert(sc->energy_bp != NULL);

  for(i = 1; i < fc->length; i++){
    for(j = i + 1; j <= fc->length; j++){
      ck_assert_int_eq(sc->energy_bp[fc->jindx[j] + i], -100 * (num - num_r) * (i + j));
    }
  }

  vrna_sc_remove(fc);
  ck_assert(fc->sc == NULL);

  /* clean up */
  vrna_fold_compound_free(fc);
  free(seq);
