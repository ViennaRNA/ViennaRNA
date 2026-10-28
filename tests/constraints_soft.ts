#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/utils/strings.h>
#include <ViennaRNA/constraints/soft.h>
#include <ViennaRNA/mfe/global.h>

#suite Constraints

#tcase  SoftConstraints

#test test_vrna_sc_fractional_up
{
  int sign;

  for (sign = -1; sign <= 1; sign += 2) {
    vrna_fold_compound_t *fc = vrna_fold_compound("AAAAAA", NULL, VRNA_OPTION_MFE | VRNA_OPTION_PF);
    double              expected;

    /* Preparing between additions must not discard fractional contributions. */
    vrna_sc_add_up(fc, 3, sign * 0.003, VRNA_OPTION_DEFAULT);
    vrna_sc_prepare(fc, VRNA_OPTION_MFE | VRNA_OPTION_PF);
    ck_assert_int_eq(fc->sc->energy_up[3][1], 0);

    vrna_sc_add_up(fc, 3, sign * 0.003, VRNA_OPTION_DEFAULT);
    vrna_sc_prepare(fc, VRNA_OPTION_MFE | VRNA_OPTION_PF);
    expected = exp(-sign * 6. / fc->exp_params->kT);
    ck_assert_int_eq(sizeof(fc->sc->up_storage[0]), sizeof(double));
    ck_assert_int_eq(fc->sc->energy_up[3][1], sign);
    ck_assert(fabs(fc->sc->exp_energy_up[3][1] - expected) < 1e-6);

    /* MFE rounds per nucleotide, rather than after summing a whole stretch. */
    vrna_sc_add_up(fc, 4, sign * 0.004, VRNA_OPTION_DEFAULT);
    vrna_sc_prepare(fc, VRNA_OPTION_MFE | VRNA_OPTION_PF);
    expected = exp(-sign * 10. / fc->exp_params->kT);
    ck_assert_int_eq(fc->sc->energy_up[3][2], sign);
    ck_assert(fabs(fc->sc->exp_energy_up[3][2] - expected) < 1e-6);

    vrna_fold_compound_free(fc);
  }
}

#test test_vrna_sc_fractional_bp
{
  int sign;

  for (sign = -1; sign <= 1; sign += 2) {
    vrna_fold_compound_t *fc = vrna_fold_compound("GCGCUUCGCCGC", NULL, VRNA_OPTION_MFE | VRNA_OPTION_PF);
    int                 idx = fc->jindx[11] + 2;
    double              expected;

    vrna_sc_add_bp(fc, 2, 11, sign * 0.003, VRNA_OPTION_DEFAULT);
    vrna_sc_prepare(fc, VRNA_OPTION_MFE | VRNA_OPTION_PF);
    ck_assert_int_eq(fc->sc->energy_bp[idx], 0);

    vrna_sc_add_bp(fc, 2, 11, sign * 0.003, VRNA_OPTION_DEFAULT);
    vrna_sc_prepare(fc, VRNA_OPTION_MFE | VRNA_OPTION_PF);
    expected = exp(-sign * 6. / fc->exp_params->kT);
    ck_assert_int_eq(sizeof(fc->sc->bp_storage[2][0].e), sizeof(double));
    ck_assert_int_eq(fc->sc->energy_bp[idx], sign);
    ck_assert(fabs(fc->sc->exp_energy_bp[idx] - expected) < 1e-6);

    vrna_fold_compound_free(fc);
  }
}

#test test_vrna_sc_fractional_set
{
  vrna_fold_compound_t *fc = vrna_fold_compound("GCGCUUCGCCGC", NULL, VRNA_OPTION_MFE | VRNA_OPTION_PF);
  FLT_OR_DBL          up[13] = {0};
  FLT_OR_DBL          bp[13][13] = {{0}};
  const FLT_OR_DBL    *rows[13];
  unsigned int        i;
  int                 idx = fc->jindx[11] + 2;

  for (i = 0; i < 13; i++)
    rows[i] = bp[i];

  vrna_sc_add_up(fc, 3, 1., VRNA_OPTION_DEFAULT);
  vrna_sc_add_bp(fc, 2, 11, 1., VRNA_OPTION_DEFAULT);
  up[3] = 0.004;
  bp[2][11] = 0.006;
  vrna_sc_set_up(fc, up, VRNA_OPTION_DEFAULT);
  vrna_sc_set_bp(fc, rows, VRNA_OPTION_DEFAULT);
  vrna_sc_prepare(fc, VRNA_OPTION_MFE | VRNA_OPTION_PF);
  ck_assert_int_eq(fc->sc->energy_up[3][1], 0);
  ck_assert_int_eq(fc->sc->energy_bp[idx], 1);
  ck_assert(fabs(fc->sc->exp_energy_up[3][1] - exp(-4. / fc->exp_params->kT)) < 1e-6);
  ck_assert(fabs(fc->sc->exp_energy_bp[idx] - exp(-6. / fc->exp_params->kT)) < 1e-6);

  /* The add and set APIs must use the same unrounded storage. */
  vrna_sc_add_up(fc, 3, 0.003, VRNA_OPTION_DEFAULT);
  vrna_sc_add_bp(fc, 2, 11, -0.003, VRNA_OPTION_DEFAULT);
  vrna_sc_prepare(fc, VRNA_OPTION_MFE | VRNA_OPTION_PF);
  ck_assert_int_eq(fc->sc->energy_up[3][1], 1);
  ck_assert_int_eq(fc->sc->energy_bp[idx], 0);
  ck_assert(fabs(fc->sc->exp_energy_up[3][1] - exp(-7. / fc->exp_params->kT)) < 1e-6);
  ck_assert(fabs(fc->sc->exp_energy_bp[idx] - exp(-3. / fc->exp_params->kT)) < 1e-6);

  vrna_fold_compound_free(fc);
}

#test test_vrna_sc_fractional_rounding_boundary
{
  vrna_fold_compound_t *fc = vrna_fold_compound("GCGCUUCGCCGC", NULL, VRNA_OPTION_MFE);
  FLT_OR_DBL          energy = 0.0049999999;

  /* Converting the stored double to float before rounding can turn an energy
   * below 0.005 kcal/mol into an apparent half-integer in internal units. */
  ck_assert(energy * 100. < 0.5);
  vrna_sc_add_up(fc, 3, energy, VRNA_OPTION_DEFAULT);
  vrna_sc_add_bp(fc, 2, 11, -energy, VRNA_OPTION_DEFAULT);
  vrna_sc_prepare(fc, VRNA_OPTION_MFE);
  ck_assert_int_eq(fc->sc->energy_up[3][1], 0);
  ck_assert_int_eq(fc->sc->energy_bp[fc->jindx[11] + 2], 0);

  vrna_fold_compound_free(fc);
}

#test test_vrna_sc_fractional_dimer
{
  vrna_fold_compound_t *fc = vrna_fold_compound("AAAA&GCGCAAAAGCGC", NULL, VRNA_OPTION_DEFAULT);
  char                structure[17];
  double              baseline = vrna_mfe_dimer(fc, structure);

  /* The strands cannot pair with each other. The second monomer must retain
   * its fractional constraints when the dimer wrapper copies them. */
  vrna_sc_add_bp(fc, 5, 16, 0.003, VRNA_OPTION_DEFAULT);
  vrna_sc_add_bp(fc, 5, 16, 0.003, VRNA_OPTION_DEFAULT);
  ck_assert(fabs(vrna_mfe_dimer(fc, structure) - baseline - 0.01) < 1e-6);
  ck_assert_str_eq(structure, "....((((....))))");

  vrna_fold_compound_free(fc);
}

#test test_vrna_sc_add_up_simple
{
  int                   i, j;

  char                  *seq = vrna_random_string(10, "ACGU");

  vrna_fold_compound_t  *fc = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);

  for (i = 1; i <= fc->length; i++)
    vrna_sc_add_up(fc, i, -1. * i, VRNA_OPTION_DEFAULT);

  vrna_sc_prepare(fc, VRNA_OPTION_MFE);

  vrna_sc_t *sc = fc->sc;

  ck_assert(sc != NULL);
  ck_assert(sc->energy_up != NULL);

  for (i = 1; i <= fc->length; i++) {
    int counter = 0;
    for (j = 1; i + j - 1 <= fc->length; j++) {
      counter += (j + i - 1) * -100;
      ck_assert_int_eq(sc->energy_up[i][j], counter);
    }
  }

  vrna_sc_remove(fc);
  ck_assert(fc->sc == NULL);

  /* clean up */
  vrna_fold_compound_free(fc);
  free(seq);
}

#test test_vrna_sc_add_up_addition
{
  int                   i, e;

  char                  *seq = vrna_random_string(10, "ACGU");

  vrna_fold_compound_t  *fc = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);

  for (e = 0, i = 1; i <= fc->length; i++) {
    vrna_sc_add_up(fc, 1, -1. * i, VRNA_OPTION_DEFAULT);
    e += -100 * i;
  }

  vrna_sc_prepare(fc, VRNA_OPTION_MFE);

  vrna_sc_t *sc = fc->sc;

  ck_assert(sc != NULL);
  ck_assert(sc->energy_up != NULL);

  for (i = 1; i <= fc->length; i++)
    ck_assert_int_eq(sc->energy_up[1][i], e);

  vrna_sc_remove(fc);
  ck_assert(fc->sc == NULL);

  /* clean up */
  vrna_fold_compound_free(fc);
  free(seq);
}

#test test_vrna_sc_add_up_addition_extended
{
  int                   i, j, e;

  char                  *seq = vrna_random_string(10, "ACGU");

  vrna_fold_compound_t  *fc = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);

  for (i = 1, e = 0; i <= fc->length; i++) {
    e += -100 * i;
    for (j = 1; j <= fc->length; j++)
      vrna_sc_add_up(fc, j, -1. * i, VRNA_OPTION_DEFAULT);
  }

  vrna_sc_prepare(fc, VRNA_OPTION_MFE);

  vrna_sc_t *sc = fc->sc;

  ck_assert(sc != NULL);
  ck_assert(sc->energy_up != NULL);

  for (i = 1; i <= fc->length; i++)
    for (j = 1; i + j - 1 <= fc->length; j++)
      ck_assert_int_eq(sc->energy_up[i][j], j * e);

  vrna_sc_remove(fc);
  ck_assert(fc->sc == NULL);

  /* clean up */
  vrna_fold_compound_free(fc);
  free(seq);
}

#test test_vrna_sc_add_bp
{
  int                   i, j, turn;

  char                  *seq = vrna_random_string(10, "ACGU");

  vrna_fold_compound_t  *fc = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);
  turn = fc->params->model_details.min_loop_size;

  for (i = 1; i < fc->length; i++)
    for (j = i + turn + 1; j <= fc->length; j++)
      vrna_sc_add_bp(fc, i, j, -1. * (i + j), VRNA_OPTION_DEFAULT);

  vrna_sc_prepare(fc, VRNA_OPTION_MFE);

  vrna_sc_t *sc = fc->sc;

  ck_assert(sc != NULL);
  ck_assert(sc->energy_bp != NULL);

  for (i = 1; i < fc->length; i++)
    for (j = i + turn + 1; j <= fc->length; j++)
      ck_assert_int_eq(sc->energy_bp[fc->jindx[j] + i], -100 * (i + j));

  vrna_sc_remove(fc);
  ck_assert(fc->sc == NULL);

  /* clean up */
  vrna_fold_compound_free(fc);
  free(seq);
}

#test test_vrna_sc_add_bp_addition
{
  int                   i, j, c, num, turn;

  num = 10;

  char                  *seq = vrna_random_string(10, "ACGU");

  vrna_fold_compound_t  *fc = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);
  turn = fc->params->model_details.min_loop_size;

  for (c = 0; c < num; c++)
    for (i = 1; i < fc->length; i++)
      for (j = i + turn + 1; j <= fc->length; j++)
        vrna_sc_add_bp(fc, i, j, -1. * (i + j), VRNA_OPTION_DEFAULT);

  vrna_sc_prepare(fc, VRNA_OPTION_MFE);

  vrna_sc_t *sc = fc->sc;

  ck_assert(sc != NULL);
  ck_assert(sc->energy_bp != NULL);

  for (i = 1; i < fc->length; i++)
    for (j = i + turn + 1; j <= fc->length; j++)
      ck_assert_int_eq(sc->energy_bp[fc->jindx[j] + i], -100 * num * (i + j));

  vrna_sc_remove(fc);
  ck_assert(fc->sc == NULL);

  /* clean up */
  vrna_fold_compound_free(fc);
  free(seq);
}

#test test_vrna_sc_add_bp_removal
{
  int                   i, j, c, num, num_r, turn;

  num   = 10;
  num_r = 3;

  char                  *seq = vrna_random_string(10, "ACGU");

  vrna_fold_compound_t  *fc = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);
  turn = fc->params->model_details.min_loop_size;

  for (c = 0; c < num; c++)
    for (i = 1; i < fc->length; i++)
      for (j = i + turn + 1; j <= fc->length; j++)
        vrna_sc_add_bp(fc, i, j, -1. * (i + j), VRNA_OPTION_DEFAULT);

  for (c = 0; c < num_r; c++)
    for (i = 1; i < fc->length; i++)
      for (j = i + turn + 1; j <= fc->length; j++)
        vrna_sc_add_bp(fc, i, j, 1. * (i + j), VRNA_OPTION_DEFAULT);

  vrna_sc_prepare(fc, VRNA_OPTION_MFE);

  vrna_sc_t *sc = fc->sc;

  ck_assert(sc != NULL);
  ck_assert(sc->energy_bp != NULL);

  for (i = 1; i < fc->length; i++)
    for (j = i + turn + 1; j <= fc->length; j++)
      ck_assert_int_eq(sc->energy_bp[fc->jindx[j] + i], -100 * (num - num_r) * (i + j));

  vrna_sc_remove(fc);
  ck_assert(fc->sc == NULL);

  /* clean up */
  vrna_fold_compound_free(fc);
  free(seq);
}


#main-pre
    srunner_set_tap(sr, "-");
