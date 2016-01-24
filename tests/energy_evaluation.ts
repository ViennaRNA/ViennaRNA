#include <ViennaRNA/params.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/loop_energies.h>

#suite Energy_Evaluating_Functions

/*
 * check for properly working E_Hairpin() function
 */

#test eval_E_hairpin
  vrna_param_t param = {0};
  int i, j, k, l;

  for (i = 0; i <= 30; ++i)
    param.hairpin[i] = i + 100;

  for (i = 0, l = 0; i <= NBPAIRS; ++i)
    for (j = 0; j <= 4; ++j)
      for (k = 0; k <=4; ++k)
        param.mismatchH[i][j][k] = ++l;

  strcpy(param.Triloops, "CGGGG GCCCC");
  param.Triloop_E[0] = 1000;
  param.Triloop_E[1] = 2000;

  strcpy(param.Tetraloops, "CGGGGG GCCCCC");
  param.Tetraloop_E[0] = 3000;
  param.Tetraloop_E[1] = 4000;

  strcpy(param.Hexaloops, "CGGGGGGG GCCCCCCC");
  param.Hexaloop_E[0] = 5000;
  param.Hexaloop_E[1] = 6000;

  param.TerminalAU = 10000;

  ck_assert_int_eq(E_Hairpin(0, 0, 0, 0, NULL, &param), 100);
  ck_assert_int_eq(E_Hairpin(1, 0, 0, 0, NULL, &param), 101);
  ck_assert_int_eq(E_Hairpin(2, 0, 0, 0, NULL, &param), 102);
  ck_assert_int_eq(E_Hairpin(3, 0, 0, 0, NULL, &param), 104);
  ck_assert_int_eq(E_Hairpin(3, 1, 2, 3, NULL, &param), 142);
  ck_assert_int_eq(E_Hairpin(3, 0, 0, 0, "CGGGG", &param), 104);

  param.model_details.special_hp = 1;

  ck_assert_int_eq(E_Hairpin(3, 0, 0, 0, "CGGGG", &param), 1000);
  ck_assert_int_eq(E_Hairpin(3, 0, 0, 0, "GCCCC", &param), 2000);
  ck_assert_int_eq(E_Hairpin(3, 0, 0, 0, "AUUUU", &param), 103); //mismatchH is ignored
  ck_assert_int_eq(E_Hairpin(3, 1, 0, 0, "AUUUU", &param), 103);
  ck_assert_int_eq(E_Hairpin(3, 2, 0, 0, "AUUUU", &param), 103);
  ck_assert_int_eq(E_Hairpin(3, 3, 0, 0, "AUUUU", &param), 10103);
  ck_assert_int_eq(E_Hairpin(3, 4, 0, 0, "AUUUU", &param), 10103);

  ck_assert_int_eq(E_Hairpin(4, 0, 0, 0, "CGGGGG", &param), 3000);
  ck_assert_int_eq(E_Hairpin(4, 0, 0, 0, "GCCCCC", &param), 4000);
  ck_assert_int_eq(E_Hairpin(4, 0, 0, 0, "AUUUUU", &param), 105);

  ck_assert_int_eq(E_Hairpin(6, 0, 0, 0, "CGGGGGGG", &param), 5000);
  ck_assert_int_eq(E_Hairpin(6, 0, 0, 0, "GCCCCCCC", &param), 6000);
  ck_assert_int_eq(E_Hairpin(6, 0, 0, 0, "AUUUUUUU", &param), 107);

/*
 * check for properly working E_Stem() function
 */

#test eval_E_stem
  vrna_param_t param = {0};

  param.dangle3[1][1] = 3;
  param.dangle3[3][1] = 2;
  param.dangle5[1][1] = 5;
  param.dangle5[1][2] = 12;
  param.mismatchExt[1][1][1] = 100;
  param.mismatchM[1][1][1] = 200;
  param.TerminalAU = 1000;
  param.MLintern[1] = 10000;
  param.MLintern[3] = 100000;

  ck_assert_int_eq(E_Stem(1, -1, -1, 0, &param), 10000);
  ck_assert_int_eq(E_Stem(1, -1, -1, 1, &param), 0);
  ck_assert_int_eq(E_Stem(3, -1, -1, 0, &param), 101000);
  ck_assert_int_eq(E_Stem(3, -1, -1, 1, &param), 1000);
  ck_assert_int_eq(E_Stem(1, 1, -1, 0, &param), 10005);
  ck_assert_int_eq(E_Stem(1, 1, -1, 1, &param), 5);
  ck_assert_int_eq(E_Stem(1, -1, 1, 0, &param), 10003);
  ck_assert_int_eq(E_Stem(1, -1, 1, 1, &param), 3);
  ck_assert_int_eq(E_Stem(1, 1, 1, 0, &param), 10200);
  ck_assert_int_eq(E_Stem(1, 1, 1, 1, &param), 100);

/*
 * check for properly working E_MLstem() function
 */

#test eval_E_MLstem
  vrna_param_t param = {0};
  param.MLintern[1] = 1;
  param.dangle3[1][0] = 3;
  param.dangle5[1][0] = 5;
  param.mismatchM[1][0][0] = 100;

  ck_assert_int_eq(E_MLstem(1, -1, -1, &param), 1);
  ck_assert_int_eq(E_MLstem(1, 0, -1, &param), 6);
  ck_assert_int_eq(E_MLstem(1, -1, 0, &param), 4);
  ck_assert_int_eq(E_MLstem(1, 0, 0, &param), 101);

  param.TerminalAU = 10000;
  param.MLintern[3] = 10;
  param.dangle3[3][0] = 33;
  param.dangle5[3][0] = 55;
  param.mismatchM[3][0][0] = 1000;

  ck_assert_int_eq(E_MLstem(3, -1, -1, &param), 10010);
  ck_assert_int_eq(E_MLstem(3, 0, -1, &param), 10065);
  ck_assert_int_eq(E_MLstem(3, -1, 0, &param), 10043);
  ck_assert_int_eq(E_MLstem(3, 0, 0, &param), 11010);

/*
 * check for properly working E_ExtLoop() function
 */

#test eval_E_ExtLoop
  vrna_param_t param = {0};
  param.dangle3[1][0] = 3;
  param.dangle5[1][0] = 5;
  param.mismatchExt[1][0][0] = 100;

  ck_assert_int_eq(E_ExtLoop(1, -1, -1, &param), 0);
  ck_assert_int_eq(E_ExtLoop(1, 0, -1, &param), 5);
  ck_assert_int_eq(E_ExtLoop(1, -1, 0, &param), 3);
  ck_assert_int_eq(E_ExtLoop(1, 0, 0, &param), 100);

  param.TerminalAU = 10000;
  param.dangle3[3][0] = 33;
  param.dangle5[3][0] = 55;
  param.mismatchExt[3][0][0] = 1000;

  ck_assert_int_eq(E_ExtLoop(3, -1, -1, &param), 10000);
  ck_assert_int_eq(E_ExtLoop(3, 0, -1, &param), 10055);
  ck_assert_int_eq(E_ExtLoop(3, -1, 0, &param), 10033);
  ck_assert_int_eq(E_ExtLoop(3, 0, 0, &param), 11000);

/*
 * check for properly working E_IntLoop() function
 */

#test eval_E_IntLoop
  vrna_param_t param = {0};

  param.stack[1][2] = 1;
  ck_assert_int_eq(E_IntLoop(0, 0, 1, 2, -1, -1, -1, -1, &param), 1);

  param.TerminalAU = 100;
  param.bulge[1] = 10;
  param.bulge[2] = 20;
  ck_assert_int_eq(E_IntLoop(0, 1, 1, 2, -1, -1, -1, -1, &param), 11);
  ck_assert_int_eq(E_IntLoop(1, 0, 1, 2, -1, -1, -1, -1, &param), 11);
  ck_assert_int_eq(E_IntLoop(0, 2, 1, 2, -1, -1, -1, -1, &param), 20);
  ck_assert_int_eq(E_IntLoop(2, 0, 1, 2, -1, -1, -1, -1, &param), 20);
  ck_assert_int_eq(E_IntLoop(0, 2, 3, 2, -1, -1, -1, -1, &param), 120);
  ck_assert_int_eq(E_IntLoop(0, 2, 2, 3, -1, -1, -1, -1, &param), 120);
  ck_assert_int_eq(E_IntLoop(0, 2, 3, 3, -1, -1, -1, -1, &param), 220);

  param.int11[1][2][3][4] = 11;
  ck_assert_int_eq(E_IntLoop(1, 1, 1, 2, 3, 4, -1, -1, &param), 11);

  param.int21[1][2][1][2][3] = 21;
  ck_assert_int_eq(E_IntLoop(1, 2, 1, 2, 1, 3, -1, 2, &param), 21);
  ck_assert_int_eq(E_IntLoop(2, 1, 2, 1, 2, -1, 3, 1, &param), 21);

  param.internal_loop[4] = 10;
  param.ninio[2] = 100;
  param.mismatch1nI[1][1][2] = 1;
  param.mismatch1nI[2][4][3] = 2;
  ck_assert_int_eq(E_IntLoop(1, 3, 1, 2, 1, 2, 3, 4, &param), 213);
  ck_assert_int_eq(E_IntLoop(3, 1, 1, 2, 1, 2, 3, 4, &param), 213);

  param.int22[1][2][1][2][3][4] = 22;
  ck_assert_int_eq(E_IntLoop(2, 2, 1, 2, 1, 4, 2, 3, &param), 22);

  param.internal_loop[5] = 20;
  param.mismatch23I[1][1][2] = 3;
  param.mismatch23I[2][4][3] = 4;
  ck_assert_int_eq(E_IntLoop(2, 3, 1, 2, 1, 2, 3, 4, &param), 127);
  ck_assert_int_eq(E_IntLoop(3, 2, 1, 2, 1, 2, 3, 4, &param), 127);

  param.internal_loop[8] = 30;
  param.mismatchI[1][1][2] = 2;
  param.mismatchI[2][4][3] = 3;
  ck_assert_int_eq(E_IntLoop(3, 5, 1, 2, 1, 2, 3, 4, &param), 235);
  ck_assert_int_eq(E_IntLoop(5, 3, 1, 2, 1, 2, 3, 4, &param), 235);
