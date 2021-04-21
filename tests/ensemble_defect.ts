#include <stdio.h>      /* printf, scanf, NULL */
#include <stdlib.h>     /* malloc, free, rand */

#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/utils/structures.h>
#include <ViennaRNA/equilibrium_probs.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/part_func.h>

#suite Ensemble_Defect

#tcase Ensemble_Defect

#test test_ensemble_defect
{
  vrna_md_t             md;
  vrna_fold_compound_t  *vc;
  const char            *sequence = "AGGAAACCUUAAUUGGUUA";
  const char  		*str1     = ".((...))(([[..))]].";
  double ed;
  short int *ptpk;

  vrna_md_set_default(&md);

  vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_PF);

  vrna_pf(vc, NULL);

  ptpk = vrna_ptable_from_string(str1, VRNA_BRACKETS_ANY);

  ck_assert(vrna_ensemble_defect(vc, str1) == 0.6140797258673892);
  ck_assert(vrna_ensemble_defect_pt(vc, ptpk) == 0.7279171755397522);

  vrna_fold_compound_free(vc);
  free(ptpk);
}
