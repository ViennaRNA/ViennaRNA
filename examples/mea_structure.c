#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include  <ViennaRNA/data_structures.h>
#include  <ViennaRNA/params/basic.h>
#include  <ViennaRNA/utils/basic.h>
#include  <ViennaRNA/eval/structures.h>
#include  <ViennaRNA/mfe/global.h>
#include  <ViennaRNA/partfunc/basic.h>
#include  <ViennaRNA/structures/mea.h>


int
main(int  argc,
     char *argv[])
{
  char *seq = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCAGGUACGA";

  /* get a vrna_fold_compound with default settings */
  vrna_fold_compound_t *vc = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);

  /* call MFE function */
  double mfe = (double)vrna_mfe(vc, NULL);

  /* rescale parameters for Boltzmann factors */
  vrna_exp_params_rescale(vc, &mfe);

  /* call PF function */
  FLT_OR_DBL en = vrna_pf(vc, NULL);

  /* compute centroid structure */
  float mea;
  char  *mea_structure = vrna_MEA(fc, 1.0, &mea);

  /* print MEA structure, its free energy and maximum expected accuracy */
  printf("%s\n%s (%6.2f MEA=%6.2f)\n",
         seq,
         mea_structure,
         vrna_eval_structure(vc, mea_structure),
         mea);

  /* free MEA structure */
  free(mea_structure);

  /* free memory occupied by vrna_fold_compound */
  vrna_fold_compound_free(vc);

  return EXIT_SUCCESS;
}
