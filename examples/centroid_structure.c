#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include  <ViennaRNA/data_structures.h>
#include  <ViennaRNA/params/basic.h>
#include  <ViennaRNA/utils/basic.h>
#include  <ViennaRNA/eval/structures.h>
#include  <ViennaRNA/mfe/global.h>
#include  <ViennaRNA/partfunc/basic.h>
#include  <ViennaRNA/structures/centroid.h>


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
  double  dist;
  char    *cent = vrna_centroid(vc, &dist);

  /* print centroid structure, its free energy and mean distance to the ensemble */
  printf("%s\n%s (%6.2f d=%6.2f)\n",
         seq,
         cent,
         vrna_eval_structure(vc, cent),
         dist);

  /* free centroid structure */
  free(cent);

  /* free memory occupied by vrna_fold_compound */
  vrna_fold_compound_free(vc);

  return EXIT_SUCCESS;
}
