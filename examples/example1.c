#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include  <ViennaRNA/data_structures.h>
#include  <ViennaRNA/params.h>
#include  <ViennaRNA/utils.h>
#include  <ViennaRNA/eval.h>
#include  <ViennaRNA/fold.h>
#include  <ViennaRNA/part_func.h>


int main(int argc, char *argv[]){

  char  *seq = "AGACGACAAGGUUGAAUCGCACCCACAGUCUAUGAGUCGGUGACAACAUUACGAAAGGCUGUAAAAUCAAUUAUUCACCACAGGGGGCCCCCGUGUCUAG";
  char  *mfe_structure = vrna_alloc(sizeof(char) * (strlen(seq) + 1));
  char  *prob_string   = vrna_alloc(sizeof(char) * (strlen(seq) + 1));

  /* get a vrna_fold_compound with default settings */
  vrna_fold_compound_t *vc = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);

  /* call MFE function */
  double mfe = (double)vrna_mfe(vc, mfe_structure);

  printf("%s\n%s (%6.2f)\n", seq, mfe_structure, mfe);

  /* rescale parameters for Boltzmann factors */
  vrna_exp_params_rescale(vc, &mfe);

  /* call PF function */
  FLT_OR_DBL en = vrna_pf(vc, prob_string);

  /* print probability string and free energy of ensemble */
  printf("%s (%6.2f)\n", prob_string, en);

  /* compute centroid structure */
  double dist;
  char *cent = vrna_centroid(vc, &dist);

  /* print centroid structure, its free energy and mean distance to the ensemble */
  printf("%s (%6.2f d=%6.2f)\n", cent, vrna_eval_structure(vc, cent), dist);

  /* free centroid structure */
  free(cent);

  /* free pseudo dot-bracket probability string */
  free(prob_string);

  /* free mfe structure */
  free(mfe_structure);

  /* free memory occupied by vrna_fold_compound */
  vrna_fold_compound_free(vc);

  return EXIT_SUCCESS;
}
