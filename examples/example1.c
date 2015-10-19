#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include  <ViennaRNA/data_structures.h>
#include  <ViennaRNA/utils.h>
#include  <ViennaRNA/eval.h>
#include  <ViennaRNA/fold.h>
#include  <ViennaRNA/part_func.h>


int main(int argc, char *argv[]){

  char  *seq = "AGACGACAAGGUUGAAUCGCACCCACAGUCUAUGAGUCGGUGACAACAUUACGAAAGGCUGUAAAAUCAAUUAUUCACCACAGGGGGCCCCCGUGUCUAG";
  char  *mfe_structure = vrna_alloc(sizeof(char) * (strlen(seq) + 1));
  char  *prob_string   = vrna_alloc(sizeof(char) * (strlen(seq) + 1));

  /* get a vrna_fold_compound with MFE and PF DP matrices and default model details */
  vrna_fold_compound *vc = vrna_get_fold_compound(seq, NULL, VRNA_OPTION_MFE | VRNA_OPTION_PF);

  /* call MFE function */
  double mfe = (double)vrna_fold(vc, mfe_structure);

  printf("%s\n%s (%6.2f)\n", seq, mfe_structure, mfe);

  /* rescale parameters for Boltzmann factors */
  vrna_exp_params_update(vc, &mfe);

  /* call PF function */
  FLT_OR_DBL en = vrna_pf_fold(vc, prob_string);

  /* print probability string and free energy of ensemble */
  printf("%s (%6.2f)\n", prob_string, en);

  /* compute centroid structure */
  double dist;
  char *cent = vrna_get_centroid_struct(vc, &dist);

  /* print centroid structure, its free energy and mean distance to the ensemble */
  printf("%s (%6.2f d=%6.2f)\n", cent, vrna_eval_structure(vc, cent), dist);

  /* free centroid structure */
  free(cent);

  /* free pseudo dot-bracket probability string */
  free(prob_string);

  /* free mfe structure */
  free(mfe_structure);

  /* free memory occupied by vrna_fold_compound */
  vrna_free_fold_compound(vc);

  return EXIT_SUCCESS;
}
