#include <stdlib.h>
#include <stdio.h>

#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/utils/strings.h>
#include <ViennaRNA/mfe.h>


int
main()
{
  /* initialize random number generator */
  vrna_init_rand();

  /* Generate a random sequence of 50 nucleotides */
  char                  *seq = vrna_random_string(50, "ACGU");

  /* Create a fold compound for the sequence */
  vrna_fold_compound_t  *fc = vrna_fold_compound(seq, NULL, VRNA_OPTION_DEFAULT);

  /* allocate memory for MFE structure (length + 1) */
  char                  *structure = (char *)vrna_alloc(sizeof(char) * (strlen(seq) + 1));

  /* predict Minmum Free Energy and corresponding secondary structure */
  float                 mfe = vrna_mfe(fc, structure);

  /* print sequence, structure and MFE */
  printf("%s\n%s [ %6.2f ]\n", seq, structure, mfe);

  /* cleanup memory */
  free(seq);
  free(structure);
  vrna_fold_compound_free(fc);

  return 0;
}
