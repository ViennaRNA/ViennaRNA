#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <ViennaRNA/model.h>
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
  char      *seq = vrna_random_string(50, "ACGU");

  /* allocate memory for MFE structure (length + 1) */
  char      *structure = (char *)vrna_alloc(sizeof(char) * (strlen(seq) + 1));

  /* create a new model details structure to store the Model Settings */
  vrna_md_t md;

  /* ALWAYS set default model settings first! */
  vrna_md_set_default(&md);

  /* change temperature and activate G-Quadruplex prediction */
  md.temperature  = 25.0; /* 25 Deg Celcius */
  md.gquad        = 1;    /* Turn-on G-Quadruples support */

  /* create a fold compound */
  vrna_fold_compound_t  *fc = vrna_fold_compound(seq, &md, VRNA_OPTION_DEFAULT);

  /* predict Minmum Free Energy and corresponding secondary structure */
  float                 mfe = vrna_mfe(fc, structure);

  /* print sequence, structure and MFE */
  printf("%s\n%s [ %6.2f ]\n", seq, structure, mfe);

  /* cleanup memory */
  free(structure);
  vrna_fold_compound_free(fc);

  return 0;
}
