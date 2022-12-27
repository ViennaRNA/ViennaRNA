#ifndef   _RNAXPLORER_REPELLANT_SAMPLING_H_
#define   _RNAXPLORER_REPELLANT_SAMPLING_H_

void
repellant_sampling(vrna_fold_compound_t *fc);

int
rnax_add_repulsion(vrna_fold_compound_t *fc,
                   const char *structure,
                   double     strength);

int
rnax_change_repulsion(vrna_fold_compound_t *fc,
                      int                  id,
                      double     strength);

#endif
