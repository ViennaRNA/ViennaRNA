#ifndef VIENNA_RNA_PACKAGE_ZSCORE_H
#define VIENNA_RNA_PACKAGE_ZSCORE_H

typedef struct vrna_zsc_dat_s *vrna_zsc_dat_t;

#define VRNA_ZSCORE_OPTIONS_NONE      1U
#define VRNA_ZSCORE_FILTER_ON         2U
#define VRNA_ZSCORE_PRE_FILTER        4U
#define VRNA_ZSCORE_REPORT_SUBSUMED   8U
#define VRNA_ZSCORE_MODEL_DEFAULT     16U
#define VRNA_ZSCORE_SETTINGS_DEFAULT  (VRNA_ZSCORE_FILTER_ON | VRNA_ZSCORE_MODEL_DEFAULT)

int
vrna_zsc_filter_init(vrna_fold_compound_t *fc,
                     double               min_z,
                     unsigned int         options);


int
vrna_zsc_filter_update(vrna_fold_compound_t *fc,
                       double               min_z,
                       unsigned int         options);


void
vrna_zsc_filter_free(vrna_fold_compound_t *fc);


int
vrna_zsc_filter_on(vrna_fold_compound_t *fc);


double
vrna_zsc_filter_threshold(vrna_fold_compound_t *fc);


double
vrna_zsc_compute(vrna_fold_compound_t *fc,
                 unsigned int         i,
                 unsigned int         j,
                 int                  e);


double
vrna_zsc_compute_raw(vrna_fold_compound_t *fc,
                     unsigned int         i,
                     unsigned int         j,
                     int                  e,
                     double               *avg,
                     double               *sd);


#endif
