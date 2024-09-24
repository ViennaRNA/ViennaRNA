struct vrna_zsc_dat_s {
  struct svm_model  *avg_model;
  struct svm_model  *sd_model;
  double            min_z;
  unsigned char     filter_on;
  double            *current_z;
  unsigned int      current_i;
  unsigned char     pre_filter;
  unsigned char     report_subsumed;
};
