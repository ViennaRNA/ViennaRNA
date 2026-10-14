typedef struct helper_arrays *helper_arrays_t;

typedef void  (*component_f)(vrna_fold_compound_t *fc,
                             unsigned int         l,
                             helper_arrays_t      eqm_prob_helpers,
                             void                 *data);

typedef void (*post_f)(vrna_fold_compound_t *fc,
                       unsigned int         l,
                       helper_arrays_t      eqm_prob_helpers,
                       void                 *data);

typedef void (*pre_f)(vrna_fold_compound_t  *fc,
                      helper_arrays_t       eqm_prob_helpers,
                      void                  *data);

struct vrna_eqm_probs_component_s {
  component_f          cb;
  void *               cb_data;
  vrna_auxdata_free_f  cb_data_free;
  pre_f                cb_init;
  post_f               cb_post;
};


struct helper_arrays {
  unsigned int  n;
  unsigned int  strands;

  /* overflow counter */
  FLT_OR_DBL    Qmax;
  int           ov;

  /* soft constraints wrappers */
  struct sc_ext_exp_dat     sc_wrapper_ext;
  struct sc_hp_exp_dat      sc_wrapper_hp;
  struct sc_int_exp_dat     sc_wrapper_int;
  struct sc_mb_exp_dat      sc_wrapper_mb;

  FLT_OR_DBL    *probs_ext;
  FLT_OR_DBL    *probs_hp;
  FLT_OR_DBL    *probs_int;
  FLT_OR_DBL    *probs_mb;

  FLT_OR_DBL    *qbh;
  FLT_OR_DBL    *access_ext;
  FLT_OR_DBL    *access_hp;
  FLT_OR_DBL    *access_int;
  FLT_OR_DBL    *access_mb;

  FLT_OR_DBL    *hxp;

  FLT_OR_DBL    *qm2_aux;

  vrna_eqm_probs_f       cb;
  vrna_eqm_probs_data_t  cb_data;
  vrna_eqm_probs_opt_t   options;

  vrna_array(vrna_ep_t)  bp_correction;

  vrna_array(component_f)         component_cbs;
  vrna_array(void *)              component_cbs_data;
  vrna_array(vrna_auxdata_free_f) component_cbs_data_free;
  vrna_array(pre_f)               component_cbs_init;
  vrna_array(post_f)              component_cbs_post;
};

