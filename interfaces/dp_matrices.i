/**********************************************/
/* BEGIN interface for fold compound          */
/**********************************************/

%immutable;
%rename (mx_mfe)  vrna_mx_mfe_t;

%nodefaultctor vrna_mx_mfe_t;
%nodefaultdtor vrna_mx_mfe_t;

/* hide most attributes in vrna_mx_mfe_t */
typedef struct {
  const vrna_mx_type_e  type;
  unsigned int          length;
  unsigned int          strands;
} vrna_mx_mfe_t;

%mutable;

%extend vrna_mx_mfe_t {
  /* expose DP matrices */
  var_array<int> *const f5;
  var_array<int> *const f3;
  var_array<int> *const c;
  var_array<int> *const fML;
  var_array<int> *const fM1;
  var_array<int> *const fM2;
  var_array<int> *const ggg;
  const int             Fc;
  const int             FcH;
  const int             FcI;
  const int             FcM;

  vrna_mx_mfe_t() { return NULL; }
  ~vrna_mx_mfe_t() {}
};


%{
  var_array<int> *
  vrna_mx_mfe_t_f5_get(vrna_mx_mfe_t *mx)
  {
    return var_array_lin_int_new(mx->length, mx->f5);
  }

  var_array<int> *
  vrna_mx_mfe_t_f3_get(vrna_mx_mfe_t *mx)
  {
    return var_array_lin_int_new(mx->length, mx->f3);
  }

  var_array<int> *
  vrna_mx_mfe_t_c_get(vrna_mx_mfe_t *mx)
  {
    return var_array_tri_int_new(mx->length, mx->c);
  }

  var_array<int> *
  vrna_mx_mfe_t_fML_get(vrna_mx_mfe_t *mx)
  {
    return var_array_tri_int_new(mx->length, mx->fML);
  }

  var_array<int> *
  vrna_mx_mfe_t_fM1_get(vrna_mx_mfe_t *mx)
  {
    return var_array_tri_int_new(mx->length, mx->fM1);
  }

  var_array<int> *
  vrna_mx_mfe_t_ggg_get(vrna_mx_mfe_t *mx)
  {
    return var_array_tri_int_new(mx->length, mx->ggg);
  }

  var_array<int> *
  vrna_mx_mfe_t_fM2_get(vrna_mx_mfe_t *mx)
  {
    return var_array_lin_int_new(mx->length, mx->fM2);
  }

  const int
  vrna_mx_mfe_t_Fc_get(vrna_mx_mfe_t *mx)
  {
    return mx->Fc;
  }

  const int
  vrna_mx_mfe_t_FcH_get(vrna_mx_mfe_t *mx)
  {
    return mx->FcH;
  }

  const int
  vrna_mx_mfe_t_FcI_get(vrna_mx_mfe_t *mx)
  {
    return mx->FcI;
  }

  const int
  vrna_mx_mfe_t_FcM_get(vrna_mx_mfe_t *mx)
  {
    return mx->FcM;
  }
%}

%immutable;
%rename (mx_pf)  vrna_mx_pf_t;

%nodefaultctor vrna_mx_pf_t;
%nodefaultdtor vrna_mx_pf_t;

/* hide most attributes in vrna_mx_pf_t */
typedef struct {
  const vrna_mx_type_e  type;
  unsigned int          length;
} vrna_mx_pf_t;

%mutable;


%extend vrna_mx_pf_t {
  /* expose DP matrices */
  var_array<FLT_OR_DBL> *const  scale;
  var_array<FLT_OR_DBL> *const  expMLbase;
  var_array<FLT_OR_DBL> *const  q;
  var_array<FLT_OR_DBL> *const  qb;
  var_array<FLT_OR_DBL> *const  qm;
  var_array<FLT_OR_DBL> *const  qm1;
  var_array<FLT_OR_DBL> *const  probs;
  var_array<FLT_OR_DBL> *const  q1k;
  var_array<FLT_OR_DBL> *const  qln;
  var_array<FLT_OR_DBL> *const  G;
  const FLT_OR_DBL              qo;
  var_array<FLT_OR_DBL> *const  qm2;
  const FLT_OR_DBL              qho;
  const FLT_OR_DBL              qio;
  const FLT_OR_DBL              qmo;

  vrna_mx_pf_t() { return NULL; }
  ~vrna_mx_pf_t() {}
};


%{
  var_array<FLT_OR_DBL> *
  vrna_mx_pf_t_scale_get(vrna_mx_pf_t *mx)
  {
    return var_array_lin_dbl_new(mx->length, mx->scale);
  }

  var_array<FLT_OR_DBL> *
  vrna_mx_pf_t_expMLbase_get(vrna_mx_pf_t *mx)
  {
    return var_array_lin_dbl_new(mx->length, mx->expMLbase);
  }

  var_array<FLT_OR_DBL> *
  vrna_mx_pf_t_qm1_get(vrna_mx_pf_t *mx)
  {
    return var_array_tri_dbl_new(mx->length, mx->qm1);
  }

  var_array<FLT_OR_DBL> *
  vrna_mx_pf_t_q_get(vrna_mx_pf_t *mx)
  {
    return var_array_tri_dbl_new(mx->length, mx->q);
  }

  var_array<FLT_OR_DBL> *
  vrna_mx_pf_t_qb_get(vrna_mx_pf_t *mx)
  {
    return var_array_tri_dbl_new(mx->length, mx->qb);
  }

  var_array<FLT_OR_DBL> *
  vrna_mx_pf_t_qm_get(vrna_mx_pf_t *mx)
  {
    return var_array_tri_dbl_new(mx->length, mx->qm);
  }

  var_array<FLT_OR_DBL> *
  vrna_mx_pf_t_G_get(vrna_mx_pf_t *mx)
  {
    return var_array_tri_dbl_new(mx->length, mx->G);
  }

  var_array<FLT_OR_DBL> *
  vrna_mx_pf_t_probs_get(vrna_mx_pf_t *mx)
  {
    return var_array_tri_dbl_new(mx->length, mx->probs);
  }

  var_array<FLT_OR_DBL> *
  vrna_mx_pf_t_q1k_get(vrna_mx_pf_t *mx)
  {
    return var_array_lin_dbl_new(mx->length, mx->q1k);
  }

  var_array<FLT_OR_DBL> *
  vrna_mx_pf_t_qln_get(vrna_mx_pf_t *mx)
  {
    return var_array_lin_dbl_new(mx->length, mx->qln);
  }

  var_array<FLT_OR_DBL> *
  vrna_mx_pf_t_qm2_get(vrna_mx_pf_t *mx)
  {
    return var_array_lin_dbl_new(mx->length, mx->qm2);
  }

  const FLT_OR_DBL
  vrna_mx_pf_t_qo_get(vrna_mx_pf_t *mx)
  {
    return mx->qo;
  }

  const FLT_OR_DBL
  vrna_mx_pf_t_qho_get(vrna_mx_pf_t *mx)
  {
    return mx->qho;
  }

  const FLT_OR_DBL
  vrna_mx_pf_t_qio_get(vrna_mx_pf_t *mx)
  {
    return mx->qio;
  }

  const FLT_OR_DBL
  vrna_mx_pf_t_qmo_get(vrna_mx_pf_t *mx)
  {
    return mx->qmo;
  }
%}


%include <ViennaRNA/dp_matrices.h>
