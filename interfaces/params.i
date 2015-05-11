/**********************************************/
/* BEGIN interface for energy parameters      */
/**********************************************/

/* do not create default constructor and hide data fields of vrna_param_t from SWIG */
%ignore paramT;
%ignore pf_paramT;

/* scripting language access through 'param' instead of 'vrna_param_t' */
%rename (param) vrna_param_t;

/* scripting language access through 'exp_param' instead of 'vrna_exp_param_t' */
%rename (exp_param) vrna_exp_param_t;

%nodefaultctor vrna_param_t;
typedef struct {
  vrna_md_t model_details;
} vrna_param_t;

/* do not create default constructor and hide data fields of vrna_param_t from SWIG */
%nodefaultctor vrna_exp_param_t;
typedef struct {
  double kT;
  vrna_md_t model_details;
} vrna_exp_param_t;

/* make a nice object oriented interface to vrna_param_t */
%extend vrna_param_t {
  vrna_param_t(){
    vrna_param_t *P = vrna_params(NULL);
    return P;
  }
  vrna_param_t(vrna_md_t *md){
    vrna_param_t *P = vrna_params(md);
    return P;
  }

  double  get_temperature(){
    return $self->temperature;
  }
}

/* make a nice object oriented interface to vrna_exp_param_t */
%extend vrna_exp_param_t {
  vrna_exp_param_t(){
    vrna_exp_param_t *P = vrna_exp_params(NULL);
    return P;
  }
  vrna_exp_param_t(vrna_md_t *md){
    vrna_exp_param_t *P = vrna_exp_params(md);
    return P;
  }

  double  get_temperature(){
    return $self->temperature;
  }
}

%ignore get_parameter_copy;
%ignore get_scaled_pf_parameters;
%ignore get_boltzmann_factors;
%ignore get_boltzmann_factor_copy;
%ignore get_scaled_alipf_parameters;
%ignore get_boltzmann_factors_ali;
%ignore scale_parameters;
%ignore get_scaled_parameters;
%ignore copy_parameters;
%ignore set_parameters;
%ignore scale_pf_parameters;
%ignore copy_pf_param;
%ignore set_pf_param;

%include <ViennaRNA/params.h>



/**********************************************/
/* BEGIN interface for parameter file I/O     */
/**********************************************/

%include <ViennaRNA/read_epars.h>

/**********************************************/
/* BEGIN interface for energy constants       */
/**********************************************/

%include  <ViennaRNA/energy_const.h>
