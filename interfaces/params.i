/**********************************************/
/* BEGIN interface for energy parameters      */
/**********************************************/

/* do not create default constructor and hide data fields of vrna_param_t from SWIG */
%ignore paramT;
%ignore pf_paramT;
%ignore vrna_param_s;
%ignore vrna_exp_param_s;

%nodefaultctor vrna_param_t;
typedef struct {} vrna_param_t;
/* do not create default constructor and hide data fields of vrna_param_t from SWIG */
%nodefaultctor vrna_exp_param_t;
typedef struct {} vrna_exp_param_t;

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

%ignore vrna_params_get;
%ignore vrna_exp_params_get;
%ignore scale_parameters;
%ignore vrna_get_energy_contributions;
%ignore get_scaled_parameters;
%ignore get_parameter_copy;
%ignore get_scaled_pf_parameters;
%ignore vrna_get_boltzmann_factors;
%ignore get_boltzmann_factors;
%ignore get_boltzmann_factor_copy;
/*
%ignore get_scaled_alipf_parameters;
%ignore get_boltzmann_factors_ali;
*/
%ignore copy_parameters;
%ignore set_parameters;
%ignore scale_pf_parameters;
%ignore copy_pf_param;
%ignore set_pf_param;

%ignore vrna_params_update;
%ignore vrna_exp_params_update;
%ignore vrna_exp_params_rescale;

%include "../src/ViennaRNA/params.h"



/**********************************************/
/* BEGIN interface for parameter file I/O     */
/**********************************************/

%include "../src/ViennaRNA/read_epars.h"

/**********************************************/
/* BEGIN interface for energy constants       */
/**********************************************/

%include  "../src/ViennaRNA/energy_const.h"
