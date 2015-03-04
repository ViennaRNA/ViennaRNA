/**********************************************/
/* BEGIN interface for energy parameters      */
/**********************************************/

/* do not create default constructor and hide data fields of vrna_param_t from SWIG */
%nodefaultctor vrna_param_t;
typedef struct {} vrna_param_t;
/* do not create default constructor and hide data fields of vrna_param_t from SWIG */
%nodefaultctor vrna_exp_param_t;
typedef struct {} vrna_exp_param_t;

/* make a nice object oriented interface to vrna_param_t */
%extend vrna_param_t {
  vrna_param_t(){
    vrna_param_t *P = vrna_params_get(NULL);
    return P;
  }
  vrna_param_t(vrna_md_t *md){
    vrna_param_t *P = vrna_params_get(md);
    return P;
  }
}

/* make a nice object oriented interface to vrna_exp_param_t */
%extend vrna_exp_param_t {
  vrna_exp_param_t(){
    vrna_exp_param_t *P = vrna_exp_params_get(NULL);
    return P;
  }
  vrna_exp_param_t(vrna_md_t *md){
    vrna_exp_param_t *P = vrna_exp_params_get(md);
    return P;
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

%include "../src/ViennaRNA/params.h"



/**********************************************/
/* BEGIN interface for parameter file I/O     */
/**********************************************/

%include "../src/ViennaRNA/read_epars.h"

/**********************************************/
/* BEGIN interface for model details          */
/**********************************************/

/* hide data fields of vrna_md_t from SWIG */
%nodefaultctor vrna_md_t;
// typedef struct {} vrna_md_t;

/* make a nice object oriented interface to vrna_md_t */
%extend vrna_md_t {

  /* the default constructor */
  vrna_md_t(){
    vrna_md_t *md = (vrna_md_t *)vrna_alloc(sizeof(vrna_md_t));
    vrna_md_set_default(md);
    return md;
  }
  /* a constructor that provides backward compatibility (for now) */
  vrna_md_t(char *type){
    vrna_md_t *md = (vrna_md_t *)vrna_alloc(sizeof(vrna_md_t));
    if(!strcmp(type, "global"))
      vrna_md_set_globals(md);
    else
      vrna_md_set_default(md);
    return md;
  }

  int get_dangles(){
    return $self->dangles;
  }
  void set_dangles(int d){
    $self->dangles = d;
  }
  void set_default(){
    vrna_md_set_default($self);
  }
  void set_from_globals(){
    vrna_md_set_globals($self);
  }
  void print(){
    printf( "temperature:     %g\n"
            "betaScale:       %g\n"
            "sfact:           %g\n"
            "dangles:         %d\n"
            "special_hp:      %d\n"
            "noLP:            %d\n"
            "noGU:            %d\n"
            "noGUclosure:     %d\n"
            "logML:           %d\n"
            "circ:            %d\n"
            "gquad:           %d\n"
            "canonicalBPonly: %d\n"
            "uniq_ML:         %d\n"
            "energy_set:      %d\n"
            "backtrack:       %d\n"
            "backtrack_type:  %c\n"
            "compute_bpp:     %d\n"
            "nonstandards:    %s\n"
            "max_bp_span:     %d\n"
            "min_loop_size:   %d\n"
            "oldAliEn:        %d\n"
            "ribo:            %d\n",
            $self->temperature,
            $self->betaScale,
            $self->sfact,
            $self->dangles,
            $self->special_hp,
            $self->noLP,
            $self->noGU,
            $self->noGUclosure,
            $self->logML,
            $self->circ,
            $self->gquad,
            $self->canonicalBPonly,
            $self->uniq_ML,
            $self->energy_set,
            $self->backtrack,
            $self->backtrack_type,
            $self->compute_bpp,
            $self->nonstandards,
            $self->max_bp_span,
            $self->min_loop_size,
            $self->oldAliEn,
            $self->ribo);
  }
}

%ignore vrna_md_set_default;
%ignore vrna_md_set_globals;
%ignore vrna_md_set_nonstandards;
%ignore vrna_md_set_dangles;
%ignore vrna_md_get_dangles;
%ignore vrna_md_set_temperature;
%ignore vrna_md_get_temperature;
%ignore vrna_md_set_special_hp;
%ignore vrna_md_get_special_hp;
%ignore vrna_md_set_gquad;
%ignore vrna_md_get_gquad;

%include "../src/ViennaRNA/model.h"

