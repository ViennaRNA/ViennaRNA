/**********************************************/
/* BEGIN interface for energy parameters      */
/**********************************************/

/* do not create default constructor and hide data fields of paramT from SWIG */
%nodefaultctor paramT;
typedef struct {} paramT;
/* do not create default constructor and hide data fields of paramT from SWIG */
%nodefaultctor pf_paramT;
typedef struct {} pf_paramT;

/* make a nice object oriented interface to paramT */
%extend paramT {
  paramT(){
    model_detailsT md;
    vrna_md_set_default(&md);
    paramT *P = vrna_get_energy_contributions(md);
    return P;
  }
  paramT(model_detailsT *md){
    paramT *P = vrna_get_energy_contributions(*md);
    return P;
  }
}

/* make a nice object oriented interface to pf_paramT */
%extend pf_paramT {
  pf_paramT(){
    model_detailsT md;
    vrna_md_set_default(&md);
    pf_paramT *P = vrna_get_boltzmann_factors(md);
    return P;
  }
  pf_paramT(model_detailsT *md){
    pf_paramT *P = vrna_get_boltzmann_factors(*md);
    return P;
  }
}

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

/* hide data fields of model_detailsT from SWIG */
%nodefaultctor model_detailsT;
// typedef struct {} model_detailsT;

/* make a nice object oriented interface to model_detailsT */
%extend model_detailsT {

  /* the default constructor */
  model_detailsT(){
    model_detailsT *md = (model_detailsT *)space(sizeof(model_detailsT));
    vrna_md_set_default(md);
    return md;
  }
  /* a constructor that provides backward compatibility (for now) */
  model_detailsT(char *type){
    model_detailsT *md = (model_detailsT *)space(sizeof(model_detailsT));
    if(!strcmp(type, "global"))
      set_model_details(md);
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
    set_model_details($self);
  }
  void print(){
    printf( "temperature:     %g\n"
            "betaScale:       %g\n"
            "pf_scale:        %g\n"
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
            $self->pf_scale,
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
%ignore set_model_details;
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

