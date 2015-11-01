/**********************************************/
/* BEGIN interface for model details          */
/**********************************************/

/* hide data fields of vrna_md_t from SWIG */
%ignore vrna_md_s;

%nodefaultctor vrna_md_t;
/* hide all attributes, except for some trivial ones */
typedef struct {
  double  temperature;
  int     dangles;
  int     noLP;
  int     noGU;
  int     noGUclosure;
  int     special_hp;
  int     max_bp_span;
} vrna_md_t;

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
  char *option_string(){
    return vrna_md_option_string($self);
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

/*
 * Hide the entire interface
 */
%ignore vrna_md_set_default;
%ignore vrna_md_update;
%ignore vrna_md_option_string;
%ignore vrna_md_set_nonstandards;
%ignore set_model_details;
%ignore option_string;
%ignore vrna_md_defaults_temperature;
%ignore vrna_md_defaults_temperature_get;
%ignore vrna_md_defaults_betaScale;
%ignore vrna_md_defaults_betaScale_get;
%ignore vrna_md_defaults_dangles;
%ignore vrna_md_defaults_dangles_get;
%ignore vrna_md_defaults_special_hp;
%ignore vrna_md_defaults_special_hp_get;
%ignore vrna_md_defaults_noLP;
%ignore vrna_md_defaults_noLP_get;
%ignore vrna_md_defaults_noGU;
%ignore vrna_md_defaults_noGU_get;
%ignore vrna_md_defaults_noGUclosure;
%ignore vrna_md_defaults_noGUclosure_get;
%ignore vrna_md_defaults_logML;
%ignore vrna_md_defaults_logML_get;
%ignore vrna_md_defaults_circ;
%ignore vrna_md_defaults_circ_get;
%ignore vrna_md_defaults_gquad;
%ignore vrna_md_defaults_gquad_get;
%ignore vrna_md_defaults_uniq_ML;
%ignore vrna_md_defaults_uniq_ML_get;
%ignore vrna_md_defaults_energy_set;
%ignore vrna_md_defaults_energy_set_get;
%ignore vrna_md_defaults_backtrack;
%ignore vrna_md_defaults_backtrack_get;
%ignore vrna_md_defaults_backtrack_type;
%ignore vrna_md_defaults_backtrack_type_get;
%ignore vrna_md_defaults_compute_bpp;
%ignore vrna_md_defaults_compute_bpp_get;
%ignore vrna_md_defaults_max_bp_span;
%ignore vrna_md_defaults_max_bp_span_get;
%ignore vrna_md_defaults_min_loop_size;
%ignore vrna_md_defaults_min_loop_size_get;
%ignore vrna_md_defaults_window_size;
%ignore vrna_md_defaults_window_size_get;
%ignore vrna_md_defaults_oldAliEn;
%ignore vrna_md_defaults_oldAliEn_get;
%ignore vrna_md_defaults_ribo;
%ignore vrna_md_defaults_ribo_get;
%ignore vrna_md_defaults_cv_fact;
%ignore vrna_md_defaults_cv_fact_get;
%ignore vrna_md_defaults_nc_fact;
%ignore vrna_md_defaults_nc_fact_get;
%ignore vrna_md_defaults_sfact;
%ignore vrna_md_defaults_sfact_get;

/*
 * Include typemaps for augmenting reading and assignment of
 * global variables with their corresponding getter/setter
 * functions in the library
 */
%include md_globals_tmaps.i

/*
 * The list of available global variables in the scripting
 * interface. Access to these variables is defined in the
 * typemap definition file above
 */
extern double temperature;
extern int    dangles;
extern double betaScale;
extern int    tetra_loop;     /* this is an alias of special_hp */
extern int    special_hp;
extern int    noLonelyPairs;  /* this is an alias of noLP */
extern int    noLP;
extern int    noGU;
extern int    no_closingGU;    /* this is an alias of noGUclosure */
extern int    noGUclosure;
extern int    logML;
extern int    circ;
extern int    gquad;
extern int    uniq_ML;
extern int    energy_set;
extern int    backtrack;
extern char   backtrack_type; /* still needs implementation */
extern int    do_backtrack;   /* this is an alias of compute_bpp */
extern int    compute_bpp;
extern int    max_bp_span;
extern int    min_loop_size;
extern int    window_size;
extern int    oldAliEn;
extern int    ribo;
extern double cv_fact;
extern double nc_fact;
extern double sfact;

%include "../src/ViennaRNA/model.h"

