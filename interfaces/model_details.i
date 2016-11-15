/**********************************************/
/* BEGIN interface for model details          */
/**********************************************/

/* scripting language access through 'md' instead of 'vrna_md_t' */
%rename (md) vrna_md_t;

%nodefaultctor vrna_md_t;
/* hide all attributes, except for some trivial ones */
typedef struct {
  double  temperature;
  double  betaScale;
  int     dangles;
  int     special_hp;
  int     noLP;
  int     noGU;
  int     noGUclosure;
  int     logML;
  int     circ;
  int     gquad;
  int     canonicalBPonly;
  int     uniq_ML;
  int     energy_set;
  int     backtrack;
  char    backtrack_type;
  int     compute_bpp;
  char    nonstandards[64];
  int     max_bp_span;
  int     min_loop_size;
  int     window_size;
  int     oldAliEn;
  int     ribo;
  double  cv_fact;
  double  nc_fact;
  double  sfact;
  int     rtype[8];
  short   alias[MAXALPHA+1];
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

  void reset(){
    vrna_md_set_default($self);
  }
  void set_from_globals(){
    set_model_details($self);
  }
  char *option_string(){
    return vrna_md_option_string($self);
  }
}

/*
 * Hide the entire interface
 */
%ignore set_model_details;
%ignore option_string;

/*
 * Include typemaps for augmenting reading and assignment of
 * global variables with their corresponding getter/setter
 * functions in the library
 */
%include globals-md.i

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

%include <ViennaRNA/model.h>

