/**********************************************/
/* BEGIN interface for fold compound          */
/**********************************************/

/* ignore all data structures we handle elsewhere */
%ignore PAIR;
%ignore plist;
%ignore cpair;
%ignore sect;
%ignore bondT;
%ignore pu_contrib;
%ignore interact;
%ignore pu_out;
%ignore constrain;
//%ignore duplexT;
%ignore folden;
%ignore snoopT;
%ignore dupVar;

/* start constructing a sane interface to vrna_fold_compound_t */

%rename(fc_type) vrna_fc_type_e;

/* scripting language access through 'fold_compound' instead of 'vrna_fold_compound_t' */
%rename(fold_compound) vrna_fold_compound_t;

/* no default constructor / destructor */
%nodefaultctor vrna_fold_compound_t;
%nodefaultdtor vrna_fold_compound_t;

/* hide all attributes in vrna_fold_compound_t */
typedef struct {} vrna_fold_compound_t;

/* scripting language takes ownership of objects returned by mfe() method */
%newobject vrna_fold_compound_t::mfe;

/* create object oriented interface for vrna_fold_compount_t */
%extend vrna_fold_compound_t {

  /* the default constructor */
  vrna_fold_compound_t(char *sequence, vrna_md_t *md=NULL, unsigned int options=VRNA_OPTION_MFE){
    return vrna_fold_compound(sequence, md, options);
  }
  ~vrna_fold_compound_t(){
    vrna_fold_compound_free($self);
  }

  vrna_fc_type_e type(){
    return $self->type;
  }

  char *mfe(float *OUTPUT){
    char *structure = vrna_alloc(sizeof(char) * ($self->length + 1));
    *OUTPUT = vrna_mfe($self, structure);
    return structure;
  }
}

%include "../src/ViennaRNA/data_structures.h"
