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

/* add callback binding methods for fold_compound */
%include callbacks-fc.i
%include callbacks-sc.i

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
    char *structure = (char *)vrna_alloc(sizeof(char) * ($self->length + 1));
    *OUTPUT = vrna_mfe($self, structure);
    return structure;
  }

  void sc_remove(){
    vrna_sc_remove($self);
  }

  void sc_add_up(const double *constraints, unsigned int options=VRNA_OPTION_MFE){
    vrna_sc_add_up($self, constraints, options);
  }

  void sc_add_bp(const double **constraints, unsigned int options=VRNA_OPTION_MFE){
    vrna_sc_add_bp($self, constraints, options);
  }

  int sc_add_hi_motif(const char *seq,
                      const char *structure,
                      double energy,
                      unsigned int options=VRNA_OPTION_MFE){

    return vrna_sc_add_hi_motif($self, seq, structure, energy, options);
  }

}


/*
 *  Rename all the preprocessor macros defined in data_structures.h
 *  (wrapped as constants)
 */
%constant unsigned char STATUS_MFE_PRE  = VRNA_STATUS_MFE_PRE;
%constant unsigned char STATUS_MFE_POST = VRNA_STATUS_MFE_POST;
%constant unsigned char STATUS_PF_PRE   = VRNA_STATUS_PF_PRE;
%constant unsigned char STATUS_PF_POST  = VRNA_STATUS_PF_POST;

%constant unsigned int OPTION_DEFAULT   = VRNA_OPTION_DEFAULT;
%constant unsigned int OPTION_MFE       = VRNA_OPTION_MFE;
%constant unsigned int OPTION_PF        = VRNA_OPTION_PF;
%constant unsigned int OPTION_HYBRID    = VRNA_OPTION_HYBRID;
%constant unsigned int OPTION_EVAL_ONLY = VRNA_OPTION_EVAL_ONLY;
%constant unsigned int OPTION_WINDOW    = VRNA_OPTION_WINDOW;

%rename(basepair) vrna_basepair_t;

typedef struct {
  int i;
  int j;
} vrna_basepair_t;

%include <ViennaRNA/data_structures.h>
