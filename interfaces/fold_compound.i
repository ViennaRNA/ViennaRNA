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
%include callbacks-ud.i
%include callbacks-subopt.i
%include callbacks-mfe-window.i

/* start constructing a sane interface to vrna_fold_compound_t */

/* first we remap the fold_compound type enum entries */
%rename(fc_type) my_fc_type_e;
%inline %{
  typedef enum {
    FC_TYPE_SINGLE      = VRNA_FC_TYPE_SINGLE,
    FC_TYPE_COMPARATIVE = VRNA_FC_TYPE_COMPARATIVE
  } my_fc_type_e;
%}

/* scripting language access through 'fold_compound' instead of 'vrna_fold_compound_t' */
%rename(fold_compound) vrna_fold_compound_t;

/* no default constructor / destructor */
%nodefaultctor vrna_fold_compound_t;
%nodefaultdtor vrna_fold_compound_t;

/* hide all attributes in vrna_fold_compound_t */
typedef struct {} vrna_fold_compound_t;

/* create object oriented interface for vrna_fold_compount_t */
%extend vrna_fold_compound_t {

  /* the default constructor, *md and option are optional, for single sequences*/
  vrna_fold_compound_t( const char *sequence,
                        vrna_md_t *md=NULL,
                        unsigned int options=VRNA_OPTION_DEFAULT){

    return vrna_fold_compound(sequence, md, options);
  }

  /*the constructor for alignments, *md and options are optional  */
  vrna_fold_compound_t( std::vector<std::string> alignment,
                        vrna_md_t *md=NULL,
                        unsigned int options=VRNA_OPTION_DEFAULT){

    std::vector<const char*>  vc;
    transform(alignment.begin(), alignment.end(), back_inserter(vc), convert_vecstring2veccharcp);
    vc.push_back(NULL); /* mark end of sequences */
    return vrna_fold_compound_comparative((const char **)&vc[0], md, options);
  }

  /* constructor for distance class partitioning, *md and options are, for single sequences*/
  vrna_fold_compound_t( const char *sequence,
                        char *s1,
                        char *s2,
                        vrna_md_t *md=NULL,
                        unsigned int options=VRNA_OPTION_DEFAULT){

    return vrna_fold_compound_TwoD(sequence,s1,s2, md, options);
  }

  ~vrna_fold_compound_t(){
    vrna_fold_compound_free($self);
  }
  

  vrna_fc_type_e type(){
    return $self->type;
  }
  
  unsigned int length(){
    return $self->length;
  }

  /*  ################
      # in centroid.h
      ################
  */

  char *centroid(double *OUTPUT){
    return vrna_centroid($self,OUTPUT);
  }

}

%newobject vrna_fold_compound_t::centroid;

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
