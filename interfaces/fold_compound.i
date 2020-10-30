/**********************************************/
/* BEGIN interface for fold compound          */
/**********************************************/

/* add callback binding methods for fold_compound */
%include callbacks-fc.i
%include callbacks-sc.i
%include callbacks-ud.i
%include callbacks-subopt.i
%include callbacks-boltzmann-sampling.i
%include callbacks-mfe-window.i
%include callbacks-pf-window.i
%include callbacks-melting.i

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

/*
  hide almost all attributes in vrna_fold_compound_t
  non-hidden attributes are made read-only using 'const' 
*/
typedef struct {
  const vrna_fc_type_e      type;
  const unsigned int        length;
  const unsigned int        strands;
  vrna_param_t      *const  params;
  vrna_exp_param_t  *const  exp_params;
} vrna_fold_compound_t;

/* create object oriented interface for vrna_fold_compount_t */
%extend vrna_fold_compound_t {

#ifdef SWIGPYTHON
%feature("autodoc")vrna_fold_compound_t::vrna_fold_compound_t;
%feature("kwargs")vrna_fold_compound_t::vrna_fold_compound_t;
#endif
  /* the default constructor, *md and option are optional, for single sequences*/
  vrna_fold_compound_t( const char    *sequence,
                        vrna_md_t     *md = NULL,
                        unsigned int  options = VRNA_OPTION_DEFAULT)
  {
    return vrna_fold_compound(sequence, md, options);
  }

  /*the constructor for alignments, *md and options are optional  */
  vrna_fold_compound_t( std::vector<std::string>  alignment,
                        vrna_md_t                 *md = NULL,
                        unsigned int              options = VRNA_OPTION_DEFAULT)
  {
    std::vector<const char*>  vc;
    transform(alignment.begin(), alignment.end(), back_inserter(vc), convert_vecstring2veccharcp);
    vc.push_back(NULL); /* mark end of sequences */
    return vrna_fold_compound_comparative((const char **)&vc[0], md, options);
  }

  /* constructor for distance class partitioning, *md and options are, for single sequences*/
  vrna_fold_compound_t( const char    *sequence,
                        char          *s1,
                        char          *s2,
                        vrna_md_t     *md = NULL,
                        unsigned int  options=VRNA_OPTION_DEFAULT)
  {
    return vrna_fold_compound_TwoD(sequence,s1,s2, md, options);
  }

  ~vrna_fold_compound_t()
  {
    vrna_fold_compound_free($self);
  }

#ifdef SWIGPYTHON
 std::string
  __str__()
  {
    std::ostringstream out;

    out << "{ ";

    if ($self->type == VRNA_FC_TYPE_SINGLE) {
      out << "sequence: \"" << $self->sequence << "\"";
    } else {
      out << "sequences: (" << "\"" << $self->sequences[0] << "\"";
      for (size_t i = 1; i < $self->n_seq; i++)
        out << ", \"" << $self->sequences[i] << "\"";
      out << ")";
    }
    out << ", length: " << $self->length;
    out << ", strands: " << $self->strands;
    out << " }";

    return std::string(out.str());
  }

%pythoncode %{
def __repr__(self):
    # reformat string representation (self.__str__()) to something
    # that looks like a constructor argument list
    strthis = self.__str__().replace(": ", "=").replace("{ ", "").replace(" }", "")
    return  "%s.%s(%s)" % (self.__class__.__module__, self.__class__.__name__, strthis) 
%}
#endif
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

%include <ViennaRNA/fold_compound.h>
