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
%immutable;
typedef struct {
  const vrna_fc_type_e  type;
  char                  *sequence;
  unsigned int          length;
  unsigned int          strands;
  vrna_param_t          *params;
  vrna_exp_param_t      *exp_params;
  vrna_mx_mfe_t         *matrices;
  vrna_mx_pf_t          *exp_matrices;
  vrna_hc_t             *hc;
} vrna_fold_compound_t;
%mutable;

/* create object oriented interface for vrna_fold_compount_t */
%extend vrna_fold_compound_t {
  var_array<unsigned int> *const strand_number;
  var_array<unsigned int> *const strand_order;
  var_array<unsigned int> *const strand_start;
  var_array<unsigned int> *const strand_end;
  var_array<int>          *const iindx;
  var_array<int>          *const jindx;

  var_array<short>        *const sequence_encoding;
  var_array<short>        *const sequence_encoding2;

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

%{
  var_array<unsigned int> *
  vrna_fold_compound_t_strand_number_get(vrna_fold_compound_t *fc)
  {
    return var_array_new(fc->length,
                         fc->strand_number,
                         VAR_ARRAY_LINEAR | VAR_ARRAY_ONE_BASED);
  }

  std::string
  vrna_fold_compound_t_sequence_get(vrna_fold_compound_t *fc)
  {
    return std::string(fc->sequence);
  }

  var_array<unsigned int> *
  vrna_fold_compound_t_strand_order_get(vrna_fold_compound_t *fc)
  {
    return var_array_new(fc->strands,
                         fc->strand_order,
                         VAR_ARRAY_LINEAR);
  }

  var_array<unsigned int> *
  vrna_fold_compound_t_strand_start_get(vrna_fold_compound_t *fc)
  {
    return var_array_new(fc->strands,
                         fc->strand_start,
                         VAR_ARRAY_LINEAR);
  }

  var_array<unsigned int> *
  vrna_fold_compound_t_strand_end_get(vrna_fold_compound_t *fc)
  {
    return var_array_new(fc->strands,
                         fc->strand_end,
                         VAR_ARRAY_LINEAR);
  }

  var_array<int> *
  vrna_fold_compound_t_iindx_get(vrna_fold_compound_t *fc)
  {
    return var_array_new(fc->length,
                         fc->iindx,
                         VAR_ARRAY_LINEAR | VAR_ARRAY_ONE_BASED);
  }

  var_array<int> *
  vrna_fold_compound_t_jindx_get(vrna_fold_compound_t *fc)
  {
    if (fc->type == VRNA_FC_TYPE_SINGLE)
      return var_array_new(fc->length,
                           fc->jindx,
                           VAR_ARRAY_LINEAR | VAR_ARRAY_ONE_BASED);

    return NULL;
  }

  var_array<short> *
  vrna_fold_compound_t_sequence_encoding_get(vrna_fold_compound_t *fc)
  {
    if (fc->type == VRNA_FC_TYPE_SINGLE)
      return var_array_new(fc->length + 1,
                           fc->sequence_encoding,
                           VAR_ARRAY_LINEAR | VAR_ARRAY_ONE_BASED);

    return NULL;
  }

  var_array<short> *
  vrna_fold_compound_t_sequence_encoding2_get(vrna_fold_compound_t *fc)
  {
    return var_array_new(fc->length + 1,
                         fc->sequence_encoding2,
                         VAR_ARRAY_LINEAR | VAR_ARRAY_ONE_BASED);
  }
%}

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
