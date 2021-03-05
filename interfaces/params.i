/**********************************************/
/* BEGIN interface for energy parameters      */
/**********************************************/
%{
#include <sstream>
%}

/* do not create default constructor */
%ignore paramT;
%ignore pf_paramT;

/* scripting language access through 'param' instead of 'vrna_param_t' */
%rename (param) vrna_param_t;

/* scripting language access through 'exp_param' instead of 'vrna_exp_param_t' */
%rename (exp_param) vrna_exp_param_t;

%nodefaultctor vrna_param_t;
typedef struct {
  const int       id;
  const int       stack[NBPAIRS + 1][NBPAIRS + 1];
  const int       hairpin[31];
  const int       bulge[MAXLOOP + 1];
  const int       internal_loop[MAXLOOP + 1];
  const int       mismatchExt[NBPAIRS + 1][5][5];
  const int       mismatchI[NBPAIRS + 1][5][5];
  const int       mismatch1nI[NBPAIRS + 1][5][5];
  const int       mismatch23I[NBPAIRS + 1][5][5];
  const int       mismatchH[NBPAIRS + 1][5][5];
  const int       mismatchM[NBPAIRS + 1][5][5];
  const int       dangle5[NBPAIRS + 1][5];
  const int       dangle3[NBPAIRS + 1][5];
  const int       int11[NBPAIRS + 1][NBPAIRS + 1][5][5];
  const int       int21[NBPAIRS + 1][NBPAIRS + 1][5][5][5];
  const int       int22[NBPAIRS + 1][NBPAIRS + 1][5][5][5][5];
  const int       ninio[5];
  const double    lxc;
  const int       MLbase;
  const int       MLintern[NBPAIRS + 1];
  const int       MLclosing;
  const int       TerminalAU;
  const int       DuplexInit;
  const int       Tetraloop_E[200];
  const char      Tetraloops[1401];
  const int       Triloop_E[40];
  const char      Triloops[241];
  const int       Hexaloop_E[40];
  const char      Hexaloops[1801];
  const int       TripleC;
  const int       MultipleCA;
  const int       MultipleCB;
  const int       gquad[VRNA_GQUAD_MAX_STACK_SIZE + 1][3 * VRNA_GQUAD_MAX_LINKER_LENGTH + 1];
  const int       gquadLayerMismatch;
  const int       gquadLayerMismatchMax;

  const double    temperature;

  const vrna_md_t model_details;
  const char      param_file[256];
} vrna_param_t;

/* do not create default constructor */
%nodefaultctor vrna_exp_param_t;

typedef struct {
  const int     id;
  const double  expstack[NBPAIRS + 1][NBPAIRS + 1];
  const double  exphairpin[31];
  const double  expbulge[MAXLOOP + 1];
  const double  expinternal[MAXLOOP + 1];
  const double  expmismatchExt[NBPAIRS + 1][5][5];
  const double  expmismatchI[NBPAIRS + 1][5][5];
  const double  expmismatch23I[NBPAIRS + 1][5][5];
  const double  expmismatch1nI[NBPAIRS + 1][5][5];
  const double  expmismatchH[NBPAIRS + 1][5][5];
  const double  expmismatchM[NBPAIRS + 1][5][5];
  const double  expdangle5[NBPAIRS + 1][5];
  const double  expdangle3[NBPAIRS + 1][5];
  const double  expint11[NBPAIRS + 1][NBPAIRS + 1][5][5];
  const double  expint21[NBPAIRS + 1][NBPAIRS + 1][5][5][5];
  const double  expint22[NBPAIRS + 1][NBPAIRS + 1][5][5][5][5];
  const double  expninio[5][MAXLOOP + 1];
  const double  lxc;
  const double  expMLbase;
  const double  expMLintern[NBPAIRS + 1];
  const double  expMLclosing;
  const double  expTermAU;
  const double  expDuplexInit;
  const double  exptetra[40];
  const double  exptri[40];
  const double  exphex[40];
  const char    Tetraloops[1401];
  const double  expTriloop[40];
  const char    Triloops[241];
  const char    Hexaloops[1801];
  const double  expTripleC;
  const double  expMultipleCA;
  const double  expMultipleCB;
  const double  expgquad[VRNA_GQUAD_MAX_STACK_SIZE + 1][3 * VRNA_GQUAD_MAX_LINKER_LENGTH + 1];
  const double  expgquadLayerMismatch;
  const int     gquadLayerMismatchMax;

  const double  kT;
  const double  pf_scale;

  const double  temperature;
  const double  alpha;

  const vrna_md_t model_details;
  const char      param_file[256];
} vrna_exp_param_t;

/* make a nice object oriented interface to vrna_param_t */
%extend vrna_param_t {
  vrna_param_t(vrna_md_t *model_details = NULL)
  {
    return vrna_params(model_details);
  }

#ifdef SWIGPYTHON
  std::string
  __str__()
  {
    std::ostringstream out;
    out << "{ model_details: RNA.md()";
    out << ", id: " << $self->id;
    out << ", param_file: \"" << $self->param_file << "\"";
    out << ", temperature: " << $self->temperature;
    out << ", TerminalAU: " << $self->TerminalAU;
    out << ", DuplexInit: " << $self->DuplexInit;

    out << ", MLclosing: " << $self->MLclosing;
    out << ", MLbase: " << $self->MLbase;
    out << ", MLintern: [" << $self->MLintern[0];
    for (size_t i = 1; i < NBPAIRS + 1; i++)
      out << ", " << $self->MLintern[i];
    out << "]";

    out << ", hairpin: [" << $self->hairpin[0];
    for (size_t i = 1; i < 31; i++)
      out << ", " << $self->hairpin[i];
    out << "]";

    out << ", bulge: [" << $self->bulge[0];
    for (size_t i = 1; i < MAXLOOP + 1; i++)
      out << ", " << $self->bulge[i];
    out << "]";

    out << ", internal_loop: [" << $self->internal_loop[0];
    for (size_t i = 1; i < 31; i++)
      out << ", " << $self->internal_loop[i];
    out << "]";

    out << ", stack: [[" << $self->stack[0][0];
    for (size_t a = 1; a < NBPAIRS + 1; a++)
      out << ", " << $self->stack[0][a];
    out << "]";
    for (size_t a = 1; a < NBPAIRS + 1; a++) {
      out << ", [" << $self->stack[a][0];
      for (size_t b = 1; b < NBPAIRS + 1; b++)
        out << ", " << $self->stack[a][b];
      out << "]";
    }
    out << "]";

    out << ", dangle5: [[" << $self->dangle5[0][0];
    for (size_t a = 1; a < 5; a++)
      out << ", " << $self->dangle5[0][a];
    out << "]";
    for (size_t a = 1; a < NBPAIRS + 1; a++) {
      out << ", [" << $self->dangle5[a][0];
      for (size_t b = 1; b < 5; b++)
        out << ", " << $self->dangle5[a][b];
      out << "]";
    }
    out << "]";

    out << ", dangle3: [[" << $self->dangle3[0][0];
    for (size_t a = 1; a < 5; a++)
      out << ", " << $self->dangle3[0][a];
    out << "]";
    for (size_t a = 1; a < NBPAIRS + 1; a++) {
      out << ", [" << $self->dangle3[a][0];
      for (size_t b = 1; b < 5; b++)
        out << ", " << $self->dangle3[a][b];
      out << "]";
    }
    out << "]";

    out << ", ninio: [" << $self->ninio[0];
    for (size_t i = 1; i < 5; i++)
      out << ", " << $self->ninio[i];
    out << "]";

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

/* make a nice object oriented interface to vrna_exp_param_t */
%extend vrna_exp_param_t {
  vrna_exp_param_t(vrna_md_t *model_details = NULL)
  {
    vrna_exp_param_t *P = vrna_exp_params(model_details);
    return P;
  }

#ifdef SWIGPYTHON
  std::string
  __str__()
  {
    std::ostringstream out;
    out << "{ model_details: RNA.md()";
    out << ", id: " << $self->id;
    out << ", temperature: " << $self->temperature;
    out << ", kT: " << $self->kT;
    out << ", alpha: " << $self->alpha;
    out << ", pf_scale: " << $self->alpha;
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

%extend vrna_fold_compound_t {

  void
  params_reset(vrna_md_t *md = NULL)
  {
    vrna_params_reset($self, md);
  }

  void
  params_subst(vrna_param_t *par = NULL)
  {
    vrna_params_subst($self, par);
  }

  void
  exp_params_rescale(void)
  {
    vrna_exp_params_rescale($self, NULL);
  }

  void
  exp_params_rescale(double fe)
  {
    vrna_exp_params_rescale($self, &fe);
  }

  void
  exp_params_reset(vrna_md_t *md = NULL)
  {
    vrna_exp_params_reset($self, md);
  }

  void
  exp_params_subst(vrna_exp_param_t *par)
  {
    vrna_exp_params_subst($self, par);
  }
}



%ignore get_parameter_copy;
%ignore get_scaled_pf_parameters;
%ignore get_boltzmann_factors;
%ignore get_boltzmann_factor_copy;
%ignore get_scaled_alipf_parameters;
%ignore get_boltzmann_factors_ali;
%ignore scale_parameters;
%ignore get_scaled_parameters;
%ignore copy_parameters;
%ignore set_parameters;
%ignore scale_pf_parameters;
%ignore copy_pf_param;
%ignore set_pf_param;

%include <ViennaRNA/params/basic.h>



/**********************************************/
/* BEGIN interface for parameter file I/O     */
/**********************************************/

%rename(params_load)                            my_params_load;
%rename(params_save)                            my_params_save;
%rename(params_load_from_string)                my_params_load_from_string;
%rename(params_load_RNA_Turner2004)             vrna_params_load_RNA_Turner2004;
%rename(params_load_RNA_Turner1999)             vrna_params_load_RNA_Turner1999;
%rename(params_load_RNA_Andronescu2007)         vrna_params_load_RNA_Andronescu2007;
%rename(params_load_RNA_Langdon2018)            vrna_params_load_RNA_Langdon2018;
%rename(params_load_RNA_misc_special_hairpins)  vrna_params_load_RNA_misc_special_hairpins;
%rename(params_load_DNA_Mathews2004)            vrna_params_load_DNA_Mathews2004;
%rename(params_load_DNA_Mathews1999)            vrna_params_load_DNA_Mathews1999;

#ifdef SWIGPYTHON
%feature("autodoc")my_params_load;
%feature("kwargs")my_params_load;
%feature("autodoc")my_params_save;
%feature("kwargs")my_params_save;
%feature("autodoc")my_params_load_from_string;
%feature("kwargs")my_params_load_from_string;
#endif
%{
  int
  my_params_load(std::string  filename = "",
                 unsigned int options = VRNA_PARAMETER_FORMAT_DEFAULT)
  {
    if (!filename.compare(""))
      return vrna_params_load_defaults();

    return vrna_params_load(filename.c_str(), options);
  }

  int
  my_params_save(std::string filename,
                 unsigned int options = VRNA_PARAMETER_FORMAT_DEFAULT)
  {
    return vrna_params_save(filename.c_str(), options);
  }

  int
  my_params_load_from_string(std::string parameters,
                             std::string name = "",
                             unsigned int options = VRNA_PARAMETER_FORMAT_DEFAULT)
  {
    return vrna_params_load_from_string(parameters.c_str(),
                                        name.c_str(),
                                        options);
  }
%}

int
my_params_load(std::string  filename = "",
               unsigned int options = VRNA_PARAMETER_FORMAT_DEFAULT);

int
my_params_save(std::string filename,
               unsigned int options = VRNA_PARAMETER_FORMAT_DEFAULT);

int
my_params_load_from_string(std::string parameters,
                           std::string name = "",
                           unsigned int options = VRNA_PARAMETER_FORMAT_DEFAULT);


%constant unsigned int PARAMETER_FORMAT_DEFAULT = VRNA_PARAMETER_FORMAT_DEFAULT;


%include <ViennaRNA/params/io.h>

/**********************************************/
/* BEGIN interface for energy constants       */
/**********************************************/

%immutable;

%include  <ViennaRNA/params/constants.h>

%include  <ViennaRNA/params/default.h>

%mutable;
