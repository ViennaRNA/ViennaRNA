/**********************************************/
/* BEGIN interface for energy parameters      */
/**********************************************/

/* do not create default constructor and hide data fields of vrna_param_t from SWIG */
%ignore paramT;
%ignore pf_paramT;

/* scripting language access through 'param' instead of 'vrna_param_t' */
%rename (param) vrna_param_t;

/* scripting language access through 'exp_param' instead of 'vrna_exp_param_t' */
%rename (exp_param) vrna_exp_param_t;

%nodefaultctor vrna_param_t;
typedef struct {
  int       hairpin[31];
  int       bulge[MAXLOOP+1];
  int       internal_loop[MAXLOOP+1];
  int       ninio[5];
  double    lxc;
  int       MLbase;
  int       MLintern[NBPAIRS+1];
  int       MLclosing;
  int       TerminalAU;
  int       DuplexInit;
  int       Tetraloop_E[200];
  char      Tetraloops[1401];
  int       Triloop_E[40];
  char      Triloops[241];
  int       Hexaloop_E[40];
  char      Hexaloops[1801];
  int       TripleC;
  int       MultipleCA;
  int       MultipleCB;
  double    temperature;
  vrna_md_t model_details;
} vrna_param_t;

/* do not create default constructor and hide data fields of vrna_param_t from SWIG */
%nodefaultctor vrna_exp_param_t;
typedef struct {
  double  exphairpin[31];
  double  expbulge[MAXLOOP+1];
  double  expinternal[MAXLOOP+1];
  double  lxc;
  double  expMLbase;
  double  expMLintern[NBPAIRS+1];
  double  expMLclosing;
  double  expTermAU;
  double  expDuplexInit;
  double  exptetra[40];
  double  exptri[40];
  double  exphex[40];
  char    Tetraloops[1401];
  double  expTriloop[40];
  char    Triloops[241];
  char    Hexaloops[1801];
  double  expTripleC;
  double  expMultipleCA;
  double  expMultipleCB;
  double  kT;
  double  pf_scale;
  double  temperature;
  double  alpha;
  vrna_md_t model_details;
} vrna_exp_param_t;

/* make a nice object oriented interface to vrna_param_t */
%extend vrna_param_t {
  vrna_param_t(){
    vrna_param_t *P = vrna_params(NULL);
    return P;
  }
  vrna_param_t(vrna_md_t *md){
    vrna_param_t *P = vrna_params(md);
    return P;
  }

  double  get_temperature(){
    return $self->temperature;
  }
}

/* make a nice object oriented interface to vrna_exp_param_t */
%extend vrna_exp_param_t {
  vrna_exp_param_t(){
    vrna_exp_param_t *P = vrna_exp_params(NULL);
    return P;
  }
  vrna_exp_param_t(vrna_md_t *md){
    vrna_exp_param_t *P = vrna_exp_params(md);
    return P;
  }

  double  get_temperature(){
    return $self->temperature;
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

%include <ViennaRNA/params.h>



/**********************************************/
/* BEGIN interface for parameter file I/O     */
/**********************************************/

%include <ViennaRNA/read_epars.h>

/**********************************************/
/* BEGIN interface for energy constants       */
/**********************************************/

%include  <ViennaRNA/energy_const.h>
