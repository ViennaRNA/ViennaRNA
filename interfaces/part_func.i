/**********************************************/
/* BEGIN interface for Partition Function     */
/* computation                                */
/**********************************************/

%rename (pf_fold) my_pf_fold;
%{
  char *my_pf_fold(char *string, float *energy) {
    char *struc;
    struc = (char *)calloc(strlen(string)+1,sizeof(char));
    *energy = pf_fold(string, struc);
    return(struc);
  }
  char *my_pf_fold(char *string, char *constraints, float *energy) {
    char *struc;
    struc = (char *)calloc(strlen(string)+1,sizeof(char));
    if (constraints && fold_constrained)
      strncpy(struc, constraints, strlen(string));
    *energy = pf_fold(string, struc);
    if (constraints)
      strncpy(constraints, struc, strlen(constraints));
    return(struc);
  }
%}

%newobject my_pf_fold;
char *my_pf_fold(char *string, float *OUTPUT);
char *my_pf_fold(char *string, char *constraints, float *OUTPUT);
%ignore pf_fold;

%rename (pf_circ_fold) my_pf_circ_fold;
%{
  char *my_pf_circ_fold(char *string, float *energy) {
    char *struc;
    struc = (char *)calloc(strlen(string)+1,sizeof(char));
    *energy = pf_circ_fold(string, struc);
    return(struc);
  }
  char *my_pf_circ_fold(char *string, char *constraints, float *energy) {
    char *struc;
    struc = (char *)calloc(strlen(string)+1,sizeof(char));
    if (constraints && fold_constrained)
      strncpy(struc, constraints, strlen(string));
    *energy = pf_circ_fold(string, struc);
    if (constraints)
      strncpy(constraints, struc, strlen(constraints));
    return(struc);
  }
%}

%newobject my_pf_circ_fold;
char *my_pf_circ_fold(char *string, char *constraints, float *OUTPUT);
char *my_pf_circ_fold(char *string, float *OUTPUT);

%ignore pf_circ_fold;

%newobject pbacktrack;
extern char *pbacktrack(char *sequence);

/* make the float precision identifier available through the interface */
%rename (pf_float_precision) vrna_pf_float_precision;

/* these functions remain for now due to backward compatibility reasons
%ignore pf_circ_fold;
%ignore pbacktrack;
%ignore pbacktrack5;
%ignore pbacktrack_circ;
%ignore free_pf_arrays;
%ignore update_pf_params;
%ignore mean_bp_distance;
%ignore init_pf_fold;
%ignore centroid;
*/

%ignore pf_fold_par;
%ignore update_pf_params_par;
%ignore export_bppm;
%ignore get_pf_arrays;
%ignore get_subseq_F;
%ignore mean_bp_distance_pr;
%ignore stackProb;
%ignore get_centroid_struct_gquad_pr;
%ignore mean_bp_dist;
%ignore expHairpinEnergy;
%ignore expLoopEnergy;
%ignore assign_plist_gquad_from_pr;

/* tell swig that these functions return objects that require memory management */
%newobject vrna_fold_compound_t::pf;

%extend vrna_fold_compound_t{

  char *pf(float *OUTPUT){

    char *structure = (char *)vrna_alloc(sizeof(char) * ($self->length + 1)); /*output is a structure pointer*/
    *OUTPUT= vrna_pf($self, structure);
    return structure;
  }

  double mean_bp_distance(){
    return vrna_mean_bp_distance($self);
  }

  double
  ensemble_defect(std::string structure) {
    return vrna_ensemble_defect($self, structure.c_str());
  }
}

%include  <ViennaRNA/part_func.h>
%include  <ViennaRNA/equilibrium_probs.h>

/* attach stochastic backtracking functions as method of fold_compound */
%newobject vrna_fold_compound_t::pbacktrack;

%newobject vrna_fold_compound_t::pbacktrack_nr_resume;

%ignore vrna_nr_memory_t;
%ignore vrna_nr_memory_s;

%rename (nr_memory)  vrna_nr_memory_t;

typedef struct {} vrna_nr_memory_t;

%nodefaultctor vrna_nr_memory_t;
%nodefaultdtor vrna_nr_memory_t;

%extend vrna_nr_memory_t {
  vrna_nr_memory_t() {
    vrna_nr_memory_t *m = (vrna_nr_memory_t *)vrna_alloc(sizeof(vrna_nr_memory_t));
    *m = NULL;
    return m;
  }
  ~vrna_nr_memory_t() {
    vrna_pbacktrack_nr_free(*$self);
    delete $self;
  }
}

/*
  we need this crazy piece of argout typemap only because we don't
  want the vrna_nr_memory_t object to be appended to the results(list),
  but prepended instead. Otherwise, a simple
  %append_output(SWIG_NewPointerObj(SWIG_as_voidptr(retval$argnum), $1_descriptor, 0));
  would have sufficed already
*/
%typemap(argout) vrna_nr_memory_t * {
#ifdef SWIGPYTHON
  PyObject *o, *o2, *o3;
  o = SWIG_NewPointerObj(SWIG_as_voidptr(retval$argnum), $1_descriptor, 1);
  if ((!$result) || ($result == Py_None)) {
    $result = o;
  } else {
    PyObject *o2 = $result;
    $result = PyTuple_New(1);
    PyTuple_SetItem($result,0,o2);
    o3 = PyTuple_New(1);
    PyTuple_SetItem(o3,0,o);
    o2 = $result;
    $result = PySequence_Concat(o3,o2);
    Py_DECREF(o2);
    Py_DECREF(o3);
  }
#elif defined(SWIGPERL5)
  /* increase output stack if necessary */
  if (argvi >= items) {
    EXTEND(sp,1);
  }
  /* move already existing return values to the back */
  for (int i = argvi; i > 0; i--) {
    ST(i) = ST(i - 1);
  }
  /* store result as first element in the stack */
  ST(0) = SWIG_NewPointerObj(SWIG_as_voidptr(retval$argnum), $1_descriptor, 0);
  /* increase return argument counter */
  argvi++;
#endif
}


%typemap(in) vrna_nr_memory_t *INOUT (vrna_nr_memory_t *retval)
{
#ifdef SWIGPYTHON
  if ($input == Py_None) {
#elif defined(SWIGPERL5)
  if (!SvOK($input)) {
#endif
    retval = new vrna_nr_memory_t();
    $1 = retval;
  } else {
    /* INOUT in */
    SWIG_ConvertPtr($input,SWIG_as_voidptrptr(&retval), 0, SWIG_POINTER_DISOWN);
    $1 = retval;
  }
}

%extend vrna_fold_compound_t {

  char *
  pbacktrack(void)
  {
    return vrna_pbacktrack($self);
  }

  char *
  pbacktrack(int length)
  {
    return vrna_pbacktrack5($self, length);
  }

  std::vector<std::string>
  pbacktrack_nr(unsigned int num_samples)
  {
    std::vector<std::string> str_vec;

    char **ptr, **output = vrna_pbacktrack_nr($self, num_samples);

    if (output) {
      for (ptr = output; *ptr != NULL; ptr++) {
        str_vec.push_back(std::string(*ptr));
        free(*ptr);
      }

      free(output);
    }

    return str_vec;
  }

  %apply vrna_nr_memory_t *INOUT { vrna_nr_memory_t *nr_memory };

  std::vector<std::string>
  pbacktrack_nr_resume(unsigned int      num_samples,
                       vrna_nr_memory_t  *nr_memory)
  {
    std::vector<std::string> str_vec;

    char **ptr, **output = vrna_pbacktrack_nr_resume($self,
                                                     num_samples,
                                                     nr_memory);

    if (output) {
      for (ptr = output; *ptr != NULL; ptr++) {
        str_vec.push_back(std::string(*ptr));
        free(*ptr);
      }

      free(output);
    }

    return str_vec;
  }

  %clear vrna_nr_memory_t *nr_memory;
}

%include  <ViennaRNA/boltzmann_sampling.h>

/**********************************************/
/* BEGIN interface for cofold partition       */
/* function                                   */
/**********************************************/

%ignore cofoldF;

%rename (co_pf_fold) my_co_pf_fold;
%{
  char *my_co_pf_fold(char *string, float *FA, float *FB, float *FcAB, float *FAB) {
    char *struc;
    vrna_dimer_pf_t temp;
    struc = (char *)calloc(strlen(string)+1,sizeof(char));
    temp=co_pf_fold(string, struc);
    *FAB = temp.FAB;
    *FcAB = temp.FcAB;
    *FA = temp.FA;
    *FB = temp.FB;
    return(struc);
  }
  char *my_co_pf_fold(char *string, char *constraints, float *FA, float *FB, float *FcAB, float *FAB) {
    char *struc;
    vrna_dimer_pf_t temp;
    struc = (char *)calloc(strlen(string)+1,sizeof(char));
    if (constraints && fold_constrained)
      strncpy(struc, constraints, strlen(string));
    temp=co_pf_fold(string, struc);
    *FAB = temp.FAB;
    *FcAB = temp.FcAB;
    *FA = temp.FA;
    *FB = temp.FB;
    if (constraints)
      strncpy(constraints, struc, strlen(constraints));
    return(struc);
  }
%}

%newobject my_co_pf_fold;
char *my_co_pf_fold(char *string, float *OUTPUT, float *OUTPUT, float *OUTPUT, float *OUTPUT);
char *my_co_pf_fold(char *string, char *constraints, float *OUTPUT, float *OUTPUT, float *OUTPUT, float *OUTPUT);

%ignore co_pf_fold;
%ignore co_pf_fold_par;
%ignore compute_probabilities;
%ignore get_concentrations;
%ignore export_co_bppm;
%ignore update_co_pf_params_par;
%ignore co_bppm_symbol;
%ignore init_co_pf_fold;
%ignore get_plist;
%ignore pairpro;
%ignore ConcEnt;

%rename (get_concentrations) my_get_concentrations;
%{
 void my_get_concentrations(double FcAB, double FcAA, double FcBB, double FEA, double FEB, double Ac_start, double Bc_start, double *AB, double *AA, double *BB, double *A, double *B) {
    vrna_dimer_conc_t *temp;
    double *concis;
    concis = (double *)calloc(4,sizeof(double));
    concis[0]=Ac_start;
    concis[1]=Bc_start;
    concis[2]=0;
    temp=get_concentrations(FcAB,FcAA,FcBB,FEA,FEB,concis);
    *AB=temp->ABc;
    *AA=temp->AAc;
    *BB=temp->BBc;
    *A=temp->Ac;
    *B=temp->Bc;
    free(concis);
    free(temp);
    return;
  }
%}

%newobject my_get_concentrations;
void my_get_concentrations(double FcAB, double FcAA, double FcBB, double FEA,double FEB, double A0, double BO, double *OUTPUT, double *OUTPUT, double *OUTPUT, double *OUTPUT, double *OUTPUT);

/* tell swig that these functions return objects that require memory management */
%newobject vrna_fold_compound_t::pf_dimer;

%extend vrna_fold_compound_t {

  %apply float *OUTPUT { float *FA, float *FB, float *FcAB, float *FAB };

  char *pf_dimer(float *FA, float *FB, float *FcAB, float *FAB){

    char *structure = (char *)vrna_alloc(sizeof(char) * ($self->length + 1)); /*output is a structure pointer*/
    vrna_dimer_pf_t temp = vrna_pf_dimer($self, structure);
    *FAB  = (float)temp.FAB;
    *FcAB = (float)temp.FcAB;
    *FA   = (float)temp.FA;
    *FB   = (float)temp.FB;
    return structure;
  }

  %clear float *FA, float *FB, float *FcAB, float *FAB;

}

%include  <ViennaRNA/part_func_co.h>
%include  <ViennaRNA/concentrations.h>


%extend vrna_fold_compound_t{

  std::vector<std::vector<double> > bpp(void){
    std::vector<std::vector<double> > probabilities;
    vrna_fold_compound_t *vc = $self;
    if(vc->exp_matrices && vc->exp_matrices->probs){
      int turn, i, j, *idx, n;
      FLT_OR_DBL *probs;

      n     = vc->length;
      idx   = vc->iindx;
      turn  = vc->exp_params->model_details.min_loop_size;
      probs = vc->exp_matrices->probs;

      probabilities.push_back(std::vector<double>(n+1, 0.));
      for(i=1; i <= n; i++){
        int u = MIN2(i + turn + 1, n);
        probabilities.push_back(std::vector<double>(u, 0.));
        for(j = u; j <= n; j++)
          probabilities[i].push_back((double)probs[idx[i] - j]);
      }
    }
    return probabilities;
  }
}

%{
double get_pr(int i, int j) {
  int ii;
  if (i>j) {ii=i; i=j; j=ii;}
  return pr[iindx[i]-j];
}
%}
double get_pr(int i, int j);
/* Get probability of pair i.j from the pr array */
