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

%extend vrna_fold_compound_t{

  char *pf(float *OUTPUT){

    char *structure = (char *)vrna_alloc(sizeof(char) * ($self->length + 1)); /*output is a structure pointer*/
    *OUTPUT= vrna_pf($self, structure);
    return structure;
  }

  double mean_bp_distance(){
    return vrna_mean_bp_distance($self);
  }
}

/* tell swig that these functions return objects that require memory management */
%newobject vrna_fold_compound_t::pf;

%include  <ViennaRNA/part_func.h>
%include  <ViennaRNA/equilibrium_probs.h>
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

/* tell swig that these functions return objects that require memory management */
%newobject vrna_fold_compound_t::pf_dimer;
char *vrna_fold_compound_t::pf_dimer(float *OUTPUT, float *OUTPUT, float *OUTPUT, float *OUTPUT);

%include  <ViennaRNA/part_func_co.h>


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
