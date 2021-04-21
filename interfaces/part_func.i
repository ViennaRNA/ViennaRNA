/**********************************************/
/* BEGIN interface for Partition Function     */
/* computation                                */
/**********************************************/

%rename (pf_fold) my_pf_fold;
%{
  char *
  my_pf_fold(char   *string,
             float  *energy)
  {
    char *struc;
    struc = (char *)calloc(strlen(string)+1,sizeof(char));
    *energy = pf_fold(string, struc);
    return(struc);
  }

  char *
  my_pf_fold(char *string,
             char *constraints,
             float *energy)
  {
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
  char *
  my_pf_circ_fold(char  *string,
                  float *energy)
  {
    char *struc;
    struc = (char *)calloc(strlen(string)+1,sizeof(char));
    *energy = pf_circ_fold(string, struc);
    return(struc);
  }

  char *
  my_pf_circ_fold(char  *string,
                  char  *constraints,
                  float *energy)
  {
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

/* make the float precision identifier available through the interface */
%rename (pf_float_precision) vrna_pf_float_precision;

/* these functions remain for now due to backward compatibility reasons
%ignore pf_circ_fold;
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

  char *
  pf(float *OUTPUT)
  {
    char *structure = (char *)vrna_alloc(sizeof(char) * ($self->length + 1)); /*output is a structure pointer*/
    *OUTPUT= vrna_pf($self, structure);
    return structure;
  }

  double
  mean_bp_distance()
  {
    return vrna_mean_bp_distance($self);
  }

  double
  ensemble_defect(std::string structure,
                  unsigned int options = VRNA_BRACKETS_RND)
  {
    double ed;
    short int         *pt;

    pt = vrna_ptable_from_string(structure.c_str(), options);

    ed = vrna_ensemble_defect_pt($self, pt);

    free(pt);

    return ed;
  }

  double
  ensemble_defect(std::vector<int> pair_table)
  {
    double ed;

    std::vector<short> pt_v_short;
    transform(pair_table.begin(), pair_table.end(), back_inserter(pt_v_short), convert_vecint2vecshort);
    return vrna_ensemble_defect_pt($self, (short*)&pt_v_short[0]);
  }

  std::vector<double>
  positional_entropy(void)
  {
    unsigned int        n;
    double              *pos_ent;
    std::vector<double> dv;

    n       = $self->length;
    pos_ent = vrna_positional_entropy($self);

    if (pos_ent)
      dv.assign(pos_ent, pos_ent + (n + 1));

    return dv;
  }

  double
  pr_structure(std::string structure)
  {
    return vrna_pr_structure($self, structure.c_str());
  }

  double
  pr_energy(double e)
  {
    return vrna_pr_energy($self, e);
  }
}

%include  <ViennaRNA/part_func.h>
%include  <ViennaRNA/equilibrium_probs.h>

/**********************************************/
/* BEGIN interface for cofold partition       */
/* function                                   */
/**********************************************/

%ignore cofoldF;

%rename (co_pf_fold) my_co_pf_fold;
%{
  char *
  my_co_pf_fold(char  *string,
                float *FA,
                float *FB,
                float *FcAB,
                float *FAB)
  {
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
  char *
  my_co_pf_fold(char  *string,
                char  *constraints,
                float *FA,
                float *FB,
                float *FcAB,
                float *FAB)
  {
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
 void
 my_get_concentrations(double FcAB,
                       double FcAA,
                       double FcBB,
                       double FEA,
                       double FEB,
                       double Ac_start,
                       double Bc_start,
                       double *AB,
                       double *AA,
                       double *BB,
                       double *A,
                       double *B)
  {
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

  char *
  pf_dimer(float *FA,
           float *FB,
           float *FcAB,
           float *FAB)
  {
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

  std::vector<std::vector<double> >
  bpp(void)
  {
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
double
get_pr(int i,
       int j)
{
  int ii;
  if (i>j) {ii=i; i=j; j=ii;}
  return pr[iindx[i]-j];
}
%}
double get_pr(int i, int j);
/* Get probability of pair i.j from the pr array */


/**********************************************/
/* BEGIN interface for centroid structure     */
/**********************************************/

%newobject vrna_fold_compound_t::centroid;

#ifdef SWIGPYTHON
%feature("autodoc") vrna_fold_compound_t::centroid;
%feature("kwargs") vrna_fold_compound_t::centroid;
#endif

%extend vrna_fold_compound_t{

  char *
  centroid(double *OUTPUT)
  {
    return vrna_centroid($self, OUTPUT);
  }

}

%include  <ViennaRNA/centroid.h>

/**********************************************/
/* BEGIN interface for MEA structure          */
/**********************************************/

%newobject vrna_fold_compound_t::MEA;

#ifdef SWIGPYTHON
%feature("autodoc") vrna_fold_compound_t::MEA;
%feature("kwargs") vrna_fold_compound_t::MEA;
#endif

%extend vrna_fold_compound_t{

  char *
  MEA(float *OUTPUT)
  {
    return vrna_MEA($self, 1., OUTPUT);
  }

  char *
  MEA(double  gamma,
      float   *OUTPUT)
  {
    return vrna_MEA($self, gamma, OUTPUT);
  }
}

%rename (MEA_from_plist)  my_MEA_from_plist;

%{
  char *
  my_MEA_from_plist(std::vector<vrna_ep_t> plist,
                    std::string            sequence,
                    double                 gamma,
                    vrna_md_t              *md,
                    float                  *OUTPUT)
  {
    vrna_ep_t               pp;
    std::vector<vrna_ep_t>  pl = plist;

    pp.i = pp.j = 0;
    pp.p = 0.;
    pp.type = VRNA_PLIST_TYPE_BASEPAIR;
    pl.push_back(pp);

    return vrna_MEA_from_plist(&pl[0],
                               sequence.c_str(),
                               gamma,
                               md,
                               OUTPUT);
  }

  char *
  my_MEA_from_plist(std::vector<vrna_ep_t> plist,
                    std::string            sequence,
                    vrna_md_t              *md,
                    float                  *OUTPUT)
  {
    vrna_ep_t               pp;
    std::vector<vrna_ep_t>  pl = plist;

    pp.i = pp.j = 0;
    pp.p = 0.;
    pp.type = VRNA_PLIST_TYPE_BASEPAIR;
    pl.push_back(pp);

    return vrna_MEA_from_plist(&pl[0],
                               sequence.c_str(),
                               1.,
                               md,
                               OUTPUT);
  }

  char *
  my_MEA_from_plist(std::vector<vrna_ep_t> plist,
                    std::string            sequence,
                    double                 gamma,
                    float                  *OUTPUT)
  {
    vrna_ep_t               pp;
    std::vector<vrna_ep_t>  pl = plist;

    pp.i = pp.j = 0;
    pp.p = 0.;
    pp.type = VRNA_PLIST_TYPE_BASEPAIR;
    pl.push_back(pp);

    return vrna_MEA_from_plist(&pl[0],
                               sequence.c_str(),
                               gamma,
                               NULL,
                               OUTPUT);
  }

  char *
  my_MEA_from_plist(std::vector<vrna_ep_t> plist,
                    std::string            sequence,
                    float                  *OUTPUT)
  {
    vrna_ep_t               pp;
    std::vector<vrna_ep_t>  pl = plist;

    pp.i = pp.j = 0;
    pp.p = 0.;
    pp.type = VRNA_PLIST_TYPE_BASEPAIR;
    pl.push_back(pp);

    return vrna_MEA_from_plist(&pl[0],
                               sequence.c_str(),
                               1.,
                               NULL,
                               OUTPUT);
  }

%}

%newobject MEA_from_plist;

char *my_MEA_from_plist(std::vector<vrna_ep_t> plist,
                        std::string            sequence,
                        double                 gamma,
                        vrna_md_t              *md,
                        float                  *OUTPUT);
char *my_MEA_from_plist(std::vector<vrna_ep_t> plist,
                        std::string            sequence,
                        vrna_md_t              *md,
                        float                  *OUTPUT);
char *my_MEA_from_plist(std::vector<vrna_ep_t> plist,
                        std::string            sequence,
                        double                 gamma,
                        float                  *OUTPUT);
char *my_MEA_from_plist(std::vector<vrna_ep_t> plist,
                        std::string            sequence,
                        float                  *OUTPUT);

%ignore MEA;
%ignore MEA_seq;

%include  <ViennaRNA/MEA.h>
