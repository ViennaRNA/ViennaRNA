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

%rename (pf_add) my_pf_add;
%{
  double
  my_pf_add(double dG1,
            double dG2,
            double kT = 0) {
    double result = 0;

    if (kT == 0)
      kT = (37. + K0) * GASCONST / 1000.;

    return (double)vrna_pf_add((FLT_OR_DBL)(dG1), (FLT_OR_DBL)(dG2), kT);
  }
%}

double my_pf_add(double dG1, double dG2, double kT = 0);

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
%ignore get_pf_arrays;
%ignore get_subseq_F;
%ignore get_centroid_struct_gquad_pr;
%ignore expHairpinEnergy;
%ignore expLoopEnergy;

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
}

%include  <ViennaRNA/partfunc/global.h>

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
%ignore export_co_bppm;
%ignore update_co_pf_params_par;
%ignore co_bppm_symbol;
%ignore init_co_pf_fold;
%ignore get_plist;
%ignore pairpro;

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

%include  <ViennaRNA/partfunc/exterior.h>
%include  <ViennaRNA/partfunc/internal.h>
%include  <ViennaRNA/partfunc/multibranch.h>
%include  <ViennaRNA/partfunc/gquad.h>


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

%include  <ViennaRNA/structures/centroid.h>

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

%include  <ViennaRNA/structures/mea.h>
