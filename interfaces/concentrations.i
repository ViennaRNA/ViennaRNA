/**********************************************/
/* BEGIN interface for multimer concentration */
/* stuff                                      */
/**********************************************/



%ignore get_concentrations;
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


%ignore vrna_equilibrium_constants;

%rename (equilibrium_constants) my_equilibrium_constants;

%{

  std::vector<double>
  my_equilibrium_constants(std::vector<double>                      dG_complexes,
                           std::vector<double>                      dG_strands,
                           std::vector<std::vector<unsigned int> >  A,
                           double                                   kT,
                           size_t                                   strands = 0,
                           size_t                                   complexes = 0)
  {
    double              *constants = NULL;
    unsigned int        **assoc_mx;
    std::vector<double> constants_vector;

    if (strands == 0)
      strands = (size_t)dG_strands.size();

    if (complexes == 0)
      complexes = (size_t)dG_complexes.size();

    assoc_mx = (unsigned int **)vrna_alloc(sizeof(unsigned int *) * A.size());
    for (size_t k = 0; k < A.size(); k++) {
      assoc_mx[k] = (unsigned int *)vrna_alloc(sizeof(unsigned int) * A[k].size());
      memcpy(assoc_mx[k], A[k].data(), sizeof(unsigned int) * A[k].size());
    }

    constants = vrna_equilibrium_constants(&(dG_complexes[0]),
                                           &(dG_strands[0]),
                                           (const unsigned int **)assoc_mx,
                                           kT,
                                           strands,
                                           complexes);

    /* reformat output array to actual vector */
    for (size_t k = 0; k < complexes; k++)
      constants_vector.push_back(constants[k]);

    for (size_t k = 0; k < A[k].size(); k++)
      free(assoc_mx[k]);

    free(assoc_mx);

    free(constants);
    
    return constants_vector;
  }
%}


std::vector<double>
my_equilibrium_constants(std::vector<double>                      dG_complexes,
                         std::vector<double>                      dG_strands,
                         std::vector<std::vector<unsigned int> >  A,
                         double                                   kT,
                         size_t                                   strands = 0,
                         size_t                                   complexes = 0);


%include  <ViennaRNA/concentrations.h>


%include  <ViennaRNA/wrap_dlib.h>
