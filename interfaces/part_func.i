/**********************************************/
/* BEGIN interface for Partition Function     */
/* computation                                */
/**********************************************/

%rename (pf_fold) my_pf_fold;
%{
  char *my_pf_fold(char *string, char *constraints, float *energy) {
    char *struc;
    struc = calloc(strlen(string)+1,sizeof(char));
    if (constraints && fold_constrained)
      strncpy(struc, constraints, strlen(string));
    *energy = pf_fold(string, struc);
    if (constraints)
      strncpy(constraints, struc, strlen(constraints));
    return(struc);
  }
%}

%newobject my_pf_fold;
char *my_pf_fold(char *string, char *constraints = NULL, float *OUTPUT);
%ignore pf_fold;

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

%include  <ViennaRNA/part_func.h>

/**********************************************/
/* BEGIN interface for cofold partition       */
/* function                                   */
/**********************************************/

%ignore cofoldF;

%rename (co_pf_fold) my_co_pf_fold;
%{
  char *my_co_pf_fold(char *string, char *constraints, float *FA, float *FB, float *FcAB, float *FAB) {
    char *struc;
    float en;
    vrna_dimer_pf_t temp;
    struc = calloc(strlen(string)+1,sizeof(char));
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
char *my_co_pf_fold(char *string, char *constraints = NULL, float *OUTPUT, float *OUTPUT, float *OUTPUT, float *OUTPUT);

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
 void my_get_concentrations(double FcAB, double FcAA, double FcBB, double FEA, double FEB, double A0, double B0, double *AB, double *AA, double *BB, double *A, double *B) {
    vrna_dimer_conc_t *temp;
    double *concis;
    concis = (double *)calloc(4,sizeof(double));
    concis[0]=A0;
    concis[1]=B0;
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

%include  <ViennaRNA/part_func_co.h>




%{
double get_pr(int i, int j) {
  int ii;
  if (i>j) {ii=i; i=j; j=ii;}
  return pr[iindx[i]-j];
}
%}
double get_pr(int i, int j);
/* Get probability of pair i.j from the pr array */

