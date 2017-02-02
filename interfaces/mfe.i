/**********************************************/
/* BEGIN interface for Minimum Free Energy    */
/* prediction                                 */
/**********************************************/

//%section "Folding Routines"
//%subsection "Minimum free Energy Folding"

%rename (fold) my_fold;

%{
  char *my_fold(char *string, float *energy) {
    char *struc;
    struc = (char *)calloc(strlen(string)+1,sizeof(char));
    *energy = fold(string, struc);
    return(struc);
  }

  char *my_fold(char *string, char *constraints, float *energy) {
    char *struc;
    struc = (char *)calloc(strlen(string)+1,sizeof(char));
    if (constraints && fold_constrained)
      strncpy(struc, constraints, strlen(string));
    *energy = fold(string, struc);
    if (constraints)
      strncpy(constraints, struc, strlen(constraints));
    return(struc);
  }
%}

%newobject my_fold;
char *my_fold(char *string, float *OUTPUT);
char *my_fold(char *string, char *constraints, float *OUTPUT);
%ignore fold;

/* these functions remain for now due to backward compatibility reasons
%ignore update_fold_params
%ignore free_arrays
%ignore circfold
*/

%rename (circfold) my_circfold;

%{
  char *my_circfold(char *string, float *energy) {
    char *struc;
    struc = (char *)calloc(strlen(string)+1,sizeof(char));
    *energy = circfold(string, struc);
    return(struc);
  }
  char *my_circfold(char *string, char *constraints, float *energy) {
    char *struc;
    struc = (char *)calloc(strlen(string)+1,sizeof(char));
    if (constraints && fold_constrained)
      strncpy(struc, constraints, strlen(string));
    *energy = circfold(string, struc);
    if (constraints)
      strncpy(constraints, struc, strlen(constraints));
    return(struc);
  }
%}

%newobject my_circfold;
char *my_circfold(char *string, float *OUTPUT);
char *my_circfold(char *string, char *constraints, float *OUTPUT);
%ignore circfold;




%rename (cofold) my_cofold;

%ignore fold_par;
%ignore update_fold_params_par;
%ignore initialize_fold;
%ignore HairpinE;
%ignore LoopEnergy;
%ignore export_circfold_arrays_par;
%ignore export_circfold_arrays;
%ignore export_fold_arrays;
%ignore export_fold_arrays_par;
%ignore backtrack_fold_from_pair;

%include  <ViennaRNA/fold.h>

/**********************************************/
/* BEGIN interface for advance MFE prediction */
/**********************************************/

%include  <ViennaRNA/mfe.h>

/**********************************************/
/* BEGIN interface for cofold                 */
/**********************************************/

%rename (cofold) my_cofold;

%{
  char *my_cofold(char *string, float *energy) {
    char *struc;
    struc = (char *)calloc(strlen(string)+1,sizeof(char));
    if (cut_point > (int)strlen(string)) {
       cut_point = -1;
    } 
    *energy = cofold(string, struc);
    return(struc);
  }

  char *my_cofold(char *string, char *constraints, float *energy) {
    char *struc;
    struc = (char *)calloc(strlen(string)+1,sizeof(char));
    if (constraints && fold_constrained)
      strncpy(struc, constraints, strlen(string));
    if (cut_point > (int)strlen(string)) {
       cut_point = -1;
    } 
    *energy = cofold(string, struc);
    if (constraints)
      strncpy(constraints, struc, strlen(constraints));
    return(struc);
  }
%}

%newobject my_cofold;
char *my_cofold(char *string, float *OUTPUT);
char *my_cofold(char *string, char *constraints, float *OUTPUT);
%ignore cofold;

/* these functions remain for now due to backward compatibility reasons
%ignore free_co_arrays;
%ignore update_cofold_params;
%ignore zukersubopt;
%ignore initialize_cofold;
*/

%ignore cofold_par;
%ignore update_cofold_params_par;
%ignore export_cofold_arrays_gq;
%ignore export_cofold_arrays;
%ignore zukersubopt_par;
%ignore get_monomere_mfes;

%include  <ViennaRNA/cofold.h>

/**********************************************/
/* BEGIN interface for alifold                */
/**********************************************/

%rename (alifold) my_alifold;

%{
#include <string>
#include <cstring>
#include <vector>

  char *my_alifold(std::vector<std::string> alignment, float *energy) {
    char *struc;
    /* convert std::vector<std::string> to vector<const char *> */
    std::vector<const char*>  vc;
    std::transform(alignment.begin(), alignment.end(), std::back_inserter(vc), convert_vecstring2veccharcp);
    vc.push_back(NULL); /* mark end of sequences */

    struc = (char *)calloc(strlen(vc[0])+1,sizeof(char));
    *energy = alifold((const char **)&vc[0], struc);
    return(struc);
  }
  char *my_alifold(std::vector<std::string> alignment, char *constraints, float *energy) {
    char *struc;
    /* convert std::vector<std::string> to vector<const char *> */
    std::vector<const char*>  vc;
    std::transform(alignment.begin(), alignment.end(), std::back_inserter(vc), convert_vecstring2veccharcp);
    vc.push_back(NULL); /* mark end of sequences */

    struc = (char *)calloc(strlen(vc[0])+1,sizeof(char));
    if (constraints && fold_constrained)
      strncpy(struc, constraints, strlen(vc[0]));
    *energy = alifold((const char **)&vc[0], struc);
    if (constraints)
      strncpy(constraints, struc, strlen(constraints));
    return(struc);
  }
%}

%newobject my_alifold;
char *my_alifold(std::vector<std::string> alignment, float *OUTPUT);
char *my_alifold(std::vector<std::string> alignment, char *constraints, float *OUTPUT);
%ignore alifold;

%ignore energy_of_alistruct;
%ignore energy_of_ali_gquad_structure;

%ignore alipf_fold_par;
%ignore alipf_fold;
%ignore alipf_circ_fold;
%ignore export_ali_bppm;
%ignore free_alipf_arrays;
%ignore alipbacktrack;
%ignore get_alipf_arrays;
%ignore update_alifold_params;


%extend vrna_fold_compound_t {

  char *mfe(float *OUTPUT){

    char *structure = (char *)vrna_alloc(sizeof(char) * ($self->length + 1));
    *OUTPUT = vrna_mfe($self, structure);
    return structure;
  }

  /* MFE for 2 RNA strands */
  char *mfe_dimer(float *OUTPUT){

    char *structure = (char*)vrna_alloc(sizeof(char) * ($self->length + 1));
    *OUTPUT = vrna_mfe_dimer($self, structure);
    return structure;
  }

  float mfe_window(FILE *file=NULL){

    return vrna_mfe_window($self,file);
  }

  /* ONLY possible if VRNA_WITH_SVM is set
  float mfe_window_zscore(double min_z,FILE *file=NULL){

    return vrna_mfe_window_zscore($self,min_z,file);
  }*/
  
}

/* tell swig that these functions return objects that require memory management */
%newobject vrna_fold_compound_t::mfe;
%newobject vrna_fold_compound_t::mfe_dimer;

%include <ViennaRNA/alifold.h>
