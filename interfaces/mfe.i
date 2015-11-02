/**********************************************/
/* BEGIN interface for Minimum Free Energy    */
/* prediction                                 */
/**********************************************/

//%section "Folding Routines"
//%subsection "Minimum free Energy Folding"

%rename (fold) my_fold;

%{
  char *my_fold(char *string, char *constraints, float *energy) {
    char *struc;
    struc = calloc(strlen(string)+1,sizeof(char));
    if (constraints && fold_constrained)
      strncpy(struc, constraints, strlen(string));
    *energy = fold(string, struc);
    if (constraints)
      strncpy(constraints, struc, strlen(constraints));
    return(struc);
  }
%}

%newobject my_fold;
char *my_fold(char *string, char *constraints = NULL, float *OUTPUT);
%ignore fold;

/* these functions remain for now due to backward compatibility reasons
%ignore update_fold_params
%ignore free_arrays
%ignore circfold
*/

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

%include  "../src/ViennaRNA/fold.h"

/**********************************************/
/* BEGIN interface for advance MFE prediction */
/**********************************************/

%include  "../src/ViennaRNA/mfe.h"

/**********************************************/
/* BEGIN interface for cofold                 */
/**********************************************/

%rename (cofold) my_cofold;

%{
  char *my_cofold(char *string, char *constraints, float *energy) {
    char *struc;
    struc = calloc(strlen(string)+1,sizeof(char));
    if (constraints && fold_constrained)
      strncpy(struc, constraints, strlen(string));
    if (cut_point > strlen(string)) {
       cut_point = -1;
    } 
    *energy = cofold(string, struc);
    if (constraints)
      strncpy(constraints, struc, strlen(constraints));
    return(struc);
  }
%}

%newobject my_cofold;
char *my_cofold(char *string, char *constraints = NULL, float *OUTPUT);
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

%include  "../src/ViennaRNA/cofold.h"

/**********************************************/
/* BEGIN interface for alifold                */
/**********************************************/

%ignore alifold;
%rename (alifold) my_alifold;

%{
  char *my_alifold(char **strings, char *constraints, float *energy) {
    char *struc;
    struc = calloc(strlen(strings[0])+1,sizeof(char));
    if (constraints && fold_constrained)
      strncpy(struc, constraints, strlen(strings[0]));
    *energy = alifold(strings, struc);
    if (constraints)
      strncpy(constraints, struc, strlen(constraints));
    return(struc);
  }
%}

%newobject my_alifold;
char *my_alifold(char **strings, char *constraints = NULL, float *OUTPUT);

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


%include "../src/ViennaRNA/alifold.h"
