/**********************************************/
/* BEGIN interface for Minimum Free Energy    */
/* prediction                                 */
/**********************************************/


/**********************************************/
/* BEGIN interface for MFE prediction         */
/**********************************************/


%rename (fold)      my_fold;
%rename (alifold) my_alifold;
%rename (cofold) my_cofold;
%rename (circfold)  my_circfold;
%rename (circalifold)  my_circalifold;


%apply  float *OUTPUT { float *energy };

%{
#include <string>
#include <cstring>
#include <vector>

  char *
  my_fold(char *string,
          float *energy)
  {
    char *struc;

    struc = (char *)calloc(strlen(string)+1,sizeof(char));
    *energy = vrna_fold(string, struc);

    return struc;
  }

  char *
  my_fold(char *string,
          char *constraints,
          float *energy)
  {
    char                  *struc;
    vrna_fold_compound_t  *fc;

    struc = (char *)calloc(strlen(string)+1,sizeof(char));
    fc    = vrna_fold_compound(string, NULL, VRNA_OPTION_DEFAULT);

    if (constraints && fold_constrained)
      vrna_hc_add_from_db(fc, constraints, VRNA_CONSTRAINT_DB_DEFAULT);

    *energy = vrna_mfe(fc, struc);

    vrna_fold_compound_free(fc);

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
    if (constraints && (!fold_constrained))
      strncpy(constraints, struc, strlen(constraints));
#endif

    return struc;
  }

  char *
  my_alifold(std::vector<std::string> alignment,
             float                    *energy)
  {
    char *struc;
    /* convert std::vector<std::string> to vector<const char *> */
    std::vector<const char*>  vc;
    std::transform(alignment.begin(), alignment.end(), std::back_inserter(vc), convert_vecstring2veccharcp);
    vc.push_back(NULL); /* mark end of sequences */

    struc = (char *)calloc(strlen(vc[0])+1,sizeof(char));

    *energy = vrna_alifold((const char **)&vc[0], struc);

    return struc;
  }

  char *
  my_alifold(std::vector<std::string> alignment,
             char                     *constraints,
             float                    *energy)
  {
    char                      *struc;
    vrna_fold_compound_t      *fc;
    std::vector<const char*>  vc;

    /* convert std::vector<std::string> to vector<const char *> */
    std::transform(alignment.begin(),
                   alignment.end(),
                   std::back_inserter(vc),
                   convert_vecstring2veccharcp);
    vc.push_back(NULL); /* mark end of sequences */

    struc = (char *)calloc(strlen(vc[0])+1,sizeof(char));

    fc = vrna_fold_compound_comparative((const char **)&vc[0], NULL, VRNA_OPTION_DEFAULT);

    if (constraints && fold_constrained)
      vrna_hc_add_from_db(fc, constraints, VRNA_CONSTRAINT_DB_DEFAULT);

    *energy = vrna_mfe(fc, struc);

    vrna_fold_compound_free(fc);

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
    if (constraints && (!fold_constrained))
      strncpy(constraints, struc, strlen(constraints));
#endif

    return struc;
  }

  char *
  my_cofold(char  *string,
            float *energy)
  {
    char *s, **tok, **ptr, *struc, *sequence;

    sequence = string;
    struc    = (char *)calloc(strlen(string)+1,sizeof(char));

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
    /* first, tokenize the input at delimiter '&' */
    tok = vrna_strsplit(string, "&");

    /*
        now, check whether there is only a single sequence.
        This may be a hint that someone is still using the
        'old' API where the split point had to be spliced out
        and explicitly specified through the global variable
        cut_point
     */
    if ((tok) && (tok[0])) {
      if (!tok[1]) {
        if (cut_point > (int)strlen(string)) {
          cut_point = -1;
        } else {
          /* we need to re-insert the delimiter now */
          sequence = vrna_cut_point_insert(string, cut_point);
        }
      }
    }
#endif

    *energy = vrna_cofold(sequence, struc);

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
    /* clean up */
    if (tok) {
      for (ptr = tok; *ptr; ptr++)
        free(*ptr);

      free(tok);
    }

    if (sequence != string)
      free(sequence);
#endif

    return struc;
  }

  char *
  my_cofold(char  *string,
            char  *constraints,
            float *energy)
  {
    char *s, **tok, **ptr, *struc, *sequence;
    vrna_fold_compound_t      *fc;

    sequence = string;
    struc    = (char *)calloc(strlen(string)+1,sizeof(char));

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
    /* first, tokenize the input at delimiter '&' */
    tok = vrna_strsplit(string, "&");

    /*
        now, check whether there is only a single sequence.
        This may be a hint that someone is still using the
        'old' API where the split point had to be spliced out
        and explicitly specified through the global variable
        cut_point
     */
    if ((tok) && (tok[0])) {
      if (!tok[1]) {
        if (cut_point > (int)strlen(string)) {
          cut_point = -1;
        } else {
          /* we need to re-insert the delimiter now */
          sequence = vrna_cut_point_insert(string, cut_point);
        }
      }
    }
#endif

    fc = vrna_fold_compound(sequence, NULL, VRNA_OPTION_DEFAULT);

    if (constraints && fold_constrained)
      vrna_hc_add_from_db(fc, constraints, VRNA_CONSTRAINT_DB_DEFAULT);

    *energy = vrna_mfe_dimer(fc, struc);

    /* clean up */
#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
    if (tok) {
      for (ptr = tok; *ptr; ptr++)
        free(*ptr);

      free(tok);
    }

    if (sequence != string)
      free(sequence);
#endif

    vrna_fold_compound_free(fc);

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
    if (constraints && (!fold_constrained))
      strncpy(constraints, struc, strlen(constraints));
#endif

    return struc;
  }

  char *
  my_circfold(char *string,
              float *energy)
  {
    char *struc;

    struc   = (char *)calloc(strlen(string)+1,sizeof(char));
    *energy = vrna_circfold(string, struc);

    return struc;
  }

  char *
  my_circfold(char *string,
              char *constraints,
              float *energy)
  {
    char                  *struc;
    vrna_md_t             md;
    vrna_fold_compound_t  *fc;

    vrna_md_set_default(&md);
    md.circ = 1;

    struc = (char *)calloc(strlen(string)+1,sizeof(char));

    fc  = vrna_fold_compound(string, &md, VRNA_OPTION_DEFAULT);

    if (constraints && fold_constrained)
      vrna_hc_add_from_db(fc, constraints, VRNA_CONSTRAINT_DB_DEFAULT);

    *energy = vrna_mfe(fc, struc);

    vrna_fold_compound_free(fc);

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
    if (constraints && (!fold_constrained))
      strncpy(constraints, struc, strlen(constraints));
#endif

    return struc;
  }

  char *
  my_circalifold(std::vector<std::string> alignment,
                 float                    *energy)
  {
    char *struc;
    /* convert std::vector<std::string> to vector<const char *> */
    std::vector<const char*>  vc;
    std::transform(alignment.begin(), alignment.end(), std::back_inserter(vc), convert_vecstring2veccharcp);
    vc.push_back(NULL); /* mark end of sequences */

    struc = (char *)calloc(strlen(vc[0])+1,sizeof(char));

    *energy = vrna_circalifold((const char **)&vc[0], struc);

    return struc;
  }

  char *
  my_circalifold(std::vector<std::string> alignment,
                 char                     *constraints,
                 float                    *energy)
  {
    char                      *struc;
    vrna_fold_compound_t      *fc;
    std::vector<const char*>  vc;
    vrna_md_t                 md;

    vrna_md_set_default(&md);
    md.circ = 1;

    /* convert std::vector<std::string> to vector<const char *> */
    std::transform(alignment.begin(),
                   alignment.end(),
                   std::back_inserter(vc),
                   convert_vecstring2veccharcp);
    vc.push_back(NULL); /* mark end of sequences */

    struc = (char *)calloc(strlen(vc[0])+1,sizeof(char));

    fc = vrna_fold_compound_comparative((const char **)&vc[0],
                                        &md,
                                        VRNA_OPTION_DEFAULT);

    if (constraints && fold_constrained)
      vrna_hc_add_from_db(fc, constraints, VRNA_CONSTRAINT_DB_DEFAULT);

    *energy = vrna_mfe(fc, struc);

    vrna_fold_compound_free(fc);

    return struc;
  }

%}

%newobject my_fold;
%newobject my_alifold;
%newobject my_cofold;
%newobject my_circfold;
%newobject my_circalifold;

char *my_fold(char *string, float *energy);
char *my_fold(char *string, char *constraints, float *energy);
char *my_alifold(std::vector<std::string> alignment, float *energy);
char *my_alifold(std::vector<std::string> alignment, char *constraints, float *energy);
char *my_cofold(char *string, float *energy);
char *my_cofold(char *string, char *constraints, float *energy);
char *my_circfold(char *string, float *energy);
char *my_circfold(char *string, char *constraints, float *energy);
char *my_circalifold(std::vector<std::string> alignment, float *energy);
char *my_circalifold(std::vector<std::string> alignment, char *constraints, float *energy);


/* tell swig that these functions return objects that require memory management */
%newobject vrna_fold_compound_t::mfe;
%newobject vrna_fold_compound_t::mfe_dimer;
%newobject vrna_fold_compound_t::backtrack;

%extend vrna_fold_compound_t {

#ifdef SWIGPYTHON
%feature("autodoc") mfe;
%feature("kwargs") mfe;
%feature("autodoc") mfe_dimer;
%feature("kwargs") mfe_dimer;
%feature("autodoc") backtrack;
#endif

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

  char *backtrack(unsigned int length,
                  float *OUTPUT) {
    char *structure = (char *)vrna_alloc(sizeof(char) * (length + 1));
    *OUTPUT = vrna_backtrack5($self, length, structure);
    return structure;
  }

  char *backtrack(float *OUTPUT) {
    char *structure = (char *)vrna_alloc(sizeof(char) * ($self->length + 1));
    *OUTPUT = vrna_backtrack5($self, $self->length, structure);
    return structure;
  }
}


%clear  float *energy;

%include  <ViennaRNA/mfe.h>


/**********************************************/
/* BEGIN backward compatibility               */
/**********************************************/

/* these functions remain for now due to backward compatibility reasons
%ignore update_fold_params
%ignore free_arrays
*/

%ignore fold;
%ignore circfold;
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


%ignore energy_of_alistruct;
%ignore energy_of_ali_gquad_structure;

%ignore alifold;
%ignore alipf_fold_par;
%ignore alipf_fold;
%ignore alipf_circ_fold;
%ignore export_ali_bppm;
%ignore free_alipf_arrays;
%ignore alipbacktrack;
%ignore get_alipf_arrays;
%ignore update_alifold_params;


%include <ViennaRNA/alifold.h>


/* these functions remain for now due to backward compatibility reasons
%ignore free_co_arrays;
%ignore update_cofold_params;
%ignore zukersubopt;
%ignore initialize_cofold;
*/

%ignore cofold;
%ignore cofold_par;
%ignore update_cofold_params_par;
%ignore export_cofold_arrays_gq;
%ignore export_cofold_arrays;
%ignore zukersubopt_par;
%ignore get_monomere_mfes;

%include  <ViennaRNA/cofold.h>
