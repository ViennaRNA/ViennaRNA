/**********************************************/
/* BEGIN interface for sliding-window MFE     */
/* prediction                                 */
/**********************************************/

%extend vrna_fold_compound_t {

#ifdef SWIGPYTHON
%feature("autodoc") mfe_window;
%feature("kwargs") mfe_window;
#endif

  float mfe_window(FILE *nullfile = NULL){

    return vrna_mfe_window($self, nullfile);
  }

#ifdef VRNA_WITH_SVM
  float mfe_window_zscore(double  min_z,
                          FILE    *nullfile = NULL) {

    return vrna_mfe_window_zscore($self, min_z, nullfile);
  }
#endif
}

/* tell swig that these functions return objects that require memory management */
%ignore Lfoldz;
%ignore Lfold;
%ignore aliLfold;
%ignore aliLfold_cb;
%ignore vrna_Lfoldz;
%ignore vrna_Lfold;
%ignore vrna_aliLfold;
%ignore vrna_aliLfold_cb;

#ifdef VRNA_WITH_SVM
%rename (Lfoldz)    my_Lfoldz;
#endif
%rename (Lfold)     my_Lfold;
%rename (aliLfold)  my_aliLfold;

%{

#ifdef VRNA_WITH_SVM
  float
  my_Lfoldz(std::string sequence,
            int         window_size,
            double      min_z,
            FILE        *nullfile = NULL) {

    return vrna_Lfoldz(sequence.c_str(),
                       window_size,
                       min_z,
                       nullfile);
  }
#endif

  float
  my_Lfold(std::string sequence,
           int        window_size,
           FILE       *nullfile = NULL) {

    return vrna_Lfold(sequence.c_str(), window_size, nullfile);
  }

  float
  my_aliLfold(std::vector<std::string> alignment,
              int                      window_size,
              FILE                     *nullfile = NULL) {

    /* convert std::vector<std::string> to vector<const char *> */
    std::vector<const char*>  aln;
    std::transform(alignment.begin(), alignment.end(), std::back_inserter(aln), convert_vecstring2veccharcp);
    aln.push_back(NULL); /* mark end of sequences */

    return vrna_aliLfold((const char **)&aln[0],
                         window_size,
                         nullfile);
  }

%}

#ifdef VRNA_WITH_SVM
float
my_Lfoldz(std::string sequence,
          int         window_size,
          double      min_z,
          FILE        *nullfile = NULL);
#endif

float
my_Lfold(std::string sequence,
         int         window_size,
         FILE        *nullfile = NULL);

float
my_aliLfold(std::vector<std::string> alignment,
            int                      window_size,
            FILE                     *nullfile = NULL);


%include <ViennaRNA/mfe_window.h>
