/**********************************************/
/* BEGIN interface for sliding-window MFE     */
/* prediction                                 */
/**********************************************/

%extend vrna_fold_compound_t {

#ifdef SWIGPYTHON
%feature("autodoc") mfe_window;
%feature("kwargs") mfe_window;
#endif

  float
  mfe_window(FILE *nullfile = NULL)
  {
    return vrna_mfe_window($self, nullfile);
  }

#ifdef VRNA_WITH_SVM
  float
  mfe_window_zscore(double  min_z,
                    FILE    *nullfile = NULL)
  {
    return vrna_mfe_window_zscore($self, min_z, nullfile);
  }

  int
  zsc_filter_init(double        min_z   = -2.0,
                  unsigned int  options = VRNA_ZSCORE_SETTINGS_DEFAULT)
  {
    return vrna_zsc_filter_init($self, min_z, options);
  }

  int
  zsc_filter_update(double        min_z,
                    unsigned int  options = VRNA_ZSCORE_OPTIONS_NONE)
  {
    return vrna_zsc_filter_update($self, min_z, options);
  }

  void
  zsc_filter_free(void)
  {
    vrna_zsc_filter_free($self);
  }

  int
  zsc_filter_on(void)
  {
    return vrna_zsc_filter_on($self);
  }

  double
  zsc_filter_threshold(void)
  {
    return vrna_zsc_filter_threshold($self);
  }

  double
  zsc_compute(unsigned int  i,
              unsigned int  j,
              int           e)
  {
    return vrna_zsc_compute($self, i, j, e);
  }

#endif
}

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
            FILE        *nullfile = NULL)
  {
    return vrna_Lfoldz(sequence.c_str(),
                       window_size,
                       min_z,
                       nullfile);
  }
#endif

  float
  my_Lfold(std::string sequence,
           int        window_size,
           FILE       *nullfile = NULL)
  {
    return vrna_Lfold(sequence.c_str(), window_size, nullfile);
  }

  float
  my_aliLfold(std::vector<std::string> alignment,
              int                      window_size,
              FILE                     *nullfile = NULL)
  {
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

#ifdef VRNA_WITH_SVM
%constant unsigned int ZSCORE_OPTIONS_NONE      =  VRNA_ZSCORE_OPTIONS_NONE;
%constant unsigned int ZSCORE_FILTER_ON         =  VRNA_ZSCORE_FILTER_ON;
%constant unsigned int ZSCORE_PRE_FILTER        = VRNA_ZSCORE_PRE_FILTER;
%constant unsigned int ZSCORE_REPORT_SUBSUMED   = VRNA_ZSCORE_REPORT_SUBSUMED;
%constant unsigned int ZSCORE_MODEL_DEFAULT     = VRNA_ZSCORE_MODEL_DEFAULT;
%constant unsigned int ZSCORE_SETTINGS_DEFAULT  = VRNA_ZSCORE_SETTINGS_DEFAULT;


%include <ViennaRNA/zscore.h>
#endif
