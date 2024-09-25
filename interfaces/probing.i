/**********************************************/
/* BEGIN interface for structure probing      */
/* derived constraints                        */
/**********************************************/

%rename(probing_data) vrna_probing_data_s;

typedef struct {} vrna_probing_data_s;

/* no default constructor / destructor */
%nodefaultctor vrna_probing_data_s;
%nodefaultdtor vrna_probing_data_s;

#ifdef SWIGPYTHON
%define SWG_PROBING_DATA_CONSTRUCTOR_DS
"The :py:func:`probing_data` constructor can be invoked in different ways.
Depending on the input parameters, the object returned will invoke
particular probing data integration methods either for single
sequences or multiple sequence alignments. Hence, it subsumes and
provides an easy interface for the functions
:py:func:`probing_data_Deigan2009`, :py:func:`probing_data_Deigan2009_comparative`,
etc.


When multiple sets of probing data are supplied, the constructor
assumes preparations for multiple sequence alignments (MSAs). As
a consequence, the parameters for the conversion methods must be
supplied as lists of parameters, one for each sequence. 


Parameters
----------
reactivities: list(double) or list(list(double)
    single sequence probing data (1-based) or multiple sequence probing data (0-based list for each sequence of 1-based data)
m: double
    slope for the Deigan et al. 2009 method
b: double
    intercept for the Deigan et al. 2009 method
ms: list(double)
    multiple slopes for the Deigan et al. 2009 method (0-based, one for each sequence)
bs: list(double)
    multiple intercepts for the Deigan et al. 2009 method (0-based, one for each sequence)
beta: double
    scaling factor for Zarringhalam et al. 2012 method
betas: double
    multiple scaling factors for Zarringhalam et al. 2012 method (0-based, one for each sequence)
pr_conversion: string
    probing data to conversion strategy 
pr_conversions: list(string)
    multiple probing data to conversion strategies (0-based, one for each sequence)
pr_default: double
    default probability for a nucleotide where reactivity data is missing for
pr_defaults: list(double)
    list of default probabilities for a nucleotide where reactivity data is missing for (0-based, one for each sequence)
"
%enddef

%feature("docstring", SWG_PROBING_DATA_CONSTRUCTOR_DS) vrna_probing_data_s::vrna_probing_data_s;
%feature("autodoc") vrna_probing_data_s::vrna_probing_data_s;
%feature("kwargs") vrna_probing_data_s::vrna_probing_data_s;
#endif

%extend vrna_probing_data_s {

  /* constructor for Deigan method single sequence */
  vrna_probing_data_s(std::vector<double> reactivities,
                      double              m,
                      double              b)
  {
    vrna_probing_data_s *obj = vrna_probing_data_Deigan2009(&(reactivities[0]),
                                                            reactivities.size(),
                                                            m,
                                                            b);
    return obj;
  }

  /* constructor for Deigan method MSA with multiple parameters */
  vrna_probing_data_s(std::vector< std::vector<double> > reactivities,
                      std::vector<double>              ms,
                      std::vector<double>              bs)
  {
    unsigned int              n_seq = reactivities.size();
    std::vector<unsigned int> ns;
    vrna_probing_data_s       *obj;
    double                    **d;

    for (unsigned int i = 0; i < reactivities.size(); i++)
      ns.push_back(reactivities[i].size());

    unsigned int options = VRNA_PROBING_METHOD_MULTI_PARAMS_0;
    if (ms.size() > 1)
      options |= VRNA_PROBING_METHOD_MULTI_PARAMS_1;

    if (bs.size() > 1)
      options |= VRNA_PROBING_METHOD_MULTI_PARAMS_2;

    d = (double **)vrna_alloc(sizeof(double *) * reactivities.size());
    for (unsigned int i = 0; i < reactivities.size(); i++)
      if (reactivities[i].size() > 0) {
        d[i] = (double *)vrna_alloc(sizeof(double) * reactivities[i].size());
        d[i] = (double *)memcpy(d[i], &(reactivities[i][0]), sizeof(double) * reactivities[i].size());
      }

    /* try interpreting the json input as file name */
    obj = vrna_probing_data_Deigan2009_comparative((const double **)d,
                                                   &(ns[0]),
                                                   reactivities.size(),
                                                   &(ms[0]),
                                                   &(bs[0]),
                                                   options);

    for (unsigned int i = 0; i < reactivities.size(); i++)
      free(d[i]);
    free(d);

    return obj;
  }


  /* constructor for Zarringhalam method single sequence */
  vrna_probing_data_s(std::vector<double> reactivities,
                      double              beta,
                      std::string         pr_conversion = VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_conversion,
                      double              pr_default    = VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_probability)
  {
    vrna_probing_data_s *obj = vrna_probing_data_Zarringhalam2012(&(reactivities[0]),
                                                                  reactivities.size(),
                                                                  beta,
                                                                  pr_conversion.c_str(),
                                                                  pr_default);
    return obj;
  }


  /* constructor for Zarringhalam method MSA with multiple parameters */
  vrna_probing_data_s(std::vector< std::vector<double> > reactivities,
                      std::vector<double>              betas,
                      std::vector<std::string>         pr_conversions = std::vector<std::string>(),
                      std::vector<double>              pr_defaults = std::vector<double>() )
  {
    unsigned int              n_seq = reactivities.size();
    std::vector<unsigned int> ns;
    vrna_probing_data_s       *obj;
    double                    **d;
    char                      **c;

    if (pr_conversions.size() == 0)
      pr_conversions.push_back(VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_conversion);

    if (pr_defaults.size() == 0)
      pr_defaults.push_back(VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_probability);

    for (unsigned int i = 0; i < reactivities.size(); i++)
      ns.push_back(reactivities[i].size());

    unsigned int options = VRNA_PROBING_METHOD_MULTI_PARAMS_0;
    if (betas.size() > 1)
      options |= VRNA_PROBING_METHOD_MULTI_PARAMS_1;

    if (pr_conversions.size() > 1)
      options |= VRNA_PROBING_METHOD_MULTI_PARAMS_2;

    if (pr_defaults.size() > 1)
      options |= VRNA_PROBING_METHOD_MULTI_PARAMS_3;

    c = (char **)vrna_alloc(sizeof(char *) * pr_conversions.size());
    for (unsigned int i = 0; i < pr_conversions.size(); i++)
      if (pr_conversions[i].size() > 0) {
        c[i] = (char *)vrna_alloc(sizeof(char) * (pr_conversions[i].size() + 1));
        c[i] = (char *)memcpy(c[i], pr_conversions[i].c_str(), sizeof(char *) * (pr_conversions[i].size() + 1));
      }

    d = (double **)vrna_alloc(sizeof(double *) * reactivities.size());
    for (unsigned int i = 0; i < reactivities.size(); i++)
      if (reactivities[i].size() > 0) {
        d[i] = (double *)vrna_alloc(sizeof(double) * reactivities[i].size());
        d[i] = (double *)memcpy(d[i], &(reactivities[i][0]), sizeof(double) * reactivities[i].size());
      }

    /* try interpreting the json input as file name */
    obj = vrna_probing_data_Zarringhalam2012_comparative((const double **)d,
                                                         &(ns[0]),
                                                         reactivities.size(),
                                                         &(betas[0]),
                                                         (const char **)c,
                                                         &(pr_defaults[0]),
                                                         options);

    for (unsigned int i = 0; i < reactivities.size(); i++) {
      free(d[i]);
      free(c[i]);
    }
    free(d);
    free(c);

    return obj;
  }


  ~vrna_probing_data_s()
  {
    vrna_probing_data_free($self);
  }
}

%rename (probing_data_free) vrna_probing_data_free;

%{
#include <vector>

  vrna_probing_data_s *
  probing_data_Deigan2009(std::vector<double> reactivities,
                          double              m,
                          double              b)
  {
    vrna_probing_data_s *obj = vrna_probing_data_Deigan2009(&(reactivities[0]),
                                                            reactivities.size(),
                                                            m,
                                                            b);
    return obj;
  }


  vrna_probing_data_s *
  probing_data_Deigan2009_comparative(std::vector<std::vector<double> > reactivities,
                                      std::vector<double>               ms,
                                      std::vector<double>               bs,
                                      unsigned int                      multi_params = VRNA_PROBING_METHOD_MULTI_PARAMS_0)
  {
    unsigned int              n_seq = reactivities.size();
    std::vector<unsigned int> ns;
    vrna_probing_data_s       *obj;
    double                    **d;

    for (unsigned int i = 0; i < reactivities.size(); i++)
      ns.push_back(reactivities[i].size());

    d = (double **)vrna_alloc(sizeof(double *) * reactivities.size());
    for (unsigned int i = 0; i < reactivities.size(); i++)
      if (reactivities[i].size() > 0) {
        d[i] = (double *)vrna_alloc(sizeof(double) * reactivities[i].size());
        d[i] = (double *)memcpy(d[i], &(reactivities[i][0]), sizeof(double) * reactivities[i].size());
      }

    /* try interpreting the json input as file name */
    obj = vrna_probing_data_Deigan2009_comparative((const double **)d,
                                                   &(ns[0]),
                                                   n_seq,
                                                   &(ms[0]),
                                                   &(bs[0]),
                                                   multi_params);

    for (unsigned int i = 0; i < reactivities.size(); i++)
      free(d[i]);
    free(d);

    return obj;
  }


  vrna_probing_data_s *
  probing_data_Zarringhalam2012(std::vector<double> reactivities,
                                double              beta,
                                std::string         pr_conversion = VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_conversion,
                                double              pr_default    = VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_probability)
  {
    vrna_probing_data_s *obj = vrna_probing_data_Zarringhalam2012(&(reactivities[0]),
                                                                  reactivities.size(),
                                                                  beta,
                                                                  pr_conversion.c_str(),
                                                                  pr_default);
    return obj;
  }


  vrna_probing_data_s *
  probing_data_Zarringhalam2012_comparative(std::vector< std::vector<double> > reactivities,
                                            std::vector<double>              betas,
                                            std::vector<std::string>         pr_conversions = std::vector<std::string>(),
                                            std::vector<double>              pr_defaults = std::vector<double>(),
                                            unsigned int                     multi_params = VRNA_PROBING_METHOD_MULTI_PARAMS_0)
  {
    unsigned int              n_seq = reactivities.size();
    std::vector<unsigned int> ns;
    vrna_probing_data_s       *obj;
    double                    **d;
    char                      **c;

    if (pr_conversions.size() == 0)
      pr_conversions.push_back(VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_conversion);

    if (pr_defaults.size() == 0)
      pr_defaults.push_back(VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_probability);

    for (unsigned int i = 0; i < reactivities.size(); i++)
      ns.push_back(reactivities[i].size());

    c = (char **)vrna_alloc(sizeof(char *) * pr_conversions.size());
    for (unsigned int i = 0; i < pr_conversions.size(); i++)
      if (pr_conversions[i].size() > 0) {
        c[i] = (char *)vrna_alloc(sizeof(char) * (pr_conversions[i].size() + 1));
        c[i] = (char *)memcpy(c[i], pr_conversions[i].c_str(), sizeof(char *) * (pr_conversions[i].size() + 1));
      }

    d = (double **)vrna_alloc(sizeof(double *) * reactivities.size());
    for (unsigned int i = 0; i < reactivities.size(); i++)
      if (reactivities[i].size() > 0) {
        d[i] = (double *)vrna_alloc(sizeof(double) * reactivities[i].size());
        d[i] = (double *)memcpy(d[i], &(reactivities[i][0]), sizeof(double) * reactivities[i].size());
      }

    /* try interpreting the json input as file name */
    obj = vrna_probing_data_Zarringhalam2012_comparative((const double **)d,
                                                         &(ns[0]),
                                                         reactivities.size(),
                                                         &(betas[0]),
                                                         (const char **)c,
                                                         &(pr_defaults[0]),
                                                         multi_params);

    for (unsigned int i = 0; i < reactivities.size(); i++) {
      free(d[i]);
      free(c[i]);
    }
    free(d);
    free(c);

    return obj;
  }


  vrna_probing_data_s *
  probing_data_Eddy2014_2(std::vector<double> reactivities,
                          std::vector<double> unpaired_data,
                          std::vector<double> paired_data)
  {
    vrna_probing_data_s *obj = vrna_probing_data_Eddy2014_2(&(reactivities[0]),
                                                            reactivities.size(),
                                                            &(unpaired_data[0]),
                                                            unpaired_data.size(),
                                                            &(paired_data[0]),
                                                            paired_data.size());
    return obj;
  }


  vrna_probing_data_s *
  probing_data_Eddy2014_2_comparative(std::vector< std::vector<double> >  reactivities,
                                      std::vector< std::vector<double> >  unpaired_data,
                                      std::vector< std::vector<double> >  paired_data,
                                      unsigned int                        multi_params = VRNA_PROBING_METHOD_MULTI_PARAMS_0)
  {
    unsigned int              n_seq = reactivities.size();
    std::vector<unsigned int> ns;
    std::vector<unsigned int> us;
    std::vector<unsigned int> bs;
    vrna_probing_data_s       *obj;
    double                    **d;
    double                    **up;
    double                    **bp;

    for (unsigned int i = 0; i < reactivities.size(); i++) {
      ns.push_back(reactivities[i].size());
      us.push_back(unpaired_data[i].size());
      bs.push_back(paired_data[i].size());
    }

    up = (double **)vrna_alloc(sizeof(double *) * unpaired_data.size());
    for (unsigned int i = 0; i < unpaired_data.size(); i++)
      if (unpaired_data[i].size() > 0) {
        up[i] = (double *)vrna_alloc(sizeof(double) * (unpaired_data[i].size() + 1));
        up[i] = (double *)memcpy(up[i], &(unpaired_data[i][0]), sizeof(double *) * (unpaired_data[i].size() + 1));
      }

    bp = (double **)vrna_alloc(sizeof(double *) * paired_data.size());
    for (unsigned int i = 0; i < paired_data.size(); i++)
      if (paired_data[i].size() > 0) {
        bp[i] = (double *)vrna_alloc(sizeof(double) * (paired_data[i].size() + 1));
        bp[i] = (double *)memcpy(bp[i], &(paired_data[i][0]), sizeof(double *) * (paired_data[i].size() + 1));
      }

    d = (double **)vrna_alloc(sizeof(double *) * reactivities.size());
    for (unsigned int i = 0; i < reactivities.size(); i++)
      if (reactivities[i].size() > 0) {
        d[i] = (double *)vrna_alloc(sizeof(double) * reactivities[i].size());
        d[i] = (double *)memcpy(d[i], &(reactivities[i][0]), sizeof(double) * reactivities[i].size());
      }

    /* try interpreting the json input as file name */
    obj = vrna_probing_data_Eddy2014_2_comparative((const double **)d,
                                                   &(ns[0]),
                                                   n_seq,
                                                   (const double **)up,
                                                   &(us[0]),
                                                   (const double **)bp,
                                                   &(bs[0]),
                                                   multi_params);

    for (unsigned int i = 0; i < reactivities.size(); i++) {
      free(d[i]);
      free(up[i]);
      free(bp[i]);
    }
    free(d);
    free(up);
    free(bp);

    return obj;
  }

%}

#ifdef SWIGPYTHON
%feature("autodoc") probing_data_Deigan2009;
%feature("kwargs") probing_data_Deigan2009;
%feature("autodoc") probing_data_Deigan2009_comparative;
%feature("kwargs") probing_data_Deigan2009_comparative;
%feature("autodoc") probing_data_Zarringhalam2012;
%feature("kwargs") probing_data_Zarringhalam2012;
%feature("autodoc") probing_data_Zarringhalam2012_comparative;
%feature("kwargs") probing_data_Zarringhalam2012_comparative;
%feature("autodoc") probing_data_Eddy2014_2;
%feature("kwargs") probing_data_Eddy2014_2;
%feature("autodoc") probing_data_Eddy2014_2_comparative;
%feature("kwargs") probing_data_Eddy2014_2_comparative;
#endif

vrna_probing_data_s *
probing_data_Deigan2009(std::vector<double> reactivities,
                        double              m,
                        double              b);


vrna_probing_data_s *
probing_data_Deigan2009_comparative(std::vector<std::vector<double> > reactivities,
                                    std::vector<double>               ms,
                                    std::vector<double>               bs,
                                    unsigned int                      multi_params = VRNA_PROBING_METHOD_MULTI_PARAMS_0);


vrna_probing_data_s *
probing_data_Zarringhalam2012(std::vector<double> reactivities,
                              double              beta,
                              std::string         pr_conversion = VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_conversion,
                              double              pr_default    = VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_probability);


vrna_probing_data_s *
probing_data_Zarringhalam2012_comparative(std::vector< std::vector<double> > reactivities,
                                          std::vector<double>              betas,
                                          std::vector<std::string>         pr_conversions = std::vector<std::string>(),
                                          std::vector<double>              pr_defaults = std::vector<double>(),
                                          unsigned int                     multi_params = VRNA_PROBING_METHOD_MULTI_PARAMS_0);


vrna_probing_data_s *
probing_data_Eddy2014_2(std::vector<double> reactivities,
                        std::vector<double> unpaired_data,
                        std::vector<double> paired_data);


vrna_probing_data_s *
probing_data_Eddy2014_2_comparative(std::vector< std::vector<double> >  reactivities,
                                    std::vector< std::vector<double> >  unpaired_data,
                                    std::vector< std::vector<double> >  paired_data,
                                    unsigned int                        multi_params = VRNA_PROBING_METHOD_MULTI_PARAMS_0);


%extend vrna_fold_compound_t {

#ifdef SWIGPYTHON
%feature("autodoc") sc_probing;
%feature("kwargs") sc_probing;
#endif

  int
  sc_probing(vrna_probing_data_t data)
  {
    return vrna_sc_probing($self, data);
  }
}



%constant unsigned int  PROBING_METHOD_DEIGAN2009                           = VRNA_PROBING_METHOD_DEIGAN2009;
%constant double        PROBING_METHOD_DEIGAN2009_DEFAULT_m                 = VRNA_PROBING_METHOD_DEIGAN2009_DEFAULT_m;
%constant double        PROBING_METHOD_DEIGAN2009_DEFAULT_b                 = VRNA_PROBING_METHOD_DEIGAN2009_DEFAULT_b;
%constant unsigned int  PROBING_METHOD_ZARRINGHALAM2012                     = VRNA_PROBING_METHOD_ZARRINGHALAM2012;
%constant double        PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_beta        = VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_beta;
%constant char *        PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_conversion  = VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_conversion;
%constant double        PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_probability = VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_probability;
%constant unsigned int  PROBING_METHOD_WASHIETL2012                         = VRNA_PROBING_METHOD_WASHIETL2012;
%constant unsigned int  PROBING_METHOD_EDDY2014_2                           = VRNA_PROBING_METHOD_EDDY2014_2;
%constant unsigned int  PROBING_METHOD_MULTI_PARAMS_0                       = VRNA_PROBING_METHOD_MULTI_PARAMS_0;
%constant unsigned int  PROBING_METHOD_MULTI_PARAMS_1                       = VRNA_PROBING_METHOD_MULTI_PARAMS_1;
%constant unsigned int  PROBING_METHOD_MULTI_PARAMS_2                       = VRNA_PROBING_METHOD_MULTI_PARAMS_2;
%constant unsigned int  PROBING_METHOD_MULTI_PARAMS_3                       = VRNA_PROBING_METHOD_MULTI_PARAMS_3;
%constant unsigned int  PROBING_METHOD_MULTI_PARAMS_DEFAULT                 = VRNA_PROBING_METHOD_MULTI_PARAMS_DEFAULT;
%constant unsigned int  PROBING_DATA_CHECK_SEQUENCE                         = VRNA_PROBING_DATA_CHECK_SEQUENCE;

%include  <ViennaRNA/probing/basic.h>
