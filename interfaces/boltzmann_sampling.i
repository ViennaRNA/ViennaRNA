/**********************************************/
/* BEGIN interface for Boltzmann Sampling     */
/**********************************************/

%newobject pbacktrack;
extern char *pbacktrack(char *sequence);

/* these functions remain for now due to backward compatibility reasons
%ignore pbacktrack;
%ignore pbacktrack5;
%ignore pbacktrack_circ;
*/

/* attach stochastic backtracking functions as method of fold_compound */
%newobject vrna_fold_compound_t::pbacktrack;

%newobject vrna_fold_compound_t::pbacktrack_nr_resume;

%ignore vrna_pbacktrack_mem_t;
%ignore vrna_nr_memory_s;

%rename (pbacktrack_mem)  vrna_pbacktrack_mem_t;

typedef struct {} vrna_pbacktrack_mem_t;

%nodefaultctor vrna_pbacktrack_mem_t;
%nodefaultdtor vrna_pbacktrack_mem_t;

%extend vrna_pbacktrack_mem_t {
  vrna_pbacktrack_mem_t() {
    vrna_pbacktrack_mem_t *m = (vrna_pbacktrack_mem_t *)vrna_alloc(sizeof(vrna_pbacktrack_mem_t));
    *m = NULL;
    return m;
  }
  ~vrna_pbacktrack_mem_t() {
    vrna_pbacktrack_mem_free(*$self);
    delete $self;
  }
}

#ifdef SWIGPYTHON
%feature("autodoc")vrna_fold_compound_t::pbacktrack5;
%feature("kwargs") vrna_fold_compound::pbacktrack5;
%feature("autodoc")vrna_fold_compound_t::pbacktrack;
%feature("kwargs") vrna_fold_compound::pbacktrack;
#endif

%extend vrna_fold_compound_t {

  char *
  pbacktrack(void)
  {
    return vrna_pbacktrack($self);
  }

  char *
  pbacktrack5(unsigned int length)
  {
    return vrna_pbacktrack5($self, length);
  }

  std::vector<std::string>
  pbacktrack(unsigned int num_samples,
             unsigned int options = VRNA_PBACKTRACK_DEFAULT)
  {
    std::vector<std::string> str_vec;
    char  **ptr, **output;

    output = vrna_pbacktrack_num($self, num_samples, options);

    if (output) {
      for (ptr = output; *ptr != NULL; ptr++) {
        str_vec.push_back(std::string(*ptr));
        free(*ptr);
      }

      free(output);
    }

    return str_vec;
  }

  std::vector<std::string>
  pbacktrack5(unsigned int num_samples,
              unsigned int length,
              unsigned int options = VRNA_PBACKTRACK_DEFAULT)
  {
    std::vector<std::string> str_vec;
    char  **ptr, **output;

    output = vrna_pbacktrack5_num($self, num_samples, length, options);

    if (output) {
      for (ptr = output; *ptr != NULL; ptr++) {
        str_vec.push_back(std::string(*ptr));
        free(*ptr);
      }

      free(output);
    }

    return str_vec;
  }

  %apply vrna_pbacktrack_mem_t *INOUT { vrna_pbacktrack_mem_t *nr_memory };

  std::vector<std::string>
  pbacktrack(unsigned int          num_samples,
             vrna_pbacktrack_mem_t *nr_memory,
             unsigned int          options = VRNA_PBACKTRACK_DEFAULT)
  {
    std::vector<std::string> str_vec;

    char **ptr, **output = vrna_pbacktrack_resume($self,
                                                  num_samples,
                                                  nr_memory,
                                                  options);

    if (output) {
      for (ptr = output; *ptr != NULL; ptr++) {
        str_vec.push_back(std::string(*ptr));
        free(*ptr);
      }

      free(output);
    }

    return str_vec;
  }

  std::vector<std::string>
  pbacktrack5(unsigned int          num_samples,
              unsigned int          length,
              vrna_pbacktrack_mem_t *nr_memory,
              unsigned int          options = VRNA_PBACKTRACK_DEFAULT)
  {
    std::vector<std::string> str_vec;

    char **ptr, **output;
    
    output = vrna_pbacktrack5_resume($self,
                                     num_samples,
                                     length,
                                     nr_memory,
                                     options);

    if (output) {
      for (ptr = output; *ptr != NULL; ptr++) {
        str_vec.push_back(std::string(*ptr));
        free(*ptr);
      }

      free(output);
    }

    return str_vec;
  }

  %clear vrna_pbacktrack_mem_t *nr_memory;
}

%constant unsigned int PBACKTRACK_DEFAULT       = VRNA_PBACKTRACK_DEFAULT;
%constant unsigned int PBACKTRACK_NON_REDUNDANT = VRNA_PBACKTRACK_NON_REDUNDANT;

%include  <ViennaRNA/boltzmann_sampling.h>
