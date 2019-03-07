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

%ignore vrna_nr_memory_t;
%ignore vrna_nr_memory_s;

%rename (nr_memory)  vrna_nr_memory_t;

typedef struct {} vrna_nr_memory_t;

%nodefaultctor vrna_nr_memory_t;
%nodefaultdtor vrna_nr_memory_t;

%extend vrna_nr_memory_t {
  vrna_nr_memory_t() {
    vrna_nr_memory_t *m = (vrna_nr_memory_t *)vrna_alloc(sizeof(vrna_nr_memory_t));
    *m = NULL;
    return m;
  }
  ~vrna_nr_memory_t() {
    vrna_pbacktrack_nr_free(*$self);
    delete $self;
  }
}

#ifdef SWIGPYTHON
%feature("autodoc")vrna_fold_compound_t::pbacktrack;
%feature("kwargs") vrna_fold_compound::pbacktrack;
%feature("autodoc")vrna_fold_compound_t::pbacktrack_nr;
%feature("kwargs") vrna_fold_compound::pbacktrack_nr;
%feature("autodoc")vrna_fold_compound_t::pbacktrack_nr_resume;
%feature("kwargs") vrna_fold_compound::pbacktrack_nr_resume;
#endif

%extend vrna_fold_compound_t {

  char *
  pbacktrack(void)
  {
    return vrna_pbacktrack($self);
  }

  char *
  pbacktrack(unsigned int length)
  {
    return vrna_pbacktrack5($self, length);
  }

  std::vector<std::string>
  pbacktrack(unsigned int length, unsigned int num_samples)
  {
    std::vector<std::string> str_vec;
    char  **ptr, **output;

    if (length == 0)
      output = vrna_pbacktrack_num($self, num_samples);
    else
      output = vrna_pbacktrack5_num($self, length, num_samples);

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
  pbacktrack_nr(unsigned int num_samples)
  {
    std::vector<std::string> str_vec;

    char **ptr, **output = vrna_pbacktrack_nr($self, num_samples);

    if (output) {
      for (ptr = output; *ptr != NULL; ptr++) {
        str_vec.push_back(std::string(*ptr));
        free(*ptr);
      }

      free(output);
    }

    return str_vec;
  }

  %apply vrna_nr_memory_t *INOUT { vrna_nr_memory_t *nr_memory };

  std::vector<std::string>
  pbacktrack_nr_resume(unsigned int      num_samples,
                       vrna_nr_memory_t  *nr_memory)
  {
    std::vector<std::string> str_vec;

    char **ptr, **output = vrna_pbacktrack_nr_resume($self,
                                                     num_samples,
                                                     nr_memory);

    if (output) {
      for (ptr = output; *ptr != NULL; ptr++) {
        str_vec.push_back(std::string(*ptr));
        free(*ptr);
      }

      free(output);
    }

    return str_vec;
  }

  %clear vrna_nr_memory_t *nr_memory;
}

%include  <ViennaRNA/boltzmann_sampling.h>
