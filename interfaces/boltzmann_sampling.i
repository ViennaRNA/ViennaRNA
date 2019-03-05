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

/*
  we need this crazy piece of argout typemap only because we don't
  want the vrna_nr_memory_t object to be appended to the results(list),
  but prepended instead. Otherwise, a simple
  %append_output(SWIG_NewPointerObj(SWIG_as_voidptr(retval$argnum), $1_descriptor, 0));
  would have sufficed already
*/
%typemap(argout) vrna_nr_memory_t * {
#ifdef SWIGPYTHON
  PyObject *o, *o2, *o3;
  o = SWIG_NewPointerObj(SWIG_as_voidptr(retval$argnum), $1_descriptor, 1);
  if ((!$result) || ($result == Py_None)) {
    $result = o;
  } else {
    PyObject *o2 = $result;
    $result = PyTuple_New(1);
    PyTuple_SetItem($result,0,o2);
    o3 = PyTuple_New(1);
    PyTuple_SetItem(o3,0,o);
    o2 = $result;
    $result = PySequence_Concat(o3,o2);
    Py_DECREF(o2);
    Py_DECREF(o3);
  }
#elif defined(SWIGPERL5)
  /* increase output stack if necessary */
  if (argvi >= items) {
    EXTEND(sp,1);
  }
  /* move already existing return values to the back */
  for (int i = argvi; i > 0; i--) {
    ST(i) = ST(i - 1);
  }
  /* store result as first element in the stack */
  ST(0) = SWIG_NewPointerObj(SWIG_as_voidptr(retval$argnum), $1_descriptor, 0);
  /* increase return argument counter */
  argvi++;
#endif
}


%typemap(in) vrna_nr_memory_t *INOUT (vrna_nr_memory_t *retval)
{
#ifdef SWIGPYTHON
  if ($input == Py_None) {
#elif defined(SWIGPERL5)
  if (!SvOK($input)) {
#endif
    retval = new vrna_nr_memory_t();
    $1 = retval;
  } else {
    /* INOUT in */
    SWIG_ConvertPtr($input,SWIG_as_voidptrptr(&retval), 0, SWIG_POINTER_DISOWN);
    $1 = retval;
  }
}

%extend vrna_fold_compound_t {

  char *
  pbacktrack(void)
  {
    return vrna_pbacktrack($self);
  }

  char *
  pbacktrack(int length)
  {
    return vrna_pbacktrack5($self, length);
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
