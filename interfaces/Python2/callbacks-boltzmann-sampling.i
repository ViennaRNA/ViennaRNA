/**********************************************/
/* BEGIN interface for subopt callback        */
/**********************************************/

#ifdef SWIGPYTHON
%{

#include <stdexcept>

typedef struct {
  PyObject *cb;
  PyObject *data;
} python_bs_callback_t;

static python_bs_callback_t *
bind_bs_callback(PyObject *PyFunc,
                 PyObject *data);

static void
python_wrap_bs_cb(const char *structure,
                  void       *data);

static python_bs_callback_t *
bind_bs_callback(PyObject *PyFunc,
                 PyObject *data)
{
  python_bs_callback_t *cb = (python_bs_callback_t *)vrna_alloc(sizeof(python_bs_callback_t));

  Py_INCREF(PyFunc);
  Py_INCREF(data);

  cb->cb    = PyFunc;  /* store callback */
  cb->data  = data;    /* bind data */

  return cb;
}


static void
release_bs_callback(python_bs_callback_t *cb)
{
  Py_DECREF(cb->cb);
  Py_DECREF(cb->data);
  free(cb); 
}


static void
python_wrap_bs_cb(const char *structure,
                  void       *data)
{
  PyObject *func, *arglist, *result, *err;
  python_bs_callback_t *cb = (python_bs_callback_t *)data;

  func = cb->cb;
  /* compose argument list */
  PyObject *py_structure, *py_energy;
  py_structure = (structure) ? PyString_FromString(structure) : Py_None;
  result       = PyObject_CallFunctionObjArgs(func,
                                              py_structure,
                                              (cb->data) ? cb->data : Py_None,
                                              NULL);

  if (py_structure != Py_None)
    Py_DECREF(py_structure);

  /* BEGIN recognizing errors in callback execution */
  if (result == NULL) {
    if ((err = PyErr_Occurred())) {
      /* print error message */
      PyErr_Print();
      /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
      if (PyErr_GivenExceptionMatches(err, PyExc_TypeError)) {
        throw std::runtime_error( "Boltzmann sampling callback must take exactly 2 arguments" );
      } else {
        throw std::runtime_error( "Some error occurred while executing Boltzmann sampling callback" );
      }
    }
    PyErr_Clear();
  }
  /* END recognizing errors in callback execution */

  Py_XDECREF(result);

  return /*void*/;
}

%}

%feature("autodoc", "2") vrna_fold_compound_t::pbacktrack;

/* now we bind the above functions as methods to the fold_compound object */
%extend vrna_fold_compound_t {

  unsigned int
  pbacktrack5(unsigned int num_samples,
              unsigned int length,
              PyObject     *PyFunc,
              PyObject     *data   = Py_None,
              unsigned int options = VRNA_PBACKTRACK_DEFAULT)
  {
    unsigned int i;
    python_bs_callback_t *cb = bind_bs_callback(PyFunc, data);

    i = vrna_pbacktrack5_cb($self,
                            num_samples,
                            length,
                            &python_wrap_bs_cb,
                            (void *)cb,
                            options);

    release_bs_callback(cb);

    return i;
  }

  unsigned int
  pbacktrack(unsigned int num_samples,
             PyObject     *PyFunc,
             PyObject     *data   = Py_None,
             unsigned int options = VRNA_PBACKTRACK_DEFAULT)
  {
    unsigned int i;
    python_bs_callback_t *cb = bind_bs_callback(PyFunc, data);

    i = vrna_pbacktrack_cb($self,
                           num_samples,
                           &python_wrap_bs_cb,
                           (void *)cb,
                           options);

    release_bs_callback(cb);

    return i;
  }

  %apply vrna_pbacktrack_mem_t *INOUT { vrna_pbacktrack_mem_t *nr_memory };

  unsigned int
  pbacktrack(unsigned int          num_samples,
                          PyObject              *PyFunc,
                          PyObject              *data,
                          vrna_pbacktrack_mem_t *nr_memory,
                          unsigned int          options = VRNA_PBACKTRACK_DEFAULT)
  {
    unsigned int i;
    python_bs_callback_t *cb = bind_bs_callback(PyFunc, data);

    i = vrna_pbacktrack_resume_cb($self,
                                  num_samples,
                                  &python_wrap_bs_cb,
                                  (void *)cb,
                                  nr_memory,
                                  options);

    release_bs_callback(cb);

    return i;
  }

  unsigned int
  pbacktrack5(unsigned int          num_samples,
              unsigned int          length,
              PyObject              *PyFunc,
              PyObject              *data,
              vrna_pbacktrack_mem_t *nr_memory,
              unsigned int          options = VRNA_PBACKTRACK_DEFAULT)
  {
    unsigned int i;
    python_bs_callback_t *cb = bind_bs_callback(PyFunc, data);

    i = vrna_pbacktrack5_resume_cb($self,
                                   num_samples,
                                   length,
                                   &python_wrap_bs_cb,
                                   (void *)cb,
                                   nr_memory,
                                   options);

    release_bs_callback(cb);

    return i;
  }

  %clear vrna_pbacktrack_mem_t *nr_memory;
}

#endif
