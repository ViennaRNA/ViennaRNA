/**********************************************/
/* BEGIN interface for PF window callback     */
/**********************************************/

#ifdef SWIGPYTHON
%{

#include <stdexcept>

typedef struct {
  PyObject *cb;
  PyObject *data;
} python_pf_window_callback_t;


static python_pf_window_callback_t *
bind_pf_window_callback(PyObject *PyFunc,
                        PyObject *data);

static void
python_wrap_pf_window_cb(FLT_OR_DBL   *pr,
                         int          pr_size,
                         int          i,
                         int          max,
                         unsigned int type,
                         void         *data);


static python_pf_window_callback_t *
bind_pf_window_callback(PyObject *PyFunc,
                        PyObject *data)
{
  python_pf_window_callback_t *cb = (python_pf_window_callback_t *)vrna_alloc(sizeof(python_pf_window_callback_t));

  Py_INCREF(PyFunc);
  Py_INCREF(data);

  cb->cb    = PyFunc;  /* store callback */
  cb->data  = data;    /* bind data */

  return cb;
}

static void
release_pf_window_callback(python_pf_window_callback_t *cb)
{
  Py_DECREF(cb->cb);
  Py_DECREF(cb->data);
  free(cb); 
}

static void
python_wrap_pf_window_cb(FLT_OR_DBL   *pr,
                         int          pr_size,
                         int          i,
                         int          max,
                         unsigned int type,
                         void *data)
{
  PyObject *func, *arglist, *result, *pr_list, *err;
  python_pf_window_callback_t *cb = (python_pf_window_callback_t *)data;

  func = cb->cb;

  if (type & VRNA_PROBS_WINDOW_UP) { /* We distinguish output for unpaired probabilities */

    /* create PYTHON list for unpaired probabilities */
    pr_list = PyList_New((Py_ssize_t) max + 1);

    /* 0th element */
    Py_INCREF(Py_None);
    PyList_SET_ITEM(pr_list, 0, Py_None);

    /* actual values in range [1, MIN(i, max)] */
    for (int cnt = 1; cnt <= pr_size; cnt++)
      PyList_SET_ITEM(pr_list, (Py_ssize_t) cnt, PyFloat_FromDouble(pr[cnt]));

    /* empty values in range [0,i - 1] */
    for (int cnt = pr_size + 1; cnt <= max; cnt++) {
      Py_INCREF(Py_None);
      PyList_SET_ITEM(pr_list, (Py_ssize_t) cnt, Py_None);
    }
  } else { /* and pairing/stacking probabilities for pair (i, j) or ensemble free energies for subsegment [i, j] */

    /* create PYTHON list for pr values */
    pr_list = PyList_New((Py_ssize_t) (pr_size + 1));

    /* empty values in range [0, i] */
    for (int cnt = 0; cnt <= i; cnt++) {
      Py_INCREF(Py_None);
      PyList_SET_ITEM(pr_list, (Py_ssize_t) cnt, Py_None);
    }

    /* actual values in range [i + 1, pr_size] */
    for (int cnt = i + 1; cnt <= pr_size; cnt++)
      PyList_SET_ITEM(pr_list, (Py_ssize_t) cnt, PyFloat_FromDouble(pr[cnt]));
  }

  /* compose argument list */
  PyObject *py_size, *py_i, *py_max, *py_type;
  py_size = PyLong_FromLong(pr_size);
  py_i    = PyLong_FromLong(i);
  py_max  = PyLong_FromLong(max);
  py_type = PyLong_FromLong(type);
  result  = PyObject_CallFunctionObjArgs(func,
                                         pr_list,
                                         py_size,
                                         py_i,
                                         py_max,
                                         py_type,
                                         (cb->data) ? cb->data : Py_None,
                                         NULL);

  Py_DECREF(py_size);
  Py_DECREF(py_i);
  Py_DECREF(py_max);
  Py_DECREF(py_type);

  /* BEGIN recognizing errors in callback execution */
  if (result == NULL) {
    if ((err = PyErr_Occurred())) {
      /* print error message */
      PyErr_Print();
      /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
      if (PyErr_GivenExceptionMatches(err, PyExc_TypeError)) {
        throw std::runtime_error( "Sliding window partition function callback must take exactly 6 arguments" );
      } else {
        throw std::runtime_error( "Some error occurred while executing sliding window partition function callback" );
      }
    }
    PyErr_Clear();
  }
  /* END recognizing errors in callback execution */

  Py_XDECREF(result);

  return /*void*/;
}

%}

/* now we bind vrna_probs_window() as method to the fold_compound object using the above callback wrapper */
%extend vrna_fold_compound_t {

%feature("autodoc") probs_window;
%feature("kwargs") probs_window;

  int
  probs_window(int          ulength,
               unsigned int options,
               PyObject     *PyFunc,
               PyObject     *data = Py_None)
  {
    python_pf_window_callback_t *cb = bind_pf_window_callback(PyFunc, data);
    int r = vrna_probs_window($self, ulength, options, &python_wrap_pf_window_cb, (void *)cb);
    release_pf_window_callback(cb);
    return r;
  }
}


/* Also add wrappers for the 'simple' callback interface of pfl_fold_*_cb() functions */
%{

  int
  pfl_fold_cb(std::string sequence,
              int         window_size,
              int         max_bp_span,
              PyObject    *PyFunc,
              PyObject    *data = Py_None)
  {
    python_pf_window_callback_t *cb = bind_pf_window_callback(PyFunc, data);
    int r = vrna_pfl_fold_cb(sequence.c_str(), window_size, max_bp_span, &python_wrap_pf_window_cb, (void *)cb);
    release_pf_window_callback(cb);
    return r;
  }

  int
  pfl_fold_up_cb(std::string  sequence,
                 int          ulength,
                 int          window_size,
                 int          max_bp_span,
                 PyObject     *PyFunc,
                 PyObject     *data = Py_None)
  {
    python_pf_window_callback_t *cb = bind_pf_window_callback(PyFunc, data);
    int r = vrna_pfl_fold_up_cb(sequence.c_str(), ulength, window_size, max_bp_span, &python_wrap_pf_window_cb, (void *)cb);
    release_pf_window_callback(cb);
    return r;
  }

%}

%feature("autodoc") pfl_fold_cb;
%feature("kwargs") pfl_fold_cb;
%feature("autodoc") pfl_fold_up_cb;
%feature("kwargs") pfl_fold_up_cb;

int
pfl_fold_cb(std::string sequence,
            int         window_size,
            int         max_bp_span,
            PyObject    *PyFunc,
            PyObject    *data = Py_None);

int
pfl_fold_up_cb(std::string  sequence,
               int          ulength,
               int          window_size,
               int          max_bp_span,
               PyObject     *PyFunc,
               PyObject     *data = Py_None);

#endif
