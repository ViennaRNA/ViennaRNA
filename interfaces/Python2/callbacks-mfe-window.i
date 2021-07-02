/**********************************************/
/* BEGIN interface for MFE window callback    */
/**********************************************/

#ifdef SWIGPYTHON
%{

#include <stdexcept>

typedef struct {
  PyObject *cb;
  PyObject *data;
} python_mfe_window_callback_t;

static python_mfe_window_callback_t *
bind_mfe_window_callback(PyObject *PyFunc,
                         PyObject *data);


static void
python_wrap_mfe_window_cb(int         start,
                          int         end,
                          const char  *structure,
                          float       energy,
                          void        *data);

#ifdef VRNA_WITH_SVM
static void
python_wrap_mfe_window_zscore_cb(int        start,
                                 int        end,
                                 const char *structure,
                                 float      energy,
                                 float      zscore,
                                 void       *data);
#endif

static python_mfe_window_callback_t *
bind_mfe_window_callback(PyObject *PyFunc,
                         PyObject *data)
{
  python_mfe_window_callback_t *cb = (python_mfe_window_callback_t *)vrna_alloc(sizeof(python_mfe_window_callback_t));

  Py_INCREF(PyFunc);
  Py_INCREF(data);

  cb->cb    = PyFunc;  /* store callback */
  cb->data  = data;    /* bind data */

  return cb;
}


static void
release_mfe_window_callback(python_mfe_window_callback_t *cb)
{
  Py_DECREF(cb->cb);
  Py_DECREF(cb->data);
  free(cb); 
}


static void
python_wrap_mfe_window_cb(int         start,
                          int         end,
                          const char  *structure,
                          float       energy,
                          void        *data)
{
  PyObject *func, *arglist, *result, *err;
  python_mfe_window_callback_t *cb = (python_mfe_window_callback_t *)data;

  func = cb->cb;
  /* compose argument list */
  PyObject *py_start, *py_end, *py_structure, *py_energy;
  py_start      = PyInt_FromLong(start);
  py_end        = PyInt_FromLong(end);
  py_structure  = PyString_FromString(structure);
  py_energy     = PyFloat_FromDouble((double)energy);
  result        = PyObject_CallFunctionObjArgs(func,
                                               py_start,
                                               py_end,
                                               py_structure,
                                               py_energy,
                                               (cb->data) ? cb->data : Py_None,
                                               NULL);

  Py_DECREF(py_start);
  Py_DECREF(py_end);
  Py_DECREF(py_structure);
  Py_DECREF(py_energy);

  /* BEGIN recognizing errors in callback execution */
  if (result == NULL) {
    if ((err = PyErr_Occurred())) {
      /* print error message */
      PyErr_Print();
      /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
      if (PyErr_GivenExceptionMatches(err, PyExc_TypeError)) {
        throw std::runtime_error( "Sliding window MFE callback must take exactly 5 arguments" );
      } else {
        throw std::runtime_error( "Some error occurred while executing sliding window MFE callback" );
      }
    }
    PyErr_Clear();
  }
  /* END recognizing errors in callback execution */

  Py_XDECREF(result);

  return /*void*/;
}

#ifdef VRNA_WITH_SVM
static void
python_wrap_mfe_window_zscore_cb(int        start,
                                 int        end,
                                 const char *structure,
                                 float      energy,
                                 float      zscore,
                                 void       *data)
{
  PyObject *func, *arglist, *result, *err;
  python_mfe_window_callback_t *cb = (python_mfe_window_callback_t *)data;

  func = cb->cb;
  /* compose argument list */
  PyObject *py_start, *py_end, *py_structure, *py_energy, *py_zscore;
  py_start      = PyInt_FromLong(start);
  py_end        = PyInt_FromLong(end);
  py_structure  = PyString_FromString(structure);
  py_energy     = PyFloat_FromDouble((double)energy);
  py_zscore     = PyFloat_FromDouble((double)zscore);
  result        = PyObject_CallFunctionObjArgs(func,
                                               py_start,
                                               py_end,
                                               py_structure,
                                               py_energy,
                                               py_zscore,
                                               (cb->data) ? cb->data : Py_None,
                                               NULL);

  Py_DECREF(py_start);
  Py_DECREF(py_end);
  Py_DECREF(py_structure);
  Py_DECREF(py_energy);
  Py_DECREF(py_zscore);

  /* BEGIN recognizing errors in callback execution */
  if (result == NULL) {
    if ((err = PyErr_Occurred())) {
      /* print error message */
      PyErr_Print();
      /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
      if (PyErr_GivenExceptionMatches(err, PyExc_TypeError)) {
        throw std::runtime_error( "Sliding window MFE callback (z-score) must take exactly 6 arguments" );
      } else {
        throw std::runtime_error( "Some error occurred while executing sliding window MFE callback (z-score)" );
      }
    }
    PyErr_Clear();
  }
  /* END recognizing errors in callback execution */

  Py_XDECREF(result);

  return /*void*/;
}
#endif

%}

/* now we bind the above functions as methods to the fold_compound object */
%extend vrna_fold_compound_t {

%feature("autodoc") mfe_window_cb;
%feature("kwargs") mfe_window_cb;

  float
  mfe_window_cb(PyObject *PyFunc,
                PyObject *data = Py_None)
  {
    float en;
    python_mfe_window_callback_t *cb = bind_mfe_window_callback(PyFunc, data);
    en = vrna_mfe_window_cb($self, &python_wrap_mfe_window_cb, (void *)cb);
    release_mfe_window_callback(cb);
    return en;
  }

#ifdef VRNA_WITH_SVM
%feature("autodoc") mfe_window_score_cb;
%feature("kwargs") mfe_window_score_cb;

  float
  mfe_window_score_cb(double   min_z,
                      PyObject *PyFunc,
                      PyObject *data = Py_None)
  {
    float en;
    python_mfe_window_callback_t *cb = bind_mfe_window_callback(PyFunc, data);
    en = vrna_mfe_window_zscore_cb($self, min_z, &python_wrap_mfe_window_zscore_cb, (void *)cb);
    release_mfe_window_callback(cb);
    return en;
  }
#endif

}

%rename (Lfold_cb) my_Lfold_cb;
%rename (Lfoldz_cb) my_Lfoldz_cb;
%rename (aliLfold_cb) my_aliLfold_cb;

%{
  float
  my_Lfold_cb(char      *string,
              int       window_size,
              PyObject  *PyFunc,
              PyObject  *data = Py_None)
  {
    float en;
    python_mfe_window_callback_t *cb = bind_mfe_window_callback(PyFunc, data);
    en = vrna_Lfold_cb(string, window_size, &python_wrap_mfe_window_cb, (void *)cb);
    release_mfe_window_callback(cb);
    return en;
  }

#ifdef VRNA_WITH_SVM
  float
  my_Lfoldz_cb(char *   string,
               int      window_size,
               double   min_z,
               PyObject *PyFunc,
               PyObject *data = Py_None)
  {
    float en;
    python_mfe_window_callback_t *cb = bind_mfe_window_callback(PyFunc, data);
    en = vrna_Lfoldz_cb(string, window_size, min_z, &python_wrap_mfe_window_zscore_cb, (void *)cb);
    release_mfe_window_callback(cb);
    return en;
  }
#endif

  float
  my_aliLfold_cb(std::vector<std::string> alignment,
                 int                      window_size,
                 PyObject                 *PyFunc,
                 PyObject                 *data = Py_None)
  {
    float en;
    python_mfe_window_callback_t *cb = bind_mfe_window_callback(PyFunc, data);
    std::vector<const char*>  vc;
    std::transform(alignment.begin(), alignment.end(), std::back_inserter(vc), convert_vecstring2veccharcp);
    vc.push_back(NULL); /* mark end of sequences */
    en = vrna_aliLfold_cb((const char **)&vc[0], window_size, &python_wrap_mfe_window_cb, (void *)cb);
    release_mfe_window_callback(cb);
    return en;
  }

%}

%feature("autodoc") my_Lfold_cb;
%feature("kwargs") my_Lfold_cb;

float
my_Lfold_cb(char      *string,
            int       window_size,
            PyObject  *PyFunc,
            PyObject  *data);

%ignore vrna_Lfold_cb;

#ifdef VRNA_WITH_SVM
%feature("autodoc") my_Lfoldz_cb;
%feature("kwargs") my_Lfoldz_cb;

float
my_Lfoldz_cb(char     *string,
             int      window_size,
             double   min_z,
             PyObject *PyFunc,
             PyObject *data);

%ignore vrna_Lfoldz_cb;
#endif

%feature("autodoc") my_aliLfold_cb;

float
my_aliLfold_cb(std::vector<std::string> alignment,
               int                      window_size,
               PyObject                 *PyFunc,
               PyObject                 *data);

%ignore vrna_aliLfold_cb;
%ignore aliLfold_cb;

#endif
