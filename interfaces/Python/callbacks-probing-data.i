/****************************************************/
/* BEGIN interface for probing data callback        */
/****************************************************/

#ifdef SWIGPYTHON
%{

#include <stdexcept>

/* TODO: data or cutoff */
typedef struct {
  PyObject *cb;
} python_probing_data_callback_t;

static python_probing_data_callback_t *
bind_probing_data_callback(PyObject *PyFunc);

/* TODO: change wrapper */
static double
python_wrap_probing_data_cb(double     reactivity,
                            void       *data);

static python_probing_data_callback_t *
bind_probing_data_callback(PyObject *PyFunc)
{

  python_probing_data_callback_t *cb = (python_probing_data_callback_t *)vrna_alloc(sizeof(python_probing_data_callback_t));

  Py_XINCREF(PyFunc);
  cb->cb    = PyFunc;  /* store callback */

  return cb;
}

static void
release_probing_data_callback(python_probing_data_callback_t *cb)
{
  Py_DECREF(cb->cb);
  free(cb); 
}


static double
python_wrap_probing_data_cb(double     reactivity,
                            void       *data)
{
  double ret;
  PyObject                  *func, *result, *err;
  python_probing_data_callback_t  *cb;

  cb    = (python_probing_data_callback_t *)data;
  func  = cb->cb;

  /* compose argument list */
  PyObject *py_r;

  py_r   = PyFloat_FromDouble((double)reactivity);
  result = PyObject_CallFunctionObjArgs(func,
                                        py_r,
                                        NULL);

  /* BEGIN recognizing errors in callback execution */
  if (result == NULL) {
    if ((err = PyErr_Occurred())) {
      /* print error message */
      PyErr_Print();
      /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
      if (PyErr_GivenExceptionMatches(err, PyExc_TypeError)) {
        throw std::runtime_error( "Reactivity transform callback must take exactly 1 argument" );
      } else {
        throw std::runtime_error( "Some error occurred while executing reactivity transform callback" );
      }
    }
    PyErr_Clear();
  } else if (PyFloat_Check(result)) {
    ret = (double)PyFloat_AsDouble(result);
  } else {
    throw
    std::runtime_error(
      "Reactivity transform callback must return value in double"
    );
  }
  /* END recognizing errors in callback execution */

  Py_XDECREF(result);

  return ret;
}

%}

%rename(probing_data) vrna_probing_data_s;
%nodefaultctor vrna_probing_data_s;
%nodefaultdtor vrna_probing_data_s;

/* now we overload the existing probing data function */
%extend vrna_probing_data_s {

  /* constructor for Deigan method single sequence */
  vrna_probing_data_s(std::vector<double> reactivities,
                      double              m,
                      double              b,
                      PyObject            *PyFuncOrNone = Py_None)
  {
    vrna_probing_data_s *obj;
    if (PyFuncOrNone != Py_None) {
      python_probing_data_callback_t *cb = bind_probing_data_callback(PyFuncOrNone);
      obj = vrna_probing_data_Deigan2009(&(reactivities[0]),
                                         reactivities.size(),
                                         m,
                                         b,
                                         &python_wrap_probing_data_cb,
                                         (void *)cb
                                        );
    } else {
      obj = vrna_probing_data_Deigan2009(&(reactivities[0]),
                                         reactivities.size(),
                                         m,
                                         b,
                                         NULL,
                                         NULL
                                        );
    }
    return obj;
  }

}


#endif
