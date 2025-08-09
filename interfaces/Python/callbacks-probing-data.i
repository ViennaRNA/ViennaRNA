/****************************************************/
/* BEGIN interface for probing data callback        */
/****************************************************/

#ifdef SWIGPYTHON
%{

#include <stdexcept>

/* TODO: data or cutoff */
typedef struct {
  double (*trans)(double, void*);
  PyObject *cb;
} python_probing_data_callback_t;

static python_probing_data_callback_t *
bind_probing_data_callback(PyObject *PyFunc);

/* TODO: change wrapper */
static double
python_wrap_probing_data_cb(double     reactivity,
                            void       *data);

static python_probing_data_callback_t *
bind_probing_data_callback(PyObject *PyFuncOrNone)
{
  
  python_probing_data_callback_t *cb = (python_probing_data_callback_t *)vrna_alloc(sizeof(python_probing_data_callback_t));

  if (PyFuncOrNone != Py_None) {
    cb->trans = &python_wrap_probing_data_cb;
    Py_XINCREF(PyFuncOrNone);
    cb->cb    = PyFuncOrNone;
  } else {
    cb->trans = NULL;
    cb->cb = NULL;
  }

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
    python_probing_data_callback_t * cb = bind_probing_data_callback(PyFuncOrNone);
    vrna_probing_data_s *obj = vrna_probing_data_Deigan2009(&(reactivities[0]),
                                                            reactivities.size(),
                                                            m,
                                                            b,
                                                            cb->trans,
                                                            (void*) cb
                                                           );
    return obj;
  }


  /* constructor for Zarringhalam method single sequence */
  vrna_probing_data_s(std::vector<double> reactivities,
                      double              beta,
                      std::string         pr_conversion = VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_conversion,
                      double              pr_default    = VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_probability,
                      PyObject            *PyFuncOrNone = Py_None)
  {
    python_probing_data_callback_t * cb = bind_probing_data_callback(PyFuncOrNone);
    vrna_probing_data_s *obj = vrna_probing_data_Zarringhalam2012(&(reactivities[0]),
                                                                  reactivities.size(),
                                                                  beta,
                                                                  pr_conversion.c_str(),
                                                                  pr_default,
                                                                  cb->trans,
                                                                  (void*) cb
                                                                 );
    return obj;
  }


  /* constructor for Eddy method single sequence */
  vrna_probing_data_s(std::vector<double> reactivities,
                      std::vector<double> unpaired_data,
                      std::vector<double> paired_data,
                      PyObject            *PyFuncOrNone = Py_None)
  {
    python_probing_data_callback_t * cb = bind_probing_data_callback(PyFuncOrNone);
    vrna_probing_data_s *obj = vrna_probing_data_Eddy2014_2(&(reactivities[0]),
                                                            reactivities.size(),
                                                            &(unpaired_data[0]),
                                                            unpaired_data.size(),
                                                            &(paired_data[0]),
                                                            paired_data.size(),
                                                            cb->trans,
                                                            (void*) cb
                                                           );
    return obj;
  }
}


#endif
