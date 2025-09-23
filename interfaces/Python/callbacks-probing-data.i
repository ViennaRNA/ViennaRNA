/****************************************************/
/* BEGIN interface for probing data callback        */
/****************************************************/

%apply PyObject *PyFuncOrNone { PyObject *transform_f, PyObject *data_free_f };


#ifdef SWIGPYTHON
%{


#include <stdexcept>

typedef struct {
  PyObject *cb_trans;
  PyObject *cb_trans_data;
  PyObject *cb_trans_data_free;
} python_probing_data_trans_callback_t;

static python_probing_data_trans_callback_t *
bind_probing_data_trans_callback(PyObject *PyFunc,
                                 PyObject *data,
                                 PyObject *data_free_f);

/* TODO: change wrapper */
static double
python_wrap_probing_data_trans_cb(double     reactivity,
                                  void       *data);

static void
delete_probing_data_trans(python_probing_data_trans_callback_t *cb)
{
  if(cb->cb_trans_data != Py_None){
    if(cb->cb_trans_data_free != Py_None){
      /* call user-defined data destructor */
      PyObject *func, *arglist, *result, *err;
      func = cb->cb_trans_data_free;
      arglist = Py_BuildValue("O", cb->cb_trans_data);
      result  = PyObject_CallObject(func, arglist);

      /* BEGIN recognizing errors in callback execution */
      if (result == NULL) {
        if ((err = PyErr_Occurred())) {
          /* print error message */
          PyErr_Print();
          /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
          if (PyErr_GivenExceptionMatches(err, PyExc_TypeError)) {
            throw std::runtime_error( "Probing data transform data free callback must take exactly 1 argument" );
          } else {
            throw std::runtime_error( "Some error occurred while executing probing data transform data free callback" );
          }
        }
        PyErr_Clear();
      }
      /* END recognizing errors in callback execution */

      Py_DECREF(arglist);
      Py_XDECREF(result);
    }
  }

  Py_DECREF(cb->cb_trans_data);
  Py_DECREF(cb->cb_trans_data_free);
}


static python_probing_data_trans_callback_t *
bind_probing_data_trans_callback(PyObject *transform_f,
                                 PyObject *data,
                                 PyObject *data_free_f)
{
  python_probing_data_trans_callback_t *cb = (python_probing_data_trans_callback_t *)vrna_alloc(sizeof(python_probing_data_trans_callback_t));

  if (transform_f != Py_None) {
    cb->cb_trans  = transform_f;
    Py_XINCREF(transform_f);
  } else {
    cb->cb_trans = NULL;
  }

  cb->cb_trans_data       = data;         /* remember data */
  cb->cb_trans_data_free  = data_free_f; /* remember delete data function */

  /* increase reference counter */
  Py_INCREF(data);
  Py_INCREF(data_free_f);

  return cb;
}

static void
release_probing_data_trans_callback(void *data)
{
  python_probing_data_trans_callback_t *cb = (python_probing_data_trans_callback_t *)data;
  /* first delete user data */
  delete_probing_data_trans(cb);

  /* now dispose of the callback */
  Py_DECREF(cb->cb_trans);
  
  /* finally free pycallback */
  free(cb);
}


static double
python_wrap_probing_data_cb(double     reactivity,
                            void       *data)
{
  double ret;
  PyObject                  *arglist, *func, *result, *err;
  python_probing_data_trans_callback_t  *cb;

  cb    = (python_probing_data_trans_callback_t *)data;
  func  = cb->cb_trans;

  /* compose argument list */
  arglist = Py_BuildValue("(d,O)", reactivity, (cb->cb_trans_data) ? cb->cb_trans_data : Py_None);
  result =  PyObject_CallObject(func, arglist);

  /* BEGIN recognizing errors in callback execution */
  if (result == NULL) {
    if ((err = PyErr_Occurred())) {
      /* print error message */
      PyErr_Print();
      /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
      if (PyErr_GivenExceptionMatches(err, PyExc_TypeError)) {
        throw std::runtime_error( "Reactivity transform_f callback must take exactly 1 argument" );
      } else {
        throw std::runtime_error( "Some error occurred while executing reactivity transform_f callback" );
      }
    }
    PyErr_Clear();
  } else if (PyFloat_Check(result)) {
    ret = (double)PyFloat_AsDouble(result);
  } else {
    throw
    std::runtime_error(
      "Reactivity transform_f callback must return value in double"
    );
  }
  /* END recognizing errors in callback execution */

  Py_DECREF(arglist);
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
                      PyObject            *transform_f = Py_None,
                      PyObject            *transform_data = Py_None,
                      PyObject            *data_free_f = Py_None)
  {
    if (transform_f != Py_None) {
      python_probing_data_trans_callback_t  *cb = bind_probing_data_trans_callback(transform_f,
                                                                                   transform_data,
                                                                                   data_free_f);

      return vrna_probing_data_deigan_trans(&(reactivities[0]),
                                            reactivities.size(),
                                            m,
                                            b,
                                            python_wrap_probing_data_cb,
                                            (void *)cb,
                                            release_probing_data_trans_callback);
    } else {
      return vrna_probing_data_deigan(&(reactivities[0]),
                                      reactivities.size(),
                                      m,
                                      b);
    }
  }


  /* constructor for Zarringhalam method single sequence */
  vrna_probing_data_s(std::vector<double> reactivities,
                      double              beta,
                      std::string         pr_conversion = VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_conversion,
                      double              pr_default    = VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_probability)
  {
    return vrna_probing_data_zarringhalam(&(reactivities[0]),
                                          reactivities.size(),
                                          beta,
                                          pr_conversion.c_str(),
                                          pr_default);
  }


  vrna_probing_data_s(std::vector<double> reactivities,
                      double              beta,
                      double              pr_default    = VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_probability,
                      PyObject            *transform_f = Py_None,
                      PyObject            *transform_data = Py_None,
                      PyObject            *data_free_f = Py_None)
  {
    if (transform_f != Py_None) {
      python_probing_data_trans_callback_t  *cb = bind_probing_data_trans_callback(transform_f,
                                                                                   transform_data,
                                                                                   data_free_f);
      return vrna_probing_data_zarringhalam_trans(&(reactivities[0]),
                                                  reactivities.size(),
                                                  beta,
                                                  pr_default,
                                                  python_wrap_probing_data_cb,
                                                  (void *)cb,
                                                  release_probing_data_trans_callback);
    } else {
      return vrna_probing_data_zarringhalam(&(reactivities[0]),
                                            reactivities.size(),
                                            beta,
                                            VRNA_PROBING_METHOD_ZARRINGHALAM2012_DEFAULT_conversion,
                                            pr_default);
    }
  }


  /* constructor for Eddy method single sequence */
  vrna_probing_data_s(std::vector<double> reactivities,
                      double              temperature,
                      unsigned int        options = VRNA_PROBING_STRATEGY_EDDY_OPTIONS_DEFAULT,
                      std::vector<double> unpaired_data = {},
                      std::vector<double> paired_data = {},
                      PyObject            *transform_f = Py_None,
                      PyObject            *transform_data = Py_None,
                      PyObject            *data_free_f = Py_None)
  {
    if (transform_f != Py_None) {
      python_probing_data_trans_callback_t  *cb = bind_probing_data_trans_callback(transform_f,
                                                                                   transform_data,
                                                                                   data_free_f);

      return vrna_probing_data_eddy_trans(&(reactivities[0]),
                                          reactivities.size(),
                                          temperature,
                                          options,
                                          &(unpaired_data[0]),
                                          unpaired_data.size(),
                                          &(paired_data[0]),
                                          paired_data.size(),
                                          python_wrap_probing_data_cb,
                                          (void*) cb,
                                          release_probing_data_trans_callback);
    } else {
      return vrna_probing_data_eddy(&(reactivities[0]),
                                    reactivities.size(),
                                    temperature,
                                    options,
                                    &(unpaired_data[0]),
                                    unpaired_data.size(),
                                    &(paired_data[0]),
                                    paired_data.size());
    }
  }
}


#endif
