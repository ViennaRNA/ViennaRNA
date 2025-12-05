/*
 * ********************************************
 * BEGIN interface for probing callbacks
 *********************************************
 */

#ifdef SWIGPYTHON
%{
#include <stdexcept>

  typedef struct {
    PyObject  *f;
    PyObject  *data;
    PyObject  *delete_data;
  } py_probing_strategy_t;


  double *
  py_wrap_probing_strategy_callback(vrna_fold_compound_t *fc,
                                    const double         *data,
                                    size_t               data_size,
                                    unsigned int         target,
                                    void                 *options);


  static void
  delete_py_probing_strategy_data(py_probing_strategy_t *cb)
  {
    if ((cb->data != Py_None) &&
        (cb->delete_data != Py_None)) {
      PyObject *func, *arglist, *result, *err;
      func    = cb->delete_data;
      arglist = Py_BuildValue("O", cb->data);
      result  = PyObject_CallObject(func, arglist);

      /* BEGIN recognizing errors in callback execution */
      if (result == NULL) {
        if ((err = PyErr_Occurred())) {
          /* print error message */
          PyErr_Print();
          /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
          if (PyErr_GivenExceptionMatches(err, PyExc_TypeError))
            throw
            std::runtime_error(
              "Probing strategy option feee callback must take exactly 1 argument");


          else
            throw
            std::runtime_error(
              "Some error occurred while executing probing strategy option free callback");
        }

        PyErr_Clear();
      }

      /* END recognizing errors in callback execution */

      Py_DECREF(arglist);
      Py_XDECREF(result);
    }

    Py_DECREF(cb->data);
    Py_DECREF(cb->delete_data);
  }


  static void
  delete_py_probing_strategy_callback(void *data)
  {
    py_probing_strategy_t *cb = (py_probing_strategy_t *)data;

    /* first delete user data */
    delete_py_probing_strategy_data(cb);

    /* now dispose of the registered callbacks */
    Py_DECREF(cb->f);

    /* finally free pycallback */
    free(cb);
  }


  vrna_probing_data_t
  probing_data_linear(std::vector<double> data,
                      std::vector<double> data_weights,
                      PyObject            *strategy_cb,
                      PyObject            *strategy_cb_options = Py_None,
                      PyObject            *strategy_cb_options_free = Py_None,
                      unsigned int        options = VRNA_PROBING_DATA_DEFAULT)
  {
    unsigned char     func_is_tuple, func_is_list;
    unsigned int      s;
    py_probing_strategy_t *cb;

    if (!PyCallable_Check(strategy_cb)) {
      PyErr_SetString(PyExc_TypeError, "Need a callable object for argument strategy_cb!");
    } else if ((strategy_cb_options_free != Py_None) && (!PyCallable_Check(strategy_cb_options_free))) {
      PyErr_SetString(PyExc_TypeError, "Need a callable object for argument strategy_cb_options_free!");
    } else {
      /* if we get to this point, at least cb is pointing to a callable object */
      cb = (py_probing_strategy_t *)vrna_alloc(sizeof(py_probing_strategy_t));

      cb->f = strategy_cb;                          /* remember callback */
      Py_INCREF(strategy_cb);                       /* Increase referenc counter */

      cb->delete_data = strategy_cb_options_free;   /* remember callback */
      Py_INCREF(strategy_cb_options_free);          /* Increase reference counter */

      cb->data = strategy_cb_options;               /* remember data */
      Py_INCREF(strategy_cb_options);               /* Increase reference counter */

      return vrna_probing_data_linear(&(data[0]),
                                      data.size() - 1,
                                      &(data_weights[0]),
                                      &py_wrap_probing_strategy_callback,
                                      (void *)cb,
                                      &delete_py_probing_strategy_callback,
                                      options);
    }

    return NULL;
  }


  double *
  py_wrap_probing_strategy_callback(vrna_fold_compound_t *fc,
                                    const double         *data,
                                    size_t               data_size,
                                    unsigned int         target,
                                    void                 *options)
  {
    double                  *ret;
    PyObject                *func, *arglist, *result, *err;
    py_probing_strategy_t   *cb = (py_probing_strategy_t *)options;

    ret   = NULL;
    func  = cb->f;

    /* compose argument list */
    PyObject *py_fc, *py_data, *py_data_size, *py_target;

    py_fc = SWIG_NewPointerObj(SWIG_as_voidptr(fc),
                               SWIGTYPE_p_vrna_fold_compound_t,
                               0);

    /* create PYTHON list for pr values */
    py_data = PyList_New((Py_ssize_t) data_size);

    /* 1-based data, so we add None as first element */
    PyList_SET_ITEM(py_data, (Py_ssize_t) 0, Py_None);

    /* actual values in range */
    for (size_t i = 1; i < data_size; i++)
      PyList_SET_ITEM(py_data, (Py_ssize_t) i, PyFloat_FromDouble(data[i]));

    py_data_size    = PyLong_FromLong(data_size);
    py_target       = PyLong_FromLong(target);
    result  = PyObject_CallFunctionObjArgs(func,
                                           py_fc,
                                           py_data,
                                           py_data_size,
                                           py_target,
                                           (cb->data) ? cb->data : Py_None, NULL);

    Py_DECREF(py_fc);
    Py_DECREF(py_data);
    Py_DECREF(py_data_size);
    Py_DECREF(py_target);

    if (result == Py_None) {
      ret = NULL;
    } else if (result == NULL) {
      /* BEGIN recognizing errors in callback execution */
      if ((err = PyErr_Occurred())) {
        /* print error message */
        PyErr_Print();
        /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
        if (PyErr_GivenExceptionMatches(err, PyExc_TypeError))
          throw
          std::runtime_error(
            "Probing strategy callbacks must take exactly 5 arguments");


        else
          throw
          std::runtime_error(
            "Some error occurred while executing probing strategy callback");
      }

      PyErr_Clear();
    } else if (PyList_Check(result)) {
      ret = (double *)vrna_alloc(sizeof(double) * data_size);
      ret[0] = 0.;
      for (size_t i = 1; i < data_size; i++) {
        PyObject *v = PyList_GetItem(result, (Py_ssize_t)i);
        ret[i] = (PyFloat_Check(v)) ? PyFloat_AsDouble(v) : 0.;
      }
    } else if ((PyTuple_Check(result)) &&
               (PyTuple_Size(result) == (Py_ssize_t)data_size)) {
      ret = (double *)vrna_alloc(sizeof(double) * data_size);
      ret[0] = 0.;
      for (size_t i = 1; i < data_size; i++) {
        PyObject *v = PyTuple_GetItem(result, (Py_ssize_t)i);
        ret[i] = (PyFloat_Check(v)) ? PyFloat_AsDouble(v) : 0.;
      }
    } else {
      throw
      std::runtime_error(
        "Probing strategy callback must return a list or tuple of pseudo energy contributions");
    }

    /* END recognizing errors in callback execution */

    Py_XDECREF(result);

    return ret;
  }


%}

%newobject probing_data_linear;

%feature("autodoc") probing_data_linear;
%feature("kwargs") probing_data_linear;

vrna_probing_data_t
probing_data_linear(std::vector<double>  data,
                    std::vector<double>  data_weights,
                    PyObject             *strategy_cb,
                    PyObject             *strategy_cb_options = Py_None,
                    PyObject             *strategy_cb_options_free = Py_None,
                    unsigned int         options = VRNA_PROBING_DATA_DEFAULT);


#endif
