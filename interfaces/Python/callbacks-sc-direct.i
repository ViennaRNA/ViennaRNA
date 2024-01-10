/*
 * ********************************************
 * BEGIN interface for direct soft constraint
 * callback implementation
 *********************************************
 */

#ifdef SWIGPYTHON
%{
#include <stdexcept>

  typedef struct {
    PyObject  *f;
    PyObject  *exp_f;
    PyObject  *data;
    PyObject  *prepare_data;
    PyObject  *delete_data;
  } py_sc_cb_direct_t;


  static int
  py_wrap_sc_direct_f_callback(vrna_fold_compound_t *fc,
                               int                  i,
                               int                  j,
                               int                  k,
                               int                  l,
                               void                 *data);


  static FLT_OR_DBL
  py_wrap_sc_direct_exp_f_callback(vrna_fold_compound_t *fc,
                                   int                  i,
                                   int                  j,
                                   int                  k,
                                   int                  l,
                                   void                 *data);


  static int
  py_wrap_sc_direct_data_prepare_callback(vrna_fold_compound_t  *fc,
                                          void                  *data,
                                          unsigned int          event,
                                          void                  *event_data);


  static void
  delete_py_sc_direct_data(py_sc_cb_direct_t *cb)
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
              "Generic soft constraint delete_data() callback must take exactly 1 argument");


          else
            throw
            std::runtime_error(
              "Some error occurred while executing generic soft constraint delete_data() callback");
        }

        PyErr_Clear();
      }

      /* END recognizing errors in callback execution */

      Py_DECREF(arglist);
      Py_XDECREF(result);
    }

    Py_DECREF(cb->data);
    Py_DECREF(cb->prepare_data);
    Py_DECREF(cb->delete_data);
  }


  static void
  delete_py_sc_direct_callback(void *data)
  {
    py_sc_cb_direct_t *cb = (py_sc_cb_direct_t *)data;

    /* first delete user data */
    delete_py_sc_direct_data(cb);

    /* now dispose of the registered callbacks */
    Py_DECREF(cb->f);
    Py_DECREF(cb->exp_f);

    /* finally free pycallback */
    free(cb);
  }


  static int
  sc_multi_cb_add_pycallback(vrna_fold_compound_t *fc,
                             PyObject             *f,
                             PyObject             *exp_f,
                             PyObject             *data,
                             PyObject             *data_prepare,
                             PyObject             *data_free,
                             unsigned int         decomp_type)
  {
    unsigned char     func_is_tuple, func_is_list;
    unsigned int      s;
    py_sc_cb_direct_t *cb;

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (!PyCallable_Check(f)) {
          PyErr_SetString(PyExc_TypeError, "Need a callable object for argument f!");
        } else if ((exp_f != Py_None) && (!PyCallable_Check(exp_f))) {
          PyErr_SetString(PyExc_TypeError, "Need a callable object for argument exp_f!");
        } else if ((data_prepare != Py_None) && (!PyCallable_Check(data_prepare))) {
          PyErr_SetString(PyExc_TypeError, "Need a callable object for argument data_prepare!");
        } else if ((data_free != Py_None) && (!PyCallable_Check(data_free))) {
          PyErr_SetString(PyExc_TypeError, "Need a callable object for argument data_free!");
        } else {
          /* if we get to this point, at least cb is pointing to a callable object */
          cb = (py_sc_cb_direct_t *)vrna_alloc(sizeof(py_sc_cb_direct_t));

          cb->f = f;                        /* remember callback */
          Py_INCREF(f);                     /* Increase referenc counter */

          cb->exp_f = exp_f;                /* remember callback */
          Py_INCREF(exp_f);                 /* Increase reference counter */

          cb->prepare_data = data_prepare;  /* remember callback */
          Py_INCREF(data_prepare);          /* Increase reference counter */

          cb->delete_data = data_free;      /* remember callback */
          Py_INCREF(data_free);             /* Increase reference counter */

          cb->data = data;                  /* remember data */
          Py_INCREF(data);                  /* Increase reference counter */

          return vrna_sc_multi_cb_add(fc,
                                      &py_wrap_sc_direct_f_callback,
                                      &py_wrap_sc_direct_exp_f_callback,
                                      (void *)cb,
                                      &py_wrap_sc_direct_data_prepare_callback,
                                      &delete_py_sc_direct_callback,
                                      decomp_type);
        }

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        throw
        std::runtime_error("sc_multi_cb_add(): Not implemented for comparative fold compound yet!");


        break;
    }

    return 0;
  }


  static int
  py_wrap_sc_direct_f_callback(vrna_fold_compound_t *fc,
                               int                  i,
                               int                  j,
                               int                  k,
                               int                  l,
                               void                 *data)
  {
    int               ret;
    PyObject          *func, *arglist, *result, *err;
    py_sc_callback_t  *cb = (py_sc_callback_t *)data;

    ret   = 0;
    func  = cb->cb_f;

    /* compose argument list */
    PyObject *py_fc, *py_i, *py_j, *py_k, *py_l;

    py_fc = SWIG_NewPointerObj(SWIG_as_voidptr(fc),
                               SWIGTYPE_p_vrna_fold_compound_t,
                               SWIG_POINTER_NEW);
    py_i    = PyLong_FromLong(i);
    py_j    = PyLong_FromLong(j);
    py_k    = PyLong_FromLong(k);
    py_l    = PyLong_FromLong(l);
    result  = PyObject_CallFunctionObjArgs(func,
                                           py_fc,
                                           py_i,
                                           py_j,
                                           py_k,
                                           py_l,
                                           (cb->data) ? cb->data : Py_None,
                                           NULL);

    Py_DECREF(py_fc);
    Py_DECREF(py_i);
    Py_DECREF(py_j);
    Py_DECREF(py_k);
    Py_DECREF(py_l);

    /* BEGIN recognizing errors in callback execution */
    if (result == NULL) {
      if ((err = PyErr_Occurred())) {
        /* print error message */
        PyErr_Print();
        /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
        if (PyErr_GivenExceptionMatches(err, PyExc_TypeError))
          throw
          std::runtime_error(
            "Generic direct soft constraint callbacks must take exactly 6 arguments");


        else
          throw
          std::runtime_error(
            "Some error occurred while executing generic direct soft constraint callback");
      }

      PyErr_Clear();
    } else if (PyLong_Check(result)) {
      ret = (int)PyLong_AsLong(result);
    } else {
      throw
      std::runtime_error(
        "Generic direct soft constraint callback must return pseudo energy value in 10 cal/mol");
    }

    /* END recognizing errors in callback execution */

    Py_XDECREF(result);

    return ret;
  }


  static FLT_OR_DBL
  py_wrap_sc_direct_exp_f_callback(vrna_fold_compound_t *fc,
                                   int                  i,
                                   int                  j,
                                   int                  k,
                                   int                  l,
                                   void                 *data)
  {
    FLT_OR_DBL        ret;
    PyObject          *func, *arglist, *result, *err;
    py_sc_callback_t  *cb;

    cb    = (py_sc_callback_t *)data;
    ret   = 1.;
    func  = cb->cb_exp_f;
    /* compose argument list */
    PyObject *py_fc, *py_i, *py_j, *py_k, *py_l;

    py_fc = SWIG_NewPointerObj(SWIG_as_voidptr(fc),
                               SWIGTYPE_p_vrna_fold_compound_t,
                               SWIG_POINTER_NEW);
    py_i  = PyLong_FromLong(i);
    py_j  = PyLong_FromLong(j);
    py_k  = PyLong_FromLong(k);
    py_l  = PyLong_FromLong(l);

    result = PyObject_CallFunctionObjArgs(func,
                                          py_fc,
                                          py_i,
                                          py_j,
                                          py_k,
                                          py_l,
                                          (cb->data) ? cb->data : Py_None,
                                          NULL);
    Py_DECREF(py_fc);
    Py_DECREF(py_i);
    Py_DECREF(py_j);
    Py_DECREF(py_k);
    Py_DECREF(py_l);

    /* BEGIN recognizing errors in callback execution */
    if (result == NULL) {
      if ((err = PyErr_Occurred())) {
        /* print error message */
        PyErr_Print();
        /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
        if (PyErr_GivenExceptionMatches(err, PyExc_TypeError))
          throw
          std::runtime_error(
            "Generic direct soft constraint callbacks (partition function) must take exactly 6 arguments");


        else
          throw
          std::runtime_error(
            "Some error occurred while executing generic direct soft constraint callback (partition function)");
      }

      PyErr_Clear();
    } else if (result == Py_None) {
      throw
      std::runtime_error(
        "Generic direct soft constraint callback (partition function) must return Boltzmann weighted pseudo energy value");
    } else {
      ret = (FLT_OR_DBL)PyFloat_AsDouble(result);
    }

    /* END recognizing errors in callback execution */

    Py_XDECREF(result);
    return ret;
  }


  static int
  py_wrap_sc_direct_data_prepare_callback(vrna_fold_compound_t  *fc,
                                          void                  *data,
                                          unsigned int          event,
                                          void                  *event_data)
  {
    int               ret;
    PyObject          *func, *arglist, *result, *err;
    py_sc_cb_direct_t *cb = (py_sc_cb_direct_t *)data;

    ret   = 0;
    func  = cb->prepare_data;

    /* compose argument list */
    PyObject *py_fc, *py_event, *py_event_data;

    py_fc = SWIG_NewPointerObj(SWIG_as_voidptr(fc),
                               SWIGTYPE_p_vrna_fold_compound_t,
                               SWIG_POINTER_NEW);
    py_event = PyLong_FromLong(event);
    if (event_data) {
      py_event_data = PyLong_FromLong(*((unsigned int *)event_data));
    } else {
      py_event_data = Py_None;
      Py_INCREF(Py_None);
    }

    result = PyObject_CallFunctionObjArgs(func,
                                          py_fc,
                                          (cb->data) ? cb->data : Py_None,
                                          py_event,
                                          py_event_data,
                                          NULL);

    Py_DECREF(py_fc);
    Py_DECREF(py_event);
    Py_DECREF(py_event_data);

    /* BEGIN recognizing errors in callback execution */
    if (result == NULL) {
      if ((err = PyErr_Occurred())) {
        /* print error message */
        PyErr_Print();
        /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
        if (PyErr_GivenExceptionMatches(err, PyExc_TypeError))
          throw
          std::runtime_error(
            "Generic direct soft constraint data prepapre callback must take exactly 4 arguments");


        else
          throw
          std::runtime_error(
            "Some error occurred while executing generic direct soft constraint data prepare callback");
      }

      PyErr_Clear();
    } else if (PyLong_Check(result)) {
      ret = (int)PyLong_AsLong(result);
    } else {
      throw
      std::runtime_error(
        "Generic direct soft constraint data prepare callback must return integer value");
    }

    /* END recognizing errors in callback execution */

    Py_XDECREF(result);

    return ret;
  }


%}


static unsigned int
sc_multi_cb_add_pycallback(vrna_fold_compound_t *fc,
                           PyObject             *f,
                           PyObject             *f_exp,
                           PyObject             *data,
                           PyObject             *data_prepare,
                           PyObject             *data_free,
                           unsigned int         decomp_type);


/* now we bind the above functions as methods to the fold_compound object */
%extend vrna_fold_compound_t {
  %feature("autodoc") sc_multi_cb_add;
  %feature("kwargs") sc_multi_cb_add;

  unsigned int
  sc_multi_cb_add(PyObject      *f,
                  PyObject      *f_exp           = Py_None,
                  PyObject      *data             = Py_None,
                  PyObject      *data_prepare  = Py_None,
                  PyObject      *data_free     = Py_None,
                  unsigned int  decomp_type       = 0)
  {
    return sc_multi_cb_add_pycallback($self,
                                      f,
                                      f_exp,
                                      data,
                                      data_prepare,
                                      data_free,
                                      decomp_type);
  }
}

#endif
