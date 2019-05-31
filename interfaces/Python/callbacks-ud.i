/**********************************************/
/* BEGIN interface for unstructured domains   */
/* feature callbacks                          */
/**********************************************/

#ifdef SWIGPYTHON
%{

#include <stdexcept>

typedef struct {
  PyObject  *prod_rule;
  PyObject  *exp_prod_rule;
  PyObject  *energy;
  PyObject  *exp_energy;
  PyObject  *data;
  PyObject  *delete_data;
  PyObject  *prob_add;
  PyObject  *prob_get;
} py_ud_callback_t;

static vrna_callback_ud_production      py_wrap_ud_prod_rule;
static vrna_callback_ud_exp_production  py_wrap_ud_exp_prod_rule;
static vrna_callback_ud_energy          py_wrap_ud_energy;
static vrna_callback_ud_exp_energy      py_wrap_ud_exp_energy;
static vrna_callback_ud_probs_add       py_wrap_ud_prob_add;
static vrna_callback_ud_probs_get       py_wrap_ud_prob_get;

static py_ud_callback_t *
new_py_ud_cb(void){

  py_ud_callback_t *cb = (py_ud_callback_t *)vrna_alloc(sizeof(py_ud_callback_t));

  Py_INCREF(Py_None);
  cb->prod_rule     = Py_None;
  Py_INCREF(Py_None);
  cb->exp_prod_rule = Py_None;
  Py_INCREF(Py_None);
  cb->energy        = Py_None;
  Py_INCREF(Py_None);
  cb->exp_energy    = Py_None;
  Py_INCREF(Py_None);
  cb->data          = Py_None;
  Py_INCREF(Py_None);
  cb->delete_data   = Py_None;
  Py_INCREF(Py_None);
  cb->prob_add      = Py_None;
  Py_INCREF(Py_None);
  cb->prob_get      = Py_None;

  return cb;
}


static void
delete_py_ud_data(py_ud_callback_t *cb)
{
  if ((cb->data != Py_None) &&
      (cb->delete_data != Py_None)) {
    PyObject *func, *arglist, *result, *err;
    func = cb->delete_data;
    arglist = Py_BuildValue("O", cb->data);
    result  = PyObject_CallObject(func, arglist);

    /* BEGIN recognizing errors in callback execution */
    if (result == NULL) {
      if ((err = PyErr_Occurred())) {
        /* print error message */
        PyErr_Print();
        /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
        if (PyErr_GivenExceptionMatches(err, PyExc_TypeError)) {
          throw std::runtime_error( "Unstructured domains delete_data() callback must take exactly 1 argument" );
        } else {
          throw std::runtime_error( "Some error occurred while executing unstructured domains delete_data() callback" );
        }
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
delete_py_ud_callback(void * data)
{
  py_ud_callback_t *cb = (py_ud_callback_t *)data;
  /* first delete user data */
  delete_py_ud_data(cb);

  /* now dispose of the registered callbacks */
  Py_DECREF(cb->prod_rule);
  Py_DECREF(cb->exp_prod_rule);
  Py_DECREF(cb->energy);
  Py_DECREF(cb->exp_energy);
  Py_DECREF(cb->prob_add);
  Py_DECREF(cb->prob_get);

  /* finally free pycallback */
  free(cb);
}


static void
ud_set_pydata(vrna_fold_compound_t *vc,
              PyObject             *data,
              PyObject             *PyFunc)
{
  py_ud_callback_t * cb;

  if ((vc->domains_up) &&
      (vc->domains_up->data)) {
    cb = (py_ud_callback_t *)vc->domains_up->data;
  } else {
    cb = new_py_ud_cb();
  }

  /* try to dispose of previous data */
  delete_py_ud_data(cb);

  /* increase reference counter */
  Py_INCREF(data);
  Py_INCREF(PyFunc);

  cb->data        = data;   /* remember data */
  cb->delete_data = PyFunc; /* remember delete data function */

  /* bind callback wrapper to fold compound */
  vrna_ud_set_data(vc, (void *)cb, &delete_py_ud_callback);
}


static void
ud_set_prod_cb(vrna_fold_compound_t *vc,
               PyObject             *prod_cb,
               PyObject             *eval_cb)
{
  /* try to dispose of previous callback */
  py_ud_callback_t * cb;

  if(vc->domains_up && vc->domains_up->data){
    cb = (py_ud_callback_t *)vc->domains_up->data;
  } else {
    cb = new_py_ud_cb();
    /* bind callback wrapper to fold compound */
    vrna_ud_set_data(vc, (void *)cb, &delete_py_ud_callback);
  }

  /* release previous callback */
  Py_DECREF(cb->prod_rule);
  Py_DECREF(cb->energy);

  Py_XINCREF(prod_cb);     /* Increase reference counter */
  Py_XINCREF(eval_cb);     /* Increase reference counter */

  cb->prod_rule = prod_cb; /* remember callback */
  cb->energy    = eval_cb; /* remember callback */

  vrna_ud_set_prod_rule_cb(vc, &py_wrap_ud_prod_rule, &py_wrap_ud_energy);
}


static void
ud_set_exp_prod_cb(vrna_fold_compound_t *vc,
                   PyObject             *prod_cb,
                   PyObject             *eval_cb)
{
  /* try to dispose of previous callback */
  py_ud_callback_t *cb;

  /* now bind the python function to the wrapper structure */
  if ((vc->domains_up) &&
      (vc->domains_up->data)) {
    cb = (py_ud_callback_t *)vc->domains_up->data;
  } else {
    cb = new_py_ud_cb();
    /* bind callback wrapper to fold compound */
    vrna_ud_set_data(vc, (void *)cb, &delete_py_ud_callback);
  }

  /* release previous callback */
  Py_DECREF(cb->exp_prod_rule);
  Py_DECREF(cb->exp_energy);

  Py_INCREF(prod_cb); /* Increase reference counter */
  Py_INCREF(eval_cb); /* Increase reference counter */

  cb->exp_prod_rule = prod_cb;
  cb->exp_energy    = eval_cb;

  vrna_ud_set_exp_prod_rule_cb(vc, &py_wrap_ud_exp_prod_rule, &py_wrap_ud_exp_energy);
}


static void
ud_set_prob_cb(vrna_fold_compound_t *vc,
               PyObject             *setter,
               PyObject             *getter)
{
  py_ud_callback_t *cb;

  /* now bind the python function to the wrapper structure */
  if ((vc->domains_up) &&
      (vc->domains_up->data)) {
    cb = (py_ud_callback_t *)vc->domains_up->data;
  } else {
    cb = new_py_ud_cb();
    /* bind callback wrapper to fold compound */
    vrna_ud_set_data(vc, (void *)cb, &delete_py_ud_callback);
  }

  /* release previous callbacks */
  Py_DECREF(cb->prob_add);
  Py_DECREF(cb->prob_get);

  Py_INCREF(setter); /* Increase reference counter */
  Py_INCREF(getter); /* Increase reference counter */

  cb->prob_add = setter;
  cb->prob_get = getter;

  vrna_ud_set_prob_cb(vc, &py_wrap_ud_prob_add, &py_wrap_ud_prob_get);
}


static void
py_wrap_ud_prod_rule(vrna_fold_compound_t *vc,
                     void                 *data)
{
  PyObject          *func, *arglist, *result, *err;
  py_ud_callback_t  *cb;

  cb    = (py_ud_callback_t *)data;
  func  = cb->prod_rule;

  /* compose argument list */
  arglist = Py_BuildValue("(O,O)", vc, (cb->data) ? cb->data : Py_None);
  result =  PyObject_CallObject(func, arglist);

  /* BEGIN recognizing errors in callback execution */
  if (result == NULL) {
    if ((err = PyErr_Occurred())) {
      /* print error message */
      PyErr_Print();
      /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
      if (PyErr_GivenExceptionMatches(err, PyExc_TypeError)) {
        throw std::runtime_error( "Unstructured domains production rule callback must take exactly 2 arguments" );
      } else {
        throw std::runtime_error( "Some error occurred while executing unstructured domains production rule callback" );
      }
    }
    PyErr_Clear();
  }
  /* END recognizing errors in callback execution */

  Py_DECREF(arglist);
  Py_XDECREF(result);
  return /*void*/;
}


static void
py_wrap_ud_exp_prod_rule(vrna_fold_compound_t *vc,
                         void                 *data)
{
  PyObject          *func, *arglist, *result, *err;
  py_ud_callback_t  *cb;

  cb    = (py_ud_callback_t *)data;
  func  = cb->exp_prod_rule;

  /* compose argument list */
  arglist = Py_BuildValue("(O,O)", vc, (cb->data) ? cb->data : Py_None);
  result =  PyObject_CallObject(func, arglist);

  /* BEGIN recognizing errors in callback execution */
  if (result == NULL) {
    if ((err = PyErr_Occurred())) {
      /* print error message */
      PyErr_Print();
      /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
      if (PyErr_GivenExceptionMatches(err, PyExc_TypeError)) {
        throw std::runtime_error( "Unstructured domains production rule callback (partition function) must take exactly 2 arguments" );
      } else {
        throw std::runtime_error( "Some error occurred while executing unstructured domains production rule callback (partition function)" );
      }
    }
    PyErr_Clear();
  }
  /* END recognizing errors in callback execution */

  Py_DECREF(arglist);
  Py_XDECREF(result);
  return /*void*/;
}


static int
py_wrap_ud_energy(vrna_fold_compound_t *vc,
                  int                  i,
                  int                  j,
                  unsigned int         looptype,
                  void                 *data)
{
  int               ret;
  PyObject          *func, *arglist, *result, *err, *py_vc, *py_i, *py_j, *py_looptype;
  py_ud_callback_t  *cb;

  ret   = 0;
  cb    = (py_ud_callback_t *)data;
  func  = cb->energy;

  /* compose argument list */
  py_vc       = SWIG_NewPointerObj(SWIG_as_voidptr(vc),
                                   SWIGTYPE_p_vrna_fold_compound_t,
                                   SWIG_POINTER_NEW);
  py_i        = PyInt_FromLong(i);
  py_j        = PyInt_FromLong(j);
  py_looptype = PyInt_FromLong(looptype);
  result      = PyObject_CallFunctionObjArgs(func,
                                             py_vc,
                                             py_i,
                                             py_j,
                                             py_looptype,
                                             (cb->data) ? cb->data : Py_None,
                                             NULL);

  Py_DECREF(py_vc);
  Py_DECREF(py_i);
  Py_DECREF(py_j);
  Py_DECREF(py_looptype);

  /* BEGIN recognizing errors in callback execution */
  if (result == NULL) {
    if ((err = PyErr_Occurred())) {
      /* print error message */
      PyErr_Print();
      /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
      if (PyErr_GivenExceptionMatches(err, PyExc_TypeError)) {
        throw std::runtime_error( "Unstructured domains energy callback must take exactly 5 arguments" );
      } else {
        throw std::runtime_error( "Some error occurred while executing unstructured domains energy callback" );
      }
    }
    PyErr_Clear();
  } else if (result == Py_None) {
    throw std::runtime_error( "Unstructured domains energy callback must return pseudo energy value" );
  } else {
    ret = (int)PyInt_AsLong(result);
  }
  /* END recognizing errors in callback execution */

  Py_XDECREF(result);
  return ret;
}


static FLT_OR_DBL
py_wrap_ud_exp_energy(vrna_fold_compound_t *vc,
                      int                  i,
                      int                  j,
                      unsigned int         looptype,
                      void                 *data)
{
  FLT_OR_DBL        ret;
  PyObject          *func, *arglist, *result, *err, *py_vc, *py_i, *py_j, *py_looptype;
  py_ud_callback_t  *cb;

  ret   = 1.;
  cb    = (py_ud_callback_t *)data;
  func  = cb->exp_energy;

  /* compose argument list */
  py_vc       = SWIG_NewPointerObj(SWIG_as_voidptr(vc),
                                   SWIGTYPE_p_vrna_fold_compound_t,
                                   SWIG_POINTER_NEW);
  py_i        = PyInt_FromLong(i);
  py_j        = PyInt_FromLong(j);
  py_looptype = PyInt_FromLong(looptype);
  result      = PyObject_CallFunctionObjArgs(func,
                                             py_vc,
                                             py_i,
                                             py_j,
                                             py_looptype,
                                             (cb->data) ? cb->data : Py_None,
                                             NULL);

  Py_DECREF(py_vc);
  Py_DECREF(py_i);
  Py_DECREF(py_j);
  Py_DECREF(py_looptype);

  /* BEGIN recognizing errors in callback execution */
  if (result == NULL) {
    if ((err = PyErr_Occurred())) {
      /* print error message */
      PyErr_Print();
      /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
      if (PyErr_GivenExceptionMatches(err, PyExc_TypeError)) {
        throw std::runtime_error( "Unstructured domains energy callback (partition function) must take exactly 5 arguments" );
      } else {
        throw std::runtime_error( "Some error occurred while executing unstructured domains energy callback (partition function)" );
      }
    }
    PyErr_Clear();
  } else if (result == Py_None) {
    throw std::runtime_error( "Unstructured domains energy callback (partition function) must return Boltzmann weighted pseudo energy value" );
  } else {
    ret = (FLT_OR_DBL)PyFloat_AsDouble(result);
  }
  /* END recognizing errors in callback execution */

  Py_XDECREF(result);
  return ret;
}


static void
py_wrap_ud_prob_add(vrna_fold_compound_t *vc,
                    int                  i,
                    int                  j,
                    unsigned int         looptype,
                    FLT_OR_DBL           prob,
                    void                 *data)
{
  PyObject          *func, *arglist, *result, *err, *py_vc, *py_i, *py_j, *py_looptype, *py_prob;
  py_ud_callback_t  *cb;

  cb    = (py_ud_callback_t *)data;
  func  = cb->prob_add;

  /* compose argument list */
  py_vc       = SWIG_NewPointerObj(SWIG_as_voidptr(vc),
                                   SWIGTYPE_p_vrna_fold_compound_t,
                                   SWIG_POINTER_NEW);
  py_i        = PyInt_FromLong(i);
  py_j        = PyInt_FromLong(j);
  py_looptype = PyInt_FromLong(looptype);
  py_prob     = PyFloat_FromDouble((double)prob);
  result      = PyObject_CallFunctionObjArgs(func,
                                             py_vc,
                                             py_i,
                                             py_j,
                                             py_looptype,
                                             py_prob,
                                             (cb->data) ? cb->data : Py_None,
                                             NULL);

  Py_DECREF(py_vc);
  Py_DECREF(py_i);
  Py_DECREF(py_j);
  Py_DECREF(py_looptype);
  Py_DECREF(py_prob);

  /* BEGIN recognizing errors in callback execution */
  if (result == NULL) {
    if ((err = PyErr_Occurred())) {
      /* print error message */
      PyErr_Print();
      /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
      if (PyErr_GivenExceptionMatches(err, PyExc_TypeError)) {
        throw std::runtime_error( "Unstructured domains add_probability() callback must take exactly 6 arguments" );
      } else {
        throw std::runtime_error( "Some error occurred while executing unstructured domains add_probability() callback" );
      }
    }
    PyErr_Clear();
  }
  /* END recognizing errors in callback execution */

  Py_XDECREF(result);
  return;
}


static FLT_OR_DBL
py_wrap_ud_prob_get(vrna_fold_compound_t *vc,
                    int                  i,
                    int                  j,
                    unsigned int         looptype,
                    int                  motif,
                    void                 *data)
{
  FLT_OR_DBL        ret;
  PyObject          *func, *arglist, *result, *err, *py_vc, *py_i, *py_j, *py_looptype, *py_motif;
  py_ud_callback_t  *cb;

  ret   = 1.;
  cb    = (py_ud_callback_t *)data;
  func  = cb->prob_get;

  /* compose argument list */
  py_vc       = SWIG_NewPointerObj(SWIG_as_voidptr(vc),
                                   SWIGTYPE_p_vrna_fold_compound_t,
                                   SWIG_POINTER_NEW);
  py_i        = PyInt_FromLong(i);
  py_j        = PyInt_FromLong(j);
  py_looptype = PyInt_FromLong(looptype);
  py_motif    = PyInt_FromLong(motif);
  result      = PyObject_CallFunctionObjArgs(func,
                                             py_vc,
                                             py_i,
                                             py_j,
                                             py_looptype,
                                             py_motif,
                                             (cb->data) ? cb->data : Py_None,
                                             NULL);

  Py_DECREF(py_vc);
  Py_DECREF(py_i);
  Py_DECREF(py_j);
  Py_DECREF(py_looptype);
  Py_DECREF(py_motif);

  /* BEGIN recognizing errors in callback execution */
  if (result == NULL) {
    if ((err = PyErr_Occurred())) {
      /* print error message */
      PyErr_Print();
      /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
      if (PyErr_GivenExceptionMatches(err, PyExc_TypeError)) {
        throw std::runtime_error( "Unstructured domains get_probability() callback must take exactly 6 arguments" );
      } else {
        throw std::runtime_error( "Some error occurred while executing unstructured domains get_probability() callback" );
      }
    }
    PyErr_Clear();
  } else if (result == Py_None) {
    throw std::runtime_error( "Unstructured domains get_probability() callback must return probability" );
  } else {
    ret = (FLT_OR_DBL)PyFloat_AsDouble(result);
  }
  /* END recognizing errors in callback execution */

  Py_XDECREF(result);
  return ret;
}

%}

static void
ud_set_pydata(vrna_fold_compound_t *vc,
              PyObject             *data,
              PyObject             *PyFuncOrNone);

static void
ud_set_prod_cb(vrna_fold_compound_t *vc,
               PyObject             *prod_cb,
               PyObject             *eval_cb);

static void
ud_set_exp_prod_cb(vrna_fold_compound_t *vc,
                   PyObject             *prod_cb,
                   PyObject             *eval_cb);

static void
ud_set_prob_cb(vrna_fold_compound_t *vc,
               PyObject             *setter,
               PyObject             *getter);

/* now we bind the above functions as methods to the fold_compound object */
%extend vrna_fold_compound_t {

%feature("autodoc") ud_set_data;
%feature("kwargs") ud_set_data;
%feature("autodoc") ud_set_prod_rule_cb;
%feature("kwargs") ud_set_prod_rule_cb;
%feature("autodoc") ud_set_exp_prod_rule_cb;
%feature("kwargs") ud_set_exp_prod_rule_cb;
%feature("autodoc") ud_set_prob_cb;
%feature("kwargs") ud_set_prob_cb;

  PyObject *
  ud_set_data(PyObject *data,
              PyObject *PyFuncOrNone = Py_None)
  {
    ud_set_pydata($self, data, PyFuncOrNone);
    Py_RETURN_NONE;
  }

  PyObject *
  ud_set_prod_rule_cb(PyObject *prod_cb,
                      PyObject *eval_cb)
  {
    if (!PyCallable_Check(prod_cb)) {
      PyErr_SetString(PyExc_TypeError, "Need a callable object!");
      Py_RETURN_NONE;
    } else if (!PyCallable_Check(eval_cb)) {
      PyErr_SetString(PyExc_TypeError, "Need a callable object!");
      Py_RETURN_NONE;
    }
    ud_set_prod_cb($self, prod_cb, eval_cb);
    Py_RETURN_NONE;
  }

  PyObject *
  ud_set_exp_prod_rule_cb(PyObject *prod_cb,
                          PyObject *eval_cb)
  {
    if (!PyCallable_Check(prod_cb)) {
      PyErr_SetString(PyExc_TypeError, "Need a callable object!");
      Py_RETURN_NONE;
    } else if (!PyCallable_Check(eval_cb)) {
      PyErr_SetString(PyExc_TypeError, "Need a callable object!");
      Py_RETURN_NONE;
    }
    ud_set_exp_prod_cb($self, prod_cb, eval_cb);
    Py_RETURN_NONE;
  }

  PyObject *
  ud_set_prob_cb(PyObject *setter,
                 PyObject *getter)
  {
    if (!PyCallable_Check(setter)) {
      PyErr_SetString(PyExc_TypeError, "Need a callable object!");
      Py_RETURN_NONE;
    } else if (!PyCallable_Check(getter)) {
      PyErr_SetString(PyExc_TypeError, "Need a callable object!");
      Py_RETURN_NONE;
    }
    ud_set_prob_cb($self, setter, getter);
    Py_RETURN_NONE;
  }
}

#endif
