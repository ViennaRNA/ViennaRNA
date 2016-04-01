/**********************************************/
/* BEGIN interface for fold compound status   */
/* callback                                   */
/**********************************************/

#ifdef SWIGPYTHON
%{

typedef struct {
  PyObject  *cb;
  PyObject  *data;
  PyObject  *delete_data;
} pycallback_t;

static void       py_wrap_fc_status_callback(unsigned char status, void *data);

static void
delete_pycallback(void * data){

  pycallback_t *cb = (pycallback_t *)data;
  /* first delete user data */
  if(cb->data != Py_None){
    if(cb->delete_data != Py_None){
      PyObject *func, *arglist, *result;
      func = cb->delete_data;
      arglist = Py_BuildValue("O", cb->data);
      result  = PyEval_CallObject(func, arglist);
      Py_DECREF(arglist);
      Py_XDECREF(result);
      Py_XDECREF(cb->delete_data);
    }
  }

  Py_XDECREF(cb->data);

  /* now dispose of the callback */
  Py_XDECREF(cb->cb);
  
  /* finally free pycallback */
  free(cb);
}

static void
fc_add_pycallback(vrna_fold_compound_t *vc,
                  PyObject *PyFunc){

  /* try to dispose of previous callback */
  pycallback_t * cb;
  if(vc->auxdata){
    cb = (pycallback_t *)vc->auxdata;
    /* release previous callback */
    Py_XDECREF(cb->cb);
  } else {
    cb = (pycallback_t *)vrna_alloc(sizeof(pycallback_t));
    cb->data = NULL;
    cb->delete_data = NULL;
  }
  cb->cb = PyFunc;    /* remember callback */
  Py_XINCREF(PyFunc); /* Increase referenc counter */

  /* finaly bind callback wrapper to fold compound */
  vc->auxdata = (void *)cb;
  if(!vc->free_auxdata){
    vc->free_auxdata = &delete_pycallback;
  }
  vrna_fold_compound_add_callback(vc, &py_wrap_fc_status_callback);
}

static void
fc_add_pydata(vrna_fold_compound_t *vc,
              PyObject *data,
              PyObject *PyFunc){

  pycallback_t * cb;
  /* try to dispose of previous data */
  if(vc->auxdata){
    cb = (pycallback_t *)vc->auxdata;
    if(cb->data != Py_None){
      if(cb->delete_data != Py_None){
        PyObject *func, *arglist, *result;
        func    = cb->delete_data;
        arglist = Py_BuildValue("O", cb->data);
        result  = PyEval_CallObject(func, arglist);
        Py_DECREF(arglist);
        Py_XDECREF(result);
      }
    }
    Py_XDECREF(cb->data);
    Py_XDECREF(cb->delete_data);
  } else {
    cb              = (pycallback_t *)vrna_alloc(sizeof(pycallback_t));
    cb->cb          = NULL;
    cb->data        = NULL;
    cb->delete_data = NULL;
  }
  cb->data        = data;   /* remember data */
  cb->delete_data = PyFunc; /* remember delete data function */
  /* increase reference counter */
  Py_XINCREF(data);
  Py_XINCREF(PyFunc);

  vc->auxdata = (void *)cb;
  if(!vc->free_auxdata){
    vc->free_auxdata = &delete_pycallback;
  }
}

static void
py_wrap_fc_status_callback( unsigned char status,
                            void *data){

  PyObject *func, *arglist, *result;
  pycallback_t *cb = (pycallback_t *)data;

  func = cb->cb;
  /* compose argument list */
  arglist = Py_BuildValue("(B,O)", status, (cb->data) ? cb->data : Py_None);
  result =  PyEval_CallObject(func, arglist);
  Py_DECREF(arglist);
  Py_XDECREF(result);
  return /*void*/;
}

%}

static void fc_add_pycallback(vrna_fold_compound_t *vc, PyObject *PyFunc);
static void fc_add_pydata(vrna_fold_compound_t *vc, PyObject *data, PyObject *PyFunc);

%typemap(in) PyObject *PyFunc {
  if (!PyCallable_Check($input)) {
      PyErr_SetString(PyExc_TypeError, "Need a callable object!");
      return NULL;
  }
  $1 = $input;
}

%typemap(in) PyObject * {
  $1 = $input;
}

/* now we bind the above functions as methods to the fold_compound object */
%extend vrna_fold_compound_t {
  void add_auxdata(PyObject *data, PyObject *free_data){
    fc_add_pydata($self, data, free_data);
  }

  void add_callback(PyObject *PyFunc){
    fc_add_pycallback($self, PyFunc);
  }
}

#endif
