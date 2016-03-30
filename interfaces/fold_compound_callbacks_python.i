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

typedef struct {
  PyObject  *cb_f;
  PyObject  *cb_bt;
  PyObject  *cb_exp_f;
  PyObject  *data;
  PyObject  *delete_data;
} py_sc_callback_t;

static void       py_wrap_fc_status_callback(unsigned char status, void *data);
static int        py_wrap_sc_f_callback(int i, int j, int k, int l, char d, void *data);
static vrna_basepair_t *py_wrap_sc_bt_callback( int i, int j, int k, int l, char d,                        void *data);
static FLT_OR_DBL py_wrap_sc_exp_f_callback(int i, int j, int k, int l, char d, void *data);

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
delete_py_sc_callback(void * data){

  py_sc_callback_t *cb = (py_sc_callback_t *)data;
  /* first delete user data */
  if(cb->delete_data != Py_None){
    PyObject *func, *arglist, *result;
    func = cb->delete_data;
    arglist = Py_BuildValue("O", cb->data);
    result  = PyEval_CallObject(func, arglist);
    Py_DECREF(arglist);
    Py_XDECREF(result);
  }

  /* now dispose of the registered callbacks */
  Py_XDECREF(cb->cb_f);
  Py_XDECREF(cb->cb_bt);
  Py_XDECREF(cb->cb_exp_f);

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
sc_add_f_pycallback(vrna_fold_compound_t *vc,
                    PyObject *PyFunc){

  /* try to dispose of previous callback */
  py_sc_callback_t * cb;
  vrna_sc_add_f(vc, &py_wrap_sc_f_callback); /* this also creates the soft constraint data structure inside vc */

  /* now bind the python function to the wrapper structure */
  if(vc->sc->data){
    cb = (py_sc_callback_t *)vc->sc->data;
    /* release previous callback */
    Py_XDECREF(cb->cb_f);
  } else {
    cb = (py_sc_callback_t *)vrna_alloc(sizeof(py_sc_callback_t));
    cb->cb_f        = NULL;
    cb->cb_bt       = NULL;
    cb->cb_exp_f    = NULL;
    cb->data        = NULL;
    cb->delete_data = NULL;
  }
  cb->cb_f = PyFunc;  /* remember callback */
  Py_XINCREF(PyFunc); /* Increase referenc counter */

  /* finaly bind callback wrapper to fold compound */
  vc->sc->data = (void *)cb;
  if(!vc->sc->free_data){
    vc->sc->free_data = &delete_py_sc_callback;
  }
}

static void
sc_add_exp_f_pycallback(vrna_fold_compound_t *vc,
                        PyObject *PyFunc){

  /* try to dispose of previous callback */
  py_sc_callback_t * cb;
  vrna_sc_add_exp_f(vc, &py_wrap_sc_exp_f_callback); /* this also creates the soft constraint data structure inside vc */

  /* now bind the python function to the wrapper structure */
  if(vc->sc->data){
    cb = (py_sc_callback_t *)vc->sc->data;
    /* release previous callback */
    Py_XDECREF(cb->cb_exp_f);
  } else {
    cb = (py_sc_callback_t *)vrna_alloc(sizeof(py_sc_callback_t));
    cb->cb_f        = NULL;
    cb->cb_bt       = NULL;
    cb->cb_exp_f    = NULL;
    cb->data        = NULL;
    cb->delete_data = NULL;
  }
  cb->cb_exp_f = PyFunc;  /* remember callback */
  Py_XINCREF(PyFunc); /* Increase referenc counter */

  /* finaly bind callback wrapper to fold compound */
  vc->sc->data = (void *)cb;
  if(!vc->sc->free_data){
    vc->sc->free_data = &delete_py_sc_callback;
  }
}

static vrna_basepair_t *
sc_add_bt_pycallback( vrna_fold_compound_t *vc,
                      PyObject *PyFunc){

  /* try to dispose of previous callback */
  py_sc_callback_t * cb;
  vrna_sc_add_bt(vc, &py_wrap_sc_bt_callback); /* this also creates the soft constraint data structure inside vc */

  /* now bind the python function to the wrapper structure */
  if(vc->sc->data){
    cb = (py_sc_callback_t *)vc->sc->data;
    /* release previous callback */
    Py_XDECREF(cb->cb_bt);
  } else {
    cb = (py_sc_callback_t *)vrna_alloc(sizeof(py_sc_callback_t));
    cb->cb_f        = NULL;
    cb->cb_bt       = NULL;
    cb->cb_exp_f    = NULL;
    cb->data        = NULL;
    cb->delete_data = NULL;
  }
  cb->cb_bt = PyFunc;  /* remember callback */
  Py_XINCREF(PyFunc); /* Increase referenc counter */

  /* finaly bind callback wrapper to fold compound */
  vc->sc->data = (void *)cb;
  if(!vc->sc->free_data){
    vc->sc->free_data = &delete_py_sc_callback;
  }
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
sc_add_pydata(vrna_fold_compound_t *vc,
              PyObject *data,
              PyObject *PyFunc){

  py_sc_callback_t * cb;

  /* create soft constraints data structure */
  if(!vc->sc)
    vrna_sc_init(vc);

  /* try to dispose of previous data */
  if(vc->sc->data){
    cb = (py_sc_callback_t *)vc->auxdata;
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
    cb              = (py_sc_callback_t *)vrna_alloc(sizeof(py_sc_callback_t));
    cb->cb_f        = NULL;
    cb->cb_bt       = NULL;
    cb->cb_exp_f    = NULL;
    cb->data        = NULL;
    cb->delete_data = NULL;
  }
  cb->data        = data;   /* remember data */
  cb->delete_data = PyFunc; /* remember delete data function */
  /* increase reference counter */
  Py_XINCREF(data);
  Py_XINCREF(PyFunc);

  vc->sc->data = (void *)cb;
  if(!vc->sc->free_data){
    vc->sc->free_data = &delete_py_sc_callback;
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

static int
py_wrap_sc_f_callback(int i,
                      int j,
                      int k,
                      int l,
                      char d,
                      void *data){

  int ret;
  PyObject *func, *arglist, *result;
  py_sc_callback_t *cb = (py_sc_callback_t *)data;

  func = cb->cb_f;
  /* compose argument list */
  arglist = Py_BuildValue("(i,i,i,i,i,O)", i, j, k, l, (int)d, (cb->data) ? cb->data : Py_None);
  result =  PyEval_CallObject(func, arglist);
  ret = (int)PyInt_AsLong(result);
  Py_DECREF(arglist);
  Py_XDECREF(result);
  return ret;
}

static vrna_basepair_t *
py_wrap_sc_bt_callback( int i,
                        int j,
                        int k,
                        int l,
                        char d,
                        void *data){

  int c, len, num_pairs;
  PyObject *func, *arglist, *result, *bp;
  py_sc_callback_t *cb = (py_sc_callback_t *)data;
  vrna_basepair_t *ptr, *pairs = NULL;
  func = cb->cb_bt;
  /* compose argument list */
  arglist = Py_BuildValue("(i,i,i,i,i,O)", i, j, k, l, (int)d, (cb->data) ? cb->data : Py_None);
  result =  PyEval_CallObject(func, arglist);
  if((result == NULL) || (result == Py_None))
    return NULL;

  if(PyList_Check(result)){
    len       = 10;
    num_pairs = 0;
    pairs     = (vrna_basepair_t *)vrna_alloc(sizeof(vrna_basepair_t) * len);
  
    /* result should be list of pairs */
    for(c=0; c < PyList_Size(result); c++){
      bp = PyList_GetItem(result, c);
      /* extract pair */
      if(PyTuple_Check(bp)){
        if(   (PyTuple_Size(bp) == 2)
          &&  PyInt_Check(PyTuple_GetItem(bp, 0))
          &&  PyInt_Check(PyTuple_GetItem(bp, 1))){
          pairs[num_pairs].i = (int)PyInt_AsLong(PyTuple_GetItem(bp, 0));
          pairs[num_pairs].j = (int)PyInt_AsLong(PyTuple_GetItem(bp, 1));
          num_pairs++;
        }
      }
      /* maybe the user was so kind to create a list of vrna_basepair_t? */
      else if(SWIG_ConvertPtr(bp, (void **) &ptr, SWIGTYPE_p_vrna_basepair_t, SWIG_POINTER_EXCEPTION) == 0){
        pairs[num_pairs] = *ptr;
        num_pairs++;
      } else {
        continue;
      }

      /* increase length if necessary */
      if(num_pairs == len){
        len = (int)(1.2 * len);
        pairs = (vrna_basepair_t *)vrna_realloc(pairs, sizeof(vrna_basepair_t)*len);
      }
    }
    /* put end marker in list */
    pairs[num_pairs].i = pairs[num_pairs].j = 0;
    pairs = (vrna_basepair_t *)vrna_realloc(pairs, sizeof(vrna_basepair_t)*(num_pairs+1));
  }
  Py_DECREF(arglist);
  Py_XDECREF(result);
  return pairs;
}

static FLT_OR_DBL
py_wrap_sc_exp_f_callback(int i,
                          int j,
                          int k,
                          int l,
                          char d,
                          void *data){

  FLT_OR_DBL ret;
  PyObject *func, *arglist, *result;
  py_sc_callback_t *cb = (py_sc_callback_t *)data;

  func = cb->cb_exp_f;
  /* compose argument list */
  arglist = Py_BuildValue("(i,i,i,i,i,O)", i, j, k, l, (int)d, (cb->data) ? cb->data : Py_None);
  result =  PyEval_CallObject(func, arglist);
  ret = (FLT_OR_DBL)PyFloat_AsDouble(result);
  Py_DECREF(arglist);
  Py_XDECREF(result);
  return ret;
}

%}

static void fc_add_pycallback(vrna_fold_compound_t *vc, PyObject *PyFunc);
static void fc_add_pydata(vrna_fold_compound_t *vc, PyObject *data, PyObject *PyFunc);

static void sc_add_f_pycallback(vrna_fold_compound_t *vc, PyObject *PyFunc);
static vrna_basepair_t *sc_add_bt_pycallback(vrna_fold_compound_t *vc, PyObject *PyFunc);
static void sc_add_exp_f_pycallback(vrna_fold_compound_t *vc, PyObject *PyFunc);
static void sc_add_pydata(vrna_fold_compound_t *vc, PyObject *data, PyObject *PyFunc);

%typemap(in) PyObject *PyFunc {
  if (!PyCallable_Check($input)) {
      PyErr_SetString(PyExc_TypeError, "Need a callable object!");
      return NULL;
  }
  $1 = $input;
}
#endif
