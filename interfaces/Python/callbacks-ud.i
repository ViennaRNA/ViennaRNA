/**********************************************/
/* BEGIN interface for unstructured domains   */
/* feature callbacks                          */
/**********************************************/

#ifdef SWIGPYTHON
%{

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

  cb->prod_rule     = Py_None;
  cb->exp_prod_rule = Py_None;
  cb->energy        = Py_None;
  cb->exp_energy    = Py_None;
  cb->data          = Py_None;
  cb->delete_data   = Py_None;
  cb->prob_add      = Py_None;
  cb->prob_get      = Py_None;

  return cb;
}

static void
delete_py_ud_callback(void * data){

  py_ud_callback_t *cb = (py_ud_callback_t *)data;
  /* first delete user data */
  if(cb->delete_data != Py_None){
    PyObject *func, *arglist, *result;
    func = cb->delete_data;
    arglist = Py_BuildValue("O", cb->data);
    result  = PyObject_CallObject(func, arglist);
    Py_DECREF(arglist);
    Py_XDECREF(result);
  }

  /* now dispose of the registered callbacks */
  Py_XDECREF(cb->prod_rule);
  Py_XDECREF(cb->exp_prod_rule);
  Py_XDECREF(cb->energy);
  Py_XDECREF(cb->exp_energy);
  Py_XDECREF(cb->prob_add);
  Py_XDECREF(cb->prob_get);

  /* finally free pycallback */
  free(cb);
}


static void
ud_set_pydata(vrna_fold_compound_t *vc,
              PyObject *data,
              PyObject *PyFunc){

  py_ud_callback_t * cb;

  /* try to dispose of previous data */
  if(vc->domains_up && vc->domains_up->data){
    cb = (py_ud_callback_t *)vc->domains_up->data;
    if(cb->data != Py_None){
      if(cb->delete_data != Py_None){
        PyObject *func, *arglist, *result;
        func    = cb->delete_data;
        arglist = Py_BuildValue("O", cb->data);
        result  = PyObject_CallObject(func, arglist);
        Py_DECREF(arglist);
        Py_XDECREF(result);
      }
    }
    Py_XDECREF(cb->data);
    Py_XDECREF(cb->delete_data);
  } else {
    cb = new_py_ud_cb();
  }
  cb->data        = data;   /* remember data */
  cb->delete_data = PyFunc; /* remember delete data function */

  /* increase reference counter */
  Py_XINCREF(data);
  Py_XINCREF(PyFunc);

  /* bind callback wrapper to fold compound */
  vrna_ud_set_data(vc, (void *)cb, &delete_py_ud_callback);
}


static void
ud_set_prod_cb( vrna_fold_compound_t *vc,
                PyObject *prod_cb,
                PyObject *eval_cb){

  /* try to dispose of previous callback */
  py_ud_callback_t * cb;
  if(vc->domains_up && vc->domains_up->data){
    cb = (py_ud_callback_t *)vc->domains_up->data;
    /* release previous callback */
    Py_XDECREF(cb->prod_rule);
    Py_XDECREF(cb->energy);
  } else {
    cb = new_py_ud_cb();
    /* bind callback wrapper to fold compound */
    vrna_ud_set_data(vc, (void *)cb, &delete_py_ud_callback);
  }
  cb->prod_rule = prod_cb; /* remember callback */
  cb->prod_rule = eval_cb; /* remember callback */

  Py_XINCREF(prod_cb);     /* Increase reference counter */
  Py_XINCREF(eval_cb);     /* Increase reference counter */

  vrna_ud_set_prod_rule_cb(vc, &py_wrap_ud_prod_rule, &py_wrap_ud_energy);
}


static void
ud_set_exp_prod_cb( vrna_fold_compound_t *vc,
                    PyObject *prod_cb,
                    PyObject *eval_cb){

  /* try to dispose of previous callback */
  py_ud_callback_t *cb;

  /* now bind the python function to the wrapper structure */
  if(vc->domains_up && vc->domains_up->data){
    cb = (py_ud_callback_t *)vc->domains_up->data;
    /* release previous callback */
    Py_XDECREF(cb->exp_prod_rule);
    Py_XDECREF(cb->exp_energy);
  } else {
    cb = new_py_ud_cb();
    /* bind callback wrapper to fold compound */
    vrna_ud_set_data(vc, (void *)cb, &delete_py_ud_callback);
  }
  cb->exp_prod_rule = prod_cb;
  cb->exp_energy    = eval_cb;

  Py_XINCREF(prod_cb); /* Increase reference counter */
  Py_XINCREF(eval_cb); /* Increase reference counter */

  vrna_ud_set_exp_prod_rule_cb(vc, &py_wrap_ud_exp_prod_rule, &py_wrap_ud_exp_energy);
}


static void
ud_set_prob_cb( vrna_fold_compound_t *vc,
                PyObject *setter,
                PyObject *getter){

  py_ud_callback_t *cb;

  /* now bind the python function to the wrapper structure */
  if(vc->domains_up && vc->domains_up->data){
    cb = (py_ud_callback_t *)vc->domains_up->data;
    /* release previous callbacks */
    Py_XDECREF(cb->prob_add);
    Py_XDECREF(cb->prob_get);
  } else {
    cb = new_py_ud_cb();
    /* bind callback wrapper to fold compound */
    vrna_ud_set_data(vc, (void *)cb, &delete_py_ud_callback);
  }
  cb->prob_add = setter;
  cb->prob_get = getter;
  Py_XINCREF(setter); /* Increase reference counter */
  Py_XINCREF(getter); /* Increase reference counter */

  vrna_ud_set_prob_cb(vc, &py_wrap_ud_prob_add, &py_wrap_ud_prob_get);
}


static void
py_wrap_ud_prod_rule( vrna_fold_compound_t *vc,
                      void *data){

  int ret;
  PyObject *func, *arglist, *result;
  py_ud_callback_t *cb = (py_ud_callback_t *)data;

  func = cb->prod_rule;
  /* compose argument list */
  arglist = Py_BuildValue("(O,O)", vc, (cb->data) ? cb->data : Py_None);
  result =  PyObject_CallObject(func, arglist);
  Py_DECREF(arglist);
  Py_XDECREF(result);
  return /*void*/;
}


static void
py_wrap_ud_exp_prod_rule( vrna_fold_compound_t *vc,
                          void *data){

  int ret;
  PyObject *func, *arglist, *result;
  py_ud_callback_t *cb = (py_ud_callback_t *)data;

  func = cb->exp_prod_rule;
  /* compose argument list */
  arglist = Py_BuildValue("(O,O)", vc, (cb->data) ? cb->data : Py_None);
  result =  PyObject_CallObject(func, arglist);
  Py_DECREF(arglist);
  Py_XDECREF(result);
  return /*void*/;
}


static int
py_wrap_ud_energy(vrna_fold_compound_t *vc,
                  int i,
                  int j,
                  unsigned int looptype,
                  void *data){

  int ret;
  PyObject *func, *arglist, *result;
  py_ud_callback_t *cb = (py_ud_callback_t *)data;

  func = cb->energy;
  /* compose argument list */
  arglist = Py_BuildValue("(O,i,i,I,O)", vc, i, j, looptype, (cb->data) ? cb->data : Py_None);
  result =  PyObject_CallObject(func, arglist);
  ret = (int)PyInt_AsLong(result);
  Py_DECREF(arglist);
  Py_XDECREF(result);
  return ret;
}


static FLT_OR_DBL
py_wrap_ud_exp_energy(vrna_fold_compound_t *vc,
                      int i,
                      int j,
                      unsigned int looptype,
                      void *data){

  FLT_OR_DBL ret;
  PyObject *func, *arglist, *result;
  py_ud_callback_t *cb = (py_ud_callback_t *)data;

  func = cb->exp_energy;
  /* compose argument list */
  arglist = Py_BuildValue("(O,i,i,I,O)", vc, i, j, looptype, (cb->data) ? cb->data : Py_None);
  result =  PyObject_CallObject(func, arglist);
  ret = (FLT_OR_DBL)PyFloat_AsDouble(result);
  Py_DECREF(arglist);
  Py_XDECREF(result);
  return ret;
}


static void
py_wrap_ud_prob_add(vrna_fold_compound_t *vc,
                    int i,
                    int j,
                    unsigned int looptype,
                    FLT_OR_DBL prob,
                    void *data){

  PyObject *func, *arglist, *result;
  py_ud_callback_t *cb = (py_ud_callback_t *)data;

  func = cb->prob_add;
  /* compose argument list */
  arglist = Py_BuildValue("(O,i,i,I,d,O)", vc, i, j, looptype, (double)prob, (cb->data) ? cb->data : Py_None);
  result =  PyObject_CallObject(func, arglist);
  Py_DECREF(arglist);
  Py_XDECREF(result);
  return;
}


static FLT_OR_DBL
py_wrap_ud_prob_get(vrna_fold_compound_t *vc,
                    int i,
                    int j,
                    unsigned int looptype,
                    int motif,
                    void *data){

  int ret;
  PyObject *func, *arglist, *result;
  py_ud_callback_t *cb = (py_ud_callback_t *)data;

  func = cb->prob_get;
  /* compose argument list */
  arglist = Py_BuildValue("(O,i,i,I,i,O)", vc, i, j, looptype, motif, (cb->data) ? cb->data : Py_None);
  result  =  PyObject_CallObject(func, arglist);
  ret     = (int)PyInt_AsLong(result);
  Py_DECREF(arglist);
  Py_XDECREF(result);
  return ret;
}

%}

static void ud_set_pydata(vrna_fold_compound_t *vc, PyObject *data, PyObject *PyFuncOrNone);
static void ud_set_prod_cb(vrna_fold_compound_t *vc, PyObject *prod_cb, PyObject *eval_cb);
static void ud_set_exp_prod_cb(vrna_fold_compound_t *vc, PyObject *prod_cb, PyObject *eval_cb);
static void ud_set_prob_cb( vrna_fold_compound_t *vc, PyObject *setter, PyObject *getter);

/* now we bind the above functions as methods to the fold_compound object */
%extend vrna_fold_compound_t {

  PyObject *ud_set_data(PyObject *data, PyObject *PyFuncOrNone=Py_None){
    ud_set_pydata($self, data, PyFuncOrNone);
    Py_RETURN_NONE;
  }

  PyObject *ud_set_prod_rule_cb(PyObject *prod_cb, PyObject *eval_cb){
    if(!PyCallable_Check(prod_cb)) {
      PyErr_SetString(PyExc_TypeError, "Need a callable object!");
      Py_RETURN_NONE;
    } else if(!PyCallable_Check(eval_cb)){
      PyErr_SetString(PyExc_TypeError, "Need a callable object!");
      Py_RETURN_NONE;
    }
    ud_set_prod_cb($self, prod_cb, eval_cb);
    Py_RETURN_NONE;
  }

  PyObject *ud_set_exp_prod_rule_cb(PyObject *prod_cb, PyObject *eval_cb){
    if(!PyCallable_Check(prod_cb)) {
      PyErr_SetString(PyExc_TypeError, "Need a callable object!");
      Py_RETURN_NONE;
    } else if(!PyCallable_Check(eval_cb)){
      PyErr_SetString(PyExc_TypeError, "Need a callable object!");
      Py_RETURN_NONE;
    }
    ud_set_exp_prod_cb($self, prod_cb, eval_cb);
    Py_RETURN_NONE;
  }

  PyObject *ud_set_prob_cb(PyObject *setter, PyObject *getter){
    if(!PyCallable_Check(setter)) {
      PyErr_SetString(PyExc_TypeError, "Need a callable object!");
      Py_RETURN_NONE;
    } else if(!PyCallable_Check(getter)){
      PyErr_SetString(PyExc_TypeError, "Need a callable object!");
      Py_RETURN_NONE;
    }
    ud_set_prob_cb($self, setter, getter);
    Py_RETURN_NONE;
  }
}

#endif
