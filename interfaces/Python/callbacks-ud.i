/**********************************************/
/* BEGIN interface for fold compound status   */
/* callback                                   */
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
} py_ud_callback_t;

static vrna_callback_ud_production      py_wrap_ud_prod_rule;
static vrna_callback_ud_exp_production  py_wrap_ud_exp_prod_rule;
static vrna_callback_ud_energy          py_wrap_ud_energy;
static vrna_callback_ud_exp_energy      py_wrap_ud_exp_energy;


static void
delete_py_ud_callback(void * data){

  py_ud_callback_t *cb = (py_ud_callback_t *)data;
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
  Py_XDECREF(cb->prod_rule);
  Py_XDECREF(cb->exp_prod_rule);
  Py_XDECREF(cb->energy);
  Py_XDECREF(cb->exp_energy);

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
        result  = PyEval_CallObject(func, arglist);
        Py_DECREF(arglist);
        Py_XDECREF(result);
      }
    }
    Py_XDECREF(cb->data);
    Py_XDECREF(cb->delete_data);
  } else {
    cb                = (py_ud_callback_t *)vrna_alloc(sizeof(py_ud_callback_t));
    cb->prod_rule     = Py_None;
    cb->exp_prod_rule = Py_None;
    cb->energy        = Py_None;
    cb->exp_energy    = Py_None;
    cb->data          = Py_None;
    cb->delete_data   = Py_None;
    /* bind callback wrapper to fold compound */
    vrna_ud_set_data(vc, (void *)cb, &delete_py_ud_callback);
  }
  cb->data        = data;   /* remember data */
  cb->delete_data = PyFunc; /* remember delete data function */

  /* increase reference counter */
  Py_XINCREF(data);
  Py_XINCREF(PyFunc);
}


static void
ud_set_prod_rule_cb(vrna_fold_compound_t *vc,
                    PyObject *PyFunc){

  /* try to dispose of previous callback */
  py_ud_callback_t * cb;
  if(vc->domains_up && vc->domains_up->data){
    cb = (py_ud_callback_t *)vc->domains_up->data;
    /* release previous callback */
    Py_XDECREF(cb->prod_rule);
  } else {
    cb = (py_ud_callback_t *)vrna_alloc(sizeof(py_ud_callback_t));
    cb->prod_rule     = Py_None;
    cb->exp_prod_rule = Py_None;
    cb->energy        = Py_None;
    cb->exp_energy    = Py_None;
    cb->data          = Py_None;
    cb->delete_data   = Py_None;
    /* bind callback wrapper to fold compound */
    vrna_ud_set_data(vc, (void *)cb, &delete_py_ud_callback);
  }
  cb->prod_rule = PyFunc; /* remember callback */

  Py_XINCREF(PyFunc);     /* Increase referenc counter */

  vrna_ud_set_prod_rule(vc, &py_wrap_ud_prod_rule);
}


static void
ud_set_exp_prod_rule_cb(vrna_fold_compound_t *vc,
                        PyObject *PyFunc){

  /* try to dispose of previous callback */
  py_ud_callback_t * cb;
  if(vc->domains_up && vc->domains_up->data){
    cb = (py_ud_callback_t *)vc->domains_up->data;
    /* release previous callback */
    Py_XDECREF(cb->exp_prod_rule);
  } else {
    cb = (py_ud_callback_t *)vrna_alloc(sizeof(py_ud_callback_t));
    cb->prod_rule     = Py_None;
    cb->exp_prod_rule = Py_None;
    cb->energy        = Py_None;
    cb->exp_energy    = Py_None;
    cb->data          = Py_None;
    cb->delete_data   = Py_None;
    /* bind callback wrapper to fold compound */
    vrna_ud_set_data(vc, (void *)cb, &delete_py_ud_callback);
  }
  cb->exp_prod_rule = PyFunc; /* remember callback */
  Py_XINCREF(PyFunc);         /* Increase referenc counter */

  vrna_ud_set_exp_prod_rule(vc, &py_wrap_ud_exp_prod_rule);
}


static void
ud_set_energy_cb( vrna_fold_compound_t *vc,
                  PyObject *PyFunc){

  /* try to dispose of previous callback */
  py_ud_callback_t *cb;

  /* now bind the python function to the wrapper structure */
  if(vc->domains_up && vc->domains_up->data){
    cb = (py_ud_callback_t *)vc->domains_up->data;
    /* release previous callback */
    Py_XDECREF(cb->energy);
  } else {
    cb = (py_ud_callback_t *)vrna_alloc(sizeof(py_ud_callback_t));
    cb->energy        = Py_None;
    cb->exp_energy    = Py_None;
    cb->prod_rule     = Py_None;
    cb->exp_prod_rule = Py_None;
    cb->data          = Py_None;
    cb->delete_data   = Py_None;
    /* bind callback wrapper to fold compound */
    vrna_ud_set_data(vc, (void *)cb, &delete_py_ud_callback);
  }
  cb->energy = PyFunc;  /* remember callback */
  Py_XINCREF(PyFunc); /* Increase referenc counter */

  vrna_ud_set_energy(vc, &py_wrap_ud_energy);
}


static void
ud_set_exp_energy_cb( vrna_fold_compound_t *vc,
                      PyObject *PyFunc){

  /* try to dispose of previous callback */
  py_ud_callback_t *cb;

  /* now bind the python function to the wrapper structure */
  if(vc->domains_up && vc->domains_up->data){
    cb = (py_ud_callback_t *)vc->domains_up->data;
    /* release previous callback */
    Py_XDECREF(cb->exp_energy);
  } else {
    cb = (py_ud_callback_t *)vrna_alloc(sizeof(py_ud_callback_t));
    cb->energy        = Py_None;
    cb->exp_energy    = Py_None;
    cb->prod_rule     = Py_None;
    cb->exp_prod_rule = Py_None;
    cb->data          = Py_None;
    cb->delete_data   = Py_None;
    /* bind callback wrapper to fold compound */
    vrna_ud_set_data(vc, (void *)cb, &delete_py_ud_callback);
  }
  cb->exp_energy = PyFunc;  /* remember callback */
  Py_XINCREF(PyFunc); /* Increase referenc counter */

  vrna_ud_set_exp_energy(vc, &py_wrap_ud_exp_energy);
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
  result =  PyEval_CallObject(func, arglist);
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
  result =  PyEval_CallObject(func, arglist);
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
  arglist = Py_BuildValue("(O,i,i,I,O)", i, j, looptype, (cb->data) ? cb->data : Py_None);
  result =  PyEval_CallObject(func, arglist);
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
  arglist = Py_BuildValue("(O,i,i,I,O)", i, j, looptype, (cb->data) ? cb->data : Py_None);
  result =  PyEval_CallObject(func, arglist);
  ret = (FLT_OR_DBL)PyFloat_AsDouble(result);
  Py_DECREF(arglist);
  Py_XDECREF(result);
  return ret;
}

%}

static void ud_set_pydata(vrna_fold_compound_t *vc, PyObject *data, PyObject *PyFuncOrNone);
static void ud_set_prod_rule_cb(vrna_fold_compound_t *vc, PyObject *PyFunc);
static void ud_set_exp_prod_rule_cb(vrna_fold_compound_t *vc, PyObject *PyFunc);
static void ud_set_energy_cb(vrna_fold_compound_t *vc, PyObject *PyFunc);
static void ud_set_exp_energy_cb(vrna_fold_compound_t *vc, PyObject *PyFunc);

/* now we bind the above functions as methods to the fold_compound object */
%extend vrna_fold_compound_t {

  void ud_set_data(PyObject *data, PyObject *PyFuncOrNone=Py_None){
    ud_set_pydata($self, data, PyFuncOrNone);
  }

  void ud_set_prod_rule(PyObject *PyFuncOrNone=Py_None){
    ud_set_prod_rule_cb($self, PyFuncOrNone);
  }

  void ud_set_exp_prod_rule(PyObject *PyFuncOrNone=Py_None){
    ud_set_exp_prod_rule_cb($self, PyFuncOrNone);
  }

}

#endif
