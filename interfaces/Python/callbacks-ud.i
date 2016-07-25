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
    cb->prod_rule     = NULL;
    cb->exp_prod_rule = NULL;
    cb->energy        = NULL;
    cb->exp_energy    = NULL;
    cb->data          = NULL;
    cb->delete_data   = NULL;
  }
  cb->data        = data;   /* remember data */
  cb->delete_data = PyFunc; /* remember delete data function */
  /* increase reference counter */
  Py_XINCREF(data);
  Py_XINCREF(PyFunc);

  vrna_ud_set_data(vc, (void *)cb, &delete_py_ud_callback);

}


%}

static void ud_set_pydata(vrna_fold_compound_t *vc, PyObject *data, PyObject *PyFuncOrNone);
/*
static void ud_set_prod_rule_cb(vrna_fold_compound_t *vc, PyObject *PyFunc);
static void ud_set_exp_prod_rule_cb(vrna_fold_compound_t *vc, PyObject *PyFunc);
static void ud_set_energy_cb(vrna_fold_compound_t *vc, PyObject *PyFunc);
static void ud_set_exp_energy_cb(vrna_fold_compound_t *vc, PyObject *PyFunc);
*/

/* now we bind the above functions as methods to the fold_compound object */
%extend vrna_fold_compound_t {

  void ud_set_data(PyObject *data, PyObject *PyFuncOrNone=Py_None){
    ud_set_pydata($self, data, PyFuncOrNone);
  }

}

#endif
