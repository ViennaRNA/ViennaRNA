/**********************************************/
/* BEGIN interface for subopt callback        */
/**********************************************/

#ifdef SWIGPYTHON
%{

typedef struct {
  PyObject *cb;
  PyObject *data;
} python_subopt_callback_t;

static python_subopt_callback_t * bind_subopt_callback(PyObject *PyFunc, PyObject *data);
static void python_wrap_subopt_cb(const char *structure, float energy, void *data);

static python_subopt_callback_t *
bind_subopt_callback(PyObject *PyFunc, PyObject *data){

  python_subopt_callback_t *cb = (python_subopt_callback_t *)vrna_alloc(sizeof(python_subopt_callback_t));

  cb->cb    = PyFunc;  /* store callback */
  cb->data  = data;    /* bind data */

  return cb;
}

static void
python_wrap_subopt_cb(const char *structure, float energy, void *data){

  PyObject *func, *arglist, *result;
  python_subopt_callback_t *cb = (python_subopt_callback_t *)data;

  func = cb->cb;
  /* compose argument list */
  arglist = Py_BuildValue("(z,d,O)", structure, (double)energy, (cb->data) ? cb->data : Py_None);
  result =  PyObject_CallObject(func, arglist);
  Py_DECREF(arglist);
  Py_XDECREF(result);

  return /*void*/;
}

%}

/* now we bind the above functions as methods to the fold_compound object */
%extend vrna_fold_compound_t {

  PyObject *subopt_cb(int delta, PyObject *PyFunc, PyObject *data = Py_None){

    python_subopt_callback_t *cb = bind_subopt_callback(PyFunc, data);
    vrna_subopt_cb($self, delta, &python_wrap_subopt_cb, (void *)cb);
    free(cb);
    Py_RETURN_NONE;
  }

}


#endif
