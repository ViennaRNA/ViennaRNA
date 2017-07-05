/**********************************************/
/* BEGIN interface for PF window callback     */
/**********************************************/

#ifdef SWIGPYTHON
%{

typedef struct {
  PyObject *cb;
  PyObject *data;
} python_pf_window_callback_t;

static python_pf_window_callback_t * bind_pf_window_callback(PyObject *PyFunc, PyObject *data);
static void python_wrap_pf_window_cb(FLT_OR_DBL *pr, int pr_size, int i, int max, unsigned int type, void *data);

static python_pf_window_callback_t *
bind_pf_window_callback(PyObject *PyFunc, PyObject *data){

  python_pf_window_callback_t *cb = (python_pf_window_callback_t *)vrna_alloc(sizeof(python_pf_window_callback_t));

  cb->cb    = PyFunc;  /* store callback */
  cb->data  = data;    /* bind data */

  return cb;
}

static void
python_wrap_pf_window_cb(FLT_OR_DBL *pr, int pr_size, int i, int max, unsigned int type, void *data){

  PyObject *func, *arglist, *result, *pr_list;
  python_pf_window_callback_t *cb = (python_pf_window_callback_t *)data;

  func = cb->cb;

  if (type & VRNA_PROBS_WINDOW_UP) { /* We distinguish output for unpaired probabilities */

    /* create PYTHON list for unpaired probabilities */
    pr_list = PyList_New((Py_ssize_t) max + 1);

    /* 0th element */
    PyList_SET_ITEM(pr_list, 0, Py_None);

    /* actual values in range [1, MIN(i, max)] */
    for (int cnt = 1; cnt <= pr_size; cnt++)
      PyList_SET_ITEM(pr_list, (Py_ssize_t) cnt, PyFloat_FromDouble(pr[cnt]));

    /* empty values in range [0,i - 1] */
    for (int cnt = pr_size + 1; cnt <= max; cnt++)
      PyList_SET_ITEM(pr_list, (Py_ssize_t) cnt, Py_None);

  } else { /* and pairing/stacking probabilities for pair (i, j) or ensemble free energies for subsegment [i, j] */

    /* create PYTHON list for pr values */
    pr_list = PyList_New((Py_ssize_t) (pr_size + 1));

    /* empty values in range [0, i] */
    for (int cnt = 0; cnt <= i; cnt++)
      PyList_SET_ITEM(pr_list, (Py_ssize_t) cnt, Py_None);

    /* actual values in range [i + 1, pr_size] */
    for (int cnt = i + 1; cnt <= pr_size; cnt++)
      PyList_SET_ITEM(pr_list, (Py_ssize_t) cnt, PyFloat_FromDouble(pr[cnt]));
  }

  /* compose argument list */
  arglist = Py_BuildValue("(O, i, i, i, i, O)", pr_list, pr_size, i, max, type, (cb->data) ? cb->data : Py_None);
  result =  PyObject_CallObject(func, arglist);
  Py_DECREF(arglist);
  Py_XDECREF(result);

  return /*void*/;
}

%}

/* now we bind vrna_probs_window() as method to the fold_compound object using the above callback wrapper */
%extend vrna_fold_compound_t {

  void probs_window(int ulength, unsigned int options, PyObject *PyFunc, PyObject *data = Py_None) {
    python_pf_window_callback_t *cb = bind_pf_window_callback(PyFunc, data);
    vrna_probs_window($self, ulength, options, &python_wrap_pf_window_cb, (void *)cb);
    free(cb);
  }
}


/* Also add wrappers for the 'simple' callback interface of pfl_fold_*_cb() functions */
%{

  void pfl_fold_cb(std::string sequence, int window_size, int max_bp_span, PyObject *PyFunc, PyObject *data = Py_None) {
    python_pf_window_callback_t *cb = bind_pf_window_callback(PyFunc, data);
    vrna_pfl_fold_cb(sequence.c_str(), window_size, max_bp_span, &python_wrap_pf_window_cb, (void *)cb);
    free(cb);
  }

  void pfl_fold_up_cb(std::string sequence, int ulength, int window_size, int max_bp_span, PyObject *PyFunc, PyObject *data = Py_None) {
    python_pf_window_callback_t *cb = bind_pf_window_callback(PyFunc, data);
    vrna_pfl_fold_up_cb(sequence.c_str(), ulength, window_size, max_bp_span, &python_wrap_pf_window_cb, (void *)cb);
    free(cb);
  }


%}

void pfl_fold_cb(std::string sequence, int window_size, int max_bp_span, PyObject *PyFunc, PyObject *data = Py_None);
void pfl_fold_up_cb(std::string sequence, int ulength, int window_size, int max_bp_span, PyObject *PyFunc, PyObject *data = Py_None);

#endif
