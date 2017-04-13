/**********************************************/
/* BEGIN interface for MFE window callback    */
/**********************************************/

#ifdef SWIGPYTHON
%{

typedef struct {
  PyObject *cb;
  PyObject *data;
} python_mfe_window_callback_t;

static python_mfe_window_callback_t * bind_mfe_window_callback(PyObject *PyFunc, PyObject *data);
static void python_wrap_mfe_window_cb(int start, int end, const char *structure, float energy, void *data);

#ifdef VRNA_WITH_SVM
static void python_wrap_mfe_window_zscore_cb(int start, int end, const char *structure, float energy, float zscore, void *data);
#endif

static python_mfe_window_callback_t *
bind_mfe_window_callback(PyObject *PyFunc, PyObject *data){

  python_mfe_window_callback_t *cb = (python_mfe_window_callback_t *)vrna_alloc(sizeof(python_mfe_window_callback_t));

  cb->cb    = PyFunc;  /* store callback */
  cb->data  = data;    /* bind data */

  return cb;
}

static void
python_wrap_mfe_window_cb(int start, int end, const char *structure, float energy, void *data){

  PyObject *func, *arglist, *result;
  python_mfe_window_callback_t *cb = (python_mfe_window_callback_t *)data;

  func = cb->cb;
  /* compose argument list */
  arglist = Py_BuildValue("(i, i, z, d,O)", start, end, structure, (double)energy, (cb->data) ? cb->data : Py_None);
  result =  PyObject_CallObject(func, arglist);
  Py_DECREF(arglist);
  Py_XDECREF(result);

  return /*void*/;
}

#ifdef VRNA_WITH_SVM
static void
python_wrap_mfe_window_zscore_cb(int start, int end, const char *structure, float energy, float zscore, void *data){

  PyObject *func, *arglist, *result;
  python_mfe_window_callback_t *cb = (python_mfe_window_callback_t *)data;

  func = cb->cb;
  /* compose argument list */
  arglist = Py_BuildValue("(i, i, z, d, d, O)", start, end, structure, (double)energy, (double)zscore, (cb->data) ? cb->data : Py_None);
  result =  PyObject_CallObject(func, arglist);
  Py_DECREF(arglist);
  Py_XDECREF(result);

  return /*void*/;
}
#endif

%}

/* now we bind the above functions as methods to the fold_compound object */
%extend vrna_fold_compound_t {

  float mfe_window_cb(PyObject *PyFunc, PyObject *data = Py_None) {
    float en;
    python_mfe_window_callback_t *cb = bind_mfe_window_callback(PyFunc, data);
    en = vrna_mfe_window_cb($self, &python_wrap_mfe_window_cb, (void *)cb);
    free(cb);
    return en;
  }

#ifdef VRNA_WITH_SVM
  float mfe_window_score_cb(double min_z, PyObject *PyFunc, PyObject *data = Py_None) {
    float en;
    python_mfe_window_callback_t *cb = bind_mfe_window_callback(PyFunc, data);
    en = vrna_mfe_window_zscore_cb($self, min_z, &python_wrap_mfe_window_zscore_cb, (void *)cb);
    free(cb);
    return en;
  }
#endif

}

%rename (Lfold_cb) my_Lfold_cb;
%rename (Lfoldz_cb) my_Lfoldz_cb;
%rename (aliLfold_cb) my_aliLfold_cb;

%{
  float my_Lfold_cb(char *string, int window_size, PyObject *PyFunc, PyObject *data = Py_None) {
    float en;
    python_mfe_window_callback_t *cb = bind_mfe_window_callback(PyFunc, data);
    en = vrna_Lfold_cb(string, window_size, &python_wrap_mfe_window_cb, (void *)cb);
    free(cb);
    return en;
  }

#ifdef VRNA_WITH_SVM
  float my_Lfoldz_cb(char *string, int window_size, double min_z, PyObject *PyFunc, PyObject *data = Py_None) {
    float en;
    python_mfe_window_callback_t *cb = bind_mfe_window_callback(PyFunc, data);
    en = vrna_Lfoldz_cb(string, window_size, min_z, &python_wrap_mfe_window_zscore_cb, (void *)cb);
    free(cb);
    return en;
  }
#endif

  float my_aliLfold_cb(std::vector<std::string> alignment, int window_size, PyObject *PyFunc, PyObject *data = Py_None) {
    float en;
    python_mfe_window_callback_t *cb = bind_mfe_window_callback(PyFunc, data);
    std::vector<const char*>  vc;
    std::transform(alignment.begin(), alignment.end(), std::back_inserter(vc), convert_vecstring2veccharcp);
    vc.push_back(NULL); /* mark end of sequences */
    en = vrna_aliLfold_cb((const char **)&vc[0], window_size, &python_wrap_mfe_window_cb, (void *)cb);
    free(cb);
    return en;
  }

%}

float my_Lfold_cb(char *string, int window_size, PyObject *PyFunc, PyObject *data);
%ignore vrna_Lfold_cb;

#ifdef VRNA_WITH_SVM
float my_Lfoldz_cb(char *string, int window_size, double min_z, PyObject *PyFunc, PyObject *data);
%ignore vrna_Lfoldz_cb;
#endif

float my_aliLfold_cb(std::vector<std::string> alignment, int window_size, PyObject *PyFunc, PyObject *data);
%ignore vrna_aliLfold_cb;
%ignore aliLfold_cb;



#endif
