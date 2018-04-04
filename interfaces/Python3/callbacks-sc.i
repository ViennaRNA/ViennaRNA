/**********************************************/
/* BEGIN interface for fold compound status   */
/* callback                                   */
/**********************************************/

#ifdef SWIGPYTHON
%{

#include <stdexcept>

typedef struct {
  PyObject  *cb_f;
  PyObject  *cb_bt;
  PyObject  *cb_exp_f;
  PyObject  *data;
  PyObject  *delete_data;
} py_sc_callback_t;

static int              py_wrap_sc_f_callback(int i, int j, int k, int l, unsigned char d, void *data);
static vrna_basepair_t  *py_wrap_sc_bt_callback( int i, int j, int k, int l, unsigned char d, void *data);
static FLT_OR_DBL       py_wrap_sc_exp_f_callback(int i, int j, int k, int l, unsigned char d, void *data);

static void
delete_py_sc_callback(void * data){

  py_sc_callback_t *cb = (py_sc_callback_t *)data;
  /* first delete user data */
  if(cb->delete_data != Py_None){
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
          throw std::runtime_error( "Generic soft constraint delete_data() callback must take exactly 1 argument" );
        } else {
          throw std::runtime_error( "Some error occurred while executing generic soft constraint delete_data() callback" );
        }
      }
      PyErr_Clear();
    }
    /* END recognizing errors in callback execution */

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
    cb->cb_f        = Py_None;
    cb->cb_bt       = Py_None;
    cb->cb_exp_f    = Py_None;
    cb->data        = Py_None;
    cb->delete_data = Py_None;
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
    cb->cb_f        = Py_None;
    cb->cb_bt       = Py_None;
    cb->cb_exp_f    = Py_None;
    cb->data        = Py_None;
    cb->delete_data = Py_None;
  }
  cb->cb_exp_f = PyFunc;  /* remember callback */
  Py_XINCREF(PyFunc); /* Increase referenc counter */

  /* finaly bind callback wrapper to fold compound */
  vc->sc->data = (void *)cb;
  if(!vc->sc->free_data){
    vc->sc->free_data = &delete_py_sc_callback;
  }
}

static void
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
    cb->cb_f        = Py_None;
    cb->cb_bt       = Py_None;
    cb->cb_exp_f    = Py_None;
    cb->data        = Py_None;
    cb->delete_data = Py_None;
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
sc_add_pydata(vrna_fold_compound_t *vc,
              PyObject *data,
              PyObject *PyFunc){

  py_sc_callback_t * cb;

  /* create soft constraints data structure */
  if(!vc->sc)
    vrna_sc_init(vc);

  /* try to dispose of previous data */
  if(vc->sc->data){
    cb = (py_sc_callback_t *)vc->sc->data;
    if(cb->data != Py_None){
      if(cb->delete_data != Py_None){
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
            if (PyErr_GivenExceptionMatches(err, PyExc_TypeError)) {
              throw std::runtime_error( "Generic soft constraint delete_data() callback must take exactly 1 argument" );
            } else {
              throw std::runtime_error( "Some error occurred while executing generic soft constraint delete_data() callback" );
            }
          }
          PyErr_Clear();
        }
        /* END recognizing errors in callback execution */

        Py_DECREF(arglist);
        Py_XDECREF(result);
      }
    }
    Py_XDECREF(cb->data);
    Py_XDECREF(cb->delete_data);
  } else {
    cb              = (py_sc_callback_t *)vrna_alloc(sizeof(py_sc_callback_t));
    cb->cb_f        = Py_None;
    cb->cb_bt       = Py_None;
    cb->cb_exp_f    = Py_None;
    cb->data        = Py_None;
    cb->delete_data = Py_None;
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

static int
py_wrap_sc_f_callback(int i,
                      int j,
                      int k,
                      int l,
                      unsigned char d,
                      void *data){

  int ret;
  PyObject *func, *arglist, *result, *err;
  py_sc_callback_t *cb = (py_sc_callback_t *)data;

  ret  = 0;
  func = cb->cb_f;

  /* compose argument list */
  arglist = Py_BuildValue("(i,i,i,i,i,O)", i, j, k, l, (int)d, (cb->data) ? cb->data : Py_None);
  result =  PyObject_CallObject(func, arglist);
  /* BEGIN recognizing errors in callback execution */
  if (result == NULL) {
    if ((err = PyErr_Occurred())) {
      /* print error message */
      PyErr_Print();
      /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
      if (PyErr_GivenExceptionMatches(err, PyExc_TypeError)) {
        throw std::runtime_error( "Generic soft constraint callbacks must take exactly 6 arguments" );
      } else {
        throw std::runtime_error( "Some error occurred while executing generic soft constraint callback" );
      }
    }
    PyErr_Clear();
  } else if (result == Py_None) {
    throw std::runtime_error( "Generic soft constraint callback must return pseudo energy value" );
  } else {
    ret = (int)PyInt_AsLong(result);
  }
  /* END recognizing errors in callback execution */

  Py_DECREF(arglist);
  Py_XDECREF(result);

  return ret;
}

static vrna_basepair_t *
py_wrap_sc_bt_callback( int i,
                        int j,
                        int k,
                        int l,
                        unsigned char d,
                        void *data){

  int c, len, num_pairs;
  PyObject *func, *arglist, *result, *bp, *err;
  py_sc_callback_t *cb = (py_sc_callback_t *)data;
  vrna_basepair_t *ptr, *pairs = NULL;
  func = cb->cb_bt;
  /* compose argument list */
  arglist = Py_BuildValue("(i,i,i,i,i,O)", i, j, k, l, (int)d, (cb->data) ? cb->data : Py_None);
  result =  PyObject_CallObject(func, arglist);

  /* BEGIN recognizing errors in callback execution */
  if (result == NULL) {
    if ((err = PyErr_Occurred())) {
      /* print error message */
      PyErr_Print();
      /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
      if (PyErr_GivenExceptionMatches(err, PyExc_TypeError)) {
        throw std::runtime_error( "Generic soft constraint callbacks must take exactly 6 arguments" );
      } else {
        throw std::runtime_error( "Some error occurred while executing generic soft constraint callback" );
      }
    }
    PyErr_Clear();
    return NULL;
  }
  /* END recognizing errors in callback execution */

  if(PyList_Check(result)){
    len       = 10;
    num_pairs = 0;
    pairs     = (vrna_basepair_t *)vrna_alloc(sizeof(vrna_basepair_t) * len);
  
    /* result should be list of pairs */
    for(c=0; c < PyList_Size(result); c++){
      bp = PyList_GetItem(result, c);
      /* maybe the user was so kind to create a list of vrna_basepair_t? */
      if(SWIG_ConvertPtr(bp, (void **) &ptr, SWIGTYPE_p_vrna_basepair_t, SWIG_POINTER_EXCEPTION) == 0){
        pairs[num_pairs] = *ptr;
        num_pairs++;
      }
      /* users may also specify base pairs as tuples of size 2 */
      else if(PyTuple_Check(bp)){
        if(   (PyTuple_Size(bp) == 2)
          &&  PyLong_Check(PyTuple_GetItem(bp, 0))
          &&  PyLong_Check(PyTuple_GetItem(bp, 1))){
          pairs[num_pairs].i = (int)PyLong_AsLong(PyTuple_GetItem(bp, 0));
          pairs[num_pairs].j = (int)PyLong_AsLong(PyTuple_GetItem(bp, 1));
          num_pairs++;
        }
      }
      /* or is it even a dictionary with i j keys? */
      else if(PyDict_Check(bp)){
        /* check whether the dictionary actually contains the correct keys */
        PyObject *bp_i, *bp_j;
        bp_i = PyDict_GetItemString(bp, "i");
        bp_j = PyDict_GetItemString(bp, "j");
        /* both dictionary keys must be present and the corresponding values have to be integer types */
        if(bp_i && bp_j && PyLong_Check(bp_i) && PyLong_Check(bp_j)){
          pairs[num_pairs].i = (int)PyLong_AsLong(bp_i);
          pairs[num_pairs].j = (int)PyLong_AsLong(bp_j);
          num_pairs++;
        }
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
                          unsigned char d,
                          void *data){

  FLT_OR_DBL ret;
  PyObject *func, *arglist, *result, *err;
  py_sc_callback_t *cb = (py_sc_callback_t *)data;

  func = cb->cb_exp_f;
  /* compose argument list */
  arglist = Py_BuildValue("(i,i,i,i,i,O)", i, j, k, l, (int)d, (cb->data) ? cb->data : Py_None);
  result =  PyObject_CallObject(func, arglist);

  /* BEGIN recognizing errors in callback execution */
  if (result == NULL) {
    if ((err = PyErr_Occurred())) {
      /* print error message */
      PyErr_Print();
      /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
      if (PyErr_GivenExceptionMatches(err, PyExc_TypeError)) {
        throw std::runtime_error( "Generic soft constraint callbacks (partition function) must take exactly 6 arguments" );
      } else {
        throw std::runtime_error( "Some error occurred while executing generic soft constraint callback (partition function)" );
      }
    }
    PyErr_Clear();
  } else if (result == Py_None) {
    throw std::runtime_error( "Generic soft constraint callback (partition function) must return Boltzmann weighted pseudo energy value" );
  } else {
    ret = (FLT_OR_DBL)PyFloat_AsDouble(result);
  }
  /* END recognizing errors in callback execution */

  Py_DECREF(arglist);
  Py_XDECREF(result);
  return ret;
}

%}

static void sc_add_f_pycallback(vrna_fold_compound_t *vc, PyObject *PyFunc);
static void sc_add_bt_pycallback(vrna_fold_compound_t *vc, PyObject *PyFunc);
static void sc_add_exp_f_pycallback(vrna_fold_compound_t *vc, PyObject *PyFunc);
static void sc_add_pydata(vrna_fold_compound_t *vc, PyObject *data, PyObject *PyFuncOrNone);

/* now we bind the above functions as methods to the fold_compound object */
%extend vrna_fold_compound_t {

  PyObject *sc_add_data(PyObject *data, PyObject *PyFuncOrNone=Py_None){
    sc_add_pydata($self, data, PyFuncOrNone);
    Py_RETURN_NONE;
  }
  
  PyObject *sc_add_f(PyObject *PyFunc){
    sc_add_f_pycallback($self, PyFunc);
    Py_RETURN_NONE;
  }

  PyObject *sc_add_bt(PyObject *PyFunc){
    sc_add_bt_pycallback($self, PyFunc);
    Py_RETURN_NONE;
  }

  PyObject *sc_add_exp_f(PyObject *PyFunc){
    sc_add_exp_f_pycallback($self, PyFunc);
    Py_RETURN_NONE;
  }
}

#endif
