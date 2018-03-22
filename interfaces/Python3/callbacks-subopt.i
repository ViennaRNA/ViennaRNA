/**********************************************/
/* BEGIN interface for subopt callback        */
/**********************************************/

#ifdef SWIGPYTHON
%{

#include <stdexcept>

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

  PyObject *func, *arglist, *result, *err;
  python_subopt_callback_t *cb = (python_subopt_callback_t *)data;

  func = cb->cb;
  /* compose argument list */
  arglist = Py_BuildValue("(z,d,O)", structure, (double)energy, (cb->data) ? cb->data : Py_None);
  result =  PyObject_CallObject(func, arglist);

  /* BEGIN recognizing errors in callback execution */
  if (result == NULL) {
    if ((err = PyErr_Occurred())) {
      /* print error message */
      PyErr_Print();
      /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
      if (PyErr_GivenExceptionMatches(err, PyExc_TypeError)) {
        throw std::runtime_error( "Subopt callback must take exactly 3 arguments" );
      } else {
        throw std::runtime_error( "Some error occurred while executing subopt callback" );
      }
    }
    PyErr_Clear();
  }
  /* END recognizing errors in callback execution */

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
