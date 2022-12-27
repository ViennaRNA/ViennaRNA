/************************************************/
/* BEGIN interface for landscape, neighborhood, */
/* and moves callbacks                          */
/************************************************/

#ifdef SWIGPYTHON
%{

#include <stdexcept>

typedef struct {
  PyObject  *cb;
  PyObject  *data;
} pycallback_simple_t;


static pycallback_simple_t *
bind_simple_callback(PyObject *PyFunc,
                     PyObject *data);

static void
python_wrap_move_update_cb(vrna_fold_compound_t *fc,
                           vrna_move_t          neighbor,
                           unsigned int         state,
                           void                 *data);

static void
release_simple_callback(pycallback_simple_t *cb);


static pycallback_simple_t *
bind_simple_callback(PyObject *PyFunc,
                     PyObject *data)
{

  pycallback_simple_t *cb = (pycallback_simple_t *)vrna_alloc(sizeof(pycallback_simple_t));

  Py_INCREF(PyFunc);
  Py_INCREF(data);
  cb->cb    = PyFunc;  /* store callback */
  cb->data  = data;    /* bind data */

  return cb;
}

static void
release_simple_callback(pycallback_simple_t *cb)
{
  Py_DECREF(cb->cb);
  Py_DECREF(cb->data);
  free(cb);
}


static void
python_wrap_move_update_cb(vrna_fold_compound_t *fc,
                           vrna_move_t          neighbor,
                           unsigned int         state,
                           void                 *data)
{
  PyObject            *func, *arglist, *result, *err, *py_fc, *py_neighbor, *py_state;
  pycallback_simple_t *cb;

  cb    = (pycallback_simple_t *)data;
  func  = cb->cb;

  /* compose argument list */
  py_fc       = SWIG_NewPointerObj(SWIG_as_voidptr(fc),
                                   SWIGTYPE_p_vrna_fold_compound_t,
                                   SWIG_POINTER_NEW);
  py_neighbor = SWIG_NewPointerObj(SWIG_as_voidptr(&neighbor),
                                   SWIGTYPE_p_vrna_move_t,
                                   SWIG_POINTER_NEW);
  py_state    = PyLong_FromLong((long int)state);
  result      = PyObject_CallFunctionObjArgs(func,
                                           py_fc,
                                           py_neighbor,
                                           py_state,
                                           (cb->data) ? cb->data : Py_None,
                                           NULL);

  Py_DECREF(py_fc);
  Py_DECREF(py_neighbor);
  Py_DECREF(py_state);

  /* BEGIN recognizing errors in callback execution */
  if (result == NULL) {
    if ((err = PyErr_Occurred())) {
      /* print error message */
      PyErr_Print();
      /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
      if (PyErr_GivenExceptionMatches(err, PyExc_TypeError)) {
        throw std::runtime_error( "move_neighbor_diff callback must take exactly 4 arguments" );
      } else {
        throw std::runtime_error( "Some error occurred while executing move_neighbor_diff callback" );
      }
    }
    PyErr_Clear();
  }
  /* END recognizing errors in callback execution */

  Py_XDECREF(result);

  return /*void*/;
}

%}

/* now we bind the above functions as methods to the fold_compound object */
%extend vrna_fold_compound_t {

%feature("autodoc") move_neighbor_diff;
%feature("kwargs") move_neighbor_diff;

  int
  move_neighbor_diff(var_array<short> &pt,
                     vrna_move_t      *move,
                     PyObject         *PyFunc,
                     PyObject         *data = Py_None,
                     unsigned int     options = VRNA_MOVESET_DEFAULT)
  {
    int                 ret = 0;
    pycallback_simple_t *cb = bind_simple_callback(PyFunc,
                                                   data);

    ret = vrna_move_neighbor_diff_cb($self,
                                     pt.data,
                                     *move,
                                     &python_wrap_move_update_cb,
                                     (void *)cb,
                                     options);

    release_simple_callback(cb);
    return ret;
  }
}

#endif
