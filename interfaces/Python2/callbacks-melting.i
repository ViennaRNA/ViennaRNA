/**********************************************/
/* BEGIN interface for heat_capacity callback */
/**********************************************/

#ifdef SWIGPYTHON
%{

#include <stdexcept>

typedef struct {
  PyObject *cb;
  PyObject *data;
} python_heat_capacity_callback_t;

static python_heat_capacity_callback_t *
bind_heat_capacity_callback(PyObject *PyFunc,
                            PyObject *data);

static void
python_wrap_heat_capacity_cb(float      temp,
                             float      hc,
                             void       *data);

static python_heat_capacity_callback_t *
bind_heat_capacity_callback(PyObject *PyFunc,
                            PyObject *data)
{

  python_heat_capacity_callback_t *cb = (python_heat_capacity_callback_t *)vrna_alloc(sizeof(python_heat_capacity_callback_t));

  Py_INCREF(PyFunc);
  Py_INCREF(data);
  cb->cb    = PyFunc;  /* store callback */
  cb->data  = data;    /* bind data */

  return cb;
}

static void
release_heat_capacity_callback(python_heat_capacity_callback_t *cb)
{
  Py_DECREF(cb->cb);
  Py_DECREF(cb->data);
  free(cb);
}


static void
python_wrap_heat_capacity_cb(float temp,
                             float hc,
                             void  *data)
{
  PyObject                        *func, *arglist, *result, *err, *py_temp, *py_hc;
  python_heat_capacity_callback_t *cb;

  cb    = (python_heat_capacity_callback_t *)data;
  func  = cb->cb;

  /* compose argument list */
  py_temp = PyFloat_FromDouble((double)temp);
  py_hc   = PyFloat_FromDouble((double)hc);
  result  = PyObject_CallFunctionObjArgs(func,
                                         py_temp,
                                         py_hc,
                                         (cb->data) ? cb->data : Py_None,
                                         NULL);

  Py_DECREF(py_temp);
  Py_DECREF(py_hc);

  /* BEGIN recognizing errors in callback execution */
  if (result == NULL) {
    if ((err = PyErr_Occurred())) {
      /* print error message */
      PyErr_Print();
      /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
      if (PyErr_GivenExceptionMatches(err, PyExc_TypeError)) {
        throw std::runtime_error( "heat_capacity callback must take exactly 3 arguments" );
      } else {
        throw std::runtime_error( "Some error occurred while executing heat_capacity callback" );
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

%feature("autodoc") heat_capacity_cb;
%feature("kwargs") heat_capacity_cb;

  PyObject *
  heat_capacity_cb(float        T_min,
                   float        T_max,
                   float        T_increment,
                   unsigned int mpoints,
                   PyObject     *PyFunc,
                   PyObject     *data = Py_None)
  {
    python_heat_capacity_callback_t *cb = bind_heat_capacity_callback(PyFunc, data);
    vrna_heat_capacity_cb($self, T_min, T_max, T_increment, mpoints, &python_wrap_heat_capacity_cb, (void *)cb);
    release_heat_capacity_callback(cb);
    Py_RETURN_NONE;
  }

}


#endif
