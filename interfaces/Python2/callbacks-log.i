/**********************************************/
/* BEGIN interface for logging callback       */
/**********************************************/

#ifdef SWIGPYTHON
%{

#include <stdexcept>

typedef struct {
  PyObject  *cb;
  PyObject  *data;
  PyObject  *delete_data;
} pycallback_log_t;

static void
py_wrap_log_callback(vrna_log_event_t  *event,
                     void              *log_data);

static void
delete_py_log(pycallback_log_t *cb)
{
  if(cb->data != Py_None){
    if(cb->delete_data != Py_None){
      /* call user-defined data destructor */
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
            throw std::runtime_error( "log_cb delete_data() callback must take exactly 1 argument" );
          } else {
            throw std::runtime_error( "Some error occurred while executing log_cb delete_data() callback" );
          }
        }
        PyErr_Clear();
      }
      /* END recognizing errors in callback execution */

      Py_DECREF(arglist);
      Py_XDECREF(result);
    }
  }

  Py_DECREF(cb->data);
  Py_DECREF(cb->delete_data);
}


static void
delete_pycallback_log(void * data)
{
  pycallback_log_t *cb = (pycallback_log_t *)data;
  /* first delete user data */
  delete_py_log(cb);

  /* now dispose of the callback */
  Py_DECREF(cb->cb);
  
  /* finally free pycallback */
  free(cb);
}


static unsigned int
log_cb_add_pycallback(PyObject  *PyFunc,
              PyObject          *data,
              PyObject          *data_release,
              vrna_log_levels_e level)
{
  unsigned int r;

  pycallback_log_t *cb;

  cb = (pycallback_log_t *)vrna_alloc(sizeof(pycallback_log_t));

  cb->cb          = PyFunc;    /* remember callback */
  cb->data        = data;
  cb->delete_data = data_release;

  Py_INCREF(data);
  Py_XINCREF(PyFunc);       /* Increase referenc counter */
  Py_XINCREF(data_release); /* Increase referenc counter */

  /* finaly bind to logging system */
  r = vrna_log_cb_add(&py_wrap_log_callback,
                      (void *)cb,
                      &delete_pycallback_log,
                      level);

  return r;
}


static void
py_wrap_log_callback(vrna_log_event_t  *event,
                     void              *data)
{

  PyObject *func, *arglist, *result, *err;
  pycallback_log_t *cb = (pycallback_log_t *)data;

  func = cb->cb;

  /* compose argument list */
  char      *message  = vrna_strdup_vprintf(event->format_string, event->params);
  char      *fn       = vrna_strdup_printf(event->file_name);

  arglist = Py_BuildValue("{s:s,s:i,s:i,s:s},O",
                "message", message,
                "level", event->level,
                "line_number", event->line_number,
                "file_name", fn,
                (cb->data) ? cb->data : Py_None);

  result =  PyObject_CallObject(func, arglist);

  free(message);
  free(fn);

  /* BEGIN recognizing errors in callback execution */
  if (result == NULL) {
    if ((err = PyErr_Occurred())) {
      /* print error message */
      PyErr_Print();
      /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
      if (PyErr_GivenExceptionMatches(err, PyExc_TypeError)) {
        throw std::runtime_error( "logging callback must take exactly 2 arguments" );
      } else {
        throw std::runtime_error( "Some error occurred while executing logging callback" );
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

static unsigned int
log_cb_add_pycallback(PyObject  *PyFunc,
              PyObject          *data,
              PyObject          *PyFuncOrNone,
              vrna_log_levels_e level);

%feature("autodoc") log_cb_add;
%feature("kwargs") log_cb_add;

%apply PyObject *PyFuncOrNone { PyObject *data_release };
%apply PyObject *PyFunc { PyObject *callback };

/* now we bind the above functions as methods to the fold_compound object */
%{
  unsigned int
  log_cb_add(PyObject *callback,
             PyObject *data           = Py_None,
             PyObject *data_release   = Py_None,
             vrna_log_levels_e level  = VRNA_LOG_LEVEL_WARNING)
  {
    return log_cb_add_pycallback(callback, data, data_release, level);
  }

%}

unsigned int
log_cb_add(PyObject *callback,
           PyObject *data           = Py_None,
           PyObject *data_release   = Py_None,
           vrna_log_levels_e level  = VRNA_LOG_LEVEL_WARNING);


%clear  PyObject *data_release;
%clear  PyObject *callback;

#endif
