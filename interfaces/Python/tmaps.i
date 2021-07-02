// convert between python and C file handle
%include "file_py3.i" // python 3 FILE *

%typemap(out) float [ANY] {
  int i;
  $result = PyList_New($1_dim0);
  for (i = 0; i < $1_dim0; i++) {
    PyObject *o = PyFloat_FromDouble((double) $1[i]);
    PyList_SetItem($result,i,o);
  }
}

%typemap(out) int [ANY] {
  int i;
  $result = PyList_New($1_dim0);
  for (i = 0; i < $1_dim0; i++) {
    PyObject *o = PyLong_FromLong((long) $1[i]);
    PyList_SetItem($result,i,o);
  }
}

/* for global variables */
%typemap(varout) int [ANY][ANY] {
  size_t i, j;
  //$1, $1_dim0, $1_dim1
  $result = PyList_New($1_dim0);
  for (i = 0; i < $1_dim0; i++) {
    PyObject *l = PyList_New($1_dim1);
    for (j = 0; j < $1_dim1; j++) {
      PyObject *o = PyLong_FromLong((long) $1[i][j]);
      PyList_SetItem(l, j, o);
    }
    PyList_SetItem($result,i,l);
  }
}

/* for structure members */
%typemap(out) int [ANY][ANY] {
  size_t i, j;
  //$1, $1_dim0, $1_dim1
  $result = PyList_New($1_dim0);
  for (i = 0; i < $1_dim0; i++) {
    PyObject *l = PyList_New($1_dim1);
    for (j = 0; j < $1_dim1; j++) {
      PyObject *o = PyLong_FromLong((long) $1[i][j]);
      PyList_SetItem(l, j, o);
    }
    PyList_SetItem($result,i,l);
  }
}

%apply int [ANY][ANY] { const int[ANY][ANY], long int [ANY][ANY], const long int [ANY][ANY], short [ANY][ANY], const short [ANY][ANY] };

// This tells SWIG to treat char ** as a special case
%typemap(in) char ** {
  /* Check if is a list */
  if (PyList_Check($input)) {
    int size = PyList_Size($input);
    int i = 0;
    $1 = (char **) malloc((size+1)*sizeof(char *));
    for (i = 0; i < size; i++) {
      PyObject *o = PyList_GetItem($input,i);
      if (PyUnicode_Check(o))
        $1[i] = PyString_AsString(PyUnicode_AsASCIIString(o));
      else {
        PyErr_SetString(PyExc_TypeError,"list must contain strings");
        free($1);
        return NULL;
      }
    }
    $1[i] = 0;
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}

// This tells SWIG to treat char *[], const char **, and const char *[] the same as char **
%apply char ** { char *[], const char **, const char *[] };

/*
  we need this crazy piece of argout typemap only because we don't
  want the vrna_pbacktrack_mem_t object to be appended to the results(list),
  but prepended instead. Otherwise, a simple
  %append_output(SWIG_NewPointerObj(SWIG_as_voidptr(retval$argnum), $1_descriptor, 0));
  would have sufficed already
*/
%typemap(argout) vrna_pbacktrack_mem_t *INOUT {
  PyObject *o, *o2, *o3;
  o = SWIG_NewPointerObj(SWIG_as_voidptr(retval$argnum), $1_descriptor, 1);
  if ((!$result) || ($result == Py_None)) {
    $result = o;
  } else {
    PyObject *o2 = $result;
    $result = PyTuple_New(1);
    PyTuple_SetItem($result,0,o2);
    o3 = PyTuple_New(1);
    PyTuple_SetItem(o3,0,o);
    o2 = $result;
    $result = PySequence_Concat(o3,o2);
    Py_DECREF(o2);
    Py_DECREF(o3);
  }
}


%typemap(in) vrna_pbacktrack_mem_t *INOUT (vrna_pbacktrack_mem_t *retval)
{
  if ($input == Py_None) {
    retval = new vrna_pbacktrack_mem_t();
    $1 = retval;
  } else {
    /* INOUT in */
    SWIG_ConvertPtr($input,SWIG_as_voidptrptr(&retval), 0, SWIG_POINTER_DISOWN);
    $1 = retval;
  }
}


%typemap(in) PyObject *PyFunc {
  if (!PyCallable_Check($input)) {
      PyErr_SetString(PyExc_TypeError, "Need a callable object!");
      return NULL;
  }
  $1 = $input;
}

%typemap(in) PyObject *PyFuncOrNone {
  if($input != Py_None){
    if (!PyCallable_Check($input)) {
        PyErr_SetString(PyExc_TypeError, "Need a callable object!");
        return NULL;
    }
  }
  $1 = $input;
}


%typemap(in) PyObject * {
  $1 = $input;
}
