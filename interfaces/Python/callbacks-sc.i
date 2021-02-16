/*
 * ********************************************
 * BEGIN interface for fold compound status
 * callback
 *********************************************
 */

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

  static int
  py_wrap_sc_f_callback(int           i,
                        int           j,
                        int           k,
                        int           l,
                        unsigned char d,
                        void          *data);


  static vrna_basepair_t *
  py_wrap_sc_bt_callback(int            i,
                         int            j,
                         int            k,
                         int            l,
                         unsigned char  d,
                         void           *data);


  static FLT_OR_DBL
  py_wrap_sc_exp_f_callback(int           i,
                            int           j,
                            int           k,
                            int           l,
                            unsigned char d,
                            void          *data);


  static void
  delete_py_sc_data(py_sc_callback_t *cb)
  {
    if ((cb->data != Py_None) &&
        (cb->delete_data != Py_None)) {
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
          if (PyErr_GivenExceptionMatches(err, PyExc_TypeError))
            throw
            std::runtime_error(
              "Generic soft constraint delete_data() callback must take exactly 1 argument");


          else
            throw
            std::runtime_error(
              "Some error occurred while executing generic soft constraint delete_data() callback");
        }

        PyErr_Clear();
      }

      /* END recognizing errors in callback execution */

      Py_DECREF(arglist);
      Py_XDECREF(result);
    }

    Py_DECREF(cb->data);
    Py_DECREF(cb->delete_data);
  }


  static void
  delete_py_sc_callback(void *data)
  {
    py_sc_callback_t *cb = (py_sc_callback_t *)data;

    /* first delete user data */
    delete_py_sc_data(cb);

    /* now dispose of the registered callbacks */
    Py_DECREF(cb->cb_f);
    Py_DECREF(cb->cb_bt);
    Py_DECREF(cb->cb_exp_f);

    /* finally free pycallback */
    free(cb);
  }


  static py_sc_callback_t *
  reuse_or_new_cb_f(vrna_sc_t *sc)
  {
    py_sc_callback_t *cb;

    cb =
      (sc->data) ? (py_sc_callback_t *)sc->data : (py_sc_callback_t *)vrna_alloc(
        sizeof(py_sc_callback_t));

    if (cb->cb_f) {
      /* release previous callback */
      Py_DECREF(cb->cb_f);
    } else {
      /* initialize the remaining soft constraint callbacks */
      Py_INCREF(Py_None);
      cb->cb_bt = Py_None;

      Py_INCREF(Py_None);
      cb->cb_exp_f = Py_None;

      Py_INCREF(Py_None);
      cb->data = Py_None;

      Py_INCREF(Py_None);
      cb->delete_data = Py_None;
    }

    return cb;
  }


  static py_sc_callback_t *
  reuse_or_new_cb_exp_f(vrna_sc_t *sc)
  {
    py_sc_callback_t *cb;

    cb =
      (sc->data) ? (py_sc_callback_t *)sc->data : (py_sc_callback_t *)vrna_alloc(
        sizeof(py_sc_callback_t));

    if (cb->cb_exp_f) {
      /* release previous callback */
      Py_DECREF(cb->cb_exp_f);
    } else {
      /* initialize the remaining soft constraint callbacks */
      Py_INCREF(Py_None);
      cb->cb_f = Py_None;

      Py_INCREF(Py_None);
      cb->cb_bt = Py_None;

      Py_INCREF(Py_None);
      cb->data = Py_None;

      Py_INCREF(Py_None);
      cb->delete_data = Py_None;
    }

    return cb;
  }


  static py_sc_callback_t *
  reuse_or_new_cb_data(vrna_sc_t *sc)
  {
    py_sc_callback_t *cb;

    cb =
      (sc->data) ? (py_sc_callback_t *)sc->data : (py_sc_callback_t *)vrna_alloc(
        sizeof(py_sc_callback_t));

    if (cb->data) {
      /* release previous callback */
      delete_py_sc_data(cb);
    } else {
      /* initialize the remaining soft constraint callbacks */
      Py_INCREF(Py_None);
      cb->cb_f = Py_None;

      Py_INCREF(Py_None);
      cb->cb_bt = Py_None;

      Py_INCREF(Py_None);
      cb->cb_exp_f = Py_None;
    }

    return cb;
  }


  static int
  sc_add_f_pycallback(vrna_fold_compound_t  *vc,
                      PyObject              *func)
  {
    unsigned char     func_is_tuple, func_is_list;
    unsigned int      s;
    /* try to dispose of previous callback */
    py_sc_callback_t  *cb;
    PyObject          *f, *err;

    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (!PyCallable_Check(func)) {
          PyErr_SetString(PyExc_TypeError, "Need a callable object!");
        } else if (vrna_sc_add_f(vc, &py_wrap_sc_f_callback)) {
          /*
           *  The above call returns 0 on any error.
           *  Otherwise it binds the wrapper function and
           *  prepares the soft constraint data structure
           *  inside vc
           */

          /* now bind the python function to the wrapper structure */
          cb = reuse_or_new_cb_f(vc->sc);

          Py_INCREF(func);  /* Increase referenc counter */
          cb->cb_f = func;  /* remember callback */

          /* finaly bind callback wrapper to fold compound */
          vc->sc->data = (void *)cb;

          /* also (re-)bind the free-callback-data function */
          vc->sc->free_data = &delete_py_sc_callback;

          return 1;
        }

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        /* first check, whether data and PyFunc are of correct type */
        func_is_tuple = PyTuple_Check(func);
        func_is_list  = PyList_Check(func);

        if (func_is_tuple || func_is_list) {
          if (!vc->scs)
            vrna_sc_init(vc);

          for (s = 0; s < vc->n_seq; s++) {
            f = (func_is_tuple) ? PyTuple_GetItem(func, s) : PyList_GetItem(func, s);

            if (f == NULL) {
              /* an error occurred */
              if ((err = PyErr_Occurred())) {
                /* print error message */
                PyErr_Print();
                /* we only treat IndexErrors differently here, as they indicate that the callback list is too short! */
                if (PyErr_GivenExceptionMatches(err, PyExc_IndexError))
                  throw
                  std::runtime_error(
                    "sc_add_f(): Comparative prediction callback list or tuple must have an entry for each sequence in the alignment");


                else
                  throw
                  std::runtime_error(
                    "sc_add_f(): Some error occurred while accessing generic soft constraint callback for sequence alignment");
              }

              PyErr_Clear();
            } else if (!PyCallable_Check(f)) {
              PyErr_SetString(PyExc_TypeError, "Need a callable object!");
            } else {
              cb = reuse_or_new_cb_f(vc->scs[s]);

              Py_INCREF(f); /* Increase reference counter */
              cb->cb_f = f; /* remember callback */

              /* finaly bind callback wrapper to fold compound */
              vc->scs[s]->data = (void *)cb;

              /* also (re-)bind the free-callback-data function */
              vc->scs[s]->free_data = &delete_py_sc_callback;

              /* finally, we bind the callback wrapper */
              vc->scs[s]->f = &py_wrap_sc_f_callback;
            }
          }

          return 1;
        } else {
          throw
          std::runtime_error(
            "sc_add_f(): Comparative prediction callbacks must be provided as list or tuple");
        }

        break;
    }

    return 0;
  }


  static int
  sc_add_exp_f_pycallback(vrna_fold_compound_t  *vc,
                          PyObject              *func)
  {
    unsigned char     func_is_tuple, func_is_list;
    unsigned int      s;
    /* try to dispose of previous callback */
    py_sc_callback_t  *cb;
    PyObject          *f, *err;

    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (!PyCallable_Check(func)) {
          PyErr_SetString(PyExc_TypeError, "Need a callable object!");
        } else if (vrna_sc_add_exp_f(vc, &py_wrap_sc_exp_f_callback)) {
          /*
           *  The above call returns 0 on any error.
           *  Otherwise it binds the wrapper function and
           *  prepares the soft constraint data structure
           *  inside vc
           */

          /* now bind the python function to the wrapper structure */
          cb = reuse_or_new_cb_exp_f(vc->sc);

          Py_INCREF(func);      /* Increase referenc counter */
          cb->cb_exp_f = func;  /* remember callback */

          /* finaly bind callback wrapper to fold compound */
          vc->sc->data = (void *)cb;

          /* also (re-)bind the free-callback-data function */
          vc->sc->free_data = &delete_py_sc_callback;

          return 1;
        }

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        /* first check, whether data and PyFunc are of correct type */
        func_is_tuple = PyTuple_Check(func);
        func_is_list  = PyList_Check(func);

        if (func_is_tuple || func_is_list) {
          if (!vc->scs)
            vrna_sc_init(vc);

          for (s = 0; s < vc->n_seq; s++) {
            f = (func_is_tuple) ? PyTuple_GetItem(func, s) : PyList_GetItem(func, s);

            if (f == NULL) {
              /* an error occurred */
              if ((err = PyErr_Occurred())) {
                /* print error message */
                PyErr_Print();
                /* we only treat IndexErrors differently here, as they indicate that the callback list is too short! */
                if (PyErr_GivenExceptionMatches(err, PyExc_IndexError))
                  throw
                  std::runtime_error(
                    "sc_add_exp_f(): Comparative prediction callback list or tuple must have an entry for each sequence in the alignment");


                else
                  throw
                  std::runtime_error(
                    "sc_add_exp_f(): Some error occurred while accessing generic soft constraint callback for sequence alignment");
              }

              PyErr_Clear();
            } else if (!PyCallable_Check(f)) {
              PyErr_SetString(PyExc_TypeError, "Need a callable object!");
            } else {
              cb = reuse_or_new_cb_exp_f(vc->scs[s]);

              Py_INCREF(f);     /* Increase reference counter */
              cb->cb_exp_f = f; /* remember callback */

              /* finaly bind callback wrapper to fold compound */
              vc->scs[s]->data = (void *)cb;

              /* also (re-)bind the free-callback-data function */
              vc->scs[s]->free_data = &delete_py_sc_callback;

              /* finally, we bind the callback wrapper */
              vc->scs[s]->f = &py_wrap_sc_f_callback;
            }
          }

          return 1;
        } else {
          throw
          std::runtime_error(
            "sc_add_exp_f(): Comparative prediction callbacks must be provided as list or tuple");
        }

        break;
    }

    return 0;
  }


  static int
  sc_add_bt_pycallback(vrna_fold_compound_t *vc,
                       PyObject             *PyFunc)
  {
    /* try to dispose of previous callback */
    py_sc_callback_t *cb;

    if (vrna_sc_add_bt(vc, &py_wrap_sc_bt_callback)) {
      /*
       *  The above call returns 0 on any error.
       *  Otherwise it binds the wrapper function and
       *  prepares the soft constraint data structure
       *  inside vc
       */

      /* now bind the python function to the wrapper structure */
      if (vc->sc->data) {
        /* re-use previously bound wrapper data structure */
        cb = (py_sc_callback_t *)vc->sc->data;
        /* release previous callback */
        Py_DECREF(cb->cb_bt);
      } else {
        cb = (py_sc_callback_t *)vrna_alloc(sizeof(py_sc_callback_t));
        Py_INCREF(Py_None);
        cb->cb_f = Py_None;
        Py_INCREF(Py_None);
        cb->cb_exp_f = Py_None;
        Py_INCREF(Py_None);
        cb->data = Py_None;
        Py_INCREF(Py_None);
        cb->delete_data = Py_None;
      }

      Py_XINCREF(PyFunc); /* Increase referenc counter */
      cb->cb_bt = PyFunc; /* remember callback */

      /* finaly bind callback wrapper to fold compound */
      vc->sc->data = (void *)cb;

      /* also (re-)bind the free-callback-data function */
      vc->sc->free_data = &delete_py_sc_callback;

      return 1;
    }

    return 0;
  }


  static int
  sc_add_pydata(vrna_fold_compound_t  *vc,
                PyObject              *data,
                PyObject              *free_data_cb)
  {
    unsigned char     data_is_tuple, data_is_list, func_is_tuple, func_is_list;
    unsigned int      s;
    py_sc_callback_t  *cb;
    PyObject          *err;

    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        /* create soft constraints data structure */
        if (!vc->sc)
          vrna_sc_init(vc);

        cb = reuse_or_new_cb_data(vc->sc);

        /* increase reference counters */
        Py_INCREF(data);
        cb->data = data;        /* remember data */

        if ((free_data_cb != Py_None) && (!PyCallable_Check(free_data_cb))) {
          PyErr_SetString(PyExc_TypeError, "Require a callable object for free_data_cb!");
        } else {
          Py_INCREF(free_data_cb);
          cb->delete_data = free_data_cb; /* remember delete data function */
        }

        vc->sc->data = (void *)cb;

        /* also (re-)bind the free-callback-data function */
        vc->sc->free_data = &delete_py_sc_callback;

        return 1;

      case VRNA_FC_TYPE_COMPARATIVE:
        /* first check, whether data and PyFunc are of correct type */
        data_is_tuple = PyTuple_Check(data);
        data_is_list  = PyList_Check(data);
        func_is_tuple = PyTuple_Check(free_data_cb);
        func_is_list  = PyList_Check(free_data_cb);

        if (data_is_tuple || data_is_list) {
          if (!vc->scs)
            vrna_sc_init(vc);

          /* try to dispose of previous data */
          for (s = 0; s < vc->n_seq; s++) {
            cb = reuse_or_new_cb_data(vc->scs[s]);

            /* increase reference counters */
            PyObject *d, *f;

            d = (data_is_tuple) ? PyTuple_GetItem(data, s) : PyList_GetItem(data, s);

            if (d == NULL) {
              /* an error occurred */
              if ((err = PyErr_Occurred())) {
                /* print error message */
                PyErr_Print();
                /* we only treat IndexErrors differently here, as they indicate that the callback list is too short! */
                if (PyErr_GivenExceptionMatches(err, PyExc_IndexError))
                  throw
                  std::runtime_error(
                    "sc_add_data(): Comparative prediction callback data list or tuple must have an entry for each sequence in the alignment");


                else
                  throw
                  std::runtime_error(
                    "sc_add_data(): Some error occurred while accessing generic soft constraint callback data list for sequence alignment");
              }

              PyErr_Clear();
            } else {
              Py_INCREF(d);
              cb->data = d; /* remember data */

              if (func_is_tuple || func_is_list) {
                f = (func_is_tuple) ? PyTuple_GetItem(free_data_cb, s) : PyList_GetItem(
                  free_data_cb,
                  s);
                if (f == NULL) {
                  /* an error occurred */
                  if ((err = PyErr_Occurred())) {
                    /* print error message */
                    PyErr_Print();
                    /* we only treat IndexErrors differently here, as they indicate that the callback list is too short! */
                    if (PyErr_GivenExceptionMatches(err, PyExc_IndexError))
                      throw
                      std::runtime_error(
                        "sc_add_data(): Comparative prediction callback data free() list or tuple must have an entry for each sequence in the alignment");


                    else
                      throw
                      std::runtime_error(
                        "sc_add_data(): Some error occurred while accessing generic soft constraint callback data free() list for sequence alignment");
                  }

                  PyErr_Clear();
                } else if ((f != Py_None) && (!PyCallable_Check(f))) {
                  PyErr_SetString(PyExc_TypeError, "Require a callable object for free_data_cb!");
                  cb->delete_data = Py_None;
                  Py_INCREF(Py_None);
                } else {
                  cb->delete_data = f; /* remember delete data function */
                  Py_INCREF(f);
                }
              } else {
                cb->delete_data = Py_None;
                Py_INCREF(Py_None);
              }
            }

            vc->scs[s]->data = (void *)cb;

            /* also (re-)bind the free-callback-data function */
            vc->scs[s]->free_data = &delete_py_sc_callback;
          }

          return 1;
        } else {
          throw
          std::runtime_error(
            "sc_add_data(): Comparative prediction data must be provided as list or tuple");
        }

        break;
    }

    return 0;
  }


  static int
  py_wrap_sc_f_callback(int           i,
                        int           j,
                        int           k,
                        int           l,
                        unsigned char d,
                        void          *data)
  {
    int               ret;
    PyObject          *func, *arglist, *result, *err;
    py_sc_callback_t  *cb = (py_sc_callback_t *)data;

    ret   = 0;
    func  = cb->cb_f;

    /* compose argument list */
    PyObject *py_i, *py_j, *py_k, *py_l, *py_d;

    py_i    = PyInt_FromLong(i);
    py_j    = PyInt_FromLong(j);
    py_k    = PyInt_FromLong(k);
    py_l    = PyInt_FromLong(l);
    py_d    = PyInt_FromLong(d);
    result  = PyObject_CallFunctionObjArgs(func,
                                           py_i,
                                           py_j,
                                           py_k,
                                           py_l,
                                           py_d,
                                           (cb->data) ? cb->data : Py_None,
                                           NULL);

    Py_DECREF(py_i);
    Py_DECREF(py_j);
    Py_DECREF(py_k);
    Py_DECREF(py_l);
    Py_DECREF(py_d);

    /* BEGIN recognizing errors in callback execution */
    if (result == NULL) {
      if ((err = PyErr_Occurred())) {
        /* print error message */
        PyErr_Print();
        /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
        if (PyErr_GivenExceptionMatches(err, PyExc_TypeError))
          throw
          std::runtime_error("Generic soft constraint callbacks must take exactly 6 arguments");


        else
          throw
          std::runtime_error(
            "Some error occurred while executing generic soft constraint callback");
      }

      PyErr_Clear();
    } else if (PyInt_Check(result)) {
      ret = (int)PyInt_AsLong(result);
    } else if (PyLong_Check(result)) {
      ret = (int)PyLong_AsLong(result);
    } else {
      throw
      std::runtime_error(
        "Generic soft constraint callback must return pseudo energy value in 10 cal/mol");
    }

    /* END recognizing errors in callback execution */

    Py_XDECREF(result);

    return ret;
  }


  static vrna_basepair_t *
  py_wrap_sc_bt_callback(int            i,
                         int            j,
                         int            k,
                         int            l,
                         unsigned char  d,
                         void           *data)
  {
    int               c, len, num_pairs;
    PyObject          *func, *arglist, *result, *bp, *err;
    py_sc_callback_t  *cb;
    vrna_basepair_t   *ptr, *pairs;

    cb    = (py_sc_callback_t *)data;
    pairs = NULL;
    func  = cb->cb_bt;
    /* compose argument list */
    PyObject *py_i, *py_j, *py_k, *py_l, *py_d;

    py_i    = PyInt_FromLong(i);
    py_j    = PyInt_FromLong(j);
    py_k    = PyInt_FromLong(k);
    py_l    = PyInt_FromLong(l);
    py_d    = PyInt_FromLong(d);
    result  = PyObject_CallFunctionObjArgs(func,
                                           py_i,
                                           py_j,
                                           py_k,
                                           py_l,
                                           py_d,
                                           (cb->data) ? cb->data : Py_None,
                                           NULL);

    Py_DECREF(py_i);
    Py_DECREF(py_j);
    Py_DECREF(py_k);
    Py_DECREF(py_l);
    Py_DECREF(py_d);

    /* BEGIN recognizing errors in callback execution */
    if (result == NULL) {
      if ((err = PyErr_Occurred())) {
        /* print error message */
        PyErr_Print();
        /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
        if (PyErr_GivenExceptionMatches(err, PyExc_TypeError))
          throw
          std::runtime_error("Generic soft constraint callbacks must take exactly 6 arguments");


        else
          throw
          std::runtime_error(
            "Some error occurred while executing generic soft constraint callback");
      }

      PyErr_Clear();
      return NULL;
    }

    /* END recognizing errors in callback execution */

    if (PyList_Check(result)) {
      len       = 10;
      num_pairs = 0;
      pairs     = (vrna_basepair_t *)vrna_alloc(sizeof(vrna_basepair_t) * len);

      /* result should be list of pairs */
      for (c = 0; c < PyList_Size(result); c++) {
        bp = PyList_GetItem(result, c);
        /* maybe the user was so kind to create a list of vrna_basepair_t? */
        if (SWIG_ConvertPtr(bp, (void **)&ptr, SWIGTYPE_p_vrna_basepair_t,
                            SWIG_POINTER_EXCEPTION) == 0) {
          pairs[num_pairs] = *ptr;
          num_pairs++;
        }
        /* users may also specify base pairs as tuples of size 2 */
        else if (PyTuple_Check(bp)) {
          if ((PyTuple_Size(bp) == 2)
              && PyInt_Check(PyTuple_GetItem(bp, 0))
              && PyInt_Check(PyTuple_GetItem(bp, 1))) {
            pairs[num_pairs].i  = (int)PyInt_AsLong(PyTuple_GetItem(bp, 0));
            pairs[num_pairs].j  = (int)PyInt_AsLong(PyTuple_GetItem(bp, 1));
            num_pairs++;
          }
        }
        /* or is it even a dictionary with i j keys? */
        else if (PyDict_Check(bp)) {
          /* check whether the dictionary actually contains the correct keys */
          PyObject *bp_i, *bp_j;
          bp_i  = PyDict_GetItemString(bp, "i");
          bp_j  = PyDict_GetItemString(bp, "j");
          /* both dictionary keys must be present and the corresponding values have to be integer types */
          if (bp_i && bp_j && PyInt_Check(bp_i) && PyInt_Check(bp_j)) {
            pairs[num_pairs].i  = (int)PyInt_AsLong(bp_i);
            pairs[num_pairs].j  = (int)PyInt_AsLong(bp_j);
            num_pairs++;
          }
        } else {
          continue;
        }

        /* increase length if necessary */
        if (num_pairs == len) {
          len   = (int)(1.2 * len);
          pairs = (vrna_basepair_t *)vrna_realloc(pairs, sizeof(vrna_basepair_t) * len);
        }
      }
      /* put end marker in list */
      pairs[num_pairs].i  = pairs[num_pairs].j = 0;
      pairs               =
        (vrna_basepair_t *)vrna_realloc(pairs, sizeof(vrna_basepair_t) * (num_pairs + 1));
    }

    Py_XDECREF(result);
    return pairs;
  }


  static FLT_OR_DBL
  py_wrap_sc_exp_f_callback(int           i,
                            int           j,
                            int           k,
                            int           l,
                            unsigned char d,
                            void          *data)
  {
    FLT_OR_DBL        ret;
    PyObject          *func, *arglist, *result, *err;
    py_sc_callback_t  *cb;

    cb    = (py_sc_callback_t *)data;
    ret   = 1.;
    func  = cb->cb_exp_f;
    /* compose argument list */
    PyObject *py_i, *py_j, *py_k, *py_l, *py_d;

    py_i  = PyInt_FromLong(i);
    py_j  = PyInt_FromLong(j);
    py_k  = PyInt_FromLong(k);
    py_l  = PyInt_FromLong(l);
    py_d  = PyInt_FromLong(d);

    result = PyObject_CallFunctionObjArgs(func,
                                          py_i,
                                          py_j,
                                          py_k,
                                          py_l,
                                          py_d,
                                          (cb->data) ? cb->data : Py_None,
                                          NULL);
    Py_DECREF(py_i);
    Py_DECREF(py_j);
    Py_DECREF(py_k);
    Py_DECREF(py_l);
    Py_DECREF(py_d);

    /* BEGIN recognizing errors in callback execution */
    if (result == NULL) {
      if ((err = PyErr_Occurred())) {
        /* print error message */
        PyErr_Print();
        /* we only treat TypeErrors differently here, as they indicate that the callback does not follow requirements! */
        if (PyErr_GivenExceptionMatches(err, PyExc_TypeError))
          throw
          std::runtime_error(
            "Generic soft constraint callbacks (partition function) must take exactly 6 arguments");


        else
          throw
          std::runtime_error(
            "Some error occurred while executing generic soft constraint callback (partition function)");
      }

      PyErr_Clear();
    } else if (result == Py_None) {
      throw
      std::runtime_error(
        "Generic soft constraint callback (partition function) must return Boltzmann weighted pseudo energy value");
    } else {
      ret = (FLT_OR_DBL)PyFloat_AsDouble(result);
    }

    /* END recognizing errors in callback execution */

    Py_XDECREF(result);
    return ret;
  }


%}


static int
sc_add_f_pycallback(vrna_fold_compound_t  *vc,
                    PyObject              *callback);


static int
sc_add_bt_pycallback(vrna_fold_compound_t *vc,
                     PyObject             *PyFunc);


static int
sc_add_exp_f_pycallback(vrna_fold_compound_t  *vc,
                        PyObject              *PyFunc);


static int
sc_add_pydata(vrna_fold_compound_t  *vc,
              PyObject              *data,
              PyObject              *callback);


/* now we bind the above functions as methods to the fold_compound object */
%extend vrna_fold_compound_t {

%feature("autodoc") sc_add_data;
%feature("kwargs") sc_add_data;
%feature("autodoc") sc_add_f;
%feature("kwargs") sc_add_f;
%feature("autodoc") sc_add_bt;
%feature("kwargs") sc_add_bt;
%feature("autodoc") sc_add_exp_f;
%feature("kwargs") sc_add_exp_f;

  int
  sc_add_data(PyObject  *data,
              PyObject  *callback = Py_None)
  {
    return sc_add_pydata($self, data, callback);
  }


  int
  sc_add_f(PyObject *callback)
  {
    return sc_add_f_pycallback($self, callback);
  }


  int
  sc_add_bt(PyObject *PyFunc)
  {
    return sc_add_bt_pycallback($self, PyFunc);
  }


  int
  sc_add_exp_f(PyObject *PyFunc)
  {
    return sc_add_exp_f_pycallback($self, PyFunc);
  }
}

#endif
