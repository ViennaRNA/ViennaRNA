/*
 * file_py3.i: Typemaps for FILE* for Python 3
 *
 * Copyright (c) 2011, Karel Slany (karel.slany AT nic.cz)
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 *     * Redistributions of source code must retain the above copyright notice,
 *       this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the organization nor the names of its
 *       contributors may be used to endorse or promote products derived from this
 *       software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

/*
 *  Modifications to synchronize Python IO stream and C FILE * stream
 *  inspired by matplotlib and implemented 2018 by Ronny Lorenz
 */

%{
#include <unistd.h>
#include <fcntl.h>
%}

%types(FILE *);

/* converts basic file descriptor flags onto a string */
%fragment("fdfl_to_str", "header") {
const char *
fdfl_to_str(int fdfl) {

  static const char * const file_mode[] = {"w+", "w", "r"};

  if (fdfl & O_RDWR) {
    return file_mode[0];
  } else if (fdfl & O_WRONLY) {
    return file_mode[1];
  } else {
    return file_mode[2];
  }
}
}

#if 0
#define SWIG_FILE3_DEBUG
#endif

%fragment("obj_to_file","header", fragment="fdfl_to_str") {
FILE *
obj_to_file(PyObject *obj, long int *start_position) {
%#if PY_VERSION_HEX >= 0x03000000
  int fd, fd2, fdfl;
  long int position;
  FILE *fp;
  PyObject *ret, *os;
  if (!PyLong_Check(obj) &&                                /* is not an integer */
      PyObject_HasAttrString(obj, "fileno") &&             /* has fileno method */
      (PyObject_CallMethod(obj, "flush", NULL) != NULL) && /* flush() succeeded */
      ((fd = PyObject_AsFileDescriptor(obj)) != -1) &&     /* got file descriptor */
      ((fdfl = fcntl(fd, F_GETFL)) != -1)                  /* got descriptor flags */
    ) {
    os = PyImport_ImportModule("os");
    if (os == NULL)
      return NULL;

    ret = PyObject_CallMethod(os, (char *)"dup", (char *)"i", fd);
    Py_DECREF(os);

    if (ret == NULL)
      return NULL;

    fd2 = (int)PyNumber_AsSsize_t(ret, NULL);
    Py_DECREF(ret);

    fp = fdopen(fd2, fdfl_to_str(fdfl)); /* the FILE* must be flushed
                                            and closed after being used */
    if (fp == NULL)
      PyErr_SetString(PyExc_IOError, "Failed to get FILE * from Python file object");

    *start_position = ftell(fp);

    if (*start_position == -1) {
      /* fp is stream */
      return fp;
    }

    ret = PyObject_CallMethod(obj, (char *)"tell", (char *)"");
    if (ret == NULL) {
      fclose(fp);
      return NULL;
    }

    position = PyNumber_AsSsize_t(ret, PyExc_OverflowError);
    Py_DECREF(ret);

    if (PyErr_Occurred()) {
      fclose(fp);
      return NULL;
    }

    if (fseek(fp, position, SEEK_SET) == -1) {
      PyErr_SetString(PyExc_IOError, "Failed to seek FILE * to PyObject position");
      return NULL;
    }

#ifdef SWIG_FILE3_DEBUG
    fprintf(stderr, "opening fd %d (fl %d \"%s\") as FILE %p, start_position %ld, position: %ld\n",
            fd, fdfl, fdfl_to_str(fdfl), (void *)fp, *start_position, position);
#endif
    return fp;
  }
%#endif
  return NULL;
}
}

/* returns -1 if error occurred */
%fragment("dispose_file", "header") {
int
dispose_file(FILE **fp, PyObject *pyfile, long int start_position) {
  if (*fp == NULL)
    return 0;

  int fd;
  long int position;
  PyObject *ret, *exc_type, *exc_value, *exc_tb;

  exc_type = exc_value = exc_tb = NULL;
  PyErr_Fetch(&exc_type, &exc_value, &exc_tb);

  position = ftell(*fp);

#ifdef SWIG_FILE3_DEBUG
  fprintf(stderr, "flushing FILE %p, start_position %ld, position %ld\n",
          (void *)fp,
          start_position,
          position);
#endif

  if (!((fflush(*fp) == 0) &&
        (fclose(*fp) == 0)))
    return -1;

  *fp = NULL;

  fd = PyObject_AsFileDescriptor(pyfile);
  if (fd == -1)
    goto fail_dispose_file;

  if (lseek(fd, start_position, SEEK_SET) != -1) {
    if (position == -1) {
      PyErr_SetString(PyExc_IOError, "Failed to obtain FILE * position");
      goto fail_dispose_file;
    }

    ret = PyObject_CallMethod(pyfile, (char *)"seek", (char *)"ii", position, 0);
    if (ret == NULL)
      goto fail_dispose_file;

    Py_DECREF(ret);
  }

  PyErr_Restore(exc_type, exc_value, exc_tb);
  return 0;

fail_dispose_file:
  Py_XDECREF(exc_type);
  Py_XDECREF(exc_value);
  Py_XDECREF(exc_tb);
  
  return -1;
}
}

%typemap(arginit, noblock = 1) FILE* {
  $1 = NULL;
}

/* typemap for FILE * arguments that may be NULL */
%typemap(check, noblock = 1) FILE * nullfile {
  /* pass, even if $1 == NULL */
}

/* typemap for all other FILE * arguments */
%typemap(check, noblock = 1) FILE * {
  if ($1 == NULL) {
    /* The generated wrapper function raises TypeError on mismatching types. */
    SWIG_exception_fail(SWIG_TypeError, "in method '" "$symname" "', argument "
                        "$argnum"" of type '" "$type""'");
  }
}

/* typemap for FILE * arguments that may be NULL */
%typemap(typecheck, noblock = 1) FILE * nullfile {
  int fd, fdfl;
  if (($input == Py_None) ||
      (!PyLong_Check($input) &&                                /* is not an integer */
      PyObject_HasAttrString($input, "fileno") &&             /* has fileno method */
      (PyObject_CallMethod($input, "flush", NULL) != NULL) && /* flush() succeeded */
      ((fd = PyObject_AsFileDescriptor($input)) != -1) &&     /* got file descriptor */
      ((fdfl = fcntl(fd, F_GETFL)) != -1)                  /* got descriptor flags */
    )) {
    $1 = 1;
  } else {
    $1 = 0;
  }
}

/* typemap for all other FILE * arguments */
%typemap(typecheck, noblock = 1) FILE * {
  int fd, fdfl;
  if (!PyLong_Check($input) &&                                /* is not an integer */
      PyObject_HasAttrString($input, "fileno") &&             /* has fileno method */
      (PyObject_CallMethod($input, "flush", NULL) != NULL) && /* flush() succeeded */
      ((fd = PyObject_AsFileDescriptor($input)) != -1) &&     /* got file descriptor */
      ((fdfl = fcntl(fd, F_GETFL)) != -1)                  /* got descriptor flags */
    ) {
    $1 = 1;
  } else {
    $1 = 0;
  }
}

%typemap(in, noblock = 1, fragment = "obj_to_file") FILE *(PyObject *pyfile = NULL, long int start_position = -1) {
  if($input == Py_None){
    $1 = NULL;
  } else {
    pyfile = $input;
    $1 = obj_to_file($input, &start_position);
  }
}

%typemap(freearg, noblock = 1, fragment = "dispose_file") FILE* {
  if (dispose_file(&$1, pyfile$argnum, start_position$argnum) == -1) {
    SWIG_exception_fail(SWIG_IOError, "closing file in method '" "$symname" "', argument "
                        "$argnum"" of type '" "$type""'");
  }
}
