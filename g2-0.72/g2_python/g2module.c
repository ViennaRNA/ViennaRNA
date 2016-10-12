/*****************************************************************************
**  Copyright (C) 2006  Tijs Michels
**  This file is part of the g2 library
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Lesser General Public
**  License as published by the Free Software Foundation; either
**  version 2.1 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
**  Lesser General Public License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License along with this library; if not, write to the Free Software
**  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
******************************************************************************/

#include "Python.h"
#include "structmember.h"

#include "g2.h"
#ifdef DO_X11
# include "X11/g2_X11.h"
#endif
#ifdef DO_PS
# include "PS/g2_PS.h"
#endif
#ifdef DO_GD
# include "GD/g2_gd.h"
#endif
#ifdef DO_FIG
# include "FIG/g2_FIG.h"
#endif
#ifdef DO_WIN32
# include "WIN32/g2_win32.h"
#endif

typedef struct {
   PyObject_HEAD
   int dev; /* (virtual) device (readonly: set by 'constructor') */
} G2;

/* in Python, say for instance:
 *
 * import g2
 *
 * graph = g2.g2_open_gd('g2_demo.png', 400, 400, g2.g2_gd_png)
 * graph.g2_set_coordinate_system(200, 200, 1.0, 1.0)
 * # no device argument; all G2 methods use G2.dev
 */

/* a few helper functions, inaccessible from Python */

typedef void list_f(int, int, double *);
typedef void list_i_f(int, int, double *, int);
typedef void list_d_f(int, int, double *, double);

static PyTypeObject G2_Type; /* forward declaration */

static double
helper_item_float(PyObject *item)
{
   if PyFloat_Check(item) return PyFloat_AsDouble(item);
   if PyInt_Check(item) return (double) PyInt_AsLong(item);
   return 0; /* list is assumed to contain floats and/or integers */
}

static PyObject *
helper_il(const G2 *self, const PyObject *args, list_f *f)
{
   PyObject *list;

   if (PyArg_ParseTuple((PyObject *)args, "O!", &PyList_Type, &list)) { /* cf. <listobject.h> */
      int s = PyList_Size(list);
      if (s) {
         double * const points = malloc(s * sizeof(double)); /* in one case the buffer holds dashes, not points */
         if (points) {
            const int np = (f == g2_set_dash) ? s : s >> 1; /* two co-ordinates per point */
            while (s--) points[s] = helper_item_float(PyList_GET_ITEM(list, s));
            (*f)(self->dev, np, points);
            free(points);
            Py_RETURN_NONE;
         }
         return PyErr_NoMemory();
      }
      PyErr_SetString(PyExc_ValueError, "empty list");
   }
   return NULL;
}

static PyObject *
helper_ili(const G2 *self, const PyObject *args, list_i_f *f)
{
   PyObject *list;
   int ip; /* number of interpolated points */

   if (PyArg_ParseTuple((PyObject *)args, "O!i", &PyList_Type, &list, &ip)) {
      int s = PyList_Size(list);
      if (s > 5) {
         double * const points = malloc(s * sizeof(double));
         if (points) {
            const int np = s >> 1;
            while (s--) points[s] = helper_item_float(PyList_GET_ITEM(list, s));
            (*f)(self->dev, np, points, ip);
            free(points);
            Py_RETURN_NONE;
         }
         return PyErr_NoMemory();
      }
      PyErr_SetString(PyExc_ValueError, "list must hold at least three points");
   }
   return NULL;
}

static PyObject *
helper_ild(const G2 *self, const PyObject *args, list_d_f *f)
{
   PyObject *list;
   double factor;

   if (PyArg_ParseTuple((PyObject *)args, "O!d", &PyList_Type, &list, &factor)) {
      int s = PyList_Size(list);
      if (s > 5) {
         double * const points = malloc(s * sizeof(double));
         if (points) {
            const int np = s >> 1;
            while (s--) points[s] = helper_item_float(PyList_GET_ITEM(list, s));
            (*f)(self->dev, np, points, factor);
            free(points);
            Py_RETURN_NONE;
         }
         return PyErr_NoMemory();
      }
      PyErr_SetString(PyExc_ValueError, "list must hold at least three points");
   }
   return NULL;
}

/* end of helper functions
 *
 * functions from g2.h
 */

/* three unbound functions */

PyDoc_STRVAR(doc_g2_ld,
             "int g2_ld()\n"
             "Return the latest accessed device (an int).");

static PyObject *
C_g2_ld(PyObject *self)
{
   return Py_BuildValue("i", g2_ld());
}

PyDoc_STRVAR(doc_g2_set_ld,
             "g2_set_ld(int dev)\n"
             "Set the latest accessed device.\n"
             "e.g. : g2_set_ld(PS_graph.dev)");

static PyObject *
C_g2_set_ld(PyObject *self, PyObject *args)
{
   int dev;

   if (PyArg_ParseTuple(args, "i", &dev)) {
      g2_set_ld(dev);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_device_exist,
             "bool g2_device_exist(int dev)\n"
             "Return True if device exists, False if not.");

static PyObject *
C_g2_device_exist(PyObject *self, PyObject *args)
{
   int dix;

   if (PyArg_ParseTuple(args, "i", &dix))
      return Py_BuildValue("i", g2_device_exist(dix));
   return NULL;
}

/* (device specific) 'constructors', all returning an object of class G2 */

PyDoc_STRVAR(doc_g2_open_vd,
             "G2 g2_open_vd()\n"
             "Construct an object of class G2 that plots on a virtual device.\n"
             "e.g. : graph = g2_open_vd()\n"
             "       graph.g2_attach(g2_open_X11(800, 600))\n"
             "       graph.g2_attach(g2_open_PS(\"foo.ps\", g2_A4, g2_PS_land))\n"
             "       graph.g2_line(3, 3, 9, 9)");

static PyObject *
C_g2_open_vd(PyObject *self)
{
   G2 * const obj = (G2 *) PyType_GenericAlloc(&G2_Type, 0);
   if (obj) {
      obj->dev = g2_open_vd();
      return (PyObject *)obj;
   }
   return NULL;
}

#ifdef _G2_X11_H

PyDoc_STRVAR(doc_g2_open_X11,
             "G2 g2_open_X11(int width, int height)\n"
             "Construct an object of class G2 that plots on an X11 window.\n"
             "e.g. : graph = g2_open_X11(800, 600)");

static PyObject *
C_g2_open_X11(PyObject *self, PyObject *args)
{
   int width, height;

   if (PyArg_ParseTuple(args, "ii", &width, &height)) {
      G2 * const obj = (G2 *) PyType_GenericAlloc(&G2_Type, 0);
      if (obj) {
         obj->dev = g2_open_X11(width, height);
         return (PyObject *)obj;
      }
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_open_X11X,
             "G2 g2_open_X11X(int width, int height, int x, int y,\n"
             "                str window_name, str icon_name, str icon_data,\n"
             "                int icon_width, int icon_height)\n"
             "Construct an object of class G2 that plots on an X11 window.");

static PyObject *
C_g2_open_X11X(PyObject *self, PyObject *args)
{
   int width, height, x, y, icon_width, icon_height;
   char *window_name, *icon_name, *icon_data;

   if (PyArg_ParseTuple(args, "iiiisssii", &width, &height, &x, &y,
                        &window_name, &icon_name, &icon_data,
                        &icon_width, &icon_height)) {
      G2 * const obj = (G2 *) PyType_GenericAlloc(&G2_Type, 0);
      if (obj) {
         obj->dev = g2_open_X11X(width, height, x, y,
                                 window_name, icon_name, icon_data,
                                 icon_width, icon_height);
         return (PyObject *)obj;
      }
   }
   return NULL;
}

#endif
#ifdef _G2_PS_H

PyDoc_STRVAR(doc_g2_open_PS,
             "G2 g2_open_PS(str filename, enum paper, enum orientation)\n"
             "Construct an object of class G2 that writes a PostScript file.\n"
             "e.g. : graph = g2_open_PS(\"foo.ps\", g2_A4, g2_PS_land)");

static PyObject *
C_g2_open_PS(PyObject *self, PyObject *args)
{
   const char *filename;
   int paper, orientation;

   if (PyArg_ParseTuple(args, "sii", &filename, &paper, &orientation)) {
      G2 * const obj = (G2 *) PyType_GenericAlloc(&G2_Type, 0);
      if (obj) {
         obj->dev = g2_open_PS(filename, paper, orientation);
         return (PyObject *)obj;
      }
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_open_EPSF,
             "G2 g2_open_EPSF(str filename)\n"
             "Construct an object of class G2 that writes\n"
             "an Encapsulated PostScript file.\n");

static PyObject *
C_g2_open_EPSF(PyObject *self, PyObject *args)
{
   const char *filename;

   if (PyArg_ParseTuple(args, "s", &filename)) {
      G2 * const obj = (G2 *) PyType_GenericAlloc(&G2_Type, 0);
      if (obj) {
         obj->dev = g2_open_EPSF(filename);
         return (PyObject *)obj;
      }
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_open_EPSF_CLIP,
             "G2 g2_open_EPSF_CLIP(str filename, int width, int height)\n"
             "Construct an object of class G2 that writes\n"
             "a clipped Encapsulated PostScript file.\n");

static PyObject *
C_g2_open_EPSF_CLIP(PyObject *self, PyObject *args)
{
   const char *filename;
   int width, height;

   if (PyArg_ParseTuple(args, "sii", &filename, &width, &height)) {
      G2 * const obj = (G2 *) PyType_GenericAlloc(&G2_Type, 0);
      if (obj) {
         obj->dev = g2_open_EPSF_CLIP(filename, width, height);
         return (PyObject *)obj;
      }
   }
   return NULL;
}

#endif
#ifdef _G2_GD_H

PyDoc_STRVAR(doc_g2_open_gd,
             "G2 g2_open_gd(str filename, int width, int height, enum gd_type)\n"
             "   'gd_type' must be one of g2_gd_jpeg, g2_gd_png or g2_gd_gif.\n"
             "Construct an object of class G2 that writes a jpeg/png/gif file.");

static PyObject *
C_g2_open_gd(PyObject *self, PyObject *args)
{
   const char *filename;
   int width, height, gd_type;

   if (PyArg_ParseTuple(args, "siii", &filename, &width, &height, &gd_type)) {
      G2 * const obj = (G2 *) PyType_GenericAlloc(&G2_Type, 0);
      if (obj) {
         obj->dev = g2_open_gd(filename, width, height, gd_type);
         return (PyObject *)obj;
      }
   }
   return NULL;
}

#endif
#ifdef _G2_FIG_H

PyDoc_STRVAR(doc_g2_open_FIG,
             "G2 g2_open_FIG(str filename)\n"
             "Construct an object of class G2 that plots in the FIG 3.2 format.\n"
             "See http://www.xfig.org.");

static PyObject *
C_g2_open_FIG(PyObject *self, PyObject *args)
{
   const char *filename;

   if (PyArg_ParseTuple(args, "s", &filename)) {
      G2 * const obj = (G2 *) PyType_GenericAlloc(&G2_Type, 0);
      if (obj) {
         obj->dev = g2_open_FIG(filename);
         return (PyObject *)obj;
      }
   }
   return NULL;
}

#endif
#ifdef _G2_WIN32_H

PyDoc_STRVAR(doc_g2_open_win32,
             "G2 g2_open_win32(int width, int height, str name, enum win32_type)\n"
             "   'name' is either a file name or a window title.\n"
             "   'win32_type' is either g2_win32 or g2_wmf32.\n"
             "Construct an object of class G2 that plots on an MS Windows window or metafile.");

static PyObject *
C_g2_open_win32(PyObject *self, PyObject *args)
{
   int width, height, win32_type;
   const char *name; /* file name or window title */

   if (PyArg_ParseTuple(args, "iisi", &width, &height, &name, &win32_type)) {
      G2 * const obj = (G2 *) PyType_GenericAlloc(&G2_Type, 0);
      if (obj) {
         obj->dev = g2_open_win32(width, height, name, win32_type);
         return (PyObject *)obj;
      }
   }
   return NULL;
}

#endif

static PyMethodDef module_functions[] = {
   { "g2_ld", (PyCFunction)C_g2_ld, METH_NOARGS, doc_g2_ld },
   { "g2_set_ld", C_g2_set_ld, METH_VARARGS, doc_g2_set_ld },
   { "g2_device_exist", C_g2_device_exist, METH_VARARGS, doc_g2_device_exist },
   /* nine functions that construct an instance of class G2 */
   { "g2_open_vd", (PyCFunction)C_g2_open_vd, METH_NOARGS, doc_g2_open_vd },
#ifdef _G2_X11_H
   { "g2_open_X11", C_g2_open_X11, METH_VARARGS, doc_g2_open_X11 },
   { "g2_open_X11X", C_g2_open_X11X, METH_VARARGS, doc_g2_open_X11X },
#endif
#ifdef _G2_PS_H
   { "g2_open_PS", C_g2_open_PS, METH_VARARGS, doc_g2_open_PS },
   { "g2_open_EPSF", C_g2_open_EPSF, METH_VARARGS, doc_g2_open_EPSF },
   { "g2_open_EPSF_CLIP", C_g2_open_EPSF_CLIP, METH_VARARGS, doc_g2_open_EPSF_CLIP },
#endif
#ifdef _G2_GD_H
   { "g2_open_gd", C_g2_open_gd, METH_VARARGS, doc_g2_open_gd },
#endif
#ifdef _G2_FIG_H
   { "g2_open_FIG", C_g2_open_FIG, METH_VARARGS, doc_g2_open_FIG },
#endif
#ifdef _G2_WIN32_H
   { "g2_open_win32", C_g2_open_win32, METH_VARARGS, doc_g2_open_win32 },
#endif
   { NULL }
};

/* (device specific) enums */

static void
add_enums(PyObject *m)
{
   typedef struct { const char *name; int value; } tInt;
   const tInt enums[] = {
#ifdef _G2_PS_H
      { "g2_A0",               g2_A0               },
      { "g2_A1",               g2_A1               },
      { "g2_A2",               g2_A2               },
      { "g2_A3",               g2_A3               },
      { "g2_A4",               g2_A4               },
      { "g2_A5",               g2_A5               },
      { "g2_A6",               g2_A6               },
      { "g2_A7",               g2_A7               },
      { "g2_A8",               g2_A8               },
      { "g2_A9",               g2_A9               },
      { "g2_B0",               g2_B0               },
      { "g2_B1",               g2_B1               },
      { "g2_B2",               g2_B2               },
      { "g2_B3",               g2_B3               },
      { "g2_B4",               g2_B4               },
      { "g2_B5",               g2_B5               },
      { "g2_B6",               g2_B6               },
      { "g2_B7",               g2_B7               },
      { "g2_B8",               g2_B8               },
      { "g2_B9",               g2_B9               },
      { "g2_B10",              g2_B10              },
      { "g2_Comm_10_Envelope", g2_Comm_10_Envelope },
      { "g2_C5_Envelope",      g2_C5_Envelope      },
      { "g2_DL_Envelope",      g2_DL_Envelope      },
      { "g2_Folio",            g2_Folio            },
      { "g2_Executive",        g2_Executive        },
      { "g2_Letter",           g2_Letter           },
      { "g2_Legal",            g2_Legal            },
      { "g2_Ledger",           g2_Ledger           },
      { "g2_Tabloid",          g2_Tabloid          },
      /* orientation */
      { "g2_PS_land",          g2_PS_land },
      { "g2_PS_port",          g2_PS_port },
      /* g2 Format */
      { "g2_PS_PostScript",    g2_PS_PostScript },
      { "g2_PS_EPSF",          g2_PS_EPSF       },
      { "g2_PS_EPSF_CLIP",     g2_PS_EPSF_CLIP  },
#endif
#ifdef _G2_GD_H
      { "g2_gd_jpeg", g2_gd_jpeg },
      { "g2_gd_png",  g2_gd_png  },
      { "g2_gd_gif",  g2_gd_gif  },
#endif
#ifdef _G2_WIN32_H
      { "g2_win32",   g2_win32 },
      { "g2_wmf32",   g2_wmf32 },
#endif
      { "QPrect",     QPrect },
      { "QPcirc",     QPcirc }
   };
   const tInt *c = enums + sizeof(enums)/sizeof(tInt);
   while (c-- > enums)
      PyModule_AddIntConstant(m, (char *)c->name, c->value);
   PyModule_AddStringConstant(m, "G2_VERSION", G2_VERSION);
}

/* methods of class G2, all of which use G2 member dev */

/* methods that operate on a single device */

PyDoc_STRVAR(doc_g2_attach,
             "g2_attach(G2 graph)\n"
             "Attach an object of class G2 to a virtual device.\n"
             "e.g. : graph.g2_attach(g2_open_X11(800, 600))");

static PyObject *
C_g2_attach(G2 *self, PyObject *args)
{
   PyObject *dev;

   if (PyArg_ParseTuple(args, "O!", &G2_Type, &dev)) {
      g2_attach(self->dev, ((G2 *)dev)->dev);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_detach,
             "g2_detach(G2 graph)\n"
             "Detach an object of class G2 from a virtual device.");

static PyObject *
C_g2_detach(G2 *self, PyObject *args)
{
   PyObject *dev;

   if (PyArg_ParseTuple(args, "O!", &G2_Type, &dev)) {
      g2_detach(self->dev, ((G2 *)dev)->dev);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_ink,
             "int g2_ink(float red, float green, float blue)\n"
             "Return an ink of the specified color.\n"
             "Arguments must be between 0 and 1.\n"
             "Use only on instances of class G2 that represent a physical device.");

static PyObject *
C_g2_ink(G2 *self, PyObject *args) /* only on a physical device */
{
   double red, green, blue;

   if (PyArg_ParseTuple(args, "ddd", &red, &green, &blue))
      return Py_BuildValue("i", g2_ink(self->dev, red, green, blue));
   return NULL;
}

/* methods that can operate on multiple devices */

PyDoc_STRVAR(doc_g2_close,
             "g2_close()\n"
             "Close and delete a device.\n"
             "Later, garbage collection will delete\n"
             "the now useless object of class G2.");

static PyObject *
C_g2_close(G2 *self)
{
   g2_close(self->dev);
   Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_g2_set_auto_flush,
             "g2_set_auto_flush(bool on_off)\n"
             "Flush or not after each graphical operation.\n"
             "Note : this slows g2 down considerably.");

static PyObject *
C_g2_set_auto_flush(G2 *self, PyObject *args)
{
   int on_off;

   if (PyArg_ParseTuple(args, "i", &on_off)) {
      g2_set_auto_flush(self->dev, on_off);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_flush,
             "g2_flush()\n"
             "Flush the device (or devices in case of a virtual device).");

static PyObject *
C_g2_flush(G2 *self)
{
   g2_flush(self->dev);
   Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_g2_save,
             "g2_save()\n"
             "Save the output. This is implied by g2_close.");

static PyObject *
C_g2_save(G2 *self)
{
   g2_save(self->dev);
   Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_g2_set_coordinate_system,
             "g2_set_coordinate_system(float x_origin, float y_origin,\n"
             "                         float x_multiply, float y_multiply)\n"
             "Set a user defined co-ordinate system.");

static PyObject *
C_g2_set_coordinate_system(G2 *self, PyObject *args)
{
   double x_origin, y_origin, x_mul, y_mul;

   if (PyArg_ParseTuple(args, "dddd", &x_origin, &y_origin, &x_mul, &y_mul)) {
      g2_set_coordinate_system(self->dev, x_origin, y_origin, x_mul, y_mul);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_pen,
             "g2_pen(int pen)\n"
             "Set pen color for all following operations. See also g2_ink.\n"
             "e.g. : graph.g2_pen(graph.g2_ink(.25, .6, 0))");

static PyObject *
C_g2_pen(G2 *self, PyObject *args)
{
   int color;

   if (PyArg_ParseTuple(args, "i", &color)) {
      g2_pen(self->dev, color);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_set_dash,
             "g2_set_dash(list pattern)\n"
             "   'pattern' : [length, ... length]\n"
             "               In the list ints and floats can be mixed freely.\n"
             "Set line dash. See also g2_set_solid.\n"
             "e.g. : g2_set_dash([4, 2])\n"
             "       for lines with dashes twice as long\n"
             "       as the white spaces between them.\n"
             "Note : this one argument form is Python specific.");

static PyObject *
C_g2_set_dash(G2 *self, PyObject *args)
{
   return helper_il(self, args, g2_set_dash);
}

/* python specific function : more explicit than passing an empty list */

PyDoc_STRVAR(doc_g2_set_solid,
             "g2_set_solid()\n"
             "Set the line style to solid. See also g2_set_dash.\n"
             "This method is Python specific.");

static PyObject *
C_g2_set_solid(G2 *self)
{
   g2_set_dash(self->dev, 0, NULL);
   Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_g2_set_font_size,
             "g2_set_font_size(float size)\n"
             "Set the font size.");

static PyObject *
C_g2_set_font_size(G2 *self, PyObject *args)
{
   double size;

   if (PyArg_ParseTuple(args, "d", &size)) {
      g2_set_font_size(self->dev, size);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_set_line_width,
             "g2_set_line_width(float size)\n"
             "Set the line width.");

static PyObject *
C_g2_set_line_width(G2 *self, PyObject *args)
{
   double size;

   if (PyArg_ParseTuple(args, "d", &size)) {
      g2_set_line_width(self->dev, size);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_clear_palette,
             "g2_clear_palette()\n"
             "Remove all inks.");

static PyObject *
C_g2_clear_palette(G2 *self)
{
   g2_clear_palette(self->dev);
   Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_g2_reset_palette,
             "g2_reset_palette()\n"
             "Remove all inks and reallocate the default colors.");

static PyObject *
C_g2_reset_palette(G2 *self)
{
   g2_reset_palette(self->dev);
   Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_g2_allocate_basic_colors,
             "g2_allocate_basic_colors()\n"
             "Allocate the default colors.");

static PyObject *
C_g2_allocate_basic_colors(G2 *self)
{
   g2_allocate_basic_colors(self->dev);
   Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_g2_clear,
             "g2_clear()\n"
             "Clear the device.");

static PyObject *
C_g2_clear(G2 *self)
{
   g2_clear(self->dev);
   Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_g2_set_background,
             "g2_set_background(int pen)\n"
             "Set the background color.\n"
             "e.g. : graph.g2_set_background(graph.g2_ink(.25, .6, 0))");

static PyObject *
C_g2_set_background(G2 *self, PyObject *args)
{
   int color;

   if (PyArg_ParseTuple(args, "i", &color)) {
      g2_set_background(self->dev, color);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_move,
             "g2_move(float x, float y)\n"
             "Move graphical cursor to (x, y).");

static PyObject *
C_g2_move(G2 *self, PyObject *args)
{
   double x, y;

   if (PyArg_ParseTuple(args, "dd", &x, &y)) {
      g2_move(self->dev, x, y);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_move_r,
             "g2_move_r(float x, float y)\n"
             "Move graphical cursor to (x, y) from\n"
             "the current graphical cursor position.");

static PyObject *
C_g2_move_r(G2 *self, PyObject *args)
{
   double x, y;

   if (PyArg_ParseTuple(args, "dd", &x, &y)) {
      g2_move_r(self->dev, x, y);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_plot,
             "g2_plot(float x, float y)\n"
             "Plot a dot at position (x, y).");

static PyObject *
C_g2_plot(G2 *self, PyObject *args)
{
   double x, y;

   if (PyArg_ParseTuple(args, "dd", &x, &y)) {
      g2_plot(self->dev, x, y);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_plot_r,
             "g2_plot_r(float x, float y)\n"
             "Plot a dot at (x, y) from the current graphical cursor position.");

static PyObject *
C_g2_plot_r(G2 *self, PyObject *args)
{
   double x, y;

   if (PyArg_ParseTuple(args, "dd", &x, &y)) {
      g2_plot_r(self->dev, x, y);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_line,
             "g2_line(float x1, float y1, float x2, float y2)\n"
             "Plot a line from (x1, y1) to (x2, y2).");

static PyObject *
C_g2_line(G2 *self, PyObject *args)
{
   double x1, y1, x2, y2;

   if (PyArg_ParseTuple(args, "dddd", &x1, &y1, &x2, &y2)) {
      g2_line(self->dev, x1, y1, x2, y2);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_line_r,
             "g2_line_r(float x, float y)\n"
             "Plot a line from the current graphical cursor position\n"
             "to (x, y) from there.");

static PyObject *
C_g2_line_r(G2 *self, PyObject *args)
{
   double x, y;

   if (PyArg_ParseTuple(args, "dd", &x, &y)) {
      g2_line_r(self->dev, x, y);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_line_to,
             "g2_line_to(float x, float y)\n"
             "Plot a line from the current graphical cursor position to (x, y).");

static PyObject *
C_g2_line_to(G2 *self, PyObject *args)
{
   double x, y;

   if (PyArg_ParseTuple(args, "dd", &x, &y)) {
      g2_line_to(self->dev, x, y);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_poly_line,
             "g2_poly_line(list points)\n"
             "   'points' : [x1, y1, x2, y2, ... xn, yn]\n"
             "              In the list ints and floats can be mixed freely.\n"
             "Plot a line through the points in the list.\n"
             "e.g. : graph.g2_poly_line([2, 4, 2.5, 6.25, 3, 9, 3.5, 12.25])\n"
             "Note : this one argument form is Python specific.");

static PyObject *
C_g2_poly_line(G2 *self, PyObject *args)
{
   return helper_il(self, args, g2_poly_line);
}

PyDoc_STRVAR(doc_g2_triangle,
             "g2_triangle(float x1, float y1, float x2, float y2,\n"
             "            float x3, float y3)\n"
             "Plot a triangle through the points (x1, y1), (x2, y2) and (x3, y3).");

static PyObject *
C_g2_triangle(G2 *self, PyObject *args)
{
   double x1, y1, x2, y2, x3, y3;

   if (PyArg_ParseTuple(args, "dddddd", &x1, &y1, &x2, &y2, &x3, &y3)) {
      g2_triangle(self->dev, x1, y1, x2, y2, x3, y3);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_filled_triangle,
             "g2_filled_triangle(float x1, float y1, float x2, float y2,\n"
             "                   float x3, float y3)\n"
             "Fill a triangle through the points (x1, y1), (x2, y2) and (x3, y3).");

static PyObject *
C_g2_filled_triangle(G2 *self, PyObject *args)
{
   double x1, y1, x2, y2, x3, y3;

   if (PyArg_ParseTuple(args, "dddddd", &x1, &y1, &x2, &y2, &x3, &y3)) {
      g2_filled_triangle(self->dev, x1, y1, x2, y2, x3, y3);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_rectangle,
             "g2_rectangle(float x1, float y1, float x2, float y2)\n"
             "Plot a rectangle through the points (x1, y1) and (x2, y2).");

static PyObject *
C_g2_rectangle(G2 *self, PyObject *args)
{
   double x1, y1, x2, y2;

   if (PyArg_ParseTuple(args, "dddd", &x1, &y1, &x2, &y2)) {
      g2_rectangle(self->dev, x1, y1, x2, y2);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_filled_rectangle,
             "g2_filled_rectangle(float x1, float y1, float x2, float y2)\n"
             "Fill a rectangle through the points (x1, y1) and (x2, y2).");

static PyObject *
C_g2_filled_rectangle(G2 *self, PyObject *args)
{
   double x1, y1, x2, y2;

   if (PyArg_ParseTuple(args, "dddd", &x1, &y1, &x2, &y2)) {
      g2_filled_rectangle(self->dev, x1, y1, x2, y2);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_polygon,
             "g2_polygon(list points)\n"
             "   'points' : [x1, y1, x2, y2, ... xn, yn]\n"
             "              In the list ints and floats can be mixed freely.\n"
             "Plot a polygon through the points in the list.\n"
             "e.g. : graph.g2_polygon([3, 1, 1, 3, 3, 5, 5, 3])\n"
             "Note : this one argument form is Python specific.");

static PyObject *
C_g2_polygon(G2 *self, PyObject *args)
{
   return helper_il(self, args, g2_polygon);
}

PyDoc_STRVAR(doc_g2_filled_polygon,
             "g2_filled_polygon(list points)\n"
             "   'points' : [x1, y1, x2, y2, ... xn, yn]\n"
             "              In the list ints and floats can be mixed freely.\n"
             "Fill a polygon through the points in the list.\n"
             "Note : this one argument form is Python specific.");

static PyObject *
C_g2_filled_polygon(G2 *self, PyObject *args)
{
   return helper_il(self, args, g2_filled_polygon);
}

PyDoc_STRVAR(doc_g2_circle,
             "g2_circle(float x, float y, float r)\n"
             "Plot a circle with center (x, y) and radius r.");

static PyObject *
C_g2_circle(G2 *self, PyObject *args)
{
   double x, y, r;

   if (PyArg_ParseTuple(args, "ddd", &x, &y, &r)) {
      g2_circle(self->dev, x, y, r);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_filled_circle,
             "g2_filled_circle(float x, float y, float r)\n"
             "Fill a circle with center (x, y) and radius r.");

static PyObject *
C_g2_filled_circle(G2 *self, PyObject *args)
{
   double x, y, r;

   if (PyArg_ParseTuple(args, "ddd", &x, &y, &r)) {
      g2_filled_circle(self->dev, x, y, r);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_ellipse,
             "g2_ellipse(float x, float y, float r1, float r2)\n"
             "Plot an ellipse with center (x, y) and radii r1 and r2.");

static PyObject *
C_g2_ellipse(G2 *self, PyObject *args)
{
   double x, y, r1, r2;

   if (PyArg_ParseTuple(args, "dddd", &x, &y, &r1, &r2)) {
      g2_ellipse(self->dev, x, y, r1, r2);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_filled_ellipse,
             "g2_filled_ellipse(float x, float y, float r1, float r2)\n"
             "Fill an ellipse with center (x, y) and radii r1 and r2.");

static PyObject *
C_g2_filled_ellipse(G2 *self, PyObject *args)
{
   double x, y, r1, r2;

   if (PyArg_ParseTuple(args, "dddd", &x, &y, &r1, &r2)) {
      g2_filled_ellipse(self->dev, x, y, r1, r2);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_arc,
             "g2_arc(float x, float y, float r1, float r2,\n"
             "       float a1, float a2)\n"
             "Plot an arc with center (x, y), radii r1 and r2,\n"
             "and angles a1 and a2.");

static PyObject *
C_g2_arc(G2 *self, PyObject *args)
{
   double x, y, r1, r2, a1, a2;

   if (PyArg_ParseTuple(args, "dddddd", &x, &y, &r1, &r2, &a1, &a2)) {
      g2_arc(self->dev, x, y, r1, r2, a1, a2);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_filled_arc,
             "g2_filled_arc(float x, float y, float r1, float r2,\n"
             "              float a1, float a2)\n"
             "Fill an arc with center (x, y), radii r1 and r2,\n"
             "and angles a1 and a2.");

static PyObject *
C_g2_filled_arc(G2 *self, PyObject *args)
{
   double x, y, r1, r2, a1, a2;

   if (PyArg_ParseTuple(args, "dddddd", &x, &y, &r1, &r2, &a1, &a2)) {
      g2_filled_arc(self->dev, x, y, r1, r2, a1, a2);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_string,
             "g2_string(float x, float y, str s)\n"
             "Write string s at position (x, y).");

static PyObject *
C_g2_string(G2 *self, PyObject *args)
{
   double x, y;
   const char *s;

   if (PyArg_ParseTuple(args, "dds", &x, &y, &s)) {
      g2_string(self->dev, x, y, s);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_image,
             "g2_image(float x, float y, list list_of_lists)\n"
             "   Each list in list_of_lists represents one row,\n"
             "   where a row is a list of pens.\n"
             "   The longest row determines the image width.\n"
             "   Shorter rows are filled with pen 0.\n"
             "   The number of lists in list_of_lists determines the image height.\n"
             "e.g. : g2_image(5, 7, [[2, 4, 6], [3, 6, 9], [4, 8, 12]])\n"
             "       plots a 3x3 bitmap at position (5, 7).\n"
             "Note : this three argument form is Python specific.");

static PyObject *
C_g2_image(G2 *self, PyObject *args)
{
   double x, y;
   PyObject *list; /* a list of lists, to be unpacked as a two-dimensional array */

   if (PyArg_ParseTuple(args, "ddO!", &x, &y, &PyList_Type, &list)) {
      int h = PyList_Size(list); /* the nr of lines, i.e. the height */
      if (h) {
         int i = 0;
         int j = 0;
         PyObject * po;
         PyObject ** const lines = malloc(h * sizeof(PyObject *));
         if (lines == NULL) return PyErr_NoMemory();
         do {
            po = PyList_GET_ITEM(list, i++);
            if (PyList_Check(po)) lines[j++] = po; /* only store lists */
         } while (i < h);
         i -= j; /* in case not all elements were lists, i.e. lines, */
         h -= i; /* decrement the height */
         if (h) { /* our outer list contained at least one nested list, i.e. line */
            PyObject * const * const e = lines + h;
            PyObject * const * l = e;
            int w = 0;
            do { /* set width to the longest list */
               if ((i = PyList_Size(*--l)) > w) w = i;
            } while (l > lines);
            if (w) { /* at least one line had a width > 0 */
               int * const pens = PyMem_New(int, w * h * sizeof(int));
               int * pen = pens;
               if (pens == NULL) {
                  free(lines);
                  return PyErr_NoMemory();
               }
               do {
                  PyObject * const pl = *l;
                  i = 0;
                  j = PyList_Size(pl);
                  while (i < j) {
                     po = PyList_GET_ITEM(pl, i);
                     pen[i++] = PyInt_Check(po) ? PyInt_AsLong(po) : 0;
                  }
                  while (i < w) pen[i++] = 0;
                  pen += w;
                  l++;
               } while (l < e); /* the previous loop ended with l == lines */
               g2_image(self->dev, x, y, w, h, pens);
               PyMem_Del(pens);
            }
         }
         free(lines);
         Py_RETURN_NONE;
      }
      PyErr_SetString(PyExc_ValueError, "empty list");
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_set_QP,
             "g2_set_QP(float size, enum shape)\n"
             "   'shape' must be either QPrect or QPcirc.\n"
             "Set QuasiPixel size and shape. See g2_plot_QP.");

static PyObject *
C_g2_set_QP(G2 *self, PyObject *args)
{
   double size;
   int shape;

   if (PyArg_ParseTuple(args, "di", &size, &shape)) {
      g2_set_QP(self->dev, size, shape);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_plot_QP,
             "g2_plot_QP(float x, float y)\n"
             "Plot a Quasi Pixel at position (x, y).\n"
             "Quasi Pixels make it easy to plot cellular automata and related\n"
             "images. QP is simply a big pixel as defined by g2_set_QP.\n"
             "Co-ordinates are scaled accordingly, so no recalculation is needed\n"
             "on the client side.");

static PyObject *
C_g2_plot_QP(G2 *self, PyObject *args)
{
   double x, y;

   if (PyArg_ParseTuple(args, "dd", &x, &y)) {
      g2_plot_QP(self->dev, x, y);
      Py_RETURN_NONE;
   }
   return NULL;
}

PyDoc_STRVAR(doc_g2_query_pointer,
             "[float x, float y, int button] g2_query_pointer()\n"
             "Query pointer position and button state (e.g. mouse for X11).\n"
             "The results are returned in a list of two floats and one int.\n"
             "e.g. : x, y, button = g2_query_pointer()\n"
             "Note : this no argument form is Python specific.");

static PyObject *
C_g2_query_pointer(G2 *self)
{
   double x, y;
   unsigned int button;
   g2_query_pointer(self->dev, &x, &y, &button);
   PyObject * const r = PyTuple_New(3);
   PyTuple_SET_ITEM(r, 0, Py_BuildValue("d", x));
   PyTuple_SET_ITEM(r, 1, Py_BuildValue("d", y));
   PyTuple_SET_ITEM(r, 2, Py_BuildValue("I", button)); /* capital, for unsigned */
   return r;
}

PyDoc_STRVAR(doc_g2_spline,
             "g2_spline(list points, int ppdp)\n"
             "   'points' : [x1, y1, x2, y2, ... xn, yn]\n"
             "              In the list ints and floats can be mixed freely.\n"
             "   'ppdp'   : the number of interpolated points per data point.\n"
             "              The higher 'ppdp', the rounder the spline curve.\n"
             "Plot a spline curve through the points in the list\n"
             "using Young's method of successive over-relaxation.\n"
             "Note : this two argument form is Python specific.");

static PyObject *
C_g2_spline(G2 *self, PyObject *args)
{
   return helper_ili(self, args, g2_spline);
}

PyDoc_STRVAR(doc_g2_filled_spline,
             "g2_filled_spline(list points, int ppdp)\n"
             "As g2_spline, but filled.");

static PyObject *
C_g2_filled_spline(G2 *self, PyObject *args)
{
   return helper_ili(self, args, g2_filled_spline);
}

PyDoc_STRVAR(doc_g2_b_spline,
             "g2_b_spline(list points, int ppdp)\n"
             "   'points' : [x1, y1, x2, y2, ... xn, yn]\n"
             "              In the list ints and floats can be mixed freely.\n"
             "   'ppdp'   : the number of interpolated points per data point.\n"
             "              The higher 'ppdp', the rounder the spline curve.\n"
             "Plot a b-spline curve through the points in the list.\n"
             "For most averaging purposes, this is the right spline.\n"
             "Note : this two argument form is Python specific.");

static PyObject *
C_g2_b_spline(G2 *self, PyObject *args)
{
   return helper_ili(self, args, g2_b_spline);
}

PyDoc_STRVAR(doc_g2_filled_b_spline,
             "g2_filled_b_spline(list points, int ppdp)\n"
             "As g2_b_spline, but filled.");

static PyObject *
C_g2_filled_b_spline(G2 *self, PyObject *args)
{
   return helper_ili(self, args, g2_filled_b_spline);
}

PyDoc_STRVAR(doc_g2_raspln,
             "g2_raspln(list points, float tfact)\n"
             "   'points' : [x1, y1, x2, y2, ... xn, yn]\n"
             "              In the list ints and floats can be mixed freely.\n"
             "   'tfact'  : tension factor between 0 (very rounded) and 2.\n"
             "              With tfact 2, the curve is essentially a polyline\n"
             "              through the given data points.\n"
             "Plot a cubic polynomial through the points in the list. Each\n"
             "Hermite polynomial between two data points consists of 40 lines.\n"
             "See g2_splines.c for further information.\n"
             "Note : this two argument form is Python specific.");

static PyObject *
C_g2_raspln(G2 *self, PyObject *args)
{
   return helper_ild(self, args, g2_raspln);
}

PyDoc_STRVAR(doc_g2_filled_raspln,
             "g2_filled_raspln(list points, float tfact)\n"
             "As g2_raspln, but filled.");

static PyObject *
C_g2_filled_raspln(G2 *self, PyObject *args)
{
   return helper_ild(self, args, g2_filled_raspln);
}

PyDoc_STRVAR(doc_g2_para_3,
             "g2_para_3(list points)\n"
             "   'points' : [x1, y1, x2, y2, ... xn, yn]\n"
             "              In the list ints and floats can be mixed freely.\n"
             "Plot a piecewise parametric interpolation polynomial of degree 3\n"
             "through the given points, using Newton's Divided Differences method.\n"
             "Note : this one argument form is Python specific.");

static PyObject *
C_g2_para_3(G2 *self, PyObject *args)
{
   return helper_il(self, args, g2_para_3);
}

PyDoc_STRVAR(doc_g2_filled_para_3,
             "g2_filled_para_3(list points)\n"
             "As g2_para_3, but filled.");

static PyObject *
C_g2_filled_para_3(G2 *self, PyObject *args)
{
   return helper_il(self, args, g2_filled_para_3);
}

PyDoc_STRVAR(doc_g2_para_5,
             "g2_para_5(list points)\n"
             "As g2_para_3, but with degree 5 instead of 3.");

static PyObject *
C_g2_para_5(G2 *self, PyObject *args)
{
   return helper_il(self, args, g2_para_5);
}

PyDoc_STRVAR(doc_g2_filled_para_5,
             "g2_filled_para_5(list points)\n"
             "As g2_para_5, but filled.");

static PyObject *
C_g2_filled_para_5(G2 *self, PyObject *args)
{
   return helper_il(self, args, g2_filled_para_5);
}

static PyMethodDef G2_methods[] = {
   /* methods that operate on a single device */
   { "g2_attach", (PyCFunction)C_g2_attach, METH_VARARGS, doc_g2_attach },
   { "g2_detach", (PyCFunction)C_g2_detach, METH_VARARGS, doc_g2_detach },
   { "g2_ink", (PyCFunction)C_g2_ink, METH_VARARGS, doc_g2_ink },
   /* methods that can operate on multiple devices */
   { "g2_close", (PyCFunction)C_g2_close, METH_NOARGS, doc_g2_close },
   { "g2_set_auto_flush", (PyCFunction)C_g2_set_auto_flush, METH_VARARGS, doc_g2_set_auto_flush },
   { "g2_flush", (PyCFunction)C_g2_flush, METH_NOARGS, doc_g2_flush },
   { "g2_save", (PyCFunction)C_g2_save, METH_NOARGS, doc_g2_save },
   { "g2_set_coordinate_system", (PyCFunction)C_g2_set_coordinate_system, METH_VARARGS, doc_g2_set_coordinate_system },
   { "g2_pen", (PyCFunction)C_g2_pen, METH_VARARGS, doc_g2_pen },
   { "g2_set_dash", (PyCFunction)C_g2_set_dash, METH_VARARGS, doc_g2_set_dash },
   { "g2_set_solid", (PyCFunction)C_g2_set_solid, METH_NOARGS, doc_g2_set_solid },
   { "g2_set_font_size", (PyCFunction)C_g2_set_font_size, METH_VARARGS, doc_g2_set_font_size },
   { "g2_set_line_width", (PyCFunction)C_g2_set_line_width, METH_VARARGS, doc_g2_set_line_width },
   { "g2_clear_palette", (PyCFunction)C_g2_clear_palette, METH_NOARGS, doc_g2_clear_palette },
   { "g2_reset_palette", (PyCFunction)C_g2_reset_palette, METH_NOARGS, doc_g2_reset_palette },
   { "g2_allocate_basic_colors", (PyCFunction)C_g2_allocate_basic_colors, METH_NOARGS, doc_g2_allocate_basic_colors },
   { "g2_clear", (PyCFunction)C_g2_clear, METH_NOARGS, doc_g2_clear },
   { "g2_set_background", (PyCFunction)C_g2_set_background, METH_VARARGS, doc_g2_set_background },
   { "g2_move", (PyCFunction)C_g2_move, METH_VARARGS, doc_g2_move },
   { "g2_move_r", (PyCFunction)C_g2_move_r, METH_VARARGS, doc_g2_move_r },
   { "g2_plot", (PyCFunction)C_g2_plot, METH_VARARGS, doc_g2_plot },
   { "g2_plot_r", (PyCFunction)C_g2_plot_r, METH_VARARGS, doc_g2_plot_r },
   { "g2_line", (PyCFunction)C_g2_line, METH_VARARGS, doc_g2_line },
   { "g2_line_r", (PyCFunction)C_g2_line_r, METH_VARARGS, doc_g2_line_r },
   { "g2_line_to", (PyCFunction)C_g2_line_to, METH_VARARGS, doc_g2_line_to },
   { "g2_poly_line", (PyCFunction)C_g2_poly_line, METH_VARARGS, doc_g2_poly_line },
   { "g2_triangle", (PyCFunction)C_g2_triangle, METH_VARARGS, doc_g2_triangle },
   { "g2_filled_triangle", (PyCFunction)C_g2_filled_triangle, METH_VARARGS, doc_g2_filled_triangle },
   { "g2_rectangle", (PyCFunction)C_g2_rectangle, METH_VARARGS, doc_g2_rectangle },
   { "g2_filled_rectangle", (PyCFunction)C_g2_filled_rectangle, METH_VARARGS, doc_g2_filled_rectangle },
   { "g2_polygon", (PyCFunction)C_g2_polygon, METH_VARARGS, doc_g2_polygon },
   { "g2_filled_polygon", (PyCFunction)C_g2_filled_polygon, METH_VARARGS, doc_g2_filled_polygon },
   { "g2_circle", (PyCFunction)C_g2_circle, METH_VARARGS, doc_g2_circle },
   { "g2_filled_circle", (PyCFunction)C_g2_filled_circle, METH_VARARGS, doc_g2_filled_circle },
   { "g2_ellipse", (PyCFunction)C_g2_ellipse, METH_VARARGS, doc_g2_ellipse },
   { "g2_filled_ellipse", (PyCFunction)C_g2_filled_ellipse, METH_VARARGS, doc_g2_filled_ellipse },
   { "g2_arc", (PyCFunction)C_g2_arc, METH_VARARGS, doc_g2_arc },
   { "g2_filled_arc", (PyCFunction)C_g2_filled_arc, METH_VARARGS, doc_g2_filled_arc },
   { "g2_string", (PyCFunction)C_g2_string, METH_VARARGS, doc_g2_string },
   { "g2_image", (PyCFunction)C_g2_image, METH_VARARGS, doc_g2_image },
   { "g2_set_QP", (PyCFunction)C_g2_set_QP, METH_VARARGS, doc_g2_set_QP },
   { "g2_plot_QP", (PyCFunction)C_g2_plot_QP, METH_VARARGS, doc_g2_plot_QP },
   { "g2_query_pointer", (PyCFunction)C_g2_query_pointer, METH_NOARGS, doc_g2_query_pointer },
   { "g2_spline", (PyCFunction)C_g2_spline, METH_VARARGS, doc_g2_spline },
   { "g2_filled_spline", (PyCFunction)C_g2_filled_spline, METH_VARARGS, doc_g2_filled_spline },
   { "g2_b_spline", (PyCFunction)C_g2_b_spline, METH_VARARGS, doc_g2_b_spline },
   { "g2_filled_b_spline", (PyCFunction)C_g2_filled_b_spline, METH_VARARGS, doc_g2_filled_b_spline },
   { "g2_raspln", (PyCFunction)C_g2_raspln, METH_VARARGS, doc_g2_raspln },
   { "g2_filled_raspln", (PyCFunction)C_g2_filled_raspln, METH_VARARGS, doc_g2_filled_raspln },
   { "g2_para_3", (PyCFunction)C_g2_para_3, METH_VARARGS, doc_g2_para_3 },
   { "g2_filled_para_3", (PyCFunction)C_g2_filled_para_3, METH_VARARGS, doc_g2_filled_para_3 },
   { "g2_para_5", (PyCFunction)C_g2_para_5, METH_VARARGS, doc_g2_para_5 },
   { "g2_filled_para_5", (PyCFunction)C_g2_filled_para_5, METH_VARARGS, doc_g2_filled_para_5 },
   { NULL }
};

static PyMemberDef G2_members[] = { /* cf. structmember.h */
   { "dev", T_INT, offsetof(G2, dev), READONLY, "g2 device number\n"
     "readonly : set by g2_open_..." }, /* flags 0 : no restriction */
   { NULL }
};

static void
G2_dealloc(G2 *self)
{
   self->ob_type->tp_free((PyObject *)self);
}

static PyTypeObject G2_Type = {
   PyObject_HEAD_INIT(NULL)
   0,                      /* ob_size*/
   "g2.G2",                /* tp_name*/
   sizeof(G2),             /* tp_basicsize*/
   0,                      /* tp_itemsize*/
   (destructor)G2_dealloc, /* tp_dealloc*/
   0,                      /* tp_print*/
   0,                      /* tp_getattr*/
   0,                      /* tp_setattr*/
   0,                      /* tp_compare*/
   0,                      /* tp_repr*/
   0,                      /* tp_as_number*/
   0,                      /* tp_as_sequence*/
   0,                      /* tp_as_mapping*/
   0,                      /* tp_hash */
   0,                      /* tp_call*/
   0,                      /* tp_str*/
   0,                      /* tp_getattro*/
   0,                      /* tp_setattro*/
   0,                      /* tp_as_buffer*/
   Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /* tp_flags*/
   "G2 class",             /* tp_doc */
   0,                      /* tp_traverse */
   0,                      /* tp_clear */
   0,                      /* tp_richcompare */
   0,                      /* tp_weaklistoffset */
   0,                      /* tp_iter */
   0,                      /* tp_iternext */
   G2_methods,             /* tp_methods */
   G2_members              /* tp_members */
};

#ifndef PyMODINIT_FUNC     /* DLL import/export declarations */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initg2(void)
{
   PyObject *m;

   if ((PyType_Ready(&G2_Type) == 0) &&
       (m = Py_InitModule3("g2", module_functions,
                           "Python interface to the g2 library."))) {
      Py_INCREF(&G2_Type);
      PyModule_AddObject(m, "G2", (PyObject *)&G2_Type);
      add_enums(m);
   }
}
