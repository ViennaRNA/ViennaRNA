/*****************************************************************************
**  Copyright (C) 1998-2001  Ljubomir Milanovic & Horst Wagner
**  This file is part of the g2 library
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Lesser General Public
**  License as published by the Free Software Foundation; either
**  version 2.1 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Lesser General Public License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License along with this library; if not, write to the Free Software
**  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
******************************************************************************/
#ifndef _G2_PHYSICAL_DEVICE_H
#define _G2_PHYSICAL_DEVICE_H

#include "g2.h"
#include "g2_funix.h"


typedef enum _g2_coor {		      /* coord. type */
    g2_IntCoor, g2_DoubleCoor
} g2_coor;


typedef struct _g2_funix_fun {	      /* funix--function paar */
    g2_funix  fx;		      /* function index */
    int       (*fun)();		      /* pointer to function */
} g2_funix_fun;


typedef struct _g2_physical_device {
    int           pid;		      /* physical device id */
    void          *pdp;		      /* pointer to something in phys. dev. */
    g2_coor       coor_type;	      /* coordinate type */
    g2_funix_fun  *ff;		      /* list of funix-function pairs */
    double        a11, a22;	      /* coordinate transformation (ud->pd) */
    double        b1,  b2;	      /*  Ar + B */

    double        x_origin;	      /* User coordinates specification */
    double        y_origin;
    double        x_mul;
    double        y_mul;
} g2_physical_device;



g2_physical_device *g2_create_physical_device(int pid,
					      void *pdp,
					      g2_coor ct,
					      const g2_funix_fun *ff,
					      double a11, double a22,
					      double b1,  double b2);
void g2_destroy_physical_device(g2_physical_device *pd);

#endif /* _G2_PHYSICAL_DEVICE_H */
