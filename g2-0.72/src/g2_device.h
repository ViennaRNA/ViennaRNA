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
#ifndef _G2_DEVICE_H
#define _G2_DEVICE_H

#include "g2_physical_device.h"
#include "g2_virtual_device.h"

extern int __g2_last_device;	    /* last accessed device (ld) */


typedef enum _g2_dev_type {
    g2_ILLEGAL=-1,		    /* illegal device type */
    g2_NDEV=0,			    /* no device associated */
    g2_PD,			    /* physical device */
    g2_VD			    /* virtual device */
} g2_device_type;


typedef struct _g2_device {
    g2_device_type t;		    /* device type */
    int            dix;		    /* dev index in array (g2_dev) */
    union {
	g2_physical_device *pd;
	g2_virtual_device  *vd;
    } d;			    /* device */
    double        x;		    /* graphical cursor */
    double        y;
    int           auto_flush;	    /* 1-on 0-off */
    double        QPd;		    /* Quasi pixel spec. */
    enum QPshape  QPshape;
} g2_device;


int g2_register_physical_device(int pid,
				void *pdp,
				g2_coor ct,
				const g2_funix_fun *ff,
				double a11, double a22,
				double b1,  double b2);
int g2_register_virtual_device();

g2_device *g2_get_device_pointer(int dix);
g2_device_type g2_get_device_type(int dix);
void g2_destroy_device(int dix);
			       

#endif /* _G2_DEVICE_H */
