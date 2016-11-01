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
#ifndef _G2_UTIL_H
#define _G2_UTIL_H

#include <stdlib.h>
#include "g2_physical_device.h"

int  dtoi(double x);

void g2_uc2pdc_int(g2_physical_device *pd, double x, double y,
		   int *ix, int *iy);
void g2_uc2pdc_double(g2_physical_device *pd, double x, double y,
		      double *dx, double *dy);
void g2_us2pds_int(g2_physical_device *pd, double x, double y,
		   int *ix, int *iy);
void g2_us2pds_double(g2_physical_device *pd, double x, double y,
		      double *dx, double *dy);
void g2_pdc2uc(g2_physical_device *pd, double ix, double iy,
		   double *x, double *y);
void g2_sort2_i(int *a, int *b);
void g2_sort2_d(double *a, double *b);
void *g2_malloc(size_t size);
void *g2_realloc(void *p, size_t size);
void g2_free(void *p);

double *g2_floatp2doublep(float *f, int N);

enum g2_log_level {Error=1, Warning, Verbose, Debug};
void g2_log(enum g2_log_level log_level, const char *format, ...);

#endif /* _G2_UTIL_H */
