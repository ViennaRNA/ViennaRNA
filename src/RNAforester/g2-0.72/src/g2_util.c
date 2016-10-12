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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

#include "g2_util.h"
#include "g2_physical_device.h"
#include "g2_config.h"


/*
 *
 * Double to integer
 *
 */
int dtoi(double x)
{
    return (int)(x+0.5);
}



/*
 *
 * Transform user coord. in physical device coord (int)
 *
 */
void g2_uc2pdc_int(g2_physical_device *pd, double x, double y,
		   int *ix, int *iy)
{
    double vx, vy;
    
    vx=pd->x_mul*x+pd->x_origin;
    vy=pd->y_mul*y+pd->y_origin;
    
    *ix = dtoi(pd->a11*vx + pd->b1);
    *iy = dtoi(pd->a22*vy + pd->b2);
}


/*
 *
 * Transform user coord. in physical device coord (double)
 *
 */
void g2_uc2pdc_double(g2_physical_device *pd, double x, double y,
		      double *dx, double *dy)
{
    double vx, vy;
    
    vx=pd->x_mul*x+pd->x_origin;
    vy=pd->y_mul*y+pd->y_origin;
    
    *dx = pd->a11*vx + pd->b1;
    *dy = pd->a22*vy + pd->b2;
}


/*
 *
 * Transform user size in physical device size (int)
 *
 */
void g2_us2pds_int(g2_physical_device *pd, double x, double y,
		   int *ix, int *iy)
{
    if(ix!=NULL)
	*ix=dtoi(x*fabs(pd->x_mul*pd->a11));
    if(iy!=NULL)
	*iy=dtoi(y*fabs(pd->y_mul*pd->a22));
}


/*
 *
 * Transform user size in physical device size (double)
 *
 */
void g2_us2pds_double(g2_physical_device *pd, double x, double y,
		      double *dx, double *dy)
{
    if(dx!=NULL)
	*dx=x*fabs(pd->x_mul*pd->a11);
    if(dy!=NULL)
	*dy=y*fabs(pd->y_mul*pd->a22);
}


/*
 *
 * Transform physical device coord in user coord
 *
 */
void g2_pdc2uc(g2_physical_device *pd, double ix, double iy,
		   double *x, double *y)
{
    double pcx, pcy;
    pcx=(ix-pd->b1)/pd->a11;
    pcy=(iy-pd->b2)/pd->a22;

    *x=(pcx-pd->x_origin)/pd->x_mul;
    *y=(pcy-pd->y_origin)/pd->y_mul;
}


/*
 *  return a < b
 */
void g2_sort2_i(int *a, int *b)
{
    if(*a>*b) {
	int t=*a;
	*a=*b; *b=t;
    }
}

void g2_sort2_d(double *a, double *b)
{
    if(*a>*b) {
	double t=*a;
	*a=*b; *b=t;
    }
}



/*
 *
 * g2 malloc (with error message)
 *
 */
void *g2_malloc(size_t size)
{
    void *rv;

    if((rv=malloc(size))==NULL) {
	fprintf(stderr, "g2_malloc: Can not allocate memory\n");
	exit(-1);
    }

    return rv;
}


/*
 *
 * g2 realloc (with error message)
 *
 */
void *g2_realloc(void *p, size_t size)
{
    void *rv;
    if((rv=realloc(p, size))==NULL) {
	fprintf(stderr, "g2_realloc: Can not allocate memory\n");
	exit(-1);
    }
    return rv;
}



/*
 *
 * g2 free
 *
 */
void g2_free(void *p)
{
    if(p!=NULL)
	free(p);
}


/*
 *
 * transform float* to double* for N elements
 *
 * Note: don't forget to free d
 *
 */
double *g2_floatp2doublep(float *f, int N)
{
    int i;
    double *d;
    d=(double *)g2_malloc(N*sizeof(double));
    for(i=0;i<N;++i)
	d[i]=f[i];
    return d;
}


/*
 *
 * log messages to stderr
 *
 */
void g2_log(enum g2_log_level log_level, const char *format, ...)
{
    va_list arg;
    if(log_level > g2_LogLevel) {
	return;
    }
    va_start(arg, format);
    vfprintf(stderr, format, arg);
    va_end(arg);    
}
