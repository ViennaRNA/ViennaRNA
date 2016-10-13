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

#include "g2_funix.h"
#include "g2_virtual_device.h"
#include "g2_graphic_pd.h"
#include "g2_control_pd.h"
#include "g2_util.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif /* PI */


/*
 *
 * Plot (physical device)
 *
 */
void g2_plot_pd(g2_physical_device *pd, double x, double y)
{
    int    ix, iy;
    double dx, dy;
    
    if(pd->ff[g2_Plot].fun!=NULL) {
	switch(pd->coor_type) {
	  case g2_IntCoor:
	    g2_uc2pdc_int(pd, x, y, &ix, &iy);
	    pd->ff[g2_Plot].fun(pd->pid, pd->pdp,
				ix, iy);
	    break;
	  case g2_DoubleCoor:
	    g2_uc2pdc_double(pd, x, y, &dx, &dy);
	    pd->ff[g2_Plot].fun(pd->pid, pd->pdp,
				dx, dy);
	    break;
	}
    } else {
	/* emulate ... with .... */  
    }
}



/*
 *
 * Draw line (low_level)
 *
 */
void g2_line_pd(g2_physical_device *pd,
		double x1, double y1, double x2, double y2)
{
    int    ix1, iy1, ix2, iy2;
    double dx1, dy1, dx2, dy2;
    
    if(pd->ff[g2_Line].fun!=NULL) {
	switch(pd->coor_type) {
	  case g2_IntCoor:
	    g2_uc2pdc_int(pd, x1, y1, &ix1, &iy1);
	    g2_uc2pdc_int(pd, x2, y2, &ix2, &iy2);
	    pd->ff[g2_Line].fun(pd->pid, pd->pdp,
				ix1, iy1, ix2, iy2);
	    break;
	  case g2_DoubleCoor:
	    g2_uc2pdc_double(pd, x1, y1, &dx1, &dy1);
	    g2_uc2pdc_double(pd, x2, y2, &dx2, &dy2);
	    pd->ff[g2_Line].fun(pd->pid, pd->pdp,
				dx1, dy1, dx2, dy2);
	    break;
	}		
    } else {
	/* emulate ... with .... */  
    }
}



/*
 *
 * Draw poly line (physical device)
 *
 */
void g2_poly_line_pd(g2_physical_device *pd,
		     int N, double *points)
{
    int i;
    int    *ipt;
    double *dpt;
    
    if(pd->ff[g2_PolyLine].fun!=NULL) {
	switch(pd->coor_type) {
	  case g2_IntCoor:
	    ipt=g2_malloc(2*N*sizeof(int));
	    for(i=0;i<2*N;i+=2)
		g2_uc2pdc_int(pd, points[i+0], points[i+1], ipt+i, ipt+i+1);
	    pd->ff[g2_PolyLine].fun(pd->pid, pd->pdp,
				    N, ipt);
	    free(ipt);
	    break;
	  case g2_DoubleCoor:
	    dpt=g2_malloc(2*N*sizeof(double));
	    for(i=0;i<2*N;i+=2)
		g2_uc2pdc_double(pd,
				  points[i+0], points[i+1],
				  dpt+i, dpt+i+1);
	    pd->ff[g2_PolyLine].fun(pd->pid, pd->pdp,
				    N, dpt);
	    free(dpt);
	    break;
	}
    } else {
	for (i=0;i<N-1;i++)		 /* emulate polyline with lines */  
	    g2_line_pd(pd,
		       points[2*i], points[2*i+1],
		       points[2*i+2], points[2*i+3]);
    }
}



/*
 *
 * Triangle (physical device)
 *
 */
void g2_triangle_pd(g2_physical_device *pd,
		    double x1, double y1,
		    double x2, double y2,
		    double x3, double y3)
{
    int    ix1, iy1, ix2, iy2, ix3, iy3;
    double dx1, dy1, dx2, dy2, dx3, dy3;
    
    if(pd->ff[g2_Triangle].fun!=NULL) {
	switch(pd->coor_type) {
	  case g2_IntCoor:
	    g2_uc2pdc_int(pd, x1, y1, &ix1, &iy1);
	    g2_uc2pdc_int(pd, x2, y2, &ix2, &iy2);
	    g2_uc2pdc_int(pd, x3, y3, &ix3, &iy3);
	    pd->ff[g2_Triangle].fun(pd->pid, pd->pdp,
				    ix1, iy1,
				    ix2, iy2,
				    ix3, iy3);
	    break;
	  case g2_DoubleCoor:
	    g2_uc2pdc_double(pd, x1, y1, &dx1, &dy1);
	    g2_uc2pdc_double(pd, x2, y2, &dx2, &dy2);
	    g2_uc2pdc_double(pd, x3, y3, &dx3, &dy3);
	    pd->ff[g2_Triangle].fun(pd->pid, pd->pdp,
				    dx1, dy1,
				    dx2, dy2,
				    dx3, dy3);
	}		
    } else {
	g2_line_pd(pd, x1, y1, x2, y2);      /* emulate triangle with lines */
	g2_line_pd(pd, x2, y2, x3, y3);
	g2_line_pd(pd, x3, y3, x1, y1);
    }
}



/*
 *
 * Filled triangle (physical device)
 *
 */
void g2_filled_triangle_pd(g2_physical_device *pd,
			   double x1, double y1,
			   double x2, double y2,
			   double x3, double y3)
{
    int    ix1, iy1, ix2, iy2, ix3, iy3;
    double dx1, dy1, dx2, dy2, dx3, dy3;
    
    if(pd->ff[g2_FilledTriangle].fun!=NULL) {
	switch(pd->coor_type) {
	  case g2_IntCoor:
	    g2_uc2pdc_int(pd, x1, y1, &ix1, &iy1);
	    g2_uc2pdc_int(pd, x2, y2, &ix2, &iy2);
	    g2_uc2pdc_int(pd, x3, y3, &ix3, &iy3);
	    pd->ff[g2_FilledTriangle].fun(pd->pid, pd->pdp,
					  ix1, iy1,
					  ix2, iy2,
					  ix3, iy3);
	    break;
	  case g2_DoubleCoor:
	    g2_uc2pdc_double(pd, x1, y1, &dx1, &dy1);
	    g2_uc2pdc_double(pd, x2, y2, &dx2, &dy2);
	    g2_uc2pdc_double(pd, x3, y3, &dx3, &dy3);
	    pd->ff[g2_FilledTriangle].fun(pd->pid, pd->pdp,
					  dx1, dy1,
					  dx2, dy2,
					  dx3, dy3);
	}		
    } else {
	double Triangle[6];     /* emulate FilledTriangle with FilledPolygon */
	Triangle[0] = x1;
	Triangle[1] = y1;
	Triangle[2] = x2;
	Triangle[3] = y2;
	Triangle[4] = x3;
	Triangle[5] = y3;
	g2_filled_polygon_pd(pd, 3, Triangle); 
    }
}



/*
 *
 * Rectangle (physical device)
 *
 */
void g2_rectangle_pd(g2_physical_device *pd,
		     double x1, double y1, double x2, double y2)
{
    int    ix1, iy1, ix2, iy2;
    double dx1, dy1, dx2, dy2;
    
    if(pd->ff[g2_Rectangle].fun!=NULL) {
	switch(pd->coor_type) {
	  case g2_IntCoor:
	    g2_uc2pdc_int(pd, x1, y1, &ix1, &iy1);
	    g2_uc2pdc_int(pd, x2, y2, &ix2, &iy2);
	    g2_sort2_i(&ix1, &ix2);
	    g2_sort2_i(&iy1, &iy2);
	    pd->ff[g2_Rectangle].fun(pd->pid, pd->pdp,
				     ix1, iy1, ix2, iy2);
	    break;
	  case g2_DoubleCoor:
	    g2_uc2pdc_double(pd, x1, y1, &dx1, &dy1);
	    g2_uc2pdc_double(pd, x2, y2, &dx2, &dy2);
	    g2_sort2_d(&dx1, &dx2);
	    g2_sort2_d(&dy1, &dy2);
	    pd->ff[g2_Rectangle].fun(pd->pid, pd->pdp,
				     dx1, dy1, dx2, dy2);
	    break;
	}		
    } else {
	g2_line_pd(pd, x1, y1, x1, y2);     /* emulate rectangle with lines */
	g2_line_pd(pd, x1, y1, x2, y1);
	g2_line_pd(pd, x2, y1, x2, y2);
	g2_line_pd(pd, x1, y2, x2, y2);
    }
}



/*
 *
 * Filled rectangle (physical device)
 *
 */
void g2_filled_rectangle_pd(g2_physical_device *pd,
			    double x1, double y1, double x2, double y2)
{
    int    ix1, iy1, ix2, iy2;
    double dx1, dy1, dx2, dy2;
    
    if(pd->ff[g2_FilledRectangle].fun!=NULL) {
	switch(pd->coor_type) {
	  case g2_IntCoor:
	    g2_uc2pdc_int(pd, x1, y1, &ix1, &iy1);
	    g2_uc2pdc_int(pd, x2, y2, &ix2, &iy2);
	    g2_sort2_i(&ix1, &ix2);
	    g2_sort2_i(&iy1, &iy2);
	    pd->ff[g2_FilledRectangle].fun(pd->pid, pd->pdp,
					   ix1, iy1, ix2, iy2);
	    break;
	  case g2_DoubleCoor:
	    g2_uc2pdc_double(pd, x1, y1, &dx1, &dy1);
	    g2_uc2pdc_double(pd, x2, y2, &dx2, &dy2);
	    g2_sort2_d(&dx1, &dx2);
	    g2_sort2_d(&dy1, &dy2);
	    pd->ff[g2_FilledRectangle].fun(pd->pid, pd->pdp,
					   dx1, dy1, dx2, dy2);
	    break;
	}		
    } else {
	double Rectangle[8];  /* emulate FilledRectangle with FilledPolygon */
	Rectangle[0] = x1;
	Rectangle[1] = y1;
	Rectangle[2] = x2;
	Rectangle[3] = y1;
	Rectangle[4] = x2;
	Rectangle[5] = y2;
	Rectangle[6] = x1;
	Rectangle[7] = y2;
	g2_filled_polygon_pd(pd,4,Rectangle);
    }
}



/*
 *
 * Polygon (physical device)
 *
 */
void g2_polygon_pd(g2_physical_device *pd,
		   int N, double *points)
{
    int i;
    int    *ipt;
    double *dpt;
    
    if(pd->ff[g2_Polygon].fun!=NULL) {
	switch(pd->coor_type) {
	  case g2_IntCoor:
	    ipt=g2_malloc(2*N*sizeof(int));
	    for(i=0;i<2*N;i+=2)
		g2_uc2pdc_int(pd, points[i+0], points[i+1], ipt+i, ipt+i+1);
	    pd->ff[g2_Polygon].fun(pd->pid, pd->pdp,
				   N, ipt);
	    free(ipt);
	    break;
	  case g2_DoubleCoor:
	    dpt=g2_malloc(2*N*sizeof(double));
	    for(i=0;i<2*N;i+=2)
		g2_uc2pdc_double(pd,
				  points[i+0], points[i+1], dpt+i, dpt+i+1);
	    pd->ff[g2_Polygon].fun(pd->pid, pd->pdp,
				   N, dpt);
	    free(dpt);
	    break;
	}
    } else {
	for (i=0;i<N-1;i++)       /* emulate polygon with lines */
	    g2_line_pd(pd,
		       points[2*i], points[2*i+1],
		       points[2*i+2], points[2*i+3]);
	g2_line_pd(pd,points[2*N-2], points[2*N-1],points[0], points[1]);
    }
}



/*
 *
 * Filled polygon (physical device)
 *
 */
void g2_filled_polygon_pd(g2_physical_device *pd,
			  int N, double *points)
{
    int i;
    int    *ipt;
    double *dpt;

    if(pd->ff[g2_FilledPolygon].fun!=NULL) {
	switch(pd->coor_type) {
	  case g2_IntCoor:
	    ipt=g2_malloc(2*N*sizeof(int));
	    for(i=0;i<2*N;i+=2)
		g2_uc2pdc_int(pd, points[i+0], points[i+1], ipt+i, ipt+i+1);
	    pd->ff[g2_FilledPolygon].fun(pd->pid, pd->pdp,
					 N, ipt);
	    free(ipt);
	    break;
	  case g2_DoubleCoor:
	    dpt=g2_malloc(2*N*sizeof(double));
	    for(i=0;i<2*N;i+=2)
		g2_uc2pdc_double(pd,
				  points[i+0], points[i+1], dpt+i, dpt+i+1);
	    pd->ff[g2_FilledPolygon].fun(pd->pid, pd->pdp,
					 N, dpt);
	    free(dpt);
	    break;
	}
    } else {
	/* emulate filled polygon with .... */  
    }
}



/*
 *
 * Ellipse (physical device)
 *
 */
void g2_ellipse_pd(g2_physical_device *pd,
		   double x, double y, double r1, double r2)
{
    int    ix, iy, ir1, ir2;
    double dx, dy, dr1, dr2;

    if(pd->ff[g2_Ellipse].fun!=NULL) {
	switch(pd->coor_type) {
	  case g2_IntCoor:
	    g2_uc2pdc_int(pd, x, y, &ix, &iy);
	    g2_us2pds_int(pd, r1, r2, &ir1, &ir2);
	    pd->ff[g2_Ellipse].fun(pd->pid, pd->pdp,
				   ix, iy, ir1, ir2);
	    break;
	  case g2_DoubleCoor:
	    g2_uc2pdc_double(pd, x, y, &dx, &dy);
	    g2_us2pds_double(pd, r1, r2, &dr1, &dr2);
	    pd->ff[g2_Ellipse].fun(pd->pid, pd->pdp,
				   dx, dy, dr1, dr2);
	    break;
	}
    } else {
	g2_arc_pd(pd, x, y, r1, r2, 0., 360.);   /* emulate ellipse with arc */
    }
}



/*
 *
 * Filled ellipse (physical device)
 *
 */
void g2_filled_ellipse_pd(g2_physical_device *pd,
			  double x, double y, double r1, double r2)
{
    int    ix, iy, ir1, ir2;
    double dx, dy, dr1, dr2;

    if(pd->ff[g2_FilledEllipse].fun!=NULL) {
	switch(pd->coor_type) {
	  case g2_IntCoor:
	    g2_uc2pdc_int(pd, x, y, &ix, &iy);
	    g2_us2pds_int(pd, r1, r2, &ir1, &ir2);
	    pd->ff[g2_FilledEllipse].fun(pd->pid, pd->pdp,
					 ix, iy, ir1, ir2);
	    break;
	  case g2_DoubleCoor:
	    g2_uc2pdc_double(pd, x, y, &dx, &dy);
	    g2_us2pds_double(pd, r1, r2, &dr1, &dr2);
	    pd->ff[g2_FilledEllipse].fun(pd->pid, pd->pdp,
					 dx, dy, dr1, dr2);
	    break;
	}
    } else {
	g2_filled_arc_pd(pd,        /* emulate filledellipse with filled arc */
			 x, y,
			 r1, r2,
			 0., 360.);
    }
}



/*
 *
 * Circle (physical device)
 *
 */
void g2_circle_pd(g2_physical_device *pd,
		  double x, double y, double r)
{
    int    ix, iy, ir;
    double dx, dy, dr;

    if(pd->ff[g2_Circle].fun!=NULL) {
	switch(pd->coor_type) {
	  case g2_IntCoor:
	    g2_uc2pdc_int(pd, x, y, &ix, &iy);
	    g2_us2pds_int(pd, r, 0, &ir, NULL);
	    pd->ff[g2_Circle].fun(pd->pid, pd->pdp,
				  ix, iy, ir);
	    break;
	  case g2_DoubleCoor:
	    g2_uc2pdc_double(pd, x, y, &dx, &dy);
	    g2_us2pds_double(pd, r, 0, &dr, NULL);
	    pd->ff[g2_Circle].fun(pd->pid, pd->pdp,
				  dx, dy, dr);
	    break;
	}
    } else {
	g2_ellipse_pd(pd, x, y, r, r);	     /* emulate circle with ellipse */
    }   
}



/*
 *
 * Filled circle (physical device)
 *
 */
void g2_filled_circle_pd(g2_physical_device *pd,
			 double x, double y, double r)
{
    int    ix, iy, ir;
    double dx, dy, dr;

    if(pd->ff[g2_FilledCircle].fun!=NULL) {
	switch(pd->coor_type) {
	  case g2_IntCoor:
	    g2_uc2pdc_int(pd, x, y, &ix, &iy);
	    g2_us2pds_int(pd, r, 0, &ir, NULL);
	    pd->ff[g2_FilledCircle].fun(pd->pid, pd->pdp,
					ix, iy, ir);
	    break;
	  case g2_DoubleCoor:
	    g2_uc2pdc_double(pd, x, y, &dx, &dy);
	    g2_us2pds_double(pd, r, 0, &dr, NULL);
	    pd->ff[g2_FilledCircle].fun(pd->pid, pd->pdp,
					dx, dy, dr);
	    break;
	}
    } else {
	g2_filled_ellipse_pd(pd, x, y, r, r);           /* emulate */
    }   
}



/*
 *
 * Arc (physical device)
 *
 */
void g2_arc_pd(g2_physical_device *pd,
	       double x, double y, double r1, double r2, double a1, double a2)
{
    int    ix, iy, ir1, ir2;
    double dx, dy, dr1, dr2;
	
    if(pd->ff[g2_Arc].fun!=NULL) {
	switch(pd->coor_type) {
	  case g2_IntCoor:
	    g2_uc2pdc_int(pd, x, y, &ix, &iy);
	    g2_us2pds_int(pd, r1, r2, &ir1, &ir2);
	    pd->ff[g2_Arc].fun(pd->pid, pd->pdp,
			       ix, iy, ir1, ir2, a1, a2);
	    break;
	  case g2_DoubleCoor:
	    g2_uc2pdc_double(pd, x, y, &dx, &dy);
	    g2_us2pds_double(pd, r1, r2, &dr1, &dr2);
	    pd->ff[g2_Arc].fun(pd->pid, pd->pdp,
			       dx, dy, dr1, dr2, a1, a2);
	    break;
	}
    } else {
	double a, da, *pt;			  /* emulate arc */
	int N, i;
	N=(a2==a1)?360:(int)fabs(a2-a1)+8;
	a=a1*2.*PI/360.;
	da=((a2>a1)? (a2-a1):360.-(a1-a2))*2.*PI/360./(N-1);
	pt=g2_malloc(2*N*sizeof(double));
	for(i=0;i<N;i++) {
	    pt[2*i+0]=x+r1*cos(a+i*da);
	    pt[2*i+1]=y+r2*sin(a+i*da);
	}
	g2_poly_line_pd(pd, N, pt);		  /* using poly_line */
	g2_free(pt);
    }
}



/*
 *
 * Filled arc (physical device)
 *
 */
void g2_filled_arc_pd(g2_physical_device *pd,
		      double x, double y,
		      double r1, double r2,
		      double a1, double a2)
{
    int    ix, iy, ir1, ir2;
    double dx, dy, dr1, dr2;
    	
    if(pd->ff[g2_FilledArc].fun!=NULL) {
	switch(pd->coor_type) {
	  case g2_IntCoor:
	    g2_uc2pdc_int(pd, x, y, &ix, &iy);
	    g2_us2pds_int(pd, r1, r2, &ir1, &ir2);
	    pd->ff[g2_FilledArc].fun(pd->pid, pd->pdp,
				     ix, iy, ir1, ir2, a1, a2);
	    break;
	  case g2_DoubleCoor:
	    g2_uc2pdc_double(pd, x, y, &dx, &dy);
	    g2_us2pds_double(pd, r1, r2, &dr1, &dr2);
	    pd->ff[g2_FilledArc].fun(pd->pid, pd->pdp,
				     dx, dy, dr1, dr2, a1, a2);
	    break;
	}
    } else {
	double a, da, *pt;			  /* emulate filled arc */
	int N, i;
	N=(a2==a1)?360:(int)fabs(a2-a1)+8;
	a=a1*2.*PI/360.;
	da=((a2>a1)? (a2-a1):360.-(a1-a2))*2.*PI/360./(N-1);
	pt=g2_malloc(2*(N+2)*sizeof(double));
	pt[0]=x;
	pt[1]=y;
	for(i=0;i<N;i++) {
	    pt[2*i+2]=x+r1*cos(a+i*da);
	    pt[2*i+3]=y+r2*sin(a+i*da);
	}
	pt[2*N+2]=x;
	pt[2*N+3]=y;
	g2_filled_polygon_pd(pd, N+2, pt);	  /* using filled polygon */
	g2_free(pt);
    }
}



/*
 *
 * Draw string (physical device)
 *
 */
void g2_string_pd(g2_physical_device *pd,
		  double x, double y, const char *text)
{
    int    ix, iy;
    double dx, dy;

    if(pd->ff[g2_String].fun!=NULL) {
	switch(pd->coor_type) {
	  case g2_IntCoor:
	    g2_uc2pdc_int(pd, x, y, &ix, &iy);
	    pd->ff[g2_String].fun(pd->pid, pd->pdp,
				  ix, iy, text);
	    break;
	  case g2_DoubleCoor:
	    g2_uc2pdc_double(pd, x, y, &dx, &dy);
	    pd->ff[g2_String].fun(pd->pid, pd->pdp,
				  dx, dy, text);
	    break;
	}
    } else {
	/* emulate ... with .... */  
    }
}



void g2_image_pd(g2_physical_device *pd,
	      double x, double y, int x_size, int y_size, int *pens)
{
    int    ix, iy;
    double dx, dy;

    if(pd->ff[g2_Image].fun!=NULL) {
	switch(pd->coor_type) {
	  case g2_IntCoor:
	    g2_uc2pdc_int(pd, x, y, &ix, &iy);
	    pd->ff[g2_Image].fun(pd->pid, pd->pdp,
				 ix, iy, x_size, y_size, pens);
	    break;
	  case g2_DoubleCoor:
	    g2_uc2pdc_double(pd, x, y, &dx, &dy);
	    pd->ff[g2_Image].fun(pd->pid, pd->pdp,
				 dx, dy, x_size, y_size, pens);
	    break;
	}
    } else {
	for(ix=0;ix<x_size;ix++) 
	    for(iy=0;iy<y_size;iy++) {
		g2_pen_pd(pd, pens[ix+x_size*iy]);
		g2_plot_pd(pd, ix+x, (y_size-iy-1)+y);
	    }
    }
}

