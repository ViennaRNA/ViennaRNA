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

#include "g2_device.h"
#include "g2_graphic_pd.h"
#include "g2_util.h"

/**
 * \ingroup interface
 * \defgroup graphic graphical output
 */

/**
 *
 * Move graphic cursor.
 *
 * \param dev device
 * \param x x coordinate
 * \param y y coordinate
 *
 * \ingroup graphic
 */
void g2_move(int dev, double x, double y)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_move: No such device: %d\n", dev);
	return;
    }
    
    devp->x=x;					  /* set graph. cursor */
    devp->y=y;

    switch(devp->t) {
      case g2_PD:
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_move(devp->d.vd->dix[i], x, y);
	break;
      case g2_ILLEGAL:
	break;
      case g2_NDEV:
	break;
    }
    __g2_last_device=dev;
}



/**
 *
 * Move graphic cursor relative to the currner graphical cursor position.
 *
 * \param dev device
 * \param dx x coordinate increment
 * \param dy y coordinate increment
 *
 * \ingroup graphic
 */
void g2_move_r(int dev, double dx, double dy)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_move_r: No such device: %d\n", dev);
	return;
    }

    devp->x+=dx;				  /* set graph. cursor */
    devp->y+=dy;
    
    switch(devp->t) {
      case g2_PD:
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_move_r(devp->d.vd->dix[i], dx, dy);
	break;
      case g2_ILLEGAL:
	break;
      case g2_NDEV:
	break;
    }
    __g2_last_device=dev;
}



/**
 *
 * Plot a point
 *
 * \param dev device
 * \param x x coordinate
 * \param y y coordinate
 *
 * \ingroup graphic
 */
void g2_plot(int dev, double x, double y)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_plot: No such device: %d\n", dev);
	return;
    }
    
    devp->x=x;					  /* set graph. cursor */
    devp->y=y;
    
    switch(devp->t) {
      case g2_PD:
	g2_plot_pd(devp->d.pd, x, y);    
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_plot(devp->d.vd->dix[i], x, y);
	break;
      case g2_ILLEGAL:
	break;
      case g2_NDEV:
	break;
    }
    
    if(devp->auto_flush)
	g2_flush(dev);
    
    __g2_last_device=dev;
}



/**
 *
 * Plot a point relative to graphical cursor.
 *
 * \param dev device
 * \param rx relative x coordinate 
 * \param ry relative y coordinate
 *
 * \ingroup graphic
 */
void g2_plot_r(int dev, double rx, double ry)
{
    g2_device *devp;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_plot_r: No such device: %d\n", dev);
	return;
    }
    
    g2_plot(dev, devp->x+rx, devp->y+ry);
    
    __g2_last_device=dev;
}



/**
 *
 * Draw a line from \p x1, \p y1 to \p x2, \p y2.
 *
 * \param dev device
 * \param x1 see above
 * \param y1 see above
 * \param x2 see above
 * \param y2 see above
 *
 * \ingroup graphic
 */
void g2_line(int dev, double x1, double y1, double x2, double y2)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_line: No such device: %d\n", dev);
	return;
    }
    
    devp->x=x2;
    devp->y=y2;

    switch(devp->t) {
      case g2_PD:
	g2_line_pd(devp->d.pd, x1, y1, x2, y2);
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_line(devp->d.vd->dix[i], x1, y1, x2, y2);
	break;
      case g2_ILLEGAL:
	break;
      case g2_NDEV:
	break;
    }
    
    if(devp->auto_flush)
	g2_flush(dev);
    
    __g2_last_device=dev;
}



/**
 *
 * Draw line relative to the graphic cursor.
 *
 * \param dev device
 * \param dx relative x coordinate
 * \param dy relative y coordinate
 *
 * \ingroup graphic
 */
void g2_line_r(int dev, double dx, double dy)
{
    g2_device *devp;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_line_r: No such device: %d\n", dev);
	return;
    }
    g2_line(dev, devp->x, devp->y, devp->x+dx, devp->y+dy);

    __g2_last_device=dev;
}



/**
 *
 * Draw line from graphic cursor to the point \a x, \a y
 *
 * \param dev device
 * \param x x coordinate
 * \param y y coordinate
 *
 * \ingroup graphic
 */
void g2_line_to(int dev, double x, double y)
{
    g2_device *devp;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_line_to: No such device: %d\n", dev);
	return;
    }
    g2_line(dev, devp->x, devp->y, x, y);

    __g2_last_device=dev;
}



/**
 *
 * Draw a poly line.
 *
 * \param dev device
 * \param N_pt number of points (Note: It is not size of \a points vector!)
 * \param points vector of coordinates: x1, y1, x2, y2, ...
 *
 * \ingroup graphic
 */
void g2_poly_line(int dev, int N_pt, double *points)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_poly_line: No such device: %d\n", dev);
	return;
    }
    
    devp->x=points[2*(N_pt-1)+0];
    devp->y=points[2*(N_pt-1)+1];
    
    switch(devp->t) {
      case g2_PD:
	g2_poly_line_pd(devp->d.pd, N_pt, points);
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_poly_line(devp->d.vd->dix[i], N_pt, points);
	break;
      case g2_ILLEGAL:
	break;
      case g2_NDEV:
	break;
    }
    
    if(devp->auto_flush)
	g2_flush(dev);

    __g2_last_device=dev;
}



/**
 *
 * Draw a triangle described by 3 corner points.
 *
 * \param dev device
 * \param x1 x coordinate of the 1st corner
 * \param y1 y coordinate of the 1st corner
 * \param x2 x coordinate of the 2nd corner
 * \param y2 y coordinate of the 2nd corner
 * \param x3 x coordinate of the 3rd corner
 * \param y3 y coordinate of the 3rd corner
 *
 * \ingroup graphic
 */
void g2_triangle(int dev,
		 double x1, double y1,
		 double x2, double y2,
		 double x3, double y3)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_triangle: No such device: %d\n", dev);
	return;
    }
    
    devp->x=x3;
    devp->y=y3;
    
    switch(devp->t) {
      case g2_PD:
	g2_triangle_pd(devp->d.pd,
		       x1, y1,
		       x2, y2,
		       x3, y3);
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_triangle(devp->d.vd->dix[i], x1, y1, x2, y2, x3, y3);
	break;
      case g2_ILLEGAL:
	break;
      case g2_NDEV:
	break;
    }
    
    if(devp->auto_flush)
	g2_flush(dev);

    __g2_last_device=dev;
}



/**
 *
 * Draw a filled triangle specified by the 3 corner points.
 *
 * \param dev device
 * \param x1 x coordinate of the 1st corner
 * \param y1 y coordinate of the 1st corner
 * \param x2 x coordinate of the 2nd corner
 * \param y2 y coordinate of the 2nd corner
 * \param x3 x coordinate of the 3rd corner
 * \param y3 y coordinate of the 3rd corner
 *
 * \ingroup graphic
 */
void g2_filled_triangle(int dev,
			double x1, double y1,
			double x2, double y2,
			double x3, double y3)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_filled_triangle: No such device: %d\n", dev);
	return;
    }
    
    devp->x=x3;
    devp->y=y3;
    
    switch(devp->t) {
      case g2_PD:
	g2_filled_triangle_pd(devp->d.pd,
			      x1, y1,
			      x2, y2,
			      x3, y3);
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_filled_triangle(devp->d.vd->dix[i], x1, y1, x2, y2, x3, y3);
	break;
      case g2_ILLEGAL:
	break;
      case g2_NDEV:
	break;
    }
    
    if(devp->auto_flush)
	g2_flush(dev);

    __g2_last_device=dev;
}



/**
 *
 * Draw a rectangle specified by the two opposite corner points.
 *
 * \param dev device
 * \param x1 x coordinate of the 1st corner
 * \param y1 y coordinate of the 1st corner
 * \param x2 x coordinate of the 3rd corner
 * \param y2 y coordinate of the 3rd corner
 *
 * \ingroup graphic
 */
void g2_rectangle(int dev, double x1, double y1, double x2, double y2)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_rectangle: No such device: %d\n", dev);
	return;
    }
    
    devp->x=x2;
    devp->y=y2;

    switch(devp->t) {
      case g2_PD:
	g2_rectangle_pd(devp->d.pd, x1, y1, x2, y2);
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_rectangle(devp->d.vd->dix[i], x1, y1, x2, y2);
	break;
      case g2_ILLEGAL:
	break;
      case g2_NDEV:
	break;
    }
    
    if(devp->auto_flush)
	g2_flush(dev);

    __g2_last_device=dev;
}



/**
 *
 * Draw a filled rectangle specified by the two opposite corner points.
 *
 * \param dev device
 * \param x1 x coordinate of the 1st corner
 * \param y1 y coordinate of the 1st corner
 * \param x2 x coordinate of the 3rd corner
 * \param y2 y coordinate of the 3rd corner
 *
 * \ingroup graphic
 */
void g2_filled_rectangle(int dev, double x1, double y1, double x2, double y2)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_filled_rectangle: No such device: %d\n", dev);
	return;
    }
    
    devp->x=x2;
    devp->y=y2;

    switch(devp->t) {
      case g2_PD:
	g2_filled_rectangle_pd(devp->d.pd, x1, y1, x2, y2);
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_filled_rectangle(devp->d.vd->dix[i], x1, y1, x2, y2);
	break;
      case g2_ILLEGAL:
	break;
      case g2_NDEV:
	break;
    }
    
    if(devp->auto_flush)
	g2_flush(dev);

    __g2_last_device=dev;
}



/**
 *
 * Draw a polygon.
 *
 * \param dev device
 * \param N_pt number of points (Note: It is not size of \a points vector!)
 * \param points vector of coordinates: x1, y1, x2, y2, ...
 *
 * \ingroup graphic
 */
void g2_polygon(int dev, int N_pt, double *points)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_polygon: No such device: %d\n", dev);
	return;
    }
    
    switch(devp->t) {
      case g2_PD:
	g2_polygon_pd(devp->d.pd, N_pt, points);
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_polygon(devp->d.vd->dix[i], N_pt, points);
	break;
      case g2_ILLEGAL:
	break;
      case g2_NDEV:
	break;
    }
    
    if(devp->auto_flush)
	g2_flush(dev);

    __g2_last_device=dev;
}



/**
 *
 * Draw a filled polygon.
 *
 * \param dev device
 * \param N_pt number of points (Note: It is not size of \a points vector!)
 * \param points vector of coordinates: x1, y1, x2, y2, ...
 *
 * \ingroup graphic
 */
void g2_filled_polygon(int dev, int N_pt, double *points)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_filled_polygon: No such device: %d\n", dev);
	return;
    }
    
    switch(devp->t) {
      case g2_PD:
	g2_filled_polygon_pd(devp->d.pd, N_pt, points);
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_filled_polygon(devp->d.vd->dix[i], N_pt, points);
	break;
      case g2_ILLEGAL:
	break;
      case g2_NDEV:
	break;
    }
    
    if(devp->auto_flush)
	g2_flush(dev);

    __g2_last_device=dev;
}



/**
 *
 * Draw an ellipse.
 *
 * \param dev device
 * \param x  x coordinate of the center
 * \param y  y coordinate of the center
 * \param r1 x radius
 * \param r2 y radius
 *
 * \ingroup graphic
 */
void g2_ellipse(int dev, double x, double y, double r1, double r2)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_ellipse: No such device: %d\n", dev);
	return;
    }
    
    devp->x=x;
    devp->y=y;
    
    switch(devp->t) {
      case g2_PD:
	g2_ellipse_pd(devp->d.pd, x, y, r1, r2);
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_ellipse(devp->d.vd->dix[i], x, y, r1, r2);
	break;
      case g2_ILLEGAL:
	break;
      case g2_NDEV:
	break;
    }

	
    if(devp->auto_flush)
	g2_flush(dev);

    __g2_last_device=dev;
}



/**
 *
 * Draw a filled ellipse.
 *
 * \param dev device
 * \param x  x coordinate of the center
 * \param y  y coordinate of the center
 * \param r1 x radius
 * \param r2 y radius
 *
 * \ingroup graphic
 */
void g2_filled_ellipse(int dev, double x, double y, double r1, double r2)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_filled_ellipse: No such device: %d\n", dev);
	return;
    }
    
    devp->x=x;
    devp->y=y;
    
    switch(devp->t) {
      case g2_PD:
	g2_filled_ellipse_pd(devp->d.pd, x, y, r1, r2);
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_filled_ellipse(devp->d.vd->dix[i], x, y, r1, r2);
	break;
      case g2_ILLEGAL:
	break;
      case g2_NDEV:
	break;
    }

	
    if(devp->auto_flush)
	g2_flush(dev);

    __g2_last_device=dev;
}



/**
 *
 * Draw a circle.
 *
 * \param dev device
 * \param x  x coordinate of the center
 * \param y  y coordinate of the center
 * \param r  radius
 *
 * \ingroup graphic
 */
void g2_circle(int dev, double x, double y, double r)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_circle: No such device: %d\n", dev);
	return;
    }
    
    devp->x=x;
    devp->y=y;
    
    switch(devp->t) {
      case g2_PD:
	g2_circle_pd(devp->d.pd, x, y, r);
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_circle(devp->d.vd->dix[i], x, y, r);
	break;
      case g2_ILLEGAL:
	break;
      case g2_NDEV:
	break;
    }

	
    if(devp->auto_flush)
	g2_flush(dev);

    __g2_last_device=dev;
}



/**
 *
 * Draw a filled circle.
 *
 * \param dev device
 * \param x  x coordinate of the center
 * \param y  y coordinate of the center
 * \param r  radius
 *
 * \ingroup graphic
 */
void g2_filled_circle(int dev, double x, double y, double r)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_filled_circle: No such device: %d\n", dev);
	return;
    }
    
    devp->x=x;
    devp->y=y;
    
    switch(devp->t) {
      case g2_PD:
	g2_filled_circle_pd(devp->d.pd, x, y, r);
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_filled_circle(devp->d.vd->dix[i], x, y, r);
	break;
      case g2_ILLEGAL:
	break;
      case g2_NDEV:
	break;
    }

    if(devp->auto_flush)
	g2_flush(dev);

    __g2_last_device=dev;
}



/**
 *
 * Draw an arc.
 *
 * \param dev device
 * \param x  x coordinate of the center
 * \param y  y coordinate of the center
 * \param r1 x radius
 * \param r2 y radius
 * \param a1 starting angle (in deg. 0-360)
 * \param a2 ending angle (in deg. 0-360)
 *
 * \ingroup graphic
 */
void g2_arc(int dev,
	    double x, double y,
	    double r1, double r2,
	    double a1, double a2)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_arc: No such device: %d\n", dev);
	return;
    }
    
    devp->x=x;
    devp->y=y;
    
    switch(devp->t) {
      case g2_PD:
	g2_arc_pd(devp->d.pd, x, y, r1, r2, a1, a2);
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_arc(devp->d.vd->dix[i], x, y, r1, r2, a1, a2);
	break;
      case g2_ILLEGAL:
	break;
      case g2_NDEV:
	break;
    }
	
    if(devp->auto_flush)
	g2_flush(dev);

    __g2_last_device=dev;
}



/**
 *
 * Draw a filled arc.
 *
 * \param dev device
 * \param x  x coordinate of the center
 * \param y  y coordinate of the center
 * \param r1 x radius
 * \param r2 y radius
 * \param a1 starting angle (in deg. 0-360)
 * \param a2 ending angle (in deg. 0-360)
 *
 * \ingroup graphic
 */
void g2_filled_arc(int dev,
		   double x, double y,
		   double r1, double r2,
		   double a1, double a2)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_filled_arc: No such device: %d\n", dev);
	return;
    }
    
    devp->x=x;
    devp->y=y;
    
    switch(devp->t) {
      case g2_PD:
	g2_filled_arc_pd(devp->d.pd, x, y, r1, r2, a1, a2);
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_filled_arc(devp->d.vd->dix[i], x, y, r1, r2, a1, a2);
	break;
      case g2_ILLEGAL:
	break;
      case g2_NDEV:
	break;
    }
	
    if(devp->auto_flush)
	g2_flush(dev);

    __g2_last_device=dev;
}



/**
 *
 * Draw string, see also g2_set_font_size().
 *
 * \param dev device
 * \param x  x coordinate
 * \param y  y coordinate
 * \param text null terminated string
 *
 * \ingroup graphic
 */
void g2_string(int dev, double x, double y, const char *text)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_string: No such device: %d\n", dev);
	return;
    }
    
    devp->x=x;
    devp->y=y;
    
    switch(devp->t) {
      case g2_PD:
	g2_string_pd(devp->d.pd, x, y, text);
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_string(devp->d.vd->dix[i], x, y, text);
	break;
      case g2_ILLEGAL:
	break;
      case g2_NDEV:
	break;
    }
	
    if(devp->auto_flush)
	g2_flush(dev);

    __g2_last_device=dev;
}



/**
 *
 * Draw a pen image
 *
 * \param dev device
 * \param x  x coordinate
 * \param y  y coordinate
 * \param x_size  x size
 * \param y_size  y size
 * \param pens vector of x_size*y_size pens: p11, p21, ... pxy, ...
 *
 * \ingroup graphic
 */
void g2_image(int dev, double x, double y, int x_size, int y_size, int *pens)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_image: No such device: %d\n", dev);
	return;
    }
    
    devp->x=x;
    devp->y=y;
    
    switch(devp->t) {
      case g2_PD:
	g2_image_pd(devp->d.pd, x, y, x_size, y_size, pens);
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_image(devp->d.vd->dix[i], x, y, x_size, y_size, pens);
	break;
      case g2_ILLEGAL:
	break;
      case g2_NDEV:
	break;
    }
	
    if(devp->auto_flush)
	g2_flush(dev);

    __g2_last_device=dev;
}



/**
 *
 * Quasi Pixel fake. Quasi pixel is introduced to make easier
 * plotting of cellular automata and related pictures. QP is simple a big pixel as
 * specified by g2_set_QP(). Coordinates are skaled accordingly, so no recalculation
 * is needed on client side.
 *
 * \param dev device
 * \param x  x coordinate
 * \param y  y coordinate
 *
 * \ingroup graphic
 */
void g2_plot_QP(int dev, double x, double y)
{
    g2_device *devp;
    double d;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_plot_QP: No such device: %d\n", dev);
	return;
    }
    
    x=dtoi(x);
    y=dtoi(y);
    d=devp->QPd;
    switch(devp->QPshape) {
      case QPrect:
	g2_filled_rectangle(dev, x*d-d/2, y*d-d/2, x*d+d/2, y*d+d/2);
	break;
      case QPcirc:
	g2_filled_circle(dev, x*d, y*d, d/2.0);
	break;
      default:
	fprintf(stderr, "g2: QP: unknown shape\n");
	break;
    }
    if(devp->auto_flush)
        g2_flush(dev);

    __g2_last_device=dev;
}






