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
#ifndef _G2_GRAPHIC_PD_H
#define _G2_GRAPHIC_PD_H

#include "g2_physical_device.h"
#include "g2_funix.h"


void g2_plot_pd(g2_physical_device *pd, double x, double y);
void g2_line_pd(g2_physical_device *pd,
		double x1, double y1, double x2, double y2);
void g2_poly_line_pd(g2_physical_device *pd, int N, double *points);
void g2_triangle_pd(g2_physical_device *pd,
		    double x1, double y1,
		    double x2, double y2,
		    double x3, double y3);
void g2_filled_triangle_pd(g2_physical_device *pd,
			   double x1, double y1,
			   double x2, double y2,
			   double x3, double y3);
void g2_rectangle_pd(g2_physical_device *pd,
		     double x1, double y1, double x2, double y2);
void g2_filled_rectangle_pd(g2_physical_device *pd,
			    double x1, double y1, double x2, double y2);
void g2_polygon_pd(g2_physical_device *pd, int N, double *points);
void g2_filled_polygon_pd(g2_physical_device *pd, int N, double *points);
void g2_ellipse_pd(g2_physical_device *pd,
		   double x, double y, double r1, double r2);
void g2_filled_ellipse_pd(g2_physical_device *pd,
			  double x, double y, double r1, double r2);
void g2_circle_pd(g2_physical_device *pd,
		  double x, double y, double r);
void g2_filled_circle_pd(g2_physical_device *pd,
			 double x, double y, double r);
void g2_arc_pd(g2_physical_device *pd,
	       double x, double y,
	       double r1, double r2,
	       double a1, double a2);
void g2_filled_arc_pd(g2_physical_device *pd,
		      double x, double y,
		      double r1, double r2,
		      double a1, double a2);
void g2_string_pd(g2_physical_device *pd,
		  double x, double y, const char *text);
void g2_image_pd(g2_physical_device *pd,
	         double x, double y, int x_size, int y_size, int *pens);


#endif /* _G2_GRAPHIC_PD_H */
