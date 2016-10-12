/*****************************************************************************
**  This is part of the g2 library
**  Copyright (C) 1998  Ljubomir Milanovic & Horst Wagner
**
**  This program is free software; you can redistribute it and/or modify
**  it under the terms of the GNU General Public License (version 2) as
**  published by the Free Software Foundation.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program; if not, write to the Free Software
**  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
******************************************************************************/
#ifndef _G2_PS_P_H
#define _G2_PS_P_H

#include "g2_PS.h"


typedef struct _g2_PS_inks {
    double r;					/* red   [0:1] */
    double g;					/* green [0:1] */
    double b;					/* blue  [0:1] */
} g2_PS_inks;



typedef struct _g2_PS_device {
  FILE                    *fp;		/* output file pointer */
  enum g2_PS_paper        paper;	/* paper type */
  enum g2_PS_orientation  orient;	/* page orientation */
  enum g2_PS_format	  format;	/* PS or EPSF format */
  long                    width,height;	/* width and height for EPSF */
  double                  x1,y1,x2,y2;	/* min. Bounding Box */
  int                     bbox;		/* Bounding Box empty flag */
  double                  w,size;	/* line width/font size (required for Bbox) */
  
  g2_PS_inks     *inks;			/* allocated colors */
  int            N_ink;			/* number of allocated colors*/
  int            pen;			/* current pen */
  int            page_counter;		/* page counter ( Clear() ) */ 
} g2_PS_device;



void g2_PS_bbox_add(g2_PS_device *ps,double x,double y,double size);
int g2_PS_write_file_header(g2_PS_device *psd);
int g2_PS_delete(int pid, void *pdp);
int g2_PS_ink(int pid, void *pdp,
	      double red, double green, double blue);
int g2_PS_pen(int pid, void *pdp, int color);
int g2_PS_set_background(int pid, void *pdp, int color);
int g2_PS_reset_palette(int pid, void *pdp);
int g2_PS_clear_palette(int pid, void *pdp);
int g2_PS_set_line_width(int pid, void *pdp, double w);
int g2_PS_set_dash(int pid, void *pdp, int N, double *data);
int g2_PS_set_font_size(int pid, void *pdp, double size);
int g2_PS_clear(int pid, void *pdp);
int g2_PS_flush(int pid, void *pdp);
int g2_PS_plot(int pid, void *pdp, double x, double y);
int g2_PS_line(int pid, void *pdp, double x1, double y1, double x2, double y2);
int g2_PS_poly_line(int pid, void *pdp, int N, double *points);
int g2_PS_polygon(int pid, void *pdp, int N, double *points);
int g2_PS_filled_polygon(int pid, void *pdp, int N, double *points);
int g2_PS_rectangle(int pid, void *pdp,
		    double x1, double y1, double x2, double y2);
int g2_PS_filled_rectangle(int pid, void *pdp,
			   double x1, double y1, double x2, double y2);
int g2_PS_triangle(int pid, void *pdp,
		   double x1, double y1,
		   double x2, double y2,
		   double x3, double y3);
int g2_PS_filled_triangle(int pid, void *pdp,
			 double x1, double y1,
			 double x2, double y2,
			 double x3, double y3);
int g2_PS_arc(int pid, void *pdp,
	      double x, double y,
	      double r1, double r2,
	      double a1, double a2);
int g2_PS_filled_arc(int pid, void *pdp,
		     double x, double y,
		     double r1, double r2,
		     double a1, double a2);
int g2_PS_ellipse(int pid, void *pdp,
		  double x, double y,
		  double r1, double r2);
int g2_PS_filled_ellipse(int pid, void *pdp,
			 double x, double y,
			 double r1, double r2);
int g2_PS_draw_string(int pid, void *pdp,
		      double x, double y, const char *text);

#endif /* _G2_PS_P_H */
