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
#ifndef _G2_FIG_P_H
#define _G2_FIG_P_H

#include "g2_FIG.h"

#include <stdio.h>

typedef struct _g2_FIG_inks {
    unsigned char red;
    unsigned char green;
    unsigned char blue;
} g2_FIG_inks;

typedef struct _g2_FIG_device {
    FILE           *fp;	          /* output file pointer */
    int            pen_color;	  /* current pen */
    int            thickness;     /* line thickness */
    int            font_size;     /* font size */
    int            line_style;    /* line style according to fig format */
    int            style_val;     /* line style value acc. to fig format */

    g2_FIG_inks    inks[512];      /* 512 user defined colors */
    int            N_inks;         /* number of allocated colors */
    fpos_t         color_file_pos; /* file position of colors, ftell */
} g2_FIG_device;



int g2_FIG_write_file_header(g2_FIG_device *figd);
int g2_FIG_delete(int pid, void *pdp);
int g2_FIG_ink(int pid, void *pdp,
	      double red, double green, double blue);
int g2_FIG_pen(int pid, void *pdp, int color);
int g2_FIG_set_background(int pid, void *pdp, int color);
int g2_FIG_reset_palette(int pid, void *pdp);
int g2_FIG_clear_palette(int pid, void *pdp);
int g2_FIG_set_line_width(int pid, void *pdp, int w);
int g2_FIG_set_dash(int pid, void *pdp, int N, int *data);
int g2_FIG_set_font_size(int pid, void *pdp, int size);
int g2_FIG_clear(int pid, void *pdp);
int g2_FIG_flush(int pid, void *pdp);
int g2_FIG_plot(int pid, void *pdp, int x, int y);
int g2_FIG_line(int pid, void *pdp, int x1, int y1, int x2, int y2);
int g2_FIG_poly_line(int pid, void *pdp, int N, int *points);
int g2_FIG_polygon(int pid, void *pdp, int N, int *points);
int g2_FIG_filled_polygon(int pid, void *pdp, int N, int *points);
int g2_FIG_arc(int pid, void *pdp,
	      int x, int y,
	      int r1, int r2,
	      double a1, double a2);
int g2_FIG_filled_arc(int pid, void *pdp,
		     int x, int y,
		     int r1, int r2,
		     double a1, double a2);
int g2_FIG_ellipse(int pid, void *pdp,
		  int x, int y,
		  int r1, int r2);
int g2_FIG_filled_ellipse(int pid, void *pdp,
			 int x, int y,
			 int r1, int r2);
int g2_FIG_draw_string(int pid, void *pdp,
		      int x, int y, const char *text);

#endif /* _G2_FIG_P_H */
