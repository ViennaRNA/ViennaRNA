/*****************************************************************************
**  Copyright (C) 1998-2004  Ljubomir Milanovic & Horst Wagner
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
#include <stdarg.h>
#include <limits.h>
#include <math.h>
#include <string.h>

#include "g2.h"
#include "g2_device.h"
#include "g2_util.h"
#include "g2_config.h"

#include "g2_FIG.h"
#include "g2_FIG_P.h"
#include "g2_FIG_funix.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif /* PI */

static int N_FIG=0;
static g2_FIG_device *g2_FIG_dev=NULL;

/**
 * \ingroup physdev
 * \defgroup FIG FIG
 *
 * FIG devices generate output in the FIG 3.2 format. For more details
 * about FIG format and xfig application please visit http://www.xfig.org .
 *
 * \note FIG is a vector-oriented (as oposed to pixel-oriented) format.
 *       Therefore ::g2_image function and splines are not optimally supported. 
 */


/**
 *
 * Create a FIG device. g2 uses A4 paper size (landscape orientation) as default.
 *
 * \param file_name fig file name
 *
 * \return physical device id
 *
 * \ingroup FIG
 */
G2L int g2_open_FIG(const char *file_name)
{
    g2_FIG_device *figout=NULL;
    int pid=-1, i;
    int vid;
    FILE *fp;

    if((fp=fopen(file_name, "w"))==NULL) {
	g2_log(Error, "g2_attach_PS: Error! Can not open file '%s'\n",
	       file_name);
	return -1;
    }
    
    if(g2_FIG_dev==NULL) {
	g2_FIG_dev=g2_malloc(sizeof(g2_FIG_device));
	N_FIG=1;				  /* first FIG device */
	figout=&g2_FIG_dev[N_FIG-1];
	pid=0;
    } else {
	for(i=0;i<N_FIG;i++) {			  /* find free place */
	    if(g2_FIG_dev[i].fp==NULL) {
		figout=&g2_FIG_dev[i];
		pid=i;
		break;
	    }
	}
	if(i==N_FIG) {				  /* free place not avail. */
	    N_FIG++;
	    g2_FIG_dev=g2_realloc(g2_FIG_dev, sizeof(g2_FIG_device)*N_FIG);
	    figout=&g2_FIG_dev[N_FIG-1];
	    pid=N_FIG-1;
	}
    }

    vid = g2_register_physical_device(pid, NULL,
				      g2_IntCoor, g2_FIG_funix,
				      1200./80., -1200./80.,
				      0.0, 595/72*1200-160);
                                   /* -180 is a hand-made correcture factor */

    /* init FIG structures */
    figout->fp=fp;
    figout->pen_color=1;
    figout->thickness=-1;
    figout->font_size=-1;
    figout->N_inks=0;

    /* write file header (incl. color placeholder) */
    g2_FIG_write_file_header(figout);
    
    /* g2 settings callbacks */
    g2_allocate_basic_colors(vid);
    g2_pen(vid, 1);
    g2_set_line_width(vid, 1);
    g2_set_dash(vid, 0, NULL);
    g2_set_font_size(vid, 12);

    return vid;
}



/*
 *
 *  Write header for fig file
 *
 */
int g2_FIG_write_file_header(g2_FIG_device *fig)
{
    int i;
    
    fprintf(fig->fp, "#FIG 3.2\n");
    fprintf(fig->fp, "#Creator: g2 %s\n", G2_VERSION);
    
    fprintf(fig->fp, "Landscape\n");
    fprintf(fig->fp, "Flush Left\n");
    fprintf(fig->fp, "Metric\n");
    fprintf(fig->fp, "A4\n");
    fprintf(fig->fp, "100\n");
    fprintf(fig->fp, "Single\n");
    fprintf(fig->fp, "-2\n");
    fprintf(fig->fp, "1200 2\n");

    /* write placeholder for all colors */
    fgetpos(fig->fp, &(fig->color_file_pos));
    for(i=0;i<512;i++)
	fprintf(fig->fp, "0 %3d #%02x%02x%02x\n",32+i, 0xff, 0xff, 0xff);

    return 0;
}


/*
 *
 *    Delete FIG device
 *
 */
int g2_FIG_delete(int pid, void *pdp)
{
    int i;
    g2_FIG_device *fig=&g2_FIG_dev[pid];

    /* now write all defined  colors to the reserved place */
    fsetpos(fig->fp, &(fig->color_file_pos));
    for(i=0;i<fig->N_inks;i++) {
	fprintf(fig->fp, "0 %3d #%02x%02x%02x\n",
		32+i,
		(int)fig->inks[i].red,
		(int)fig->inks[i].green,
		(int)fig->inks[i].blue);
    }
    
    fclose(fig->fp);
    fig->fp=NULL;				  /* free place */
    return 0;
}



int g2_FIG_ink(int pid, void *pdp,
	      double red, double green, double blue)
{
    /* the fig pen space is betwen 32 and 543, total of 512 colors,
       in g2 implementation from 0000 to 0777.
       On the other hand ink space is between #000000 and #ffffff. */

    g2_FIG_device *fig=&g2_FIG_dev[pid];

    if(fig->N_inks>=512)
	return -1;

    fig->inks[fig->N_inks].red   = red*0xff;
    fig->inks[fig->N_inks].green = green*0xff;
    fig->inks[fig->N_inks].blue  = blue*0xff;

    return fig->N_inks++;
}



int g2_FIG_pen(int pid, void *pdp, int color)
{
    g2_FIG_device *fig=&g2_FIG_dev[pid];

    if(color>=fig->N_inks) {
	return -1;
    }

    fig->pen_color = 32+color;

    return 0;
}



int g2_FIG_set_background(int pid, void *pdp, int color)
{
    return -1;
}



int g2_FIG_clear_palette(int pid, void *pdp)
{
    return -1;
}



int g2_FIG_set_line_width(int pid, void *pdp, int w)
{
    g2_FIG_device *fig=&g2_FIG_dev[pid];
    fig->thickness = w*80./1200.;
    return 0;
}



int g2_FIG_set_dash(int pid, void *pdp, int N, int *data)
{
    g2_FIG_device *fig=&g2_FIG_dev[pid];
    int black, white;
    
    if(N==0 || data==NULL) {
	fig->line_style=0;
	fig->style_val=-1;
	return 0;
    }
    if(N<2) {
	return -1;
    }
    black = data[0]*80./1200.;  /* FIG format has no sofistificated dash concept */
    white = data[1]*80./1200.;  /* we will do out best */
    if(black<4) {
	fig->line_style = 2;
    } else {
	fig->line_style = 1;
    }
    fig->style_val = white;
    return 0;
}



int g2_FIG_set_font_size(int pid, void *pdp, int size)
{
    g2_FIG_device *fig=&g2_FIG_dev[pid];
    fig->font_size = size*72./1200.;
    return 0;
}



int g2_FIG_clear(int pid, void *pdp)
{
    return -1;
}



int g2_FIG_flush(int pid, void *pdp)
{
    g2_FIG_device *fig=&g2_FIG_dev[pid];
    fflush(fig->fp);
    return 0;
}



int g2_FIG_plot(int pid, void *pdp, int x, int y)
{
    g2_FIG_device *fig=&g2_FIG_dev[pid];
    fprintf(fig->fp, "2 1 0 %d %d -1 1 1 -1 -1 1 0 -1 0 0 2\n",
	    1, fig->pen_color);
    fprintf(fig->fp, " %d %d\n %d %d\n", x, y, x+1, y);
    return 0;
}



int g2_FIG_line(int pid, void *pdp, int x1, int y1, int x2, int y2)
{
    g2_FIG_device *fig=&g2_FIG_dev[pid];
    fprintf(fig->fp, "2 1 %d %d %d -1 1 1 -1 %d 1 0 -1 0 0 2\n",
	    fig->line_style, fig->thickness, fig->pen_color, fig->style_val);
    fprintf(fig->fp, " %d %d\n %d %d\n", x1, y1, x2, y2);
    return 0;
}



int g2_FIG_poly_line(int pid, void *pdp, int N, int *points)
{
    g2_FIG_device *fig=&g2_FIG_dev[pid];
    int i;

    fprintf(fig->fp,"2 1 %d %d %d -1 1 1 -1 %d 1 0 -1 0 0 %d\n",
	    fig->line_style, fig->thickness, fig->pen_color, fig->style_val, N);
    for(i=0;i<2*N;i+=2) {
	fprintf(fig->fp, " %d %d\n", points[i], points[i+1]);
    }
    return 0;
}



int g2_FIG_polygon(int pid, void *pdp, int N, int *points)
{
    g2_FIG_device *fig=&g2_FIG_dev[pid];
    int i;

    if(N<2) {
	return -1;
    }
    
    fprintf(fig->fp,"2 3 %d %d %d -1 1 1 -1 %d 1 0 -1 0 0 %d\n",
	    fig->line_style, fig->thickness, fig->pen_color, fig->style_val, N+1);
    for(i=0;i<2*N;i+=2) {
	fprintf(fig->fp, " %d %d\n", points[i], points[i+1]);
    }
    fprintf(fig->fp, " %d %d\n", points[0], points[1]);
    return 0;
}



int g2_FIG_filled_polygon(int pid, void *pdp, int N, int *points)
{
    g2_FIG_device *fig=&g2_FIG_dev[pid];
    int i;

    if(N<2) {
	return -1;
    }
    
    fprintf(fig->fp,"2 3 %d %d %d %d 1 1 20 %d 1 0 -1 0 0 %d\n",
	    fig->line_style, fig->thickness, fig->pen_color, fig->pen_color,
	    fig->style_val, N+1);
    for(i=0;i<2*N;i+=2) {
	fprintf(fig->fp, " %d %d\n", points[i], points[i+1]);
    }
    fprintf(fig->fp, " %d %d\n", points[0], points[1]);
    return 0;
}


int g2_FIG_arc(int pid, void *pdp,
	       int x, int y,
	       int r1, int r2,
	       double a1, double a2)
{
    g2_FIG_device *fig=&g2_FIG_dev[pid];
    double a0, d;
    double a0_rad, da_rad;
    int N, i;

    a0=fmod(a1, 360.) + (a1<0? 360:0);   /* map a1 to [0, 360) */
    d=a2>a1? a2-a1:a2-a1+360;

    N=3+d/18;

    a0_rad = a0*2.*PI/360.;
    da_rad =  d*2.*PI/360./(N-1);

    fprintf(fig->fp, "3 2 %d %d %d -1 1 -1-1 %d 0 0 0 %d\n",
	    fig->line_style, fig->thickness, fig->pen_color, fig->style_val, N);
    for(i=0;i<N;i++) {
	fprintf(fig->fp, " %d %d\n",
		(int)(x+r1*cos(a0_rad+i*da_rad)),
		(int)(y-r2*sin(a0_rad+i*da_rad)));
    }
    fprintf(fig->fp, " -1");
    for(i=1;i<N-1;i++) {
	fprintf(fig->fp, " -1");
    }
    fprintf(fig->fp, " -1\n");
    
    return 0;
}


int g2_FIG_filled_arc(int pid, void *pdp,
		      int x, int y,
		      int r1, int r2,
		      double a1, double a2)
{
    g2_FIG_device *fig=&g2_FIG_dev[pid];
    double a0, d;
    double a0_rad, da_rad;
    int N, i;

    a0=fmod(a1, 360.) + (a1<0? 360:0);   /* map a1 to [0, 360) */
    d=a2>a1? a2-a1:a2-a1+360;

    N=3+d/18;

    a0_rad = a0*2.*PI/360.;
    da_rad =  d*2.*PI/360./(N-1);

    fprintf(fig->fp, "3 2 %d %d %d %d 1 -1 20 %d 0 0 0 %d\n",
	    fig->line_style, fig->thickness, fig->pen_color, fig->pen_color,
	    fig->style_val, N+2);
    for(i=0;i<N;i++) {
	fprintf(fig->fp, " %d %d\n",
		(int)(x+r1*cos(a0_rad+i*da_rad)), (int)(y-r2*sin(a0_rad+i*da_rad)));
    }
    fprintf(fig->fp, " %d %d\n", x, y); 
    fprintf(fig->fp, " %d %d\n",
	    (int)(x+r1*cos(a0_rad+0*da_rad)),
	    (int)(y-r2*sin(a0_rad+0*da_rad))); 
    fprintf(fig->fp, " 0");
    for(i=1;i<N-1;i++) {
	fprintf(fig->fp, " -1");
    }
    fprintf(fig->fp, " 0 0 0\n");
    
    return 0;
}


int g2_FIG_ellipse(int pid, void *pdp,
		  int x, int y,
		  int r1, int r2)
{
    g2_FIG_device *fig=&g2_FIG_dev[pid];
    fprintf(fig->fp,"1 1 %d %d %d -1  1 -1 -1 %d 1 0.0 %d %d %d %d %d %d %d %d\n",
	    fig->line_style, fig->thickness, fig->pen_color, fig->style_val,
	    x, y, r1, r2, x, y, x+r1, y+r2);

    return 0;
}


int g2_FIG_filled_ellipse(int pid, void *pdp,
			 int x, int y,
			 int r1, int r2)
{
    g2_FIG_device *fig=&g2_FIG_dev[pid];
    fprintf(fig->fp,"1 1 %d %d %d %d 1 -1 20 %d 1 0.0 %d %d %d %d %d %d %d %d\n",
	    fig->line_style, fig->thickness, fig->pen_color, fig->pen_color, fig->style_val,
	    x, y, r1, r2, x, y, x+r1, y+r2);

    return 0;
}
  

int g2_FIG_draw_string(int pid, void *pdp,
		      int x, int y, const char *text)
{
    const char *c;
    g2_FIG_device *fig=&g2_FIG_dev[pid];
    if(fig->font_size<=0)
	return 0;
    fprintf(fig->fp,"4 0 %d 1 -1 0 %d 0 5 %d %d %d %d ",
	    fig->pen_color, fig->font_size, fig->font_size, fig->font_size*strlen(text), x, y);
    for(c=text;*c;c++) {
	if(*c=='\\')
	    fputc('\\', fig->fp); /* escape \\ */
	fputc(*c, fig->fp);
    }
    fprintf(fig->fp,"\\001\n");
    return 0;
}
