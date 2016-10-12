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
#include <stdarg.h>
#include <limits.h>
#include <math.h>
#include <string.h>

#include "g2.h"
#include "g2_device.h"
#include "g2_util.h"
#include "g2_config.h"

#include "g2_PS.h"
#include "g2_PS_P.h"
#include "g2_PS_funix.h"
#include "g2_PS_definitions.h"

static int N_PS=0;
static g2_PS_device *g2_PS_dev=NULL;

/**
 * \ingroup physdev
 * \defgroup PS PostScript
 */


/*
 *
 * Attach generic PS device
 *
 */
G2L int g2_open_PS_generic(const char *file_name,
	       enum g2_PS_paper paper,
	       enum g2_PS_orientation orientation,
		   enum g2_PS_format format,
		   long width,
		   long height)
{
    g2_PS_device *psout=NULL;
    int pid=-1, i;
    int vid;
    FILE *fp;

    if((fp=fopen(file_name, "w"))==NULL) {
	fprintf(stderr, "g2_attach_PS: Error! Can not open file '%s'\n",
		file_name);
	return -1;
    }
    
    if(g2_PS_dev==NULL) {
	g2_PS_dev=g2_malloc(sizeof(g2_PS_device));
	N_PS=1;					  /* first PS device */
	psout=&g2_PS_dev[N_PS-1];
	pid=0;
    } else {
	for(i=0;i<N_PS;i++)			  /* find free place */
	    if(g2_PS_dev[i].fp==NULL) {
		psout=&g2_PS_dev[i];
		pid=i;
		break;
	    }
	if(i==N_PS) {				  /* free place not avail. */
	    N_PS++;
	    g2_PS_dev=g2_realloc(g2_PS_dev,
				 sizeof(g2_PS_device)*N_PS);
	    psout=&g2_PS_dev[N_PS-1];
	    pid=N_PS-1;
	}
    }

    vid = g2_register_physical_device(pid, NULL,
				      g2_DoubleCoor, g2_PS_funix,
				      1.0, 1.0,
				      0.0, 0.0);

    psout->fp=fp;			      /* init PostScript structures */
    psout->paper=paper;
    psout->orient=orientation;
    psout->format=format;
    psout->width=width;
    psout->height=height;
    psout->inks=NULL;
    psout->N_ink=0;
    psout->pen=0;
    psout->page_counter=0;
    psout->bbox = 0;
	
    g2_PS_write_file_header(psout);
    
    g2_PS_set_line_width(pid, NULL, 0.0);
    g2_PS_set_font_size(pid, NULL, 12.0);

					    /* g2 settings */
    g2_allocate_basic_colors(vid);
    g2_pen(vid, 1);
    
    return vid;
}


/**
 *
 * Create a PS device.
 *
 * \param file_name postscript file name
 * \param paper paper type, see ::g2_PS_paper and \ref appendix Appendix
 * \param orientation paper orientation, see ::g2_PS_orientation
 *
 * \return physical device id
 *
 * \ingroup PS
 */
G2L int g2_open_PS(const char *file_name,
	       enum g2_PS_paper paper,
	       enum g2_PS_orientation orientation)
{
    return g2_open_PS_generic(file_name,paper,orientation,g2_PS_PostScript,0,0);
}


/**
 *
 * Create an encapsulated PS device.
 *
 * \param file_name postscript file name
 *
 * \return physical device id
 *
 * \ingroup PS
 */
G2L int g2_open_EPSF(const char *file_name)
{
    return g2_open_PS_generic(file_name,0,0,g2_PS_EPSF,0,0);
}


/**
 *
 * Create an encapsulated PS device with clipping.
 *
 * \param file_name postscript file name
 * \param width clipping region width
 * \param height clipping region height
 *
 * \return physical device id
 *
 * \ingroup PS
 */
G2L int g2_open_EPSF_CLIP(const char *file_name,
			long width, long height)
{
    return g2_open_PS_generic(file_name,0,0,g2_PS_EPSF_CLIP,width,height);
}

/*
 *
 *  Write header for postscript file
 *
 */
int g2_PS_write_file_header(g2_PS_device *ps)
{
    int i;
    if (ps->format == g2_PS_PostScript)
	{
	    fprintf(ps->fp,"%%!PS-Adobe-2.0\n");
	    switch(ps->orient) 
		{
		  case g2_PS_land:
		    fprintf(ps->fp,"%%%%Orientation: Landscape\n");
		    break;
		  case g2_PS_port:
		    fprintf(ps->fp,"%%%%Orientation: Portrait\n");
		    break;
		}
	}
    else if (ps->format == g2_PS_EPSF_CLIP)
    {
	fprintf(ps->fp,"%%!PS-Adobe-3.0 EPSF-2.0\n");
	fprintf(ps->fp,"%%%%BoundingBox: 0 0 %ld %ld\n",ps->width,ps->height);
    }
    else if (ps->format == g2_PS_EPSF)
    {
	fprintf(ps->fp,"%%!PS-Adobe-3.0 EPSF-2.0\n");
	fprintf(ps->fp,"%%%%BoundingBox: (atend)\n");
    }
	
    fprintf(ps->fp,"%%%%Creator: g2 %s\n", G2_VERSION);
    fprintf(ps->fp, "%%%%EndComments\n");
    
    if (ps->format == g2_PS_EPSF_CLIP)
    {
	fprintf(ps->fp,"0 0 moveto\n");
	fprintf(ps->fp,"0 %ld rlineto\n",ps->height);
	fprintf(ps->fp,"%ld 0 rlineto\n",ps->width);
	fprintf(ps->fp,"0 %ld rlineto\n",-ps->height);
	fprintf(ps->fp,"closepath\n");
	fprintf(ps->fp,"clip\n");
    }
    
    for(i=0;g2_PS_operators[i]!=NULL;i++)
	fputs(g2_PS_operators[i], ps->fp);
    
    fprintf(ps->fp,"newpath\n");
    if((ps->orient==g2_PS_land) && (ps->format == g2_PS_PostScript))
	fprintf(ps->fp,"%d 0 translate 90 rotate\n",
		g2_PS_paper_size[ps->paper][0]);

    fputs("%%PageTrailer\n%%Page: 1 1\n", ps->fp);
    
    return 0;
}

/*
 *
 *    Add circle at (x,y) of dimension size to bounding box
 *
 */
void g2_PS_bbox_add(g2_PS_device *ps,double x,double y,double size)
	{
	if (ps->bbox == 0) /* bbox is empty */
		{
		ps->x1=x-size;
		ps->x2=x+size;
		ps->y1=y-size;
		ps->y2=y+size;
		ps->bbox=1;
		return;
		}

	if (ps->x1 > x-size)
		ps->x1=x-size;
	else if (ps->x2 < x+size)
		ps->x2=x+size;

	if (ps->y1 > y-size)
		ps->y1=y-size;
	else if (ps->y2 < y+size)
		ps->y2=y+size;
	return;
	}

/*
 *
 *    Delete PS device
 *
 */
int g2_PS_delete(int pid, void *pdp)
{
    g2_PS_device *ps=&g2_PS_dev[pid];
	fprintf(ps->fp,"\nshowpage\n");
    fprintf(ps->fp,"%%%%PageTrailer\n");
    fprintf(ps->fp,"%%%%EndPage\n");
	fprintf(ps->fp,"%%%%Trailer\n");
	if (ps->format == g2_PS_EPSF)
		{
		fprintf(ps->fp,"%%%%BoundingBox: %d %d %d %d\n",
				(int)floor(ps->x1),(int)floor(ps->y1),
				(int)ceil(ps->x2),(int)ceil(ps->y2));
		}
    fprintf(ps->fp,"%%%%EOF\n");
	fclose(ps->fp);
    free(ps->inks);
    
    ps->fp=NULL;				  /* free place */
    
    return 0;
}



int g2_PS_ink(int pid, void *pdp,
	      double red, double green, double blue)
{
    g2_PS_device *ps=&g2_PS_dev[pid];
    ps->N_ink++;
    if(ps->inks==NULL)
	ps->inks=(g2_PS_inks *)g2_malloc(ps->N_ink*
					  sizeof(g2_PS_inks));
    else
	ps->inks=(g2_PS_inks *)g2_realloc((void *)ps->inks,
					   ps->N_ink*
					   sizeof(g2_PS_inks));
    
    ps->inks[ps->N_ink-1].r=red;
    ps->inks[ps->N_ink-1].g=green;
    ps->inks[ps->N_ink-1].b=blue;
 
    return ps->N_ink-1;
}



int g2_PS_pen(int pid, void *pdp, int color)
{
    g2_PS_device *ps=&g2_PS_dev[pid];
    if(color>=ps->N_ink || color<0)
	return -1;
    fprintf(ps->fp,"%.4g %.4g %.4g setrgbcolor\n",
	    (double)ps->inks[color].r,
	    (double)ps->inks[color].g,
	    (double)ps->inks[color].b);
    ps->pen=color;
    return 0;
}



int g2_PS_set_background(int pid, void *pdp, int color)
{
    return 0;
}



int g2_PS_clear_palette(int pid, void *pdp)
{
    g2_PS_device *ps=&g2_PS_dev[pid];
    free((void *)ps->inks);
    ps->N_ink=0;
    ps->inks=NULL;
    return 0;
}



int g2_PS_set_line_width(int pid, void *pdp, double w)
{
    g2_PS_device *ps=&g2_PS_dev[pid];
    fprintf(ps->fp,"%.4g setlinewidth\n", w);
	ps->w = w;
    return 0;
}



int g2_PS_set_dash(int pid, void *pdp, int N, double *data)
{
    g2_PS_device *ps=&g2_PS_dev[pid];
    int i;
    fprintf(ps->fp, "[ ");
    for(i=0;i<N;i++)
	if(data[i]>0.0)
	    fprintf(ps->fp, "%.4g ", data[i]);
    fprintf(ps->fp, "] 0 setdash\n");
    return 0;
}



int g2_PS_set_font_size(int pid, void *pdp, double size)
{
    g2_PS_device *ps=&g2_PS_dev[pid];
    if(size<=0.0)
	return -1;
    fprintf(ps->fp,"%s findfont %.4g scalefont setfont\n", g2_PSFont, size);
    ps->size = size;
    return 0;
}



int g2_PS_clear(int pid, void *pdp)
{
    g2_PS_device *ps=&g2_PS_dev[pid];
    fprintf(ps->fp, "gsave showpage\n");
    fprintf(ps->fp, "%%%%PageTrailer\n%%%%EndPage\n");
    fprintf(ps->fp, "%%%%Page: %d %d\ngrestore newpath\n",
	    ps->page_counter+1, ps->page_counter+1);
    ps->page_counter++;
    return 0;
}



int g2_PS_flush(int pid, void *pdp)
{
    g2_PS_device *ps=&g2_PS_dev[pid];
    fflush(ps->fp);
    return 0;
}



int g2_PS_plot(int pid, void *pdp, double x, double y)
{
    g2_PS_device *ps=&g2_PS_dev[pid];
    fprintf(ps->fp,"%.4g %.4g P\n", x, y);
	g2_PS_bbox_add(ps,x,y,1);
    return 0;
}



int g2_PS_line(int pid, void *pdp, double x1, double y1, double x2, double y2)
{
    g2_PS_device *ps=&g2_PS_dev[pid];
    fprintf(ps->fp,"%.4g %.4g M %.4g %.4g L St\n",
	    x1, y1, x2, y2);
	g2_PS_bbox_add(ps,x1,y1,ps->w);
	g2_PS_bbox_add(ps,x2,y2,ps->w);
    return 0;
}



int g2_PS_poly_line(int pid, void *pdp, int N, double *points)
{
    g2_PS_device *ps=&g2_PS_dev[pid];
    int i;
    fprintf(ps->fp,"%.4g %.4g M\n", points[0], points[1]); 
	g2_PS_bbox_add(ps,points[0], points[1],ps->w);
    for(i=2;i<2*N;i+=2)
		{
		fprintf(ps->fp, "%.4g %.4g L\n", points[i], points[i+1]);
		g2_PS_bbox_add(ps,points[i], points[i+1],ps->w);
		}
    fprintf(ps->fp, "St\n");
    return 0;
}



int g2_PS_polygon(int pid, void *pdp, int N, double *points)
{
    g2_PS_device *ps=&g2_PS_dev[pid];
    int i;
    fprintf(ps->fp,"%.4g %.4g M\n",points[0], points[1]); 
	g2_PS_bbox_add(ps,points[0], points[1],ps->w);
    for(i=2;i<2*N;i+=2)
		{
		fprintf(ps->fp, "%.4g %.4g L\n", points[i], points[i+1]);
		g2_PS_bbox_add(ps,points[i], points[i+1],ps->w);
		}
    fprintf(ps->fp, "%.4g %.4g L St\n", points[0], points[1]);
    return 0;
}



int g2_PS_filled_polygon(int pid, void *pdp, int N, double *points)
{
    g2_PS_device *ps=&g2_PS_dev[pid];
    int i;
    fprintf(ps->fp,"newpath %.4g %.4g M\n",points[0], points[1]); 
	g2_PS_bbox_add(ps,points[0], points[1],ps->w);
    for(i=2;i<2*N;i+=2)
		{
		fprintf(ps->fp, "%.4g %.4g L\n", points[i], points[i+1]);
		g2_PS_bbox_add(ps,points[i], points[i+1],ps->w);
		}
    fprintf(ps->fp, "%.4g %.4g L fill St\n", points[0], points[1]);
    return 0;
}



int g2_PS_rectangle(int pid, void *pdp,
		    double x1, double y1, double x2, double y2)
{
    g2_PS_device *ps=&g2_PS_dev[pid];
    fprintf(ps->fp,"%.4g %.4g %.4g %.4g R\n",
	    x2, y2, x1, y1);
	g2_PS_bbox_add(ps,x1,y1,ps->w);
	g2_PS_bbox_add(ps,x2,y2,ps->w);
    return 0;
}



int g2_PS_filled_rectangle(int pid, void *pdp,
			   double x1, double y1, double x2, double y2)
{
    g2_PS_device *ps=&g2_PS_dev[pid];
    fprintf(ps->fp,"%.4g %.4g %.4g %.4g FR\n",
	    x2, y2, x1, y1);
	g2_PS_bbox_add(ps,x1,y1,ps->w);
	g2_PS_bbox_add(ps,x2,y2,ps->w);
    return 0;
}



int g2_PS_triangle(int pid, void *pdp,
		   double x1, double y1,
		   double x2, double y2,
		   double x3, double y3)
{
    g2_PS_device *ps=&g2_PS_dev[pid];
    fprintf(ps->fp,"%.4g %.4g %.4g %.4g %.4g %.4g T\n",
	    x1,y1,x2,y2,x3,y3);
	g2_PS_bbox_add(ps,x1,y1,ps->w);
	g2_PS_bbox_add(ps,x2,y2,ps->w);
	g2_PS_bbox_add(ps,x3,y3,ps->w);
    return 0;
}


int g2_PS_filled_triangle(int pid, void *pdp,
			 double x1, double y1,
			 double x2, double y2,
			 double x3, double y3)
{
    g2_PS_device *ps=&g2_PS_dev[pid];
    fprintf(ps->fp,"%.4g %.4g %.4g %.4g %.4g %.4g FT\n",
	    x1,y1,x2,y2,x3,y3);
	g2_PS_bbox_add(ps,x1,y1,ps->w);
	g2_PS_bbox_add(ps,x2,y2,ps->w);
	g2_PS_bbox_add(ps,x3,y3,ps->w);
    return 0;
}



int g2_PS_arc(int pid, void *pdp,
	      double x, double y,
	      double r1, double r2,
	      double a1, double a2)
{
    g2_PS_device *ps=&g2_PS_dev[pid];
    if(a1==a2) {
	a1=0;
	a2=360;
    }
    fprintf(ps->fp,"%.4g %.4g %.4g %.4g %.4g %.4g A\n",
	    a1, a2, r1, r2, x, y);
	g2_PS_bbox_add(ps,x+r1,y+r2,ps->w);
	g2_PS_bbox_add(ps,x-r1,y-r2,ps->w);
    return 0;
}



int g2_PS_filled_arc(int pid, void *pdp,
		     double x, double y,
		     double r1, double r2,
		     double a1, double a2)
{
    g2_PS_device *ps=&g2_PS_dev[pid];
    if(a1==a2) {
	a1=0;
	a2=360;
    }
    fprintf(ps->fp,"%.4g %.4g %.4g %.4g %.4g %.4g FA\n",
	    a1, a2, r1, r2, x, y);
	g2_PS_bbox_add(ps,x+r1,y+r2,ps->w);
	g2_PS_bbox_add(ps,x-r1,y-r2,ps->w);
    return 0;
}



int g2_PS_ellipse(int pid, void *pdp,
		  double x, double y,
		  double r1, double r2)
{
    g2_PS_device *ps=&g2_PS_dev[pid];
    fprintf(ps->fp,"0 360 %.4g %.4g %.4g %.4g A\n",
	    r1, r2, x, y);
	g2_PS_bbox_add(ps,x+r1,y+r2,ps->w);
	g2_PS_bbox_add(ps,x-r1,y-r2,ps->w);
    return 0;
}
 

 

int g2_PS_filled_ellipse(int pid, void *pdp,
			 double x, double y,
			 double r1, double r2)
{
    g2_PS_device *ps=&g2_PS_dev[pid];
    fprintf(ps->fp,"0 360 %.4g %.4g %.4g %.4g FA\n",
	    r1, r2, x, y);
	g2_PS_bbox_add(ps,x+r1,y+r2,ps->w);
	g2_PS_bbox_add(ps,x-r1,y-r2,ps->w);
    return 0;
}

  

int g2_PS_draw_string(int pid, void *pdp,
		      double x, double y, const char *text)
{
    g2_PS_device *ps=&g2_PS_dev[pid];
    fputc('(', ps->fp);
    for(;*text!='\0';text++)
	switch(*text) {
	  case '(':
	    fputs("\\(", ps->fp);
	    break;
	  case ')':
	    fputs("\\)", ps->fp);
	    break;
	  case '\\':
	    fputs("\\\\", ps->fp);
	    break;
	  default:
	    fputc(*text, ps->fp);
	    break;
	}
    fprintf(ps->fp,") %.4g %.4g S\n", x, y);
	g2_PS_bbox_add(ps,x,y,ps->size);
	g2_PS_bbox_add(ps,x+ps->size*strlen(text),y,ps->size);
    return 0;
}

