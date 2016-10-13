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
#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "g2.h"
#include "g2_device.h"
#include "g2_util.h"
#include "g2_X11_P.h"
#include "g2_X11.h"
#include "g2_X11_funix.h"
#include "g2_config.h"


static int N_X11=0;
static g2_X11_device *g2_X11_dev=NULL;

/**
 * \ingroup physdev
 * \defgroup X11 X11
 */


/**
 *
 * Open a simple X11 window (physical device device).
 *
 * \param width window width
 * \param height window height
 * \return physical device id
 *
 * \ingroup X11
 */
int g2_open_X11(int width, int height)
{
    return g2_open_X11X(width, height,
			10, 10,
			NULL, NULL,
			NULL, -1, -1);
}



/**
 *
 * Open a X11 window (physical device device). If \a icon_width or \a
 * icon_height is smaller than 0, the \a icon_data is interpreted as a
 * file name.
 *
 * \param width window width
 * \param height window height
 * \param x x position on screen
 * \param y y position on screen
 * \param window_name hint for window manager
 * \param icon_name hint for window manager
 * \param icon_data icon bitmap (\a icon_width * \a icon_height bits) or file name containing bitmap
 *                  (if \a icon_width <= 0 or \a icon_height <= 0) 
 * \param icon_width icon width
 * \param icon_height icon height
 * \return physical device id
 *
 * \ingroup X11
 */
int g2_open_X11X(int width, int height,
		 int x, int y,
		 char *window_name, char *icon_name,
		 char *icon_data, int icon_width, int icon_height)
{
    g2_X11_device *xout=NULL;
    int pid=-1, i;
    char name[32];
    int vid;
    
    if(g2_X11_dev==NULL) {
	g2_X11_dev=g2_malloc(sizeof(g2_X11_device));
	N_X11=1;				  /* first X11 device */
	xout=&g2_X11_dev[N_X11-1];
	pid=0;
    } else {
	for(i=0;i<N_X11;i++)			  /* find free place */
	    if(g2_X11_dev[i].display==NULL) {
		xout=&g2_X11_dev[i];
		pid=i;
		break;
	    }
	if(i==N_X11) {
	    N_X11++;
	    g2_X11_dev=g2_realloc(g2_X11_dev,
				  sizeof(g2_X11_device)*N_X11);
	    xout=&g2_X11_dev[N_X11-1];
	    pid=N_X11-1;
	}
    }
    
    xout->width=width;			  /* set window size */
    xout->height=height;
    
    xout->NofInks=0;			  /* reset inks */
    xout->inks=NULL;
    
    vid = g2_register_physical_device(pid, NULL,
				      g2_IntCoor, g2_X11_funix,
				      1.0, -1.0,
				      0.0, height-1);
    
    sprintf(name, "g2: %d", vid);	   /* set window and icon names */
    if(window_name==NULL)
	window_name=name;
    if(icon_name==NULL)
	icon_name=name;
    
    g2_X11_init_X11X(pid, width, height,
		     x, y,
      		     window_name, icon_name,
		     icon_data, icon_width, icon_height);
    
    
					    /* g2 calls */
    g2_allocate_basic_colors(vid);
    g2_set_background(vid, 0);
    g2_pen(vid, 1);
    
    return vid;
}

     


/*
 *
 * Extended version of the InitX11
 *                      
 */
int g2_X11_init_X11X(int pid, int width, int height,
		     int xposition, int yposition,
      		     char *window_name, char *icon_name,
		     char *icon_data,
		     unsigned int icon_width, unsigned int icon_height)
{
    g2_X11_device *xout=&g2_X11_dev[pid];
    Window root;
    XSetWindowAttributes wattr;
    XEvent event;
    Pixmap iconPixmap;
    XSizeHints sizehints;
    int xhot, yhot, rv;
    XColor w_scr, w_exa, r_scr, r_exa;
    XClassHint class_hint;

    if((xout->display=XOpenDisplay(NULL))==NULL) { 
	g2_log(Error, "g2: can't open X11 display - check DISPLAY environment variable\n");
	exit(-1);
    }


    xout->root=RootWindow(xout->display, DefaultScreen(xout->display));
    root=xout->root;
    
    wattr.event_mask=ExposureMask;
    xout->window=XCreateWindow(xout->display, root,
			       xposition, yposition,
			       width, height,
			       0,
			       CopyFromParent, InputOutput, CopyFromParent,
			       CWEventMask,
			       &wattr);

    xout->gc=XCreateGC(xout->display, xout->window,
		       0lu, NULL);

    xout->colormap=DefaultColormap(xout->display,
				   DefaultScreen(xout->display));
        
    XAllocNamedColor(xout->display, xout->colormap,
		     "red", &r_scr, &r_exa);
    XAllocNamedColor(xout->display, xout->colormap,
		     "white", &w_scr, &w_exa);


    
    if(icon_data!=NULL) {
	if(icon_width<=0 || icon_height<=0) {	  /* read icon from file */
	    rv=XReadBitmapFile(xout->display, xout->window,
			       icon_data, &icon_width, &icon_height,
			       &iconPixmap, &xhot, &yhot);
	} else {				  /* icon is bitmap */
	    iconPixmap=XCreatePixmapFromBitmapData(xout->display,
						   xout->window,
						   icon_data,
						   icon_width, icon_height,
						   1ul, 0ul, 1);
	    rv=BitmapSuccess;
	}
	switch(rv) {
	  case BitmapOpenFailed:
	    fputs("g2(OpenXX): bitmap file open failed\n",stderr);
	    break;
	  case BitmapFileInvalid:
	    fputs("g2(OpenXX): bitmap file invalid\n",stderr);
	    break;
	  case BitmapNoMemory:
	    fputs("g2(OpenXX): no enough memeory for bitmap\n",stderr);
	    break;
	}
    } else {					  /* no icon data avail. */
	unsigned char bitmapData[] = {
	    0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
	    0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
	    0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
	    0x00,0x00,0x00,0x00,0x00,0x80,0x03,0x00,0x00,0x00,0xe0,
	    0x0d,0x00,0x00,0x00,0x60,0x0c,0x00,0x00,0x00,0x20,0x18,
	    0x00,0x00,0x00,0x00,0x10,0x00,0xf8,0xe3,0x07,0x08,0x00,
	    0xfe,0xfa,0x07,0x0c,0x00,0xbf,0x6e,0x07,0x06,0x80,0x0f,
	    0xf5,0x00,0x01,0x80,0x05,0x34,0x80,0x09,0xc0,0x03,0x78,
	    0xe0,0x18,0x80,0x00,0x70,0xe0,0x1e,0xc0,0x01,0x70,0x70,
	    0x1b,0xc0,0x01,0x50,0x00,0x00,0xc0,0x01,0x70,0x00,0x00,
	    0xc0,0x00,0x70,0x00,0x00,0x40,0x03,0x38,0x00,0x00,0x80,
	    0x05,0x50,0x00,0x00,0x80,0x0a,0x6e,0x00,0x00,0x00,0xfe,
	    0x37,0x00,0x00,0x00,0x6a,0x59,0x00,0x00,0x00,0xbc,0x57,
	    0x00,0x00,0x00,0xe0,0x50,0x00,0x00,0x00,0x00,0x60,0x00,
	    0x00,0x00,0x00,0x78,0x00,0x00,0x00,0x00,0x7e,0x00,0x00,
	    0x00,0x7c,0x2b,0x00,0x00,0x00,0xa4,0x1d,0x00,0x00,0x00,
	    0xb8,0x06,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
	    0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
	    0x00,0x00};
	iconPixmap=XCreatePixmapFromBitmapData(xout->display, xout->window,
					       (char*)bitmapData,
					       40u, 40u,
					       w_scr.pixel, r_scr.pixel,
					       1ul);
    }

    sizehints.x = xposition;
    sizehints.y = yposition;
    sizehints.min_width = width;
    sizehints.max_width = width;
    sizehints.min_height = height;
    sizehints.max_height = height;
    sizehints.flags = PPosition | PMinSize | PMaxSize;
    XSetStandardProperties(xout->display, xout->window,
			   window_name, icon_name, iconPixmap,
			   (char **)NULL, 0, &sizehints);
    
    class_hint.res_name =  "g2";
    class_hint.res_class = "G2";
    XSetClassHint(xout->display, xout->window, &class_hint);

 
    XMapRaised(xout->display, xout->window);
    
    /*    XSetWindowBackground(xout->display, xout->window, w_scr.pixel); */
    XClearWindow(xout->display,xout->window);
    
    g2_X11_paper(pid, NULL, 0);
    g2_X11_set_font_size(pid, NULL, 12);
						  /* wait expose event */
						  /* (no back. store) */
    while(!XCheckWindowEvent(xout->display,xout->window,
			     ExposureMask,&event))
	;


    wattr.event_mask=NoEventMask;		  /* set NoEventMask */
    wattr.backing_store=Always;			  /* set backing store */
    XChangeWindowAttributes(xout->display, xout->window,
			    CWEventMask|CWBackingStore,
			    &wattr);

    xout->dest = xout->window;

    xout->backing_pixmap = None;

    if(XDoesBackingStore(XDefaultScreenOfDisplay(xout->display))!=Always) {
	if(g2_EmulateBackingStore) {
	    g2_log(Warning, "g2: Warning! Backing store is not available. Allocating pixmap instead.\n");

	    xout->backing_pixmap = XCreatePixmap(xout->display, xout->window,
						 xout->width, xout->height,
						 DefaultDepth(xout->display, DefaultScreen(xout->display)));
	    XSetWindowBackgroundPixmap(xout->display, xout->window,
				       xout->backing_pixmap);
	  
	    XSetForeground (xout->display, xout->gc, w_scr.pixel);
          
	    XFillRectangle(xout->display, xout->backing_pixmap, xout->gc,
			   0, 0, xout->width, xout->height);
	    xout->dest = xout->backing_pixmap;
	} else {
		g2_log(Warning, "g2: Warning! Backing store is not available.\n");
	}
    }

    XFlush(xout->display);
    return 0;
}




int g2_X11_delete(int pid, void *pdp)
{
    g2_X11_device *xout=&g2_X11_dev[pid];
    XUnmapWindow(xout->display, xout->window);
    if (xout->backing_pixmap != None)
	XFreePixmap(xout->display,xout->backing_pixmap);
    XDestroyWindow(xout->display, xout->window);
    XDestroyWindow(xout->display, xout->root);
    XFreeGC(xout->display, xout->gc);
    XFreeColormap(xout->display, xout->colormap);
    XCloseDisplay(xout->display);
    if(xout->inks!=NULL)
	g2_free(xout->inks);
    xout->width=xout->height=0;
    xout->display=NULL;
    return 0;
}



int g2_X11_clear(int pid, void *pdp)
{
    g2_X11_device *xout=&g2_X11_dev[pid];
    
    if (xout->backing_pixmap == None) {
	XClearWindow(xout->display,xout->window);
    } else {
	XSetForeground (xout->display, xout->gc,
			xout->background);
	XFillRectangle(xout->display, xout->dest, xout->gc,
		       0, 0, xout->width, xout->height);
    }
    g2_X11_flush(pid, pdp);
    return 0;
}



int g2_X11_flush(int pid, void *pdp)
{
    g2_X11_device *xout=&g2_X11_dev[pid];
    if( xout->backing_pixmap != None ) {
	XCopyArea(xout->display, xout->dest, xout->window, xout->gc,
		  0, 0, xout->width, xout->height, 0, 0);
    }
    XFlush(xout->display);
    return 0;
}



int g2_X11_ink(int pid, void *pdp,
	       double red, double green, double blue)
{
    g2_X11_device *xout=&g2_X11_dev[pid];
    XColor color;

    color.flags=DoRed|DoGreen|DoBlue;
    color.red   = (int)(red   * USHRT_MAX);
    color.green = (int)(green * USHRT_MAX);
    color.blue  = (int)(blue  * USHRT_MAX);
    if(XAllocColor(xout->display,xout->colormap,&color)) {
	xout->NofInks++;
	if(xout->inks==NULL)
	    xout->inks=
	      (unsigned long *)g2_malloc(xout->NofInks*sizeof(unsigned long));
	else
	    xout->inks=
	      (unsigned long *)g2_realloc((void *)xout->inks,
					  xout->NofInks*sizeof(unsigned long));
	if(xout->inks==NULL) {
	    fputs("g2: not enough memory\n",stderr);
	    return -1;
	}
	xout->inks[xout->NofInks-1]=color.pixel;
	return xout->NofInks-1;
    } else {
	fputs("g2: color is not available\n",stderr);
	return -1;
    }
}


int g2_X11_clear_palette(int pid, void *pdp)
{
    g2_X11_device *xout=&g2_X11_dev[pid];
    XFreeColors(xout->display,xout->colormap,
		xout->inks,xout->NofInks,0x0ul);
    xout->NofInks=0;
    if(xout->inks!=NULL)
        free(xout->inks);
    xout->inks=NULL;
    return 0;
}


int g2_X11_set_background(int pid, void *pdp, int color)
{
    g2_X11_device *xout=&g2_X11_dev[pid];
    if(color>=xout->NofInks || color<0)
	return -1;
    if (xout->backing_pixmap == None)
      {
        XSetWindowBackground(xout->display,xout->dest,
			 xout->inks[color]);
      }
    else
      {
	xout->background = xout->inks[color];
      }
    g2_X11_clear(pid,pdp);
    return 0;
}


int g2_X11_pen(int pid, void *pdp, int color)
{
    g2_X11_device *xout=&g2_X11_dev[pid];
    if(color>=xout->NofInks || color<0)
	return -1;
    XSetForeground(xout->display, xout->gc, xout->inks[color]);
    return 0;
}



int g2_X11_paper(int pid, void *pdp, int color)
{
    g2_X11_device *xout=&g2_X11_dev[pid];
    if(color>=xout->NofInks || color<0)
	return -1;
    XSetBackground(xout->display, xout->gc, xout->inks[color]);
    return 0;
}


int g2_X11_set_line_width(int pid, void *pdp, int w)
{
    g2_X11_device *xout=&g2_X11_dev[pid];
    XGCValues val;
    
    val.line_width=w;
    XChangeGC(xout->display, xout->gc, GCLineWidth, &val);
    return 0;
}


int g2_X11_set_dash(int pid, void *pdp, int n, int *data)
{
    g2_X11_device *xout=&g2_X11_dev[pid];
    XGCValues val;
    int i;
    
    if(n<=0 || data==NULL) {
	val.line_style=LineSolid;
	XChangeGC(xout->display, xout->gc, GCLineStyle,&val);
    } else {
	char *ch_data;
	ch_data=g2_malloc(n*sizeof(char));
	val.line_style=LineOnOffDash;
	for(i=0;i<n;i++)
	    if(data[i]>0)
		ch_data[i]=(char)data[i];
	    else
		ch_data[i]=1;
	XChangeGC(xout->display, xout->gc, GCLineStyle, &val);
	XSetDashes(xout->display, xout->gc, 0, ch_data, n);
	g2_free(ch_data);
    }
    return 0;
}


int g2_X11_set_font_size(int pid, void *pdp, int size)
{
    g2_X11_device *xout=&g2_X11_dev[pid];
    XFontStruct *fnt_str;
    char font_name[256];
    int sizei, d, n;
    
    sizei=dtoi(size);

    if(sizei<=0)
	sizei=1;				  /* set to smallest size */
    
    for(n=1;n<32;n++) {
	d=((n&0x01)? -1:1)*(n>>1);
	sprintf(font_name, g2_X11Font, sizei+d);
	fnt_str=XLoadQueryFont(xout->display, font_name);
	if(fnt_str==NULL) {
	    if(!d)
		fprintf(stderr,"g2: can not load font: '%s'\n",font_name);
	} else {
	    XSetFont(xout->display,xout->gc,fnt_str->fid);
	    if(d)
		fprintf(stderr,"g2: using '%s' instead\n",font_name);
	    return 0;
	}
    }
    fprintf(stderr, "g2: are you sure about %d point size\n", size);
    return -1;
}


int g2_X11_plot(int pid, void *pdp, int x, int y)
{
    g2_X11_device *xout=&g2_X11_dev[pid];
    XDrawPoint(xout->display, xout->dest, xout->gc,
	       x, y);
    return 0;
}


int g2_X11_line(int pid, void *pdp, int x1, int y1, int x2, int y2)
{
    g2_X11_device *xout=&g2_X11_dev[pid];
    XDrawLine(xout->display,xout->dest,xout->gc,
	      x1, y1, x2, y2);
    return 0;
}


int g2_X11_poly_line(int pid, void *pdp, int N, int *p)
{
    g2_X11_device *xout=&g2_X11_dev[pid];
    XPoint *points;
    int i;
    points=g2_malloc(N*sizeof(XPoint));
    for(i=0;i<N;i++) {
	points[i].x=(short)p[i*2];
	points[i].y=(short)p[i*2+1];
    }
    XDrawLines(xout->display,xout->dest,xout->gc,
	       points, N,
	       CoordModeOrigin);
    g2_free(points);
    return 0;
}


int g2_X11_polygon(int pid, void *pdp, int N, int *p)
{
    g2_X11_device *xout=&g2_X11_dev[pid];
    XPoint *points;
    int i;
    points=g2_malloc((N+1)*sizeof(XPoint));
    for(i=0;i<N;i++) {
	points[i].x=(short)p[i*2];
	points[i].y=(short)p[i*2+1];
    }
    points[N].x=(short)p[0];
    points[N].y=(short)p[1];
    XDrawLines(xout->display,xout->dest,xout->gc,
	       points, N+1,
	       CoordModeOrigin);
    g2_free(points);
    return 0;
}


int g2_X11_filled_polygon(int pid, void *pdp, int N, int *p)
{
    g2_X11_device *xout=&g2_X11_dev[pid];
    XPoint *points;
    int i;
    points=g2_malloc((N+1)*sizeof(XPoint));
    for(i=0;i<N;i++) {
	points[i].x=(short)p[i*2];
	points[i].y=(short)p[i*2+1];
    }
    points[N].x=(short)p[0];
    points[N].y=(short)p[1];
    XFillPolygon(xout->display,xout->dest,xout->gc,
		 points, N+1,
		 Complex, CoordModeOrigin);
    g2_free(points);
    return 0;
}


int g2_X11_triangle(int pid, void *pdp,
		    int x1, int y1,
		    int x2, int y2,
		    int x3, int y3)
{
    g2_X11_device *xout=&g2_X11_dev[pid];
    XPoint points[4];
    points[0].x=x1; points[0].y=y1; 
    points[1].x=x2; points[1].y=y2; 
    points[2].x=x3; points[2].y=y3; 
    points[3].x=x1; points[3].y=y1; 
    XDrawLines(xout->display,xout->dest,xout->gc,
	       points, 4, CoordModeOrigin);
    return 0;
}


int g2_X11_filled_triangle(int pid, void *pdp,
			   int x1, int y1,
			   int x2, int y2,
			   int x3, int y3)
{
    g2_X11_device *xout=&g2_X11_dev[pid];
    XPoint points[4];
    points[0].x=x1; points[0].y=y1; 
    points[1].x=x2; points[1].y=y2; 
    points[2].x=x3; points[2].y=y3; 
    points[3].x=x1; points[3].y=y1; 
    XFillPolygon(xout->display,xout->dest,xout->gc,
		 points, 4, Convex, CoordModeOrigin);
    return 0;
}


int g2_X11_rectangle(int pid, void *pdp, int x1, int y1, int x2, int y2)
{
    g2_X11_device *xout=&g2_X11_dev[pid];
    XDrawRectangle(xout->display,xout->dest,xout->gc,
		   x1, y1, x2-x1, y2-y1);
    return 0;
}



int g2_X11_filled_rectangle(int pid, void *pdp, int x1, int y1, int x2, int y2)
{
    g2_X11_device *xout=&g2_X11_dev[pid];
    XDrawRectangle(xout->display,xout->dest,xout->gc,
		   x1, y1, x2-x1, y2-y1);
    XFillRectangle(xout->display,xout->dest,xout->gc,
		   x1, y1, x2-x1, y2-y1);
    return 0;
}


int g2_X11_arc(int pid, void *pdp, int x, int y,
	       int r1, int r2, double a1, double a2)
{
    g2_X11_device *xout=&g2_X11_dev[pid];
    double a0, d;

    a0=fmod(a1, 360.) + (a1<0? 360:0);   /* map a1 to [0, 360) */
    d=a2>a1? a2-a1:a2-a1+360;
     
    XDrawArc(xout->display,xout->dest,xout->gc,
	     x-r1, y-r2,
	     r1*2, r2*2,
	     (int)(a0*64.), (int)(d*64.));
    return 0;
}
 

int g2_X11_filled_arc(int pid, void *pdp, int x, int y,
		      int r1, int r2, double a1, double a2)
{
    g2_X11_device *xout=&g2_X11_dev[pid];
    double a0, d;
    
    a0=fmod(a1, 360.) + (a1<0? 360:0);   /* map a1 to [0, 360) */
    d=a2>a1? a2-a1:a2-a1+360;
    
    XDrawArc(xout->display,xout->dest,xout->gc,
	     x-r1, y-r2,
	     r1*2, r2*2,
	     (int)(a0*64.), (int)(d*64.));
    XFillArc(xout->display,xout->dest,xout->gc,
	     x-r1, y-r2,
	     r1*2, r2*2,
	     (int)(a0*64.), (int)(d*64.));
    return 0;
}
 

int g2_X11_ellipse(int pid, void *pdp, int x, int y, int r1, int r2)
{
    g2_X11_device *xout=&g2_X11_dev[pid];
    XDrawArc(xout->display,xout->dest,xout->gc,
	     x-r1, y-r2,
	     r1*2, r2*2,
	     0,360*64);
    return 0;
}


int g2_X11_filled_ellipse(int pid, void *pdp, int x, int y, int r1, int r2)
{
    g2_X11_device *xout=&g2_X11_dev[pid];
    XDrawArc(xout->display,xout->dest,xout->gc,
	     x-r1, y-r2,
	     r1*2, r2*2,
	     0,360*64);
    XFillArc(xout->display,xout->dest,xout->gc,
	     x-r1, y-r2,
	     r1*2, r2*2,
	     0,360*64);
    XFlush(xout->display);
    return 0;
}
 

int g2_X11_draw_string(int pid, void *pdp, int x, int y, const char *text)
{
    g2_X11_device *xout=&g2_X11_dev[pid];
    XDrawString(xout->display,xout->dest,xout->gc,
		x, y, text, strlen(text));
    return 0;
}


int g2_X11_image(int pid, void *pdp,
		 int x, int y, int width, int height, int *pen_array)
{
    g2_X11_device *xout=&g2_X11_dev[pid];
    XImage *image=NULL;
    Screen *screen;
    unsigned long *ink_array;
    int i;

    ink_array=malloc(sizeof(unsigned long)*width*height);
    for(i=0;i<width*height;i++)
	ink_array[i]=xout->inks[pen_array[i]];

    screen=DefaultScreenOfDisplay(xout->display);
    image=XCreateImage(xout->display,
		       DefaultVisualOfScreen(screen),
		       DefaultDepthOfScreen(screen),
		       ZPixmap,
		       0,			  /* offset */
		       (char *)ink_array,
		       width, height,
		       sizeof(unsigned long)*8,	  /* bitmap pad */
		       0);			  /* bytes per line */
    /* XInitImage(image); problems with AIX ?!! */
    XPutImage(xout->display, xout->dest, xout->gc,
	      image,
	      0, 0,
	      x, y, width, height);
	
    XDestroyImage(image);
    free(ink_array);
    return 0;
}


int g2_X11_query_pointer(int pid, void *pdp,
			 int *x, int *y, unsigned int *button)
{
    Bool rv;
    g2_X11_device *xout=&g2_X11_dev[pid];
    Window root, child;
    int rx, ry;

    rv = XQueryPointer(xout->display, xout->window,
		       &root, &child, &rx, &ry,
		       x, y, button);

    if(rv)
	return 0;
    else
	return 1;
}



int g2_X11_get_pd_handles(int pid, void *pdp, void *handles[G2_PD_HANDLES_SIZE])
{
    g2_X11_device *xout=&g2_X11_dev[pid];

    handles[0]=xout->display;
    handles[1]=&xout->window;
    handles[2]=&xout->root;
    handles[3]=&xout->colormap;
    handles[4]=&xout->gc;
    handles[5]=&xout->dest;
    return 0;
}
