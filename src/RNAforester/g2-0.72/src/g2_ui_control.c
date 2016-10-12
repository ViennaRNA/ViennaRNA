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

#include "g2.h"
#include "g2_device.h"
#include "g2_control_pd.h"
#include "g2_util.h"

/**
 * \ingroup interface
 * \defgroup color color manipulations
 *
 * The color concept used in the g2 library is inspired by Sir Clive
 * Sinclair solution implemented in the ZX Spectrum computer. With the
 * g2_pen() function it is possible to choose a pen created by the
 * g2_ink() function. Note that g2_ink function is only defined for
 * physical devices. The predefined colors (see g2_test demo program)
 * have pens from 0 till 26 (inclusive).
 *
 * Some basic colors are: 
 *  - 0 white
 *  - 1 black
 *  - 3 blue
 *  - 7 green
 *  - 19 red
 *  - 25 yellow
 *
 */


/**
 * \ingroup interface
 * \defgroup control output control
 */
 
/**
 *
 * Flush output buffers.
 *
 * \param dev device id
 *
 * \ingroup control
 */
void g2_flush(int dev)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_flush: No such device: %d\n", dev);
	return;
    }
    
    switch(devp->t) {
      case g2_PD:
	g2_flush_pd(devp->d.pd);
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_flush(devp->d.vd->dix[i]);
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
 * Save output
 *
 * \param dev device id
 *
 * \ingroup control
 */
void g2_save(int dev)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_save: No such device: %d\n", dev);
	return;
    }
    
    switch(devp->t) {
      case g2_PD:
	g2_save_pd(devp->d.pd);
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_save(devp->d.vd->dix[i]);
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
 * Clear device
 *
 * \param dev device number
 *
 * \ingroup control
 */
void g2_clear(int dev)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_clear: No such device: %d\n", dev);
	return;
    }
    
    switch(devp->t) {
      case g2_PD:
	g2_clear_pd(devp->d.pd);
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_clear(devp->d.vd->dix[i]);
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
 * Set pen color for all following operations, see also g2_ink().
 *
 * \param dev device
 * \param color pen (either one of default pens 0-26, or a pen returned by g2_ink() )
 *
 * \ingroup color
 */
void g2_pen(int dev, int color)
{
    g2_device *devp;
    int i;

    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_pen: No such device: %d\n", dev);
	return;
    }
    
    switch(devp->t) {
      case g2_PD:
	g2_pen_pd(devp->d.pd, color);
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_pen(devp->d.vd->dix[i], color);
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
 * Set the background color
 *
 * \param dev device
 * \param color pen (either one of default pens 0-26, or a pen returned by g2_ink() )
 *
 * \ingroup color
 */
void g2_set_background(int dev, int color)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_set_background: No such device: %d\n", dev);
	return;
    }
    
    switch(devp->t) {
      case g2_PD:
	g2_set_background_pd(devp->d.pd, color);
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_set_background(devp->d.vd->dix[i], color);
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
 * Create an ink. To put ink into the pen use g2_pen().
 *
 * \param pd_dev physical device
 * \param red red component (0-1) according to the RGB color model
 * \param green green component (0-1) according to the RGB color model
 * \param blue blue component (0-1) according to the RGB color model
 * \return new pen, see g2_pen()
 *
 * \ingroup color
 */
int g2_ink(int pd_dev, double red, double green, double blue)
{
    g2_device *devp;
    int rv=-1;

    if((devp=g2_get_device_pointer(pd_dev))==NULL) {
	fprintf(stderr, "g2_ink: No such device: %d\n", pd_dev);
	return -1;
    }

    if(red   < 0.0) red=0.0;
    if(green < 0.0) green=0.0;
    if(blue  < 0.0) blue=0.0;
    if(red   > 1.0) red=1.0;
    if(green > 1.0) green=1.0;
    if(blue  > 1.0) blue=1.0;
    
    switch(devp->t) {
      case g2_PD:
	rv=g2_ink_pd(devp->d.pd, red, green, blue);
	break;
      case g2_VD:
	fprintf(stderr, "g2_ink: g2_ink is enabled only for phys. devices\n");
	break;
      case g2_ILLEGAL:
	break;
      case g2_NDEV:
	break;
    }
    __g2_last_device=pd_dev;
    return rv;
}



/**
 *
 * Clear collor palette (remove all inks) and reallocate basic colors.
 *
 * \param dev device
 *
 * \ingroup color
 */
void g2_reset_palette(int dev)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_reset_palette: No such device: %d\n", dev);
	return;
    }
    
    switch(devp->t) {
      case g2_PD:
	g2_clear_palette(dev);
	g2_allocate_basic_colors(dev);
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_reset_palette(devp->d.vd->dix[i]);
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
 * Remove all inks.
 *
 * \param dev device
 *
 * \ingroup color
 */
void g2_clear_palette(int dev)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_clear_palette: No such device: %d\n", dev);
	return;
    }
    
    switch(devp->t) {
      case g2_PD:
	g2_clear_palette_pd(devp->d.pd);
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_clear_palette(devp->d.vd->dix[i]);
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
 * Allocate basic colors
 *
 * \param dev device
 *
 * \ingroup color
 */
void g2_allocate_basic_colors(int dev)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_allocate_basic_colors: No such device: %d\n", dev);
	return;
    }
    
    switch(devp->t) {
      case g2_PD:
	g2_allocate_basic_colors_pd(devp->d.pd);
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_allocate_basic_colors(devp->d.vd->dix[i]);
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
 * Set font size
 *
 * \param dev device
 * \param size new font size
 *
 * \ingroup control
 */
void g2_set_font_size(int dev, double size)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_set_font_size: No such device: %d\n", dev);
	return;
    }
    
    switch(devp->t) {
      case g2_PD:
	g2_set_font_size_pd(devp->d.pd, size);
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_set_font_size(devp->d.vd->dix[i], size);
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
 * Set line width.
 *
 * \param dev device
 * \param w new line width
 *
 * \ingroup control
 */
void g2_set_line_width(int dev, double w)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_set_line_width: No such device: %d\n", dev);
	return;
    }
    
    switch(devp->t) {
      case g2_PD:
	g2_set_line_width_pd(devp->d.pd, w);
	break;
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_set_line_width(devp->d.vd->dix[i], w);
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
 * Set line dash. Set \a N to 0 and \a dashes to NULL to restore solid line.
 *
 * \param dev device
 * \param N number of dash components (0 for solid line)
 * \param dashes vector of dash lengths (black, white, black, ...)
 *
 * \ingroup control
 */
void g2_set_dash(int dev, int N, double *dashes)
{
    g2_device *devp;
    int i;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_set_dash: No such device: %d\n", dev);
	return;
    }
    
    switch(devp->t) {
      case g2_PD:
	g2_set_dash_pd(devp->d.pd, N, dashes);
	break;
	
      case g2_VD:
	for(i=0;i<devp->d.vd->N;i++)
	    g2_set_dash(devp->d.vd->dix[i], N, dashes);
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
 * Set QuasiPixel size and shape.
 *
 * \param dev device
 * \param d size
 * \param shape shape (rectangle or circle, see ::QPshape )
 *
 * \ingroup control
 */
void g2_set_QP(int dev, double d, enum QPshape shape)
{
    g2_device *devp;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_set_QP: No such device: %d\n", dev);
	return;
    }
    
    devp->QPd=d;
    devp->QPshape=shape;
    
    __g2_last_device=dev;
}


/**
 *
 * Query pointer (e.g. mouse for X11) position and button state. See
 * the demo program pointer.c for an example.
 *
 * \param dev device
 * \param x returns pointer x coordinate
 * \param y returns pointer y coordinate
 * \param button returns button state
 *
 * \ingroup control
 */
void g2_query_pointer(int dev, double *x, double *y, unsigned int *button)
{
    g2_device *devp;
    
    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_query_pointer: No such device: %d\n", dev);
	return;
    }
    
    switch(devp->t) {
      case g2_PD:
	g2_query_pointer_pd(devp->d.pd, x, y, button);
	break;
      case g2_VD:
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
 * Get pointers to physical device specific handles. This function
 * should be used only if you are familiar with the g2 source code.
 * For details see physical device source code (e.g. in src/X11/).
 * Example usage can be found in demo/handles.c.
 *
 * \param pd physical device
 * \param handles returns pointers to physical device low level handles
 *
 * \ingroup control
 */
void g2_get_pd_handles(int pd, void *handles[G2_PD_HANDLES_SIZE])
{
    g2_device *devp;
    int i;

    for(i=0;i<G2_PD_HANDLES_SIZE;i++) {
	handles[i]=NULL;
    }
    if((devp=g2_get_device_pointer(pd))==NULL) {
	g2_log(Error, "g2: Error! g2_get_pd_handles: No such device: %d\n", pd);
	return;
    }
    
    switch(devp->t) {
      case g2_PD:
	g2_get_pd_handles_pd(devp->d.pd, handles);
	break;
      case g2_VD:
	break;
      case g2_ILLEGAL:
	break;
      case g2_NDEV:
	break;
    }
}
