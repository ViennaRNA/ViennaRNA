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
#include "g2_physical_device.h"
#include "g2_util.h"
#include "g2_funix.h"



/*
 *
 * Flush output
 *
 */
void g2_flush_pd(g2_physical_device *pd)
{
    if(pd->ff[g2_Flush].fun!=NULL) {
	pd->ff[g2_Flush].fun(pd->pid, pd->pdp);
    } else {
	/* emulate ... with .... */  
    }
}



/*
 *
 * Save output
 *
 */
void g2_save_pd(g2_physical_device *pd)
{
    if(pd->ff[g2_Save].fun!=NULL) {
	pd->ff[g2_Save].fun(pd->pid, pd->pdp);
    } else {
	/* emulate ... with .... */  
    }
}



/*
 *
 * Clear device
 *
 */
void g2_clear_pd(g2_physical_device *pd)
{
    if(pd->ff[g2_Clear].fun!=NULL) {
	pd->ff[g2_Clear].fun(pd->pid, pd->pdp);
    } else {
	/* emulate ... with .... */  
    }
}



/*
 *
 * Set pen
 *
 */
void g2_pen_pd(g2_physical_device *pd, int color)
{
    if(pd->ff[g2_Pen].fun!=NULL) {
	pd->ff[g2_Pen].fun(pd->pid, pd->pdp, color);
    } else {
	/* emulate ... with .... */  
    }
}



/*
 *
 * Set background color
 *
 */
void g2_set_background_pd(g2_physical_device *pd, int color)
{
    if(pd->ff[g2_SetBackground].fun!=NULL) {
	pd->ff[g2_SetBackground].fun(pd->pid, pd->pdp, color);
    } else {
	/* emulate ... with .... */  
    }
}



/*
 *
 * Set ink
 *
 */
int g2_ink_pd(g2_physical_device *pd, double red, double green, double blue)
{
    int rv=-1;
    
    if(pd->ff[g2_Ink].fun!=NULL) {
	rv=pd->ff[g2_Ink].fun(pd->pid, pd->pdp,
			      red, green, blue);
    } else {
	/* emulate ... with .... */  
    }
    return rv;
}




/*
 *
 * Clear palette
 *
 */
void g2_clear_palette_pd(g2_physical_device *pd)
{
    if(pd->ff[g2_ClearPalette].fun!=NULL) {
	pd->ff[g2_ClearPalette].fun(pd->pid, pd->pdp);
    } else {
	/* emulate ... with .... */  
    }
}


/*
 *
 * Allocate basic colors
 *
 */
void g2_allocate_basic_colors_pd(g2_physical_device *pd)
{
    double ct[3]={0.0, 0.5, 1.0};
    int r, g, b;
    
    if(pd->ff[g2_Ink].fun!=NULL) {
	pd->ff[g2_Ink].fun(pd->pid, pd->pdp,	  /* white */
			   1.0, 1.0, 1.0);
	pd->ff[g2_Ink].fun(pd->pid, pd->pdp,	  /* black */
			   0.0, 0.0, 0.0);
	for(r=0;r<3;r++)
	    for(g=0;g<3;g++)
		for(b=0;b<3;b++)
		    if((r==2 && g==2 && b==2) ||
		       (r==0 && g==0 && b==0))
			continue;
		    else
			pd->ff[g2_Ink].fun(pd->pid, pd->pdp,
					   ct[r], ct[g], ct[b]);
    }  
}


/*
 *
 * Set font size
 *
 */
void g2_set_font_size_pd(g2_physical_device *pd, double size)
{
    int is;
    double ds;

    if(pd->ff[g2_SetFontSize].fun!=NULL) {
	switch(pd->coor_type) {
	  case g2_IntCoor:
	    is=dtoi(size*fabs(pd->a22));	  /* to pd coordinates */
	    pd->ff[g2_SetFontSize].fun(pd->pid, pd->pdp, is);
	    break;
	  case g2_DoubleCoor:
	    ds=size*fabs(pd->a22);		  /* to pd coordinates */
	    pd->ff[g2_SetFontSize].fun(pd->pid, pd->pdp, ds);
	    break;
	}
    } else {
	/* emulate ... with .... */  
    }
}



/*
 *
 * Set line width
 *
 */
void g2_set_line_width_pd(g2_physical_device *pd, double w)
{
    int iw;
    double dw;

    if(pd->ff[g2_SetLineWidth].fun!=NULL) {
	switch(pd->coor_type) {
	  case g2_IntCoor:
	    iw=dtoi(w*fabs(pd->a22));     /* to pd coordinates */
	    pd->ff[g2_SetLineWidth].fun(pd->pid, pd->pdp, iw);
	    break;
	  case g2_DoubleCoor:
	    dw=w*fabs(pd->a22);	        /* to pd coordinates */
	    pd->ff[g2_SetLineWidth].fun(pd->pid, pd->pdp, dw);
	    break;
	}
    } else {
	/* emulate ... with .... */  
    }
}



/*
 *
 * Set dash
 *
 */
void g2_set_dash_pd(g2_physical_device *pd, int N, double *dashes)
{
    int j;
    double *dd=NULL;
    int    *id=NULL;
    
    if(pd->ff[g2_SetDash].fun!=NULL) {
	switch(pd->coor_type) {
	  case g2_IntCoor:
	    if(dashes!=NULL) {
		id=g2_malloc(N*sizeof(int));
		for(j=0;j<N;j++)
		    id[j]=dtoi(dashes[j]*fabs(pd->a22));
		pd->ff[g2_SetDash].fun(pd->pid, pd->pdp, N, id);
		g2_free(id);
	    } else
		pd->ff[g2_SetDash].fun(pd->pid, pd->pdp, 0, NULL);
	    break;
	  case g2_DoubleCoor:
	    if(dashes!=NULL) {
		dd=g2_malloc(N*sizeof(double));
		for(j=0;j<N;j++)
		    dd[j]=dashes[j]*fabs(pd->a22);
		pd->ff[g2_SetDash].fun(pd->pid, pd->pdp, N, dd);
		g2_free(dd);
		break;
	    } else
		pd->ff[g2_SetDash].fun(pd->pid, pd->pdp, 0, NULL);
	}
    } else {
	/* emulate ... with .... */  
    }
}



/*
 *
 * Query pointer position and button state
 *
 */
void g2_query_pointer_pd(g2_physical_device *pd,
			 double *x, double *y, unsigned int *button)
{
    int    ix, iy;
    double dx, dy;

    if(pd->ff[g2_QueryPointer].fun!=NULL) {
	switch(pd->coor_type) {
	  case g2_IntCoor:
	    pd->ff[g2_QueryPointer].fun(pd->pid, pd->pdp,
				 &ix, &iy, button);
	    g2_pdc2uc(pd, ix, iy, x, y);
	    break;
	  case g2_DoubleCoor:
	    pd->ff[g2_QueryPointer].fun(pd->pid, pd->pdp,
				        &dx, &dy, button);
	    g2_pdc2uc(pd, dx, dy, x, y);
	    break;
	}
    } else {
	/* no emulation for query pointer */
    }
}


/*
 *
 * Get low level handles
 *
 */
void g2_get_pd_handles_pd(g2_physical_device *pd, void *handles[G2_PD_HANDLES_SIZE])
{
    if(pd->ff[g2_GetPDHandles].fun!=NULL) {
	pd->ff[g2_GetPDHandles].fun(pd->pid, pd->pdp, handles);
    } else {
	/* no emulation for low level handles */
    }
}
