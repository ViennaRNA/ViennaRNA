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

#include "g2.h"
#include "g2_device.h"
#include "g2_physical_device.h"
#include "g2_virtual_device.h"
#include "g2_util.h"

int __g2_last_device=-1;			  /* last acces. device (ld) */


static int       g2_dev_size=0;			  /* size of devices array */
static g2_device *g2_dev=NULL;			  /* devices array */

static int  g2_get_free_device();
static void g2_init_device(int dix);


/*
 *
 * Register physical device
 *
 */
int g2_register_physical_device(int pid,
				void *pdp,
				g2_coor ct,
				const g2_funix_fun *ff,
				double a11, double a22,
				double b1,  double b2)
{
    int dix;

    dix=g2_get_free_device();
    g2_init_device(dix);
    g2_dev[dix].t=g2_PD;
    g2_dev[dix].d.pd=g2_create_physical_device(pid, pdp,
					       ct, ff,
					       a11, a22,
					       b1,  b2);

    __g2_last_device=dix;
    return dix;
}



/*
 *
 * Register virtual device
 *
 */
int g2_register_virtual_device()
{
    int dix;
    
    dix=g2_get_free_device();
    g2_init_device(dix);    
    g2_dev[dix].t=g2_VD;
    g2_dev[dix].d.vd=g2_create_virtual_device();

    return dix;
}


/*
 *
 * Init device
 * 
 */
void g2_init_device(int dix)
{
    g2_dev[dix].t=g2_ILLEGAL;
    g2_dev[dix].dix=dix;    
    g2_dev[dix].x=0.0;				  /* set cursor */
    g2_dev[dix].y=0.0;
    g2_dev[dix].auto_flush=1;			  /* set auto flush */
    g2_dev[dix].QPd=1.0;			  /* Quasi pixel spec. */
    g2_dev[dix].QPshape=QPrect;
}

/*
 *
 * Return pointer to device dix
 *
 */
g2_device *g2_get_device_pointer(int dix)
{
    if(dix<0 || dix>=g2_dev_size)
        return NULL;
    if(g2_dev[dix].t==g2_NDEV)
        return NULL;
    
    return &g2_dev[dix];
}


/*
 *
 * Return device type
 *
 */
g2_device_type g2_get_device_type(int dix)
{
    if(dix<0 || dix>=g2_dev_size)
	return g2_ILLEGAL;
    return g2_dev[dix].t;
}


/*
 *
 * 1 if device exist otherwise 0
 *
 */
int g2_device_exist(int dix)
{
    if(dix<0 || dix>=g2_dev_size ||
       g2_dev[dix].t==g2_ILLEGAL || g2_dev[dix].t==g2_NDEV)
        return 0;
    return 1;
}


/*
 *
 * Destroy device
 *
 */
void g2_destroy_device(int dix)
{
    int i;
    
    for(i=0;i<g2_dev_size;i++)			  /* detach from all vd */
	if(g2_dev[i].t==g2_VD)
	    g2_detach(i, dix);
    
    switch(g2_dev[dix].t) {
      case g2_NDEV:
	break;
      case g2_PD:
	g2_destroy_physical_device(g2_dev[dix].d.pd);
	g2_dev[dix].t=g2_NDEV;
	break;
      case g2_VD:
	g2_destroy_virtual_device(g2_dev[dix].d.vd);
	g2_dev[dix].t=g2_NDEV;
	break;
      case g2_ILLEGAL:
	break;
    }
}



/*
 *
 * get free place for new device
 *
 */
int g2_get_free_device()
{
    int i, dix;

    if(g2_dev==NULL) {				  /* if NULL initialize */
	g2_dev_size=1;
	g2_dev=g2_malloc(sizeof(g2_device));
	g2_dev[0].t=g2_NDEV;			  /* set to free */
	g2_dev[0].d.pd=NULL;			  /* no device yet */
    }

    dix=-1;
    for(i=0;i<g2_dev_size;i++)			  /* find empty place */
      if(g2_dev[i].t==g2_NDEV) {
	  dix=i;
	  break;
      }

    if(dix<0) {					  /* no place for device */
	dix=g2_dev_size++;
	g2_dev=g2_realloc(g2_dev, g2_dev_size*sizeof(g2_device));
	g2_dev[dix].t=g2_NDEV;			  /* set to free */
	g2_dev[dix].d.pd=NULL;			  /* no device yet */
    }
    
    return dix;
}


