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
#include "g2_virtual_device.h"
#include "g2_util.h"

/**
 * \ingroup interface
 * \defgroup vd virtual device related functions
 *
 * Virtual device is a method to redirect g2 output to multiple devices. Here is an example:
 *
 * \code
 *    int d1 = g2_open_X11(100, 100);      create first X11 window
 *    int d2 = g2_open_X11(100, 100);      create 2nd X11 window
 *
 *    int vd = g2_open_vd();               open a new virtual device
 *
 *    g2_attach(vd, d1);                   attach d1 (1st window) to virtual device
 *    g2_attach(vd, d2);                   attach d2 (2nd window) to virtual device
 *
 *    g2_plot(d1, 11, 11);                 output to the 1st X11 window
 *    g2_plot(d2, 12, 12);                 output to the 2nd X11 window 
 *    g2_plot(vd, 13, 13);                 output to both X11 windows
 * \endcode
 *
 */



/**
 *
 * Create a new virtual device.
 *
 * \return virtual device ID
 *
 * \ingroup vd
 */
int g2_open_vd(void)
{
    int dix;
    dix=g2_register_virtual_device();
    __g2_last_device=dix;
    return dix;
}


/**
 * Attach a device to virtual device \a vd_dev.
 *
 * \param vd_dev virtual device (create virtual device by calling g2_open_vd() )
 * \param dev device
 *
 * \ingroup vd
 */
void g2_attach(int vd_dev, int dev)
{
    g2_device *vd_devp, *devp;

    if((vd_devp=g2_get_device_pointer(vd_dev))==NULL) {
	fprintf(stderr, "g2_attach: No such device: %d\n", vd_dev);
	return;
    }

    if((devp=g2_get_device_pointer(dev))==NULL) {
	fprintf(stderr, "g2_attach: No such device: %d\n", dev);
	return;
    }

    if(vd_devp->t!=g2_VD) {
	fprintf(stderr, "g2_attach: Device %d is not virtual.\n", vd_dev);
	return;
    }

    if(devp->t==g2_VD)				  /* if virtual device */
	if(g2_is_attached(dev, vd_dev)) {	  /* check recurency */
	    fprintf(stderr,
		    "g2_attach: Device %d is already attached to %d.\n",
		    dev, vd_dev);
	    return;
	}

    vd_devp->d.vd->N++;
    vd_devp->d.vd->dix=g2_realloc(vd_devp->d.vd->dix,
				  vd_devp->d.vd->N*sizeof(int));

    vd_devp->d.vd->dix[vd_devp->d.vd->N-1]=dev;

    __g2_last_device=vd_dev;
}


/**
 *
 * Dettach a device from the virtual device \a vd_dev.
 *
 * \param vd_dev virtual device
 * \param dev device
 *
 * \ingroup vd
 */
void g2_detach(int vd_dev, int dev)
{
    g2_device *vd_devp;
    int i;

    if((vd_devp=g2_get_device_pointer(vd_dev))==NULL) {
	fprintf(stderr, "g2_detach: No such device: %d\n", vd_dev);
	return;
    }

    if(vd_devp->t!=g2_VD) {
	fprintf(stderr, "g2_detach: Device %d is not virtual.\n", vd_dev);
	return;
    }

    for(i=0;i<vd_devp->d.vd->N;i++)
	if(vd_devp->d.vd->dix[i]==dev) {
	    if(vd_devp->d.vd->N>1)
		vd_devp->d.vd->dix[i]=vd_devp->d.vd->dix[vd_devp->d.vd->N-1];
	    vd_devp->d.vd->N--;
	    if(vd_devp->d.vd->N!=0)
		vd_devp->d.vd->dix=g2_realloc(vd_devp->d.vd->dix,
					      vd_devp->d.vd->N*sizeof(int));
	    return;
	}

    __g2_last_device=vd_dev;
}

