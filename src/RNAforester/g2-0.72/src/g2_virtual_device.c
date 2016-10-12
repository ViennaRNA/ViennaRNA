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
#include "g2_virtual_device.h"
#include "g2_device.h"
#include "g2_util.h"


/*
 *
 */
g2_virtual_device *g2_create_virtual_device()
{
    g2_virtual_device *rd;

    rd=g2_malloc(sizeof(g2_virtual_device));
    rd->N=0;
    rd->dix=g2_malloc(sizeof(int));		  /* later we use realloc */
    return rd;
}


/*
 *
 */
void g2_destroy_virtual_device(g2_virtual_device *vd)
{
    g2_free(vd->dix);
    g2_free(vd);
}



/*
 * Return 1 if dev is attached to vd
 */
int g2_is_attached(int vd, int dev)
{
    g2_device *vdp, *devp;
    int i;

    if(vd==dev)
	return 1;
    
    if((devp=g2_get_device_pointer(dev))==NULL)
	return 0;

    if((vdp=g2_get_device_pointer(vd))==NULL)
	return 0;

    if(devp==vdp)
	return 1;
    
    if(vdp->t!=g2_VD)
	return 0;

    for(i=0;i<vdp->d.vd->N;i++) {
	if(vdp->d.vd->dix[i]==dev)
	    return 1;
	if(g2_is_attached(vdp->d.vd->dix[i], dev))
	    return 1;
    }
    
    return 0;
}

