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
#include <stdlib.h>
#include <stdio.h>
#include "g2_physical_device.h"
#include "g2_funix.h"
#include "g2_util.h"


g2_physical_device *g2_create_physical_device(int pid,
					      void *pdp,
					      g2_coor ct,
					      const g2_funix_fun *ff,
					      double a11, double a22,
					      double b1,  double b2)
{
    g2_physical_device *rd;
    int i, j;

    rd=g2_malloc(sizeof(g2_physical_device));

    rd->pid=pid;	      /* physical device id (handled by driver) */
    rd->pdp=pdp;	      /* pointer to something */
    rd->coor_type=ct;	      /* coord. type */
    rd->a11=a11;	      /* device->physical device transformation */
    rd->a22=a22;
    rd->b1=b1;
    rd->b2=b2;

    rd->x_origin=0.0;	      /* User coordinates specification */
    rd->y_origin=0.0;
    rd->x_mul=1.0;
    rd->y_mul=1.0;

    rd->ff=g2_malloc(G2_N_FUNIX*sizeof(g2_funix_fun));

    for(i=0;i<G2_N_FUNIX;i++) {
	rd->ff[i].fx=i;
	rd->ff[i].fun=NULL;
	for(j=0;ff[j].fx!=g2_FUNIX_NULL;j++)
	    if(ff[j].fx==i) {
		rd->ff[i].fun = ff[j].fun;
		break;
	    }
    }
    
    return rd;
}


/*
 *
 * Destroy physical device
 *
 */
void g2_destroy_physical_device(g2_physical_device *pd)
{
    g2_free(pd->ff);
    g2_free(pd);
}


