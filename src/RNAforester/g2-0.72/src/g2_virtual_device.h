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
#ifndef _G2_VIRTUAL_DEVICE_H
#define _G2_VIRTUAL_DEVICE_H

typedef struct _g2_virtual_device {
    int  N;		             /* number of attached devices */
    int  *dix;			     /* index of attached devices */
} g2_virtual_device;


g2_virtual_device *g2_create_virtual_device();
void g2_destroy_virtual_device(g2_virtual_device *vd);
int g2_is_attached(int vd, int dev);

#endif /* _G2_VIRTUAL_DEVICE_H */
