/*****************************************************************************
**  Copyright (C) 1998-2005  Ljubomir Milanovic & Horst Wagner
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
#include <g2.h>
#include <g2_X11.h>
#include <g2_PS.h>

#include <X11/Xlib.h>

int main()
{
    int d1, d2;
    void *handles[G2_PD_HANDLES_SIZE];
    Display *display;
    Drawable *drawable;
    GC *gc;


    d1=g2_open_X11(100, 100);
    g2_line(d1, 10, 10, 90, 90);

    g2_get_pd_handles(d1, handles);
    display = handles[0];
    drawable = handles[5];
    gc = handles[4];
    
    printf("%p %p %p\n", display, drawable, gc);
    
    XDrawString(display, *drawable, *gc,
		10,10, "OK", 2);

    g2_flush(d1);

    d2=g2_open_PS("handles.ps", g2_A4, g2_PS_port);
    g2_get_pd_handles(d2, handles);
    printf("%p %p %p\n", handles[0], handles[1], handles[2]);

    
    getchar();
    return 0;
}

