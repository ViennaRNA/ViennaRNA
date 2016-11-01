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
#include <g2.h>
#include <stdlib.h>
#ifdef DO_X11
#include <g2_X11.h>
#endif
#ifdef DO_WIN32
#include <g2_win32.h>
#endif


#define N 300
#define W 400
#define H 400

int main()
{
    int d,dwin,dX11;
	int x[N];
	int y[N];
	int c[N];
	int i,t;

    d=g2_open_vd();				  /* open virtual device */

    printf("Random Walk demo\n");

#ifdef DO_WIN32
    printf("Adding win32..\n");
    dwin=g2_open_win32(W, H, "simple_animation", 0);
    g2_attach(d, dwin);
#endif
#ifdef DO_X11
    printf("Adding X11.. (X11 might flicker and not be usable)\n");
    dX11=g2_open_X11(W, H);
    g2_attach(d, dX11);
#endif
    g2_set_auto_flush(d,0);

	for (i=0;i<N;i++)
	{
		x[i] = rand()%(W/4)+H/2;
		y[i] = rand()%(H/4)+W/2;
		c[i] = rand()%16;
	}

	for (t=0; t<1000; t++)
	{
		g2_clear(d);
		for (i=0; i<N; i++)
			{
			g2_filled_circle(d,(double)x[i],(double)y[i],(double)5);
			g2_pen(d,c[i]);
			x[i] += (rand()%2) * 2 - 1;
			y[i] += (rand()%2) * 2 - 1;
			}
	g2_flush(d);
	}
    getchar();
    return 0;
}
