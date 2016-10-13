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

#ifdef DO_PS
#include <g2_PS.h>
#endif
#ifdef DO_FIG
#include <g2_FIG.h>
#endif
#ifdef DO_X11
#include <g2_X11.h>
#endif
#ifdef DO_GD
#include <g2_gd.h>
#endif
#ifdef DO_WIN32
#include <g2_win32.h>
#endif
#ifdef DO_WMF32
#include <g2_win32.h>
#endif

#define maxdev 10

int main()
{
    int i, j;
    int d, dev[maxdev]={-1, -1, -1, -1, -1,-1, -1, -1, -1, -1};
    int ndev=0;
    char str[256];
    double pts[10];
    double y;
#include "penguin.c"

    printf("\nG2_VERSION %s\n", G2_VERSION);
    
    d=g2_open_vd();				  /* open virtual device */

    printf("Adding..");
    
#ifdef DO_PS
    printf("..PS");
    dev[ndev]=g2_open_PS("g2_test.ps", g2_A4, g2_PS_land);
    g2_attach(d, dev[ndev]);
    ndev++;
	
    printf("..EPSF");
    dev[ndev]=g2_open_EPSF("g2_test.eps");
    g2_attach(d, dev[ndev]);
    ndev++;
    printf("..EPSF_CLIP");
    dev[ndev]=g2_open_EPSF_CLIP("g2_test_clip.eps",200,200);
    g2_attach(d, dev[ndev]);
    ndev++;
#endif
#ifdef DO_FIG
    printf("..FIG");
    dev[ndev]=g2_open_FIG("g2_test.fig");
    g2_attach(d, dev[ndev]);
    ndev++;
#endif
#ifdef DO_X11
    printf("..X11");
    dev[ndev]=g2_open_X11(775, 575);
    g2_attach(d, dev[ndev]);
    ndev++;
#endif
#ifdef DO_GD
    printf("..GD(png)");
    dev[ndev]=g2_open_gd("g2_test.png", 775, 575,g2_gd_png);
    g2_attach(d, dev[ndev]);
    ndev++;
    printf("..GD(jpeg)");
    dev[ndev]=g2_open_gd("g2_test.jpg", 775, 575,g2_gd_jpeg);
    g2_attach(d, dev[ndev]);
    ndev++;
#endif
#ifdef DO_WIN32
    printf("..WIN32");
    dev[ndev]=g2_open_win32(775, 575,"g2_test",0);
    g2_attach(d, dev[ndev]);
    ndev++;
#endif
#ifdef DO_WMF32
    printf("..WMF32");
    dev[ndev]=g2_open_win32(775, 575,"g2_test.emf",1);
    g2_attach(d, dev[ndev]);
    ndev++;
#endif
    g2_set_auto_flush(d,0);

    printf("\n");

    /* g2_set_coordinate_system(d, 775, 575, -0.75, -1.0); */
    
    for(i=0;i<27;i++) {
	g2_pen(d, i);
	g2_filled_circle(d, i*20+10, 10, 10); 
	g2_pen(d, 1);
	g2_circle(d, i*20+10, 10, 10);
	sprintf(str, "%d", i);
	g2_string(d, i*20+7, 21, str);
    }

    for(j=0;j<ndev;j++)
	if(dev[j]>=0)
	    for(i=0;i<=64;i++) {
		g2_move(dev[j], 2*i+575, 5);
		g2_pen(dev[j], g2_ink(dev[j], i/64., 0, 0));
		g2_line_r(dev[j], 0, 20);
		g2_pen(dev[j], g2_ink(dev[j], 0, i/64., 0));
		g2_line_r(dev[j], 10, 20);
		g2_pen(dev[j], g2_ink(dev[j], 0, 0, i/64.));
		g2_line_r(dev[j], -10, 20);
	    }

    g2_pen(d, 1);
    g2_line(d, 200, 50, 350, 50);
    g2_line(d, 200, 48, 350, 48);
    g2_line(d, 200, 46, 350, 46);
    g2_line(d, 200, 46, 200, 75);
    g2_line(d, 198, 46, 198, 75);
    g2_line(d, 196, 46, 196, 75);
    g2_string(d, 200, 50, "012abcABC#())(\\-+~*!$%&");

    g2_pen(d, 1);
    for(i=1;i<25;i++) {
	g2_line(d, 15, i*20+50, 15, i*20+50+i);
	g2_set_font_size(d, 12);
	sprintf(str, "%2d:", i);
	g2_string(d, 20, i*20+50, str);
	g2_set_font_size(d, i);
	g2_string(d, 40, i*20+50, "hello, world");
    }

    g2_plot(d, 150, 70);
    g2_line(d, 147, 68, 153, 68);
		
    y=100;
    g2_line(d, 120, y, 170, y+50);
    g2_triangle(d, 150, y, 250, y, 200, y+50);
    g2_rectangle(d, 300, y, 400, y+50);
    g2_circle(d, 450, y+25, 25);
    g2_ellipse(d, 550, y+25, 45, 25);
    g2_arc(d, 650, y+25, 25, 45, 90, 360);
    
    pts[0]=4;
    pts[1]=4;
    g2_set_dash(d, 2, pts);
    g2_line(d, 120+5, y, 170+5, y+50);
    g2_triangle(d, 150+10, y+4, 250-10, y+4, 200, y+50-5);
    g2_rectangle(d, 305, y+5, 395, y+50-5);
    g2_circle(d, 450, y+25, 20);
    g2_ellipse(d, 550, y+25, 40, 20);
    g2_arc(d, 650, y+25, 20, 40, 90, 360);
    g2_set_dash(d, 0, NULL);

    y=200;
    g2_filled_triangle(d, 150, y, 250, y, 200, y+50);
    g2_filled_rectangle(d, 300, y, 400, y+50);
    g2_filled_circle(d, 450, y+25, 25);
    g2_filled_ellipse(d, 550, y+25, 45, 25);
    g2_filled_arc(d, 650, y+25, 25, 45, 90, 360);


    y=300;
    pts[0]=150; pts[1]=y;
    pts[2]=175; pts[3]=y+100;
    pts[4]=200; pts[5]=y;
    pts[6]=225; pts[7]=y+100;
    pts[8]=250; pts[9]=y;
    g2_poly_line(d, 5, pts);
    g2_pen(d, 19);
    g2_b_spline(d, 5, pts, 20);
    g2_pen(d, 1);
    
    pts[0]=300; pts[1]=y;
    pts[2]=350; pts[3]=y;
    pts[4]=375; pts[5]=y+50;
    pts[6]=325; pts[7]=y+90;
    pts[8]=275; pts[9]=y+50;
    g2_polygon(d, 5, pts);

    pts[0]=450; pts[1]=y;
    pts[2]=500; pts[3]=y;
    pts[4]=525; pts[5]=y+50;
    pts[6]=475; pts[7]=y+90;
    pts[8]=425; pts[9]=y+50;
    g2_filled_polygon(d, 5, pts);

    
    g2_image(d, 55., 50., 53, 62, penguin);
    g2_image(d, 75., 130., 53, 62, penguin);
    g2_pen(d, 1);
    
    g2_line(d, 225, 448, 200+19*25, 448);
    for(i=1;i<20;i++) {
        g2_pen(d,i+1);
	g2_set_line_width(d, i);
	g2_move(d, 200+i*25, 450);
	g2_line_to(d, 200+i*25, 550);
    }
    g2_pen(d,1);
    
    g2_set_line_width(d, 5);
    for(i=1;i<10;i++) {
	pts[0]=1*i;
	pts[1]=2*i;
	pts[2]=3*i;
	g2_set_dash(d, 3, pts);
	g2_line(d, 550, 300+i*8, 750, 350+i*8); 
    }

    g2_set_dash(d, 0, NULL);
    g2_set_line_width(d, 4);
    g2_arc(d, 740, 180, 25, 100, -45+15, -45-15);
    g2_arc(d, 740, 180, 20,  95, -45-15, -45+15);
    g2_filled_arc(d, 740, 180, 12, 50, -45+15, -45-15);

    g2_set_line_width(d, 1);
    g2_circle(d, 400, 400, 20);
    g2_ellipse(d, 400, 400, 25, 25);
    g2_arc(d, 400, 400, 30, 30, 0, 360);
    
    g2_flush(d);
    printf("\nDone.\n[Enter]\n");
    getchar();
    g2_close(d);
    return 0;
}





