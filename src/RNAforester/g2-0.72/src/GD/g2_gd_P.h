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
#ifndef _G2_GD_P_H
#define _G2_GD_P_H

#include <gd.h>
#include <gdfontt.h>
#include <gdfonts.h>
#include <gdfontmb.h>
#include <gdfontl.h>
#include <gdfontg.h>
#include "g2_config.h"
#include "g2_gd.h"

 static unsigned short g2_Basic_Colors[32][3] = {
 {65535,65535,65535},{0,0,0},
 {44461,55512,59110},{0,65535,65535},{0,0,65535},
 {0,0,32896},{61166,33410,61166},
 {65535,0,65535},{41120,8224,61680},{65535,49344,52171},
 {61680,59110,35980},{32639,65535,54484},
 {39064,64507,39064},{0,65535,0},{12850,52685,12850},
 {39578,52685,12850},{44461,65535,12079},{65535,65535,0},
 {65535,55255,0},{61166,56797,33410},{56026,42405,8224},{47288,34438,2827},
 {42405,10794,10794},{53970,26985,7710},
 {52685,34181,16191},{65535,32639,20560},
 {61680,32896,32896},{64250,32896,29298},
 {65535,35980,0},{65535,42405,0},{65535,17733,0},
 {65535,0,0}};


typedef struct g2_gd_STRUCT {
	int height;
	int width;
	FILE *f;
	gdImagePtr im;
	gdImagePtr brush;
	int LineWidth;
	int NoOfInks;
	int Inks[256];
	int CurCol;
	int OldCol;
	int BackCol;
	gdFontPtr FontSize;
	enum g2_gd_type gd_type;
	} g2_gd_STRUCT;



int g2_gd_init_gd(int pid, void *pdp, int vid, int width,int height);
int g2_gd_init_gdX(int pid, void *pdp,
		     int width, int height,
		     int xposition, int yposition,
		     char *windowname, char *iconname,
		     char *icondata, int iconwidth, int iconheight);
int g2_gd_Delete(int pid, void *pdp);
int g2_gd_Clear(int pid, void *pdp);
int g2_gd_Flush(int pid, void *pdp);
int g2_gd_Ink(int pid, void *pdp,
	       double red, double green, double blue);
int g2_gd_ClearPalette(int pid, void *pdp);
int g2_gd_ResetPalette(int pid, void *pdp);
int g2_gd_SetBackground(int pid, void *pdp, int color);
int g2_gd_Pen(int pid, void *pdp, int color);
int g2_gd_Paper(int pid, void *pdp, int color);
int g2_gd_SetLineWidth(int pid, void *pdp, int w);
int g2_gd_SetDash(int pid, void *pdp, int n, char *data);
int g2_gd_SetFontSize(int pid, void *pdp, int size);
int g2_gd_Plot(int pid, void *pdp, int x, int y);
int g2_gd_Line(int pid, void *pdp, int x1, int y1, int x2, int y2);
int g2_gd_PolyLine(int pid, void *pdp, int N, int *p);
int g2_gd_Polygon(int pid, void *pdp, int N, int *p);
int g2_gd_FilledPolygon(int pid, void *pdp, int N, int *p);
int g2_gd_Triangle(int pid, void *pdp,
		    int x1, int y1,
		    int x2, int y2,
		    int x3, int y3);
int g2_gd_FilledTriangle(int pid, void *pdp, int x1, int y1,
			   int x2, int y2,
			   int x3, int y3);
int g2_gd_Rectangle(int pid, void *pdp, int x1, int y1, int x2, int y2);
int g2_gd_FilledRectangle(int pid, void *pdp,
			    int x1, int y1, int x2, int y2);
int g2_gd_Circle(int pid, void *pdp, int x, int y, int r);
int g2_gd_FilledCircle(int pid, void *pdp, int x, int y, int r);
int g2_gd_Arc(int pid, void *pdp, int x, int y,
	       int r1, int r2, double a1, double a2);
int g2_gd_FilledArc(int pid, void *pdp, int x, int y,
		      int r1, int r2, double a1, double a2);
int g2_gd_Ellipse(int pid, void *pdp, int x, int y, int r1, int r2);
int g2_gd_FilledEllipse(int pid, void *pdp, int x, int y, int r1, int r2);
int g2_gd_DrawString(int pid, void *pdp, int x, int y, const char *text);

int g2_gd_AllocateBasicColors(int pid, void *pdp);

#endif /* _G2_GD_P_H */
