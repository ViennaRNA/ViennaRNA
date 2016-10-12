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
#ifndef _G2_WIN32_P_H
#define _G2_WIN32_P_H

#include "windows.h"
#include <stdio.h>

extern HMENU hmenu;
extern HANDLE ghModule;
extern HINSTANCE g2res_DLL;

typedef struct g2_win32_STRUCT {
	HANDLE hThread;
	HWND   hwndThreadWindow;
	HANDLE hinst;
	HBITMAP hBitmap;
	HDC hMemDC;
	HPEN hPen;
	HPEN hNullPen;
	HBRUSH hBrush;
	HFONT hFont;
	int nWidth;
	int nHeight;
	int x;
	int y;
	int NoOfInks;
	int PenWidth;
	int PenStyle;
	DWORD *PenDash;
	int Pen;
	int BkColor;
	COLORREF PenColor;
	COLORREF *Inks;
	char *title;
	int type;
	int messageloop;
	} g2_win32_STRUCT;


BOOL SaveBitmapAs(HWND hWnd,struct XPGTHREADINFO *pThreadInfo);
LRESULT CALLBACK g2_WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam);
long WINAPI g2_StartThread(struct g2_win32_STRUCT *pdp);
void errhandler(LPSTR errtxt,HWND hwnd);


int g2_win32_init_win32(int pid, void *pdp, int vid, int width,int height);
int g2_win32_init_win32X(int pid, void *pdp,
		     int width, int height,
		     int xposition, int yposition,
		     char *windowname, char *iconname,
		     char *icondata, int iconwidth, int iconheight);
int g2_win32_Delete(int pid, void *pdp);
int g2_win32_Clear(int pid, void *pdp);
int g2_win32_Flush(int pid, void *pdp);
int g2_win32_Ink(int pid, void *pdp,
	       double red, double green, double blue);
int g2_win32_ClearPalette(int pid, void *pdp);
int g2_win32_ResetPalette(int pid, void *pdp);
int g2_win32_SetBackground(int pid, void *pdp, int color);
int g2_win32_Pen(int pid, void *pdp, int color);
int g2_win32_Paper(int pid, void *pdp, int color);
int g2_win32_SetLineWidth(int pid, void *pdp, int w);
int g2_win32_SetDash(int pid, void *pdp, int n, int *data);
int g2_win32_SetFontSize(int pid, void *pdp, int size);
int g2_win32_Plot(int pid, void *pdp, int x, int y);
int g2_win32_Line(int pid, void *pdp, int x1, int y1, int x2, int y2);
int g2_win32_PolyLine(int pid, void *pdp, int N, int *p);
int g2_win32_Polygon(int pid, void *pdp, int N, int *p);
int g2_win32_FilledPolygon(int pid, void *pdp, int N, int *p);
int g2_win32_Triangle(int pid, void *pdp,
		    int x1, int y1,
		    int x2, int y2,
		    int x3, int y3);
int g2_win32_FilledTriangle(int pid, void *pdp, int x1, int y1,
			   int x2, int y2,
			   int x3, int y3);
int g2_win32_Rectangle(int pid, void *pdp, int x1, int y1, int x2, int y2);
int g2_win32_FilledRectangle(int pid, void *pdp,
			    int x1, int y1, int x2, int y2);
int g2_win32_Circle(int pid, void *pdp, int x, int y, int r);
int g2_win32_FilledCircle(int pid, void *pdp, int x, int y, int r);
int g2_win32_Arc(int pid, void *pdp, int x, int y,
	       int r1, int r2, double a1, double a2);
int g2_win32_FilledArc(int pid, void *pdp, int x, int y,
		      int r1, int r2, double a1, double a2);
int g2_win32_Ellipse(int pid, void *pdp, int x, int y, int r1, int r2);
int g2_win32_FilledEllipse(int pid, void *pdp, int x, int y, int r1, int r2);
int g2_win32_DrawString(int pid, void *pdp, int x, int y, const char *text);
int g2_win32_QueryPointer(int pid, void *pdp, int *x, int *y, unsigned int *button);

int g2_win32_AllocateBasicColors(int pid, void *pdp);

int g2_win32_Cleanup(int pid, void *pdp);

#endif /* _G2_WIN32_P_H */
