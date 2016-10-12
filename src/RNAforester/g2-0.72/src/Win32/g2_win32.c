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
#include <stdlib.h>
#include <math.h>

#include "g2.h"
#include "g2_device.h"
#include "g2_util.h"
#include "g2_config.h"

#include "g2_win32.h"
#include "g2_win32_P.h"
#include "g2_win32_funix.h"

#include "resource.h"

/* Global Definitions */
int g2_win32_registered = FALSE;
HINSTANCE g2res_DLL;	/* Instance of the resource DLL */

#define PDP ((struct g2_win32_STRUCT *)pdp)

#ifndef PI
#define PI 3.14159265358979323846
#endif /* PI */

#define sgn(x) (x>0?1:x?-1:0)

/* someday their might be a DLL version of g2 */
#ifdef G2DLL
BOOL WINAPI DllMain( HANDLE hModule, DWORD fdwreason,  LPVOID lpReserved )
{
    switch(fdwreason) {
    case DLL_PROCESS_ATTACH:
    // The DLL is being mapped into process's address space
    //  Do any required initialization on a per application basis, return FALSE if failed
    MessageBox(NULL, "DLL Process Attach", "DLL Message 1", MB_OK);    
	break;
    case DLL_THREAD_ATTACH:
    // A thread is created. Do any required initialization on a per thread basis
    MessageBox(NULL, "DLL Thread Attach", "DLL Message 1", MB_OK);    
    break;
    case DLL_THREAD_DETACH:
    // Thread exits with  cleanup
    MessageBox(NULL, "DLL Thread Detach", "DLL Message 1", MB_OK);    
    break;
    case DLL_PROCESS_DETACH:
    // The DLL unmapped from process's address space. Do necessary cleanup
    MessageBox(NULL, "DLL Process Detach", "DLL Message 1", MB_OK);    
    break;
	default:
	MessageBox(NULL, "DLL default", "DLL Message 1", MB_OK);    

    }
    return TRUE;
}
#endif


g2_win32_SetPen(int pid, void *pdp)
	{
	HGDIOBJ oldpen;
    LOGBRUSH    logBrush ;
	logBrush.lbStyle = PS_SOLID;
	logBrush.lbColor = PDP->Inks[PDP->Pen];
    logBrush.lbHatch = 0 ;

	oldpen = PDP->hPen;
	PDP->hPen = ExtCreatePen (logBrush.lbStyle | PS_GEOMETRIC | 
                         PS_ENDCAP_FLAT | PS_JOIN_BEVEL | (PDP->PenStyle > 0)*PS_USERSTYLE,
                         PDP->PenWidth, &logBrush,
                         PDP->PenStyle , PDP->PenDash) ;
	if (PDP->hPen != NULL) 
		{
		SelectObject(PDP->hMemDC,PDP->hPen);
/*		if (PDP->type == g2_win32)*/
		DeleteObject(oldpen);
		}
	else
		{
		errhandler("Pen",NULL);
		PDP->hPen = oldpen;
		}

	oldpen = PDP->hNullPen;
	PDP->hNullPen = CreatePen(PS_SOLID,1,PDP->Inks[PDP->Pen]);
	if (PDP->hNullPen != NULL) 
		{
/*		if (PDP->type == g2_win32)*/
		DeleteObject(oldpen);
		}
	else
		{
		errhandler("Pen",NULL);
		PDP->hNullPen = oldpen;
		}

	return 0;
	}


int g2_win32_Cleanup(int pid, void *pdp)
	{
	struct g2_win32_STRUCT *thispdp;

	thispdp = pdp;
	g2_win32_ClearPalette(pid,pdp);
	if (PDP->hBrush != NULL) DeleteObject(PDP->hBrush);
	if (PDP->hPen != NULL) DeleteObject(PDP->hPen);
	if (PDP->hNullPen != NULL) DeleteObject(PDP->hPen);
	if (PDP->hFont != NULL) DeleteObject(PDP->hFont);
	if (PDP->hBitmap != NULL) DeleteObject(PDP->hBitmap);
	if (PDP->hMemDC != NULL) DeleteDC(PDP->hMemDC);
	if (PDP->PenDash != NULL) free(PDP->PenDash);
	free(thispdp);
	return 0;
	}



int g2_win32_Delete(int pid, void *pdp)
	{
	switch(PDP->type)
		{
		case g2_win32:
			SendMessage(PDP->hwndThreadWindow,WM_CLOSE,(WPARAM)NULL,(LPARAM)NULL);
			break;
		case g2_wmf32:
			{
			CloseEnhMetaFile(PDP->hMemDC);
			g2_win32_Cleanup(pid,pdp);
			break;
			}
		}
	return 0;
    }


int g2_win32_Clear(int pid, void *pdp)
	{
	int OldPen;

	OldPen = PDP->Pen;
	g2_win32_Pen(pid,pdp,PDP->BkColor);
	g2_win32_FilledRectangle(pid,pdp,0,0,PDP->nWidth,PDP->nHeight);
	g2_win32_Pen(pid,pdp,OldPen);
	return 0;
    }

int g2_win32_Flush(int pid, void *pdp)
	{
	InvalidateRect(PDP->hwndThreadWindow, (RECT *)NULL, TRUE); 
	return 0;
    }

int g2_win32_Pen(int pid, void *pdp, int color)
	{
	struct tagLOGBRUSH logbrush;
	HGDIOBJ oldbrush;
	
	if(color>=PDP->NoOfInks || color<0)
		{
		fprintf(stderr,"g2_WIN32: Ink %d not defined\n",color);
		return -1;
		}

	PDP->Pen = color;
	PDP->PenColor = PDP->Inks[color];
	g2_win32_SetPen(pid,pdp);

	logbrush.lbStyle = BS_SOLID;
	logbrush.lbColor = PDP->PenColor;
	oldbrush = PDP->hBrush;
	PDP->hBrush = CreateBrushIndirect(&logbrush);
	if (PDP->hBrush == NULL) 
		{
		errhandler("Pen (CreateBrush)",NULL);
		PDP->hBrush = oldbrush;
		}
	else
/*		if (PDP->type == g2_win32)*/
		DeleteObject(oldbrush);
	return 0;
    }

int g2_win32_Ink(int pid, void *pdp, double red, double green, double blue)
	{
	BYTE rc,gc,bc;

	rc = (BYTE)((int)(red*255));
	gc = (BYTE)((int)(green*255));
	bc = (BYTE)((int)(blue*255));

	PDP->NoOfInks++;

	if(PDP->Inks==NULL)
		PDP->Inks=(COLORREF *)malloc(PDP->NoOfInks*sizeof(COLORREF));
	else
		PDP->Inks=(COLORREF *)realloc((void *)PDP->Inks,PDP->NoOfInks*sizeof(COLORREF));

	if(PDP->Inks==NULL) 
		{
		fputs("g2: not enough memory\n",stderr);
		return -1;
		}

	PDP->Inks[PDP->NoOfInks-1]=RGB(rc,gc,bc);
	return PDP->NoOfInks-1;
	}	


int g2_win32_ClearPalette(int pid, void *pdp)
	{
	if (PDP->Inks != NULL)
		free(PDP->Inks);
	PDP->Inks = NULL;
	PDP->NoOfInks = 0;
	return 0;
    }


int g2_win32_SetBackground(int pid, void *pdp, int color)
	{
	PDP->BkColor = color;
	SetBkColor(PDP->hMemDC,PDP->BkColor);
	return 0;
    }

int g2_win32_SetLineWidth(int pid, void *pdp, int w)
	{
	PDP->PenWidth = w;
	g2_win32_SetPen(pid,pdp);
	return 0;
    }

int g2_win32_SetDash(int pid, void *pdp, int n, int *data)
	{
	if (PDP->PenDash != NULL)
		free(PDP->PenDash);
	PDP->PenDash = NULL;
	PDP->PenStyle = n;
	if (n > 0) 
		{
		int i;
		PDP->PenDash = (DWORD *)malloc(n*sizeof(DWORD));
		for (i=0;i<n;i++)
			PDP->PenDash[i] = data[i];
		}
	g2_win32_SetPen(pid,pdp);
	return 0;
	}


int g2_win32_SetFontSize(int pid, void *pdp, int size)
	{
	//static LOGFONT lf = {10,0,0,0,0,0,0,0,0,0,0,0,0,"Arial\0"};
    HGDIOBJ oldfont;
   
	oldfont = PDP->hFont;
	//lf.lfHeight = size;
	//PDP->hFont = CreateFontIndirect(&lf);
	PDP->hFont = CreateFont(-size, 0, 0, 0, FW_NORMAL, 0, 0, 0, 0, OUT_TT_ONLY_PRECIS , 0, PROOF_QUALITY,0, "Times New Roman\0");
	if (PDP->hFont == NULL) 
		{
		errhandler("Font (CreateFont)",NULL);
		PDP->hFont = oldfont;
		}
	else
		{
		SelectObject(PDP->hMemDC,PDP->hFont);
		if (oldfont != NULL)// && PDP->type == g2_win32)
			DeleteObject(oldfont);
		}
	return 0;
    }

int g2_win32_Plot(int pid, void *pdp, int x, int y)
	{
	return SetPixel(PDP->hMemDC,x,y,PDP->PenColor);
    }

int g2_win32_Line(int pid, void *pdp, int x1, int y1, int x2, int y2)
	{
	MoveToEx(PDP->hMemDC,x1,y1,NULL);
	LineTo(PDP->hMemDC,x2,y2);
	SetPixel(PDP->hMemDC,x1,y1,PDP->PenColor);
	SetPixel(PDP->hMemDC,x2,y2,PDP->PenColor);
	// specifically draw end points since windows does not include one endpoint
	return 0;
    }

int g2_win32_PolyLine(int pid, void *pdp, int N, int *points)
	{
	POINT *PointList;
	int i;

	PointList = (POINT *)malloc(N*sizeof(POINT));
	if (PointList == NULL)
		{
		fprintf(stderr,"g2_win32: not enough memory !\n");
		return -1;
		}
	for (i=0;i<N;i++)
		{
		PointList[i].x = points[2*i];
		PointList[i].y = points[2*i+1];
		}
	Polyline(PDP->hMemDC,PointList,N);
	free(PointList);
	return 0;
    }



int g2_win32_Rectangle(int pid, void *pdp, int x, int y, int x2, int y2)
	{
	SelectObject(PDP->hMemDC,GetStockObject(NULL_BRUSH));
	Rectangle(PDP->hMemDC,x,y,x2+1,y2+1); // add one since windows excludes lower right point
	return 0;
    }

int g2_win32_FilledRectangle(int pid, void *pdp, int x1, int y1, int x2, int y2)
	{
	SelectObject(PDP->hMemDC,PDP->hBrush);
	SelectObject(PDP->hMemDC,PDP->hNullPen);
	return Rectangle(PDP->hMemDC,x1,y1,x2+1,y2+1); // add one since windows excludes lower right point
	SelectObject(PDP->hMemDC,PDP->hPen);
	return 0;
	}

int g2_win32_Polygon(int pid, void *pdp, int N, int *points)
	{
	POINT *PointList;
	int i;
	PointList = (POINT *)malloc(N*sizeof(POINT));
	if (PointList == NULL)
		{
		fprintf(stderr,"g2: not enough memory !\n");
		return -1;
		}
	for (i=0;i<N;i++)
		{
		PointList[i].x = points[2*i];
		PointList[i].y = points[2*i+1];
		}
	SelectObject(PDP->hMemDC,GetStockObject(NULL_BRUSH));
	Polygon(PDP->hMemDC,PointList,N);
	free(PointList);
	return 0;
    }

int g2_win32_FilledPolygon(int pid, void *pdp, int N, int *points)
	{
	POINT *PointList;
	int i;
	PointList = (POINT *)malloc(N*sizeof(POINT));
	if (PointList == NULL)
		{
		fprintf(stderr,"g2: not enough memory !\n");
		return -1;
		}
	for (i=0;i<N;i++)
		{
		PointList[i].x = points[2*i];
		PointList[i].y = points[2*i+1];
		}
	SelectObject(PDP->hMemDC,PDP->hBrush);
	SelectObject(PDP->hMemDC,PDP->hNullPen);
	Polygon(PDP->hMemDC,PointList,N);
	SelectObject(PDP->hMemDC,PDP->hPen);
	free(PointList);
	return 0;
    }

int g2_win32_Ellipse(int pid, void *pdp, int x, int y, int r1, int r2)
	{
	SelectObject(PDP->hMemDC,GetStockObject(NULL_BRUSH));
	return Ellipse(PDP->hMemDC,x-r1,y-r2,x+r1+1,y+r2+1); // add one since windows is end exclusive
    }

int g2_win32_FilledEllipse(int pid, void *pdp, int x, int y, int r1, int r2)
	{
	SelectObject(PDP->hMemDC,PDP->hBrush);
	SelectObject(PDP->hMemDC,PDP->hNullPen);
	return Ellipse(PDP->hMemDC,x-r1,y-r2,x+r1+1,y+r2+1); // add one since windows is end exclusive
	SelectObject(PDP->hMemDC,PDP->hPen);
	return 0;
    }

int g2_win32_Arc(int pid, void *pdp, int x, int y, int r1, int r2, double a1, double a2)
	{
	a1 *= PI/180.;
	a2 *= PI/180.;
	SelectObject(PDP->hMemDC,GetStockObject(NULL_BRUSH));
	if (PDP->type == g2_win32)
		return Arc(PDP->hMemDC,x-r1,y-r2,x+r1,y+r2,dtoi(x+r1*cos(a1)),dtoi(y-r2*sin(a1)),dtoi(x+r1*cos(a2)),dtoi(y-r2*sin(a2)));
	else
		return Arc(PDP->hMemDC,x-r1,y-r2,x+r1,y+r2,dtoi(x+r1*cos(a2)),dtoi(y-r2*sin(a2)),dtoi(x+r1*cos(a1)),dtoi(y-r2*sin(a1)));
    }

int g2_win32_FilledArc(int pid, void *pdp, int x, int y, int r1, int r2, double a1, double a2)
	{
	a1 *= PI/180.;
	a2 *= PI/180.;
	SelectObject(PDP->hMemDC,PDP->hBrush);
	SelectObject(PDP->hMemDC,PDP->hNullPen);
	if (PDP->type == g2_win32)
		Pie(PDP->hMemDC,x-r1,y-r2,x+r1,y+r2,dtoi(x+r1*cos(a1)),dtoi(y-r2*sin(a1)),dtoi(x+r1*cos(a2)),dtoi(y-r2*sin(a2)));
	else
		Pie(PDP->hMemDC,x-r1,y-r2,x+r1,y+r2,dtoi(x+r1*cos(a2)),dtoi(y-r2*sin(a2)),dtoi(x+r1*cos(a1)),dtoi(y-r2*sin(a1)));
	SelectObject(PDP->hMemDC,PDP->hPen);
	return 0;
    }

int g2_win32_DrawString(int pid, void *pdp, int x, int y, const char *text)
	{
	SetTextColor(PDP->hMemDC,PDP->PenColor);
	SetBkMode(PDP->hMemDC,TRANSPARENT);
	return TextOut(PDP->hMemDC,x,y,text,strlen(text));
    }

int g2_win32_QueryPointer(int pid, void *pdp, int *x, int *y, unsigned int *button)
//
// Thanks to input by Martin stéphane
//
	{ 
	POINT point;

	GetCursorPos(&point);

	ScreenToClient(PDP->hwndThreadWindow,&point);

	*y=point.y;
	*x=point.x;
	*button=0;

	if (PDP->hwndThreadWindow != GetForegroundWindow())
		return; // return if our window does not have the focus

	if (GetKeyState(VK_LBUTTON)<0) 
		*button=*button+256;

	if (GetKeyState(VK_MBUTTON)<0) 
		*button=*button+512;

	if (GetKeyState(VK_RBUTTON)<0) 
		*button=*button+1024;

	return 0;
	}

void errhandler(LPSTR errtxt,HWND hwnd)
{
	LPVOID lpMessageBuffer;
	char szError[255];
	DWORD LastError;

	LastError = GetLastError();

	FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER |FORMAT_MESSAGE_FROM_SYSTEM,
	NULL,LastError,MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), //The user default language
	(LPTSTR) &lpMessageBuffer,0,NULL );

	sprintf(szError,"Error in %s\n Error code %d : %s\n" ,errtxt,LastError,lpMessageBuffer);
	MessageBox(hwnd,szError,"error",MB_OK);

	LocalFree( lpMessageBuffer );
}


int InitApplication()
{
  WNDCLASS  wc;

  // Fill in window class structure with parameters that describe the
  // main window.

    g2res_DLL = LoadLibrary("g2res.dll");

  if (g2res_DLL == NULL) 
		printf("Warning: Could not load g2 resource DLL\n Menu and Icon are disabled\n");
  

  wc.style = CS_HREDRAW | CS_VREDRAW; // Class style(s).
  wc.lpfnWndProc = g2_WndProc; // Function to retrieve messages for windows of this class.
  wc.cbClsExtra = 0;	// No per-class extra data.
  wc.cbWndExtra = 0;	// No per-window extra data.
  wc.hInstance = 0;		// Application that owns the class.
  wc.hIcon = NULL;
  if (g2res_DLL != NULL)
	  {
		wc.hIcon = LoadIcon(g2res_DLL, MAKEINTRESOURCE(IDI_ICON1));
		if (wc.hIcon == NULL) 
	  		errhandler("Error in LoadIcon",NULL);
	  }
  wc.hCursor = LoadCursor(NULL, IDC_ARROW);
  wc.hbrBackground = NULL;
  wc.lpszMenuName = "G2WIN32";	
  wc.lpszClassName = "g2Window";	// Name used in call to CreateWindow.

  /* Register the window class and return success/failure code. */


	if(!RegisterClass(&wc)) 
	  {
	  errhandler("RegisterClass",NULL);
	  return FALSE;
	  }
  g2_win32_registered = TRUE;
  return TRUE;
}


/**
 * \ingroup physdev
 * \defgroup win32 MS Windows
 */

/**
 *
 * Create a Windows device.
 *
 * \param width window width
 * \param height window height
 * \param title window title
 * \param type window type, see ::g2_win32_type
 *
 * \return physical device id
 *
 * \ingroup win32
 */
int  g2_open_win32(int width, int height, const char *title, int type)
	{
	int pid=0,vid;
	long ThreadID;

	g2_win32_STRUCT *pdp;
	
	pdp = (g2_win32_STRUCT *)malloc(sizeof(g2_win32_STRUCT));

	PDP->type = type;
	PDP->NoOfInks = 0;
	PDP->Inks = NULL;
	PDP->PenStyle = 0;
	PDP->PenWidth = 1;
	PDP->PenColor = RGB(1,1,1);
	PDP->PenDash = NULL;
	PDP->nWidth = width;
	PDP->nHeight = height;
	PDP->messageloop = 0;
	PDP->hFont = NULL;
	
	switch(type) {
		case g2_win32:
			{
			if(g2_win32_registered == FALSE)
				InitApplication();
		
			PDP->x = 1;
			PDP->y = 1;
			if (title != NULL)
 				PDP->title = title;
			else
				PDP->title = "g2 Window";
			PDP->hThread = CreateThread(NULL, 0,
                                         (LPTHREAD_START_ROUTINE)g2_StartThread,
                                         (LPVOID)(pdp),
                                         0,
           								 (LPDWORD) &ThreadID );

			if (PDP->hThread == NULL) 
				fprintf(stderr,"g2_win32: Thread could not be started\n");

		    SetThreadPriority(PDP->hThread,THREAD_PRIORITY_ABOVE_NORMAL);
			//Wait till window is created by Thread
			while( PDP->messageloop == 0)
				Sleep(10);
			break;
			}
		case g2_wmf32 :
			{
            DWORD dwInchesX;
            DWORD dwInchesY;
            DWORD dwDPI = 72;    
		    RECT   Rect = { 0, 0, 0, 0 };
			TCHAR   szDesc[] = "Created by g2\0\0";
			HDC     hScreenDC;
			float   PixelsX, PixelsY, MMX, MMY;
			
			dwInchesX = PDP->nWidth/72;
            dwInchesY = PDP->nHeight/72;
			// dwInchesX x dwInchesY in .01mm units
			SetRect( &Rect, 0, 0,dwInchesX*2540, dwInchesY*2540 );
 
			// Get a Reference DC
			hScreenDC = GetDC( NULL );
 
			// Get the physical characteristics of the reference DC
			PixelsX = (float)GetDeviceCaps( hScreenDC, HORZRES );
			PixelsY = (float)GetDeviceCaps( hScreenDC, VERTRES );
			MMX = (float)GetDeviceCaps( hScreenDC, HORZSIZE );
			MMY = (float)GetDeviceCaps( hScreenDC, VERTSIZE );
 
			// Create the Metafile
		    PDP->hMemDC = CreateEnhMetaFile(hScreenDC, title, &Rect, szDesc);
		    //tstDC = CreateEnhMetaFile(hScreenDC, "test.emf", &Rect, szDesc);

			// Release the reference DC
			ReleaseDC( NULL, hScreenDC );
			// Anisotropic mapping mode
			SetMapMode( PDP->hMemDC, MM_ANISOTROPIC );
			// Set the Windows extent
			SetWindowExtEx( PDP->hMemDC, dwInchesX*dwDPI, dwInchesY*dwDPI, NULL );
 
			// Set the viewport extent to reflect
			// dwInchesX" x dwInchesY" in device units
			SetViewportExtEx( PDP->hMemDC,
                      (int)((float)dwInchesX*25.4f*PixelsX/MMX),
                      (int)((float)dwInchesY*25.4f*PixelsY/MMY),
                      NULL );
//			printf("viewport: %d %d\n",(int)((float)dwInchesX*25.4f*PixelsX/MMX),
//										(int)((float)dwInchesY*25.4f*PixelsY/MMY));
			// create the device context
//			PDP->hMemDC = CreateMetaFile(NULL) ;
       
//			PDP->hMemDC =  CreateEnhMetaFile( (HDC)NULL, title, &WmfRect, "Created by g2"); 
//			PDP->hMemDC = CreateMetaFile(title);
//			SetMapMode(PDP->hMemDC,MM_LOMETRIC );
//			SetWindowExtEx(PDP->hMemDC,width,height,NULL);
//			SetViewportExtEx(PDP->hMemDC,width*1000,height*1000,NULL);
			if (PDP->hMemDC == NULL) errhandler("Could not Create Metafile !\n",NULL);
			SetArcDirection(PDP->hMemDC,AD_CLOCKWISE);
			break;
			}
		default:
			return height;
		}
		SetTextAlign(PDP->hMemDC,TA_BOTTOM | TA_LEFT);
		vid = g2_register_physical_device(pid, pdp,
				      g2_IntCoor, g2_win32_funix,
				      1.0, -1.0,
				      0.0, height-1);

		g2_allocate_basic_colors(vid);
		g2_pen(vid,1);
		g2_set_background(vid, 0);
		g2_set_font_size(vid, 10);
		g2_clear(vid);

    return vid;
	}

#undef PDP 
