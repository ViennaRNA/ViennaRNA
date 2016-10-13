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

#include <windows.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "g2_win32_P.h"

#include "resource.h"



int WINAPI About(HWND hDlg,UINT message, WPARAM wParam,LPARAM lParam)
{
	switch (message){
	case WM_INITDIALOG:
		return TRUE;

	case WM_COMMAND:
		if (wParam == IDOK)
			EndDialog(hDlg,wParam);
		break;
	}

	return FALSE;
}



LRESULT CALLBACK g2_WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	PAINTSTRUCT ps;
	HDC hDC;
	RECT Rect;
	struct 	g2_win32_STRUCT *pdp;

	pdp = (struct g2_win32_STRUCT *)GetWindowLong(hWnd, GWL_USERDATA);
	switch (message) { 

		case WM_PAINT:
			if (pdp == NULL) break;
			if (pdp->hBitmap == NULL) break;

//			printf("Received WM_PAINT\n");

			hDC = BeginPaint (hWnd, &ps);
			GetClientRect(hWnd,&Rect);
			BitBlt( hDC, Rect.left, Rect.top, Rect.right - Rect.left, 
			Rect.bottom - Rect.top, pdp->hMemDC, 0, 0, SRCCOPY );
			EndPaint (hWnd, &ps);
			return 0L;

		case WM_DESTROY:

//			printf("Received WM_DESTROY\n");
			g2_win32_Cleanup(0,pdp);
			ExitThread(0);
			return 0L;
			
		case WM_COMMAND:
			switch(LOWORD(wParam)){
			
			case ID_FILE_COPY:
				if (OpenClipboard(hWnd))
                   {
                    HBITMAP     hBitmap;
                   
                    EmptyClipboard();

                    if (pdp->hBitmap)
                        {
                        if (hBitmap = CopyImage(pdp->hBitmap,IMAGE_BITMAP,
												0,0,LR_COPYRETURNORG));
                            SetClipboardData(CF_BITMAP, hBitmap);
                        }
					CloseClipboard();
                    }
            return 0L;

			case ID_FILE_ABOUT:
				if(DialogBox(g2res_DLL,"ABOUTBOX",hWnd,(DLGPROC)About) == -1)
					errhandler("Failed to create Dialogbox",NULL);
				return 0L;

			case ID_FILE_CLOS:
				g2_win32_Delete(0,pdp);
				return 0L;

			default:
				return (DefWindowProc(hWnd, message, wParam, lParam));
			}

		default:
			return (DefWindowProc(hWnd, message, wParam, lParam));
	}
/* we should never get here */
return -1;
}

long WINAPI g2_StartThread(struct g2_win32_STRUCT *pdp)
{
RECT	Rect;
HWND hWnd;
MSG msg;
HDC hDC;
HMENU hmenu;
DWORD style;
RECT frame;

hmenu = NULL;
if (g2res_DLL != NULL)
	{
	hmenu = LoadMenu(g2res_DLL,"G2WIN32");
	if (hmenu == NULL) errhandler("Load menu failed",NULL);
	}

style = WS_POPUPWINDOW |WS_OVERLAPPED | WS_CAPTION |  WS_MINIMIZEBOX;
frame.left = 0;
frame.top = 0;
frame.right = pdp->nWidth;
frame.bottom = pdp->nHeight - ((hmenu==NULL)?GetSystemMetrics(SM_CYMENU):0);

AdjustWindowRect(&frame,style,1);

/* Save the instance handle in static variable, which will be used in  */
/* many subsequence calls from this application to Windows.            */

  /* Create a main window for this application instance.  */
pdp->hinst = GetModuleHandle(NULL);

hWnd = CreateWindow(
	 "g2Window",           // See RegisterClass() call.
	 pdp->title,		 // Text for window title bar.
	 style,
	 pdp->x, pdp->y,
	 frame.right - frame.left,   // width
	 frame.bottom - frame.top,   // height
	 NULL,                  // Overlapped windows have no parent.
	 hmenu,                 // Use the window class menu.
	 0,                     // This instance owns this window.
	 NULL                   // Pointer not needed.
  );
 
  // If window could not be created, return "failure" 
 
if (!hWnd)
	 {
	 errhandler("CreateWindow",NULL);
	 return (FALSE);   // return failure :((
	 }

#define WIDTH(x)	(x.right-x.left+1)	// Macro to get rect width
#define HEIGHT(x)	(x.bottom-x.top+1)	// Macro to get rect height

// How big is the window?
GetClientRect( hWnd, &Rect );

// Need a DC
hDC = GetDC( hWnd );
SetBkColor(hDC,RGB(255,255,255));
// Create a bitmap big enough to hold the window's image
pdp->hBitmap = CreateCompatibleBitmap( hDC, WIDTH(Rect), HEIGHT(Rect) );
// printf("memdc size: %d %d\n",WIDTH(Rect),HEIGHT(Rect));
// Create MemDC
pdp->hMemDC = CreateCompatibleDC(hDC);
SelectObject( pdp->hMemDC, pdp->hBitmap );
// clean up
ReleaseDC( hWnd, hDC );
#undef WIDTH
#undef HEIGHT



SetWindowLong(hWnd, GWL_USERDATA, (long)pdp); 
 
 pdp->hwndThreadWindow = hWnd;
 // Make the window visible; update its client area; and return "success" 
 ShowWindow(hWnd, SW_SHOWDEFAULT); // Show the window
 UpdateWindow(hWnd);     // Sends WM_PAINT message
 //printf("pdp->messageloop ->= 1;\n");
 pdp->messageloop = 1;
 //printf("pdp->messageloop = 1;\n");
 
 while (GetMessage(&msg, NULL, 0, 0)) {
        TranslateMessage(&msg);
		DispatchMessage(&msg);
    }

 
 return (TRUE);        // Returns success  :)
}



