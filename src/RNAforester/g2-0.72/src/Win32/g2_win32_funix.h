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
#ifndef _G2_WIN32_FUNIX_H
#define _G2_WIN32_FUNIX_H

#include "g2_win32_P.h"
#include "g2_physical_device.h"

const g2_funix_fun g2_win32_funix[] = {
    { g2_Delete,          g2_win32_Delete },
    { g2_Ink,             g2_win32_Ink    },
    { g2_Pen,             g2_win32_Pen    },
    { g2_SetBackground,   g2_win32_SetBackground },
    { g2_ClearPalette,    g2_win32_ClearPalette },
    { g2_SetLineWidth,    g2_win32_SetLineWidth },
    { g2_SetDash,         g2_win32_SetDash },
    { g2_SetFontSize,     g2_win32_SetFontSize },
    { g2_Clear,           g2_win32_Clear },		  
    { g2_Flush,           g2_win32_Flush },
    { g2_Save,            g2_win32_Flush },
    { g2_Plot,            g2_win32_Plot },
    { g2_Line,            g2_win32_Line },	  
    { g2_PolyLine,        g2_win32_PolyLine },	  
    { g2_Polygon,         g2_win32_Polygon },		  
    { g2_FilledPolygon,   g2_win32_FilledPolygon },
    { g2_Rectangle,       g2_win32_Rectangle },	  
    { g2_FilledRectangle, g2_win32_FilledRectangle },	  
    { g2_Arc,             g2_win32_Arc },
    { g2_FilledArc,       g2_win32_FilledArc },
    { g2_Ellipse,         g2_win32_Ellipse },
    { g2_FilledEllipse,   g2_win32_FilledEllipse },
    { g2_Circle,          NULL },
    { g2_FilledCircle,    NULL },
    { g2_String,          g2_win32_DrawString },
	{ g2_QueryPointer,    g2_win32_QueryPointer },
    { g2_FUNIX_NULL,      NULL } };


#endif /* _G2_WIN32_FUNIX_H */
