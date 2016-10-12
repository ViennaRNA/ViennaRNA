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
#ifndef _G2_X11_FUNIX_H
#define _G2_X11_FUNIX_H

#include "g2_X11_P.h"
#include "g2_physical_device.h"


const g2_funix_fun g2_X11_funix[] = {
    { g2_Delete,          g2_X11_delete },
    { g2_Ink,             g2_X11_ink    },
    { g2_Pen,             g2_X11_pen    },
    { g2_SetBackground,   g2_X11_set_background },
    { g2_ClearPalette,    g2_X11_clear_palette },
    { g2_SetLineWidth,    g2_X11_set_line_width },
    { g2_SetDash,         g2_X11_set_dash },
    { g2_SetFontSize,     g2_X11_set_font_size },
    { g2_Clear,           g2_X11_clear },		  
    { g2_Flush,           g2_X11_flush },
    { g2_Save,            g2_X11_flush },
    { g2_Plot,            g2_X11_plot },
    { g2_Line,            g2_X11_line },	  
    { g2_PolyLine,        g2_X11_poly_line },	  
    { g2_Polygon,         g2_X11_polygon },		  
    { g2_FilledPolygon,   g2_X11_filled_polygon },
    { g2_Rectangle,       g2_X11_rectangle },	  
    { g2_FilledRectangle, g2_X11_filled_rectangle },	  
    { g2_Triangle,        g2_X11_triangle },
    { g2_FilledTriangle,  g2_X11_filled_triangle },
    { g2_Arc,             g2_X11_arc },
    { g2_FilledArc,       g2_X11_filled_arc },
    { g2_Ellipse,         g2_X11_ellipse },
    { g2_FilledEllipse,   g2_X11_filled_ellipse },
    { g2_Circle,          NULL },
    { g2_FilledCircle,    NULL },
    { g2_String,          g2_X11_draw_string },
    { g2_Image,           NULL },
    { g2_QueryPointer,    g2_X11_query_pointer },
    { g2_GetPDHandles,    g2_X11_get_pd_handles },
    { g2_FUNIX_NULL,      NULL } };


#endif /* _G2_X11_FUNIX_H */
