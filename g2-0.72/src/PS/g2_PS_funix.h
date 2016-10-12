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
#ifndef _G2_PS_FUNIX_H
#define _G2_PS_FUNIX_H

#include "g2_PS_P.h"
#include "g2_physical_device.h"

const g2_funix_fun g2_PS_funix[] = {
    { g2_Delete,          g2_PS_delete },
    { g2_Ink,             g2_PS_ink    },
    { g2_Pen,             g2_PS_pen    },
    { g2_SetBackground,   g2_PS_set_background },
    { g2_ClearPalette,    g2_PS_clear_palette },
    { g2_SetLineWidth,    g2_PS_set_line_width },
    { g2_SetDash,         g2_PS_set_dash },
    { g2_SetFontSize,     g2_PS_set_font_size },
    { g2_Clear,           g2_PS_clear },		  
    { g2_Flush,           g2_PS_flush },
    { g2_Save,            g2_PS_flush },
    { g2_Plot,            g2_PS_plot },
    { g2_Line,            g2_PS_line },	  
    { g2_PolyLine,        g2_PS_poly_line },	  
    { g2_Polygon,         g2_PS_polygon },		  
    { g2_FilledPolygon,   g2_PS_filled_polygon },
    { g2_Rectangle,       g2_PS_rectangle },	  
    { g2_FilledRectangle, g2_PS_filled_rectangle },	  
    { g2_Triangle,        g2_PS_triangle },
    { g2_FilledTriangle,  g2_PS_filled_triangle },
    { g2_Arc,             g2_PS_arc },
    { g2_FilledArc,       g2_PS_filled_arc },
    { g2_Ellipse,         g2_PS_ellipse },
    { g2_FilledEllipse,   g2_PS_filled_ellipse },
    { g2_Circle,          NULL },
    { g2_FilledCircle,    NULL },
    { g2_String,          g2_PS_draw_string },
    { g2_FUNIX_NULL,      NULL } };


#endif /* _G2_PS_FUNIX_H */
