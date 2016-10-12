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
#ifndef _G2_GD_TOKEN_FUN_H
#define _G2_GD_TOKEN_FUN_H
		
	
#include "g2_virtual_device.h"
#include "g2_gd_P.h"

const g2_funix_fun g2_gd_funix[] = 
		{{g2_Delete,g2_gd_Delete},
		{g2_Ink, g2_gd_Ink},
		{g2_Pen,g2_gd_Pen},
		{g2_SetBackground,g2_gd_SetBackground},
		{g2_ClearPalette,g2_gd_ClearPalette},	 
		{g2_SetLineWidth,g2_gd_SetLineWidth},
		 /*	{g2_SetDash,g2_gd_SetDash},	*/
	 	{g2_SetFontSize,g2_gd_SetFontSize},
		{g2_Clear,g2_gd_Clear},
		{g2_Flush,g2_gd_Flush},
		{g2_Plot,g2_gd_Plot},
		{g2_Line,g2_gd_Line},
		 /*	{g2_PolyLine,g2_gd_PolyLine}, */
		{g2_Rectangle,g2_gd_Rectangle},
		{g2_FilledRectangle,g2_gd_FilledRectangle},
		 /*	{g2_Polygon,g2_gd_Polygon}, */
		{g2_FilledPolygon,g2_gd_FilledPolygon},
		 /*	{g2_Circle,g2_gd_Circle}, */
		 /*	{g2_FilledCircle,g2_gd_FilledCircle}, */
		 /*	{g2_Ellipse,g2_gd_Ellipse}, */
		 /*	{g2_FilledEllipse,g2_gd_FilledEllipse}, */
		 /*     {g2_Arc,g2_gd_Arc}, */
		 /* {g2_FilledArc,g2_gd_FilledArc}, */
		{g2_String,g2_gd_DrawString},
		{g2_FUNIX_NULL, NULL}};


#endif /* _G2_GD_TOKEN_FUN_H */
