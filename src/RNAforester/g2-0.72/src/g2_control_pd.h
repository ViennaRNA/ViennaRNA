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
#ifndef _G2_CONTROL_PD_H
#define _G2_CONTROL_PD_H

#include "g2_physical_device.h"

void g2_flush_pd(g2_physical_device *pd);
void g2_save_pd(g2_physical_device *pd);
void g2_clear_pd(g2_physical_device *pd);
void g2_pen_pd(g2_physical_device *pd, int color);
void g2_set_background_pd(g2_physical_device *pd, int color);
int  g2_ink_pd(g2_physical_device *pd, double red, double green, double blue);
void g2_clear_palette_pd(g2_physical_device *pd);
void g2_allocate_basic_colors_pd(g2_physical_device *pd);
void g2_set_font_size_pd(g2_physical_device *pd, double size);
void g2_set_line_width_pd(g2_physical_device *pd, double w);
void g2_set_dash_pd(g2_physical_device *pd, int N, double *dashes);
void g2_query_pointer_pd(g2_physical_device *pd,
			 double *x, double *y, unsigned int *button);
void g2_get_pd_handles_pd(g2_physical_device *pd, void *handles[G2_PD_HANDLES_SIZE]);

#endif /* _G2_CONTROL_PD_H */
