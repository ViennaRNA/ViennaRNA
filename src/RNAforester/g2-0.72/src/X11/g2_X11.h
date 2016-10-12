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
#ifndef _G2_X11_H
#define _G2_X11_H

#if defined(__cplusplus)
extern "C"
{
#endif


int g2_open_X11(int width, int height);

int g2_open_X11X(int width, int height,
		 int x, int y,
		 char *window_name, char *icon_name,
		 char *icon_data,
		 int icon_width, int icon_height);

#if defined(__cplusplus)
} /* end extern "C" */
#endif

#endif /* _G2_X11_H */
