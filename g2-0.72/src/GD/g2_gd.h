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
#ifndef _G2_GD_H
#define _G2_GD_H

#if defined(__cplusplus)
extern "C"
{
#endif

#include <stdio.h>

/* Common Library header for DLL and application */
#ifdef WIN32
#ifdef G2DLL
#ifdef MAKEDLL
/* Create DLL */
#define G2L __declspec( dllexport)
#else
/* Use DLL */
#define G2L __declspec( dllimport)
#endif
#else 
/* Use static win32 */
#define G2L
#endif
#else
/* Use non-win32 */
#define G2L
#endif

/**
 *  g2 gd bitmap types
 *
 * \ingroup GD
 */
enum g2_gd_type {
    g2_gd_jpeg = 0,	/**< jpeg */
    g2_gd_png  = 1,	/**< png */
    g2_gd_gif  = 2	/**< gif */
};


G2L int  g2_open_gd(const char *filename, int width, int height,enum g2_gd_type gd_type);

#if defined(__cplusplus)
} /* end extern "C" */
#endif


#endif /* _G2_GD_H */
