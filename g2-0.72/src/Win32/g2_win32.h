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
#ifndef _G2_WIN32_H
#define _G2_WIN32_H

#if defined(__cplusplus)
extern "C"
{
#endif

#include <stdio.h>

// Common Library header for DLL and application
#ifdef WIN32
#ifdef G2DLL
#ifdef MAKEDLL
/* Create DLL */
#define LINKDLL __declspec( dllexport)
#else
/* Use DLL */
#define LINKDLL __declspec( dllimport)
#endif
#else 
/* Use static win32 */
#define LINKDLL
#endif
#else
/* Use non-win32 */
#define LINKDLL
#endif

/**
 * Window type
 *
 * \ingroup win32
 */
enum g2_win32_type {
    g2_win32,			/**< regular window */
    g2_wmf32			/**< windows meta file */
};


LINKDLL int  g2_open_win32(int width, int height, const char *filename,int type);


#if defined(__cplusplus)
} /* end extern "C" */
#endif


#endif /* _G2_WIN32_H */
