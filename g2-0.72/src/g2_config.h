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
#ifndef _G2_CONFIG_H
#define _G2_CONFIG_H

/*
 *
 *  Some g2 defines
 *
 *  Default behaviour is usually ok for most cases
 *
 */

/* g2 warnings are printed to stderr, larger 0=quiet, 1=errors, 2=warnings, 3=verbose, 4=debug */
#define g2_LogLevel 1

/* set to 0 to disable backing store emulation */
#define g2_EmulateBackingStore  1

/* X11 font, note %d for font size */
#define g2_X11Font "-*-times-medium-r-normal--%d-*-*-*-*-*-*-*"

/* PostScript font */
#define g2_PSFont "/Times-Roman"

#endif /* _G2_CONFIG_H */
