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
#include <g2.h>
#include <g2_gd.h>

int main()
{
    int d;

#ifdef DO_GIF
    d=g2_open_gd("simple.gif", 100, 100, g2_gd_gif);
    g2_line(d, 10, 10, 90, 90);
    g2_close(d);
#endif

    d=g2_open_gd("simple.png", 100, 100, g2_gd_png);
    g2_line(d, 10, 10, 90, 90);
    g2_close(d);

    return 0;
}

