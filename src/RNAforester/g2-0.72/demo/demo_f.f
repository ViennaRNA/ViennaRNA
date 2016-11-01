c
c  Copyright (C) 1998-2001  Ljubomir Milanovic & Horst Wagner
c  This file is part of the g2 library
c
c  This library is free software; you can redistribute it and/or
c  modify it under the terms of the GNU Lesser General Public
c  License as published by the Free Software Foundation; either
c  version 2.1 of the License, or (at your option) any later version.
c
c  This library is distributed in the hope that it will be useful,
c  but WITHOUT ANY WARRANTY; without even the implied warranty of
c  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
c  Lesser General Public License for more details.
c
c  You should have received a copy of the GNU Lesser General Public
c  License along with this library; if not, write to the Free Software
c  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
c
        real demo_f
        real  a, b
        real d, d1, d2
        real color

        d=g2_open_vd()
        write (6,*) d        
        d1=g2_open_X11(100.0, 100.0)
        write (6,*) d1
        d2=g2_open_PS('demo_f.ps', 4.0, 1.0)
        write (6,*) d2

        call g2_attach(d, d1)
        call g2_attach(d, d2)

        call g2_plot(d, 50.0, 50.0)
        call g2_arc(d, 50.0, 50.0, 30.0, 20.0, 45.0, 180.0)

        color=g2_ink(d1, 1.0, 0.0, 0.0)
        call g2_pen(d1, color)
        write (6,*) color
        call g2_string(d1, 15.0, 75.0, 'TEST (Window)')

        color=g2_ink(d2, 0.0, 1.0, 0.0)
        call g2_pen(d2, color)
        write (6,*) color        
        call g2_string(d2, 15.0, 75.0, 'TEST (File)')

        call g2_pen(d, 1.0)
        call g2_circle(d, 20.0, 20.0, 10.0)
        call g2_string(d, 20.0, 20.0, 'All devices!')
        call g2_flush(d)

        call g2_close(d2)

        read (*,*) a

        stop
        end


