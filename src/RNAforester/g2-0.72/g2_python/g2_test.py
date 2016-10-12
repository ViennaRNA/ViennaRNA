##############################################################################
##  Copyright (C) 1998-2001  Ljubomir Milanovic & Horst Wagner
##  Copyright (C) 2006  Tijs Michels
##  This file is part of the g2 library
##
##  This library is free software; you can redistribute it and/or
##  modify it under the terms of the GNU Lesser General Public
##  License as published by the Free Software Foundation; either
##  version 2.1 of the License, or (at your option) any later version.
##
##  This library is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
##  Lesser General Public License for more details.
##
##  You should have received a copy of the GNU Lesser General Public
##  License along with this library; if not, write to the Free Software
##  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
##############################################################################

import sys
import g2
from g2 import *
from penguin import penguin

dev = []

print 'G2 version', G2_VERSION
print 'adding..',
if hasattr(g2, 'g2_open_PS'):
    print '..PS ..EPSF ..EPSF_CLIP',
    dev.append(g2_open_PS('g2_test.ps', g2_A4, g2_PS_land))
    dev.append(g2_open_EPSF('g2_test.eps'))
    dev.append(g2_open_EPSF_CLIP('g2_test_clip.eps', 200, 200))
if hasattr(g2, 'g2_open_FIG'):
    print '..FIG',
    dev.append(g2_open_FIG('g2_test.fig'))
if hasattr(g2, 'g2_open_X11'):
    print '..X11',
    dev.append(g2_open_X11(775, 575))
if hasattr(g2, 'g2_open_PS'):
    print '..GD(png) ..GD(jpeg)',
    dev.append(g2_open_gd('g2_test.png', 775, 575, g2_gd_png))
    dev.append(g2_open_gd('g2_test.jpg', 775, 575, g2_gd_jpeg))
if hasattr(g2, 'g2_open_win32'):
    print '..WIN32 ..WMF32'
    dev.append(g2_open_win32(775, 575, 'g2_test',     g2_win32))
    dev.append(g2_open_win32(775, 575, 'g2_test.emf', g2_wmf32))
print

graph = g2_open_vd() # open virtual device
for a in dev:
    graph.g2_attach(a)
graph.g2_set_auto_flush(False)

# graph.g2_set_coordinate_system(775, 575, -0.75, -1.0)

for i in range(0,27):
    x = i*20
    graph.g2_pen(i)
    graph.g2_filled_circle(x+10, 10, 10)
    graph.g2_pen(1)
    graph.g2_circle(x+10, 10, 10)
    graph.g2_string(x+7,  21, str(i))

for a in dev:
    # work on one physical device at a time
    for i in range(0,65):
        a.g2_move(i+i+575, 5)
        a.g2_pen(a.g2_ink(i/64., 0, 0))
        a.g2_line_r(0, 20)
        a.g2_pen(a.g2_ink(0, i/64., 0))
        a.g2_line_r(10, 20)
        a.g2_pen(a.g2_ink(0, 0, i/64.))
        a.g2_line_r(-10, 20)

graph.g2_pen(1)
lines = ((200, 50, 350, 50),
         (200, 48, 350, 48),
         (200, 46, 350, 46),
         (200, 46, 200, 75),
         (198, 46, 198, 75),
         (196, 46, 196, 75))
for x1, y1, x2, y2 in lines:
    graph.g2_line(x1, y1, x2, y2) # this is why we love Python
graph.g2_string(200, 50, '012abcABC#())(\\-+~*!$%&')

# graph.g2_pen(1)
for i in range(1,25):
    y = i*20+50
    graph.g2_line(15, y, 15, y+i)
    graph.g2_set_font_size(12)
    graph.g2_string(20, y, '%2d:' % i)
    graph.g2_set_font_size(i)
    graph.g2_string(40, y, 'hello, world')

graph.g2_plot(150, 70)
graph.g2_line(147, 68, 153, 68)

y = 100
graph.g2_line(120, y, 170, y+50)
graph.g2_triangle(150, y, 250, y, 200, y+50)
graph.g2_rectangle(300, y, 400, y+50)
graph.g2_circle(450, y+25, 25)
graph.g2_ellipse(550, y+25, 45, 25)
graph.g2_arc(650, y+25, 25, 45, 90, 360)

graph.g2_set_dash([4, 4]) # Python specific
graph.g2_line(120+5, y, 170+5, y+50)
graph.g2_triangle(150+10, y+4, 250-10, y+4, 200, y+50-5)
graph.g2_rectangle(305, y+5, 395, y+50-5)
graph.g2_circle(450, y+25, 20)
graph.g2_ellipse(550, y+25, 40, 20)
graph.g2_arc(650, y+25, 20, 40, 90, 360)
graph.g2_set_solid() # Python only

y = 200
graph.g2_filled_triangle(150, y, 250, y, 200, y+50)
graph.g2_filled_rectangle(300, y, 400, y+50)
graph.g2_filled_circle(450, y+25, 25)
graph.g2_filled_ellipse(550, y+25, 45, 25)
graph.g2_filled_arc(650, y+25, 25, 45, 90, 360)

y = 300
pts = [150, y,
       175, y+100,
       200, y,
       225, y+100,
       250, y]
graph.g2_poly_line(pts) # Python specific

graph.g2_pen(19)
graph.g2_b_spline(pts, 20) # Python specific
graph.g2_pen(1)

pts = [300, y,
       350, y,
       375, y+50,
       325, y+90,
       275, y+50]
graph.g2_polygon(pts) # Python specific

pts = [450, y,
       500, y,
       525, y+50,
       475, y+90,
       425, y+50]
graph.g2_filled_polygon(pts) # Python specific

graph.g2_image(55.,  50., penguin) # Python specific
graph.g2_image(75., 130., penguin) # Python specific
graph.g2_pen(1)

graph.g2_line(225, 448, 200+19*25, 448)
for i in range(1,20):
    graph.g2_pen(i+1)
    graph.g2_set_line_width(i)
    graph.g2_move(200+i*25, 450)
    graph.g2_line_to(200+i*25, 550)
graph.g2_pen(1)

graph.g2_set_line_width(5)
for i in range(1,10):
    graph.g2_set_dash([i, i*2, i*3]) # Python specific
    graph.g2_line(550, 300+i*8, 750, 350+i*8)

graph.g2_set_solid() # Python only
graph.g2_set_line_width(4)
graph.g2_arc(       740, 180, 25, 100, -45+15, -45-15)
graph.g2_arc(       740, 180, 20,  95, -45-15, -45+15)
graph.g2_filled_arc(740, 180, 12,  50, -45+15, -45-15)

graph.g2_set_line_width(1)
graph.g2_circle( 400, 400, 20)
graph.g2_ellipse(400, 400, 25, 25)
graph.g2_arc(    400, 400, 30, 30, 0, 360)

graph.g2_flush()
print '\nDone.\n[Enter]\n'
sys.stdin.read(1)
graph.g2_close()
