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

c
c  This file will be processed by cpp, the outpu g2testf.f is feeded into fortran compiler
c  See Makefile for details.
c

c	program main
c	implicit none
	integer  ndev, i, j
	real d, c, dev(0:10)
	character str*256
	real pts(0:10)
	real y
	include 'penguin.inc'
c	integer image(64)
c	data image/
c     & 0, 0, 2, 2, 2, 2, 0, 0,
c     & 0, 2, 0, 0, 0, 0, 2, 0,
c     & 2, 0, 3, 0, 0, 3, 0, 2,
c     & 2, 0, 0, 0, 0, 0, 0, 2,
c     & 2, 3, 0, 0, 0, 0, 3, 2,
c     & 2, 0, 3, 3, 3, 3, 0, 2,
c     & 0, 2, 0, 0, 0, 0, 2, 0,
c     & 0, 0, 2, 2, 2, 2, 0, 0 /
c /*open virtual device */
	d=g2_open_vd()
	print*,"Adding.. VD=", d
	ndev = 1

#ifdef DO_PS
	dev(ndev)=g2_open_ps('g2testf.ps', 4., 0.)
	print*,"..PS=",dev(ndev)
	call g2_attach(d, dev(ndev))
	ndev = ndev + 1
#endif
#ifdef DO_X11
	dev(ndev)=g2_open_x11(775., 575.)
	print*,"..X11=",dev(ndev)
	call g2_attach(d,dev(ndev))
	ndev = ndev + 1
#endif
#ifdef DO_GD
	str='g2testf.png'
	str(12:12)=char(0)
	dev(ndev)=g2_open_gd(str, 775., 575., 0.)
	print*,"..GD=",dev(ndev)
	call g2_attach(d,dev(ndev))
	ndev = ndev + 1
#endif
	call g2_set_auto_flush(d,0.)

c        call g2_set_coordinate_system(d, 775., 575., -0.75, -1.0)

	do i=0,27
	 call g2_pen(d, float(i))
	 call g2_filled_circle(d, float(i*20+10), 10., 10.) 
	 call g2_pen(d, 1.)
	 call g2_circle(d, float(i*20+10), 10., 10.)
	 write(str(1:4),'(i3,a1)') i,char(0)
	 call g2_string(d, float(i*20+7), 21., str(1:4))
	enddo
    
    	do j=1,ndev
	 if(dev(j).gt.0) then
	  do i=0,64
	   call g2_move(dev(j), float(2*i+575), 5.)
	   c = g2_ink(dev(j), float(i)/64., 0., 0.)
	   call g2_pen(dev(j), c)
	   call g2_line_r(dev(j), 0., 20.)
	   c = g2_ink(dev(j), 0., float(i)/64., 0.)
	   call g2_pen(dev(j), c )
	   call g2_line_r(dev(j), 10., 20.)
	   c = g2_ink(dev(j), 0., 0., float(i)/64.)
	   call g2_pen(dev(j), c)
	   call g2_line_r(dev(j), -10., 20.)
	  enddo
	 endif
	enddo
	call g2_pen(d, 1.)
	call g2_line(d, 200., 50., 350., 50.)
	call g2_line(d, 200., 48., 350., 48.)
	call g2_line(d, 200., 46., 350., 46.)
	call g2_line(d, 200., 46., 200., 75.)
	call g2_line(d, 198., 46., 198., 75.)
	call g2_line(d, 196., 46., 196., 75.)
	str="012abcABC#())(\\-+~*!$%&"//char(0)
	call g2_string(d, 200., 50., str)
    
    	call g2_pen(d, 1.)
    	do i=1,25
	 call g2_line(d, 15., float(i*20+50), 15., float(i*20+50+i))
	 call g2_set_font_size(d, 12.)
	 write(str(1:3),'(i2,a1)') i,char(0)
	 call g2_string(d, 20., float(i*20+50), str)
	 call g2_set_font_size(d, float(i))
	 str='Hello, world!'//char(0)
	 call g2_string(d, 40., float(i*20+50), str)
	enddo


	call g2_plot(d, 150., 70.)
	call g2_line(d, 147., 68., 153., 68.)
		
       y=100.
	call g2_line(d, 100., y, 150., y+50.)
	call g2_triangle(d, 150., y, 250., y, 200., y+50.)
	call g2_rectangle(d, 300., y, 400., y+50.)
	call g2_circle(d, 450., y+25., 25.)
	call g2_ellipse(d, 550., y+25., 45., 25.)
	call g2_arc(d, 650., y+25., 25., 45., 90., 360.)
    
       y=200.
	call g2_filled_triangle(d, 150., y, 250., y, 200., y+50.)
	call g2_filled_rectangle(d, 300., y, 400., y+50.)
	call g2_filled_circle(d, 450., y+25., 25.)
	call g2_filled_ellipse(d, 550., y+25., 45., 25.)
	call g2_filled_arc(d, 650., y+25., 25., 45., 90., 360.)


	y=300.
	pts(0)=150.
	pts(1)=y
	pts(2)=175.
	pts(3)=y+100.
	pts(4)=200.
	pts(5)=y
	pts(6)=225.
	pts(7)=y+100.
	pts(8)=250.
	pts(9)=y
	call g2_poly_line(d, 5., pts)
	call g2_pen(d, 19.)
	call g2_b_spline(d, 5., pts, 20.)
	call g2_pen(d, 1.)
	
	pts(0)=300.
	pts(1)=y
	pts(2)=350.
	pts(3)=y
	pts(4)=375.
	pts(5)=y+50.
	pts(6)=325.
	pts(7)=y+90.
	pts(8)=275.
	pts(9)=y+50.
	call g2_polygon(d, 5., pts)

	pts(0)=450.
	pts(1)=y
	pts(2)=500.
	pts(3)=y
	pts(4)=525.
	pts(5)=y+50.
	pts(6)=475.
	pts(7)=y+90.
	pts(8)=425.
	pts(9)=y+50.
	call g2_filled_polygon(d, 5, pts)
	
	call g2_image(d, 55., 50., 53., 62., penguin)
	call g2_image(d, 75., 130., 53., 62., penguin)
	call g2_pen(d, 1.)
    
	call g2_line(d, 225., 448., float(200+19*25), 448.)
    	do i=1,20
    	 call g2_pen(d,float(i+1))
	 call g2_set_line_width(d, float(i))
	 call g2_move(d, float(200+i*25), 450.)
	 call g2_line_to(d, float(200+i*25), 550.)
       enddo
	call g2_pen(d,1.)
    
	call g2_set_line_width(d, 5.)
      do i=1,10
	pts(0)=float(1*i)
	pts(1)=float(2*i)
	pts(2)=float(3*i)
	call g2_set_dash(d, 3., pts)
	call g2_line(d, 550., float(300+i*8), 750., float(350+i*8)) 
      enddo

	call g2_set_dash(d, 0., pts)
	call g2_set_line_width(d, 5.)
	call g2_arc(d, 740., 180., 25., 100., -45.+15., -45.-15.)
	call g2_filled_arc(d, 740., 180., 12., 50., -45.+15., -45.-15.)
c----------------------------------------------------------
	call g2_flush(d)
	print*,"Done...(Enter)"
	read(*,'(a)') str
	call g2_close(d)
	end
c #eof#
