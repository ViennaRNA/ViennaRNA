!
! VMS descrip.mms for g2 demos
!

!
! Change these pointers to the correct location of libGD,libpng, etc
!
GD_LIBS = [.-.-.-.GD-1_8_4]libgd.olb/lib,[.-.-.-.lpng120]libpng.olb/lib,[.-.-.-.ZLIB]libz.olb/lib,[.-.-.-.JPEG-6B]libjpeg.olb/lib

.IFDEF debug
CC_FLAGS        = /noopt/debug/nolist/warnings \
                  /include_directory=([-.src],[-.src.X11],[-.src.PS],[-.src.GD]) \
                  /define=(DO_PS, DO_X11, DO_GD)
LINK_FLAGS      = 
.ELSE
CC_FLAGS        = /nolist/warnings \
                  /include_directory=([-.src],[-.src.X11],[-.src.PS],[-.src.GD]) \
                  /define=(DO_PS, DO_X11, DO_GD)
LINK_FLAGS      = 
.ENDIF

!
! Rules
!
.c.obj	: ; cc $(CC_FLAGS)/obj=$(mms$target) $(mms$source)


!
! exe files
!
all	: g2_test.exe, simple_X11.exe, simple_PS.exe
	@ write sys$output "Demo files compiled."


g2_test.exe	: g2_test.obj
	link $(LINK_FLAGS)/exe=g2_test.exe g2_test.obj, \
             [-]g2.olb/lib, []xlink.opt/opt, $(GD_LIBS)

simple_X11.exe	: simple_X11.obj
	link $(LINK_FLAGS)/exe=simple_X11.exe simple_X11.obj, \
             [-]g2.olb/lib, []xlink.opt/opt

simple_PS.exe	: simple_PS.obj
	link $(LINK_FLAGS)/exe=simple_PS.exe simple_PS.obj, [-]g2.olb/lib



clean	:
	delete *.exe;* , *.obj;* , *.olb;*


!
! Dependencies
!
g2_test.obj	: g2_test.c
simple_X11.obj	: simple_X11.c
simple_PS.obj	: simple_PS.c
