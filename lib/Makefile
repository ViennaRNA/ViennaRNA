# Edit the options in CFLAGS and ARFLAGS to suit your system.
# make install will copy libRNA.a to LIBDIR
SHELL	= /bin/sh

RANLIB	= ranlib
AR	= ar
ARFLAGS = rv

RM	= rm -f

CC	= cc -mpentiumpro -Wall -O -g1  # insert your favorite compiler here

# add -DSUN4 if compiling for SunOS 4.1*  or -DHP9  for hpux-9.?? for faster 
# handling of floting point exceptions
# add -DLARGE_PF if you want pf_fold() to use doubles instead of floats
CFLAGS  = -I../H

LIBDIR	= /usr/local/lib

OBJ	= \
	fold.o \
	part_func.o \
	fold_vars.o \
	PS_dot.o \
	naview.o \
	RNAstruct.o \
	treedist.o \
	stringdist.o \
	ProfileDist.o \
	dist_vars.o \
	inverse.o \
	energy_par.o \
	read_epars.o \
	utils.o

all:	libRNA.a
	
energy_par.o: intloops.h

libRNA.a: $(OBJ)
	$(RM) $@
	$(AR) $(ARFLAGS) $@ $(OBJ)
	-$(RANLIB) $@

install:
	cp libRNA.a $(LIBDIR) 

clean:	
	$(RM) libRNA.a $(OBJ)

