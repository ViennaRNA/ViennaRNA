# Edit CFLAGS to suit your compiler; make install will copy executables
# and man pages to BINDIR and LIBDIR

NAME    =  AnalyseDists
SNAME   =  AnalyseSeqs

#CC	= gcc  # add your favorite compiler here
CFLAGS  = -O

BINDIR  = ../bin
MANDIR  = ../man/man1

LIBS    = -lm

OBJS    = split.o \
	  distance_matrix.o \
          cluster.o \
          treeplot.o \
	  AD_main.o \
	  utils.o

SOBJS   = statgeom.o \
	  distance_matrix.o \
          cluster.o \
	  AS_main.o \
	  split.o \
          treeplot.o \
	  PS3D.o \
	  utils.o

all: $(NAME) $(SNAME)

$(NAME):    $(OBJS)
	    $(CC) -o $(NAME) $(CFLAGS) $(OBJS) $(LIBS)

$(SNAME):   $(SOBJS)
	    $(CC) -o $(SNAME) $(CFLAGS) $(SOBJS) $(LIBS)

utils.o:    ../lib/utils.c
	    $(CC) -c ../lib/utils.c $(CFLAGS)

clean:
	rm -f $(OBJS) $(SOBJS) $(NAME) $(SNAME)

install:    all
	    cp $(NAME) $(BINDIR)
	    cp $(SNAME) $(BINDIR)
	    cp ../man/$(NAME).man $(MANDIR)/$(NAME).1
	    cp ../man/$(SNAME).man $(MANDIR)/$(SNAME).1
