# makefile

CFLAGS	= -I../H -g

LDFLAGS	= -L../lib -lRNA -lm

all:	RNAplot b2ct

RNAplot: RNAplot.o
	$(CC) -o $@ $@.o $(LDFLAGS)

clean:	
	rm -f *.o
	rm -f RNAplot b2ct

# End of file
