NAME1    = RNAsubopt

VIENNA   = ..

CFLAGS   = -g -Wall -mpentiumpro -DSUBOPT -I$(VIENNA)/H  -DDMALLOC

LDFLAGS  = -O -g1 -mpentiumpro -L$(VIENNA)/lib -lRNA -lm  -ldmalloc


FILES   = RNAsubopt.c\
          fold.c\
          list.c\
          subopt.c
          

OBJ     = RNAsubopt.o\
          fold.o\
          list.o\
          subopt.o
          

all:	$(NAME1) popt

$(NAME1):   $(OBJ)
	    $(CC) -o $(NAME1) $(OBJ) $(LDFLAGS)

popt:	    popt.o
	    $(CC) -o popt popt.o $(LDFLAGS)

fold.o:	$(VIENNA)/lib/fold.c
	   $(CC) $(CFLAGS) -c $?

veryclean:      
	   rm -f *.o $(NAME1)

clean:      
	   rm -f *.o
