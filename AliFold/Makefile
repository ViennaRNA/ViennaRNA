NAME    = RNAalifold

VIENNA   = ..

CFLAGS   = -g -O -Wall -mpentiumpro -I$(VIENNA)/H  # -DDMALLOC

LDFLAGS  = -g -O -mpentiumpro -L$(VIENNA)/lib -lRNA -lm # -ldmalloc


FILES   = RNAalifold.c\
          alifold.c 

OBJ     = RNAalifold.o\
          alifold.o \
	  alipfold.o

all:	$(NAME)

RNAalifold:   $(OBJ)
	$(CC) -o $(NAME) $(OBJ) $(LDFLAGS)

veryclean:      
	rm -f *.o $(NAME)

clean:      
	rm -f *.o
