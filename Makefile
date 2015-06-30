#***************************************************************************
#*
#* Copyright 2015 Ian Liu Rodrigues and Edson Borin
#*
#***************************************************************************

# C Compiler
CC=icc
# Compiler flags
CFLAGS=-O3 -std=gnu99 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -fopenmp
LINK_FLAGS=-fopenmp -lm 

INCLUDES=-I ./
OBJECTS=cmp.o gather.o interpol.o log.o semblance.o su.o

cmp.x: $(OBJECTS)
	$(CC) $(OBJECTS) $(LINK_FLAGS) -o cmp.x

%.o : %.c
	$(CC) $(INCLUDES) -c $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS) cmp.x c.su cmp.coher.su cmp.stack.su
