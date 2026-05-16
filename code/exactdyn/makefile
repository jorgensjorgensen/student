# Target program name
TARGET = exact

# Compiler (cc or gcc)
CC = gcc

# Compiler flags
#CFLAGS = -Aa -c -O -D_INCLUDE_AES_SOURCE
#CFLAGS = -ansi -pedantic -c -O3 -DALLOC_FULL
CFLAGS = -Wall -c -Ofast

# linker flags
# -lm -lefence
LDFLAGS = -lm

# Binary object files
OBJS = exact.o calgebra.o gauss.o setup.o tools.o construct.o\
       io.o mfeglobals.o functions.o discretize.o utils.o alloc.o

# Debugging - add the following:
# -g -y to CFLAGS
# /usr/lib/end.o to OBJS

$(TARGET): $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) -o $(TARGET)
exact.o: exact.c calgebra.h gauss.h setup.h tools.h construct.h\
         io.h mfeglobals.h functions.h discretize.h utils.h\
         mfetypes.h alloc.h
	$(CC) $(CFLAGS) exact.c
calgebra.o: calgebra.c calgebra.h
	$(CC) $(CFLAGS) calgebra.c
gauss.o: gauss.c gauss.h calgebra.h mfetypes.h utils.h alloc.h
	$(CC) $(CFLAGS) gauss.c
setup.o: setup.c setup.h alloc.h mfetypes.h calgebra.h
	$(CC) $(CFLAGS) setup.c
tools.o: tools.c tools.h mfetypes.h calgebra.h
	$(CC) $(CFLAGS) tools.c
construct.o: construct.c construct.h tools.h utils.h mfetypes.h\
             calgebra.h
	$(CC) $(CFLAGS) construct.c
io.o: io.c io.h mfeglobals.h alloc.h mfetypes.h calgebra.h
	$(CC) $(CFLAGS) io.c
mfeglobals.o: mfeglobals.c mfeglobals.h	mfetypes.h calgebra.h
	$(CC) $(CFLAGS) mfeglobals.c
functions.o: functions.c functions.h mfetypes.h mfeglobals.h\
             calgebra.h
	$(CC) $(CFLAGS) functions.c
discretize.o: discretize.c discretize.h functions.h mfetypes.h\
              mfeglobals.h calgebra.h
	$(CC) $(CFLAGS) discretize.c
utils.o: utils.c utils.h mfetypes.h calgebra.h
	$(CC) $(CFLAGS) utils.c
alloc.o: alloc.c alloc.h
	$(CC) $(CFLAGS) alloc.c

clean:
	rm -f *.o $(TARGET)

