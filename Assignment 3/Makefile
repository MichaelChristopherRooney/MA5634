OBJECTS = metropolis.o ranlxd.o
CC = gcc
CFLAGS = -Wall -lm
EXEC_NAME = prog

all: metropolis.o ranlxd.o
	$(CC) -o $(EXEC_NAME) metropolis.o ranlxd.o $(CFLAGS)

metropolis.o: metropolis.c ranlxd.h
	$(CC) -c metropolis.c $(CFLAGS)

ranlxd.o: ranlxd.c ranlxd.h
	$(CC) -c ranlxd.c $(CFLAGS)

clean:
	rm -f $(OBJECTS) $(EXEC_NAME)

