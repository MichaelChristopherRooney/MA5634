OBJECTS = metropolis.o ranlxd.o variance.o
CC = gcc
CFLAGS = -Wall -lm
EXEC_NAME = prog

all: metropolis.o ranlxd.o variance.o
	$(CC) -o $(EXEC_NAME) metropolis.o ranlxd.o variance.o $(CFLAGS)

metropolis.o: metropolis.c ranlxd.h
	$(CC) -c metropolis.c $(CFLAGS)

ranlxd.o: ranlxd.c ranlxd.h
	$(CC) -c ranlxd.c $(CFLAGS)

variance.o: variance.c
	$(CC) -c variance.c $(CFLAGS)

clean:
	rm -f $(OBJECTS) $(EXEC_NAME)

