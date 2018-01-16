OBJECTS = main.o metropolis.o ranlxd.o variance.o
CC = gcc
CFLAGS = -Wall -lm
EXEC_NAME = prog

all: main.o metropolis.o ranlxd.o variance.o
	$(CC) -o $(EXEC_NAME) main.o metropolis.o ranlxd.o variance.o $(CFLAGS)

main.o:	main.c common.h ranlxd.h
	$(CC) -c main.c $(CFLAGS)

metropolis.o: metropolis.c ranlxd.h common.h
	$(CC) -c metropolis.c $(CFLAGS)

ranlxd.o: ranlxd.c ranlxd.h common.h
	$(CC) -c ranlxd.c $(CFLAGS)

variance.o: variance.c common.h
	$(CC) -c variance.c $(CFLAGS)

clean:
	rm -f $(OBJECTS) $(EXEC_NAME)

