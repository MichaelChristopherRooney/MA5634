OBJECTS = main.o metropolis.o ranlxd.o variance.o gaussian_rv.o
CC = gcc
CFLAGS = -Wall -lm
EXEC_NAME = prog

all: main.o metropolis.o ranlxd.o variance.o gaussian_rv.o
	$(CC) -o $(EXEC_NAME) main.o metropolis.o ranlxd.o gaussian_rv.o variance.o $(CFLAGS)

main.o:	main.c common.h ranlxd.h
	$(CC) -c main.c $(CFLAGS)

metropolis.o: metropolis.c ranlxd.h common.h
	$(CC) -c metropolis.c $(CFLAGS)

ranlxd.o: ranlxd.c ranlxd.h common.h
	$(CC) -c ranlxd.c $(CFLAGS)

variance.o: variance.c common.h
	$(CC) -c variance.c $(CFLAGS)

gaussian_rv.o: gaussian_rv.c common.h
	$(CC) -c gaussian_rv.c $(CFLAGS)

clean:
	rm -f $(OBJECTS) $(EXEC_NAME)

