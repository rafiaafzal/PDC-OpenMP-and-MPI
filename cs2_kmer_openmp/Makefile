CC     = gcc
CFLAGS = -O3 -march=native -std=gnu11
all: serial parallel
serial: serial.c common.h
	$(CC) $(CFLAGS) -o $@ $<
parallel: parallel.c common.h
	$(CC) $(CFLAGS) -fopenmp -o $@ $<
clean:
	rm -f serial parallel
