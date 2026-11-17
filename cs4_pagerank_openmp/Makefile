CC     = gcc
CFLAGS = -O3 -march=native -std=gnu11
all: serial parallel
serial: serial.c pr_common.h
	$(CC) $(CFLAGS) -o $@ $< -lm
parallel: parallel.c pr_common.h
	$(CC) $(CFLAGS) -fopenmp -o $@ $< -lm
clean:
	rm -f serial parallel
