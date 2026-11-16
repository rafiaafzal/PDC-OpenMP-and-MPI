CC     = gcc
MPICC  = mpicc
CFLAGS = -O3 -march=native -std=gnu11
all: serial parallel
serial: serial.c
	$(CC) $(CFLAGS) -o $@ $< -lm
parallel: parallel.c
	$(MPICC) $(CFLAGS) -o $@ $< -lm
clean:
	rm -f serial parallel
