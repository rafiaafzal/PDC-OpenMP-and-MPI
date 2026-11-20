CC     = gcc
MPICC  = mpicc
CFLAGS = -O3 -march=native -std=gnu11
all: serial parallel parallel_mpi
serial: serial.c rng.h
	$(CC) $(CFLAGS) -o $@ $< -lm
parallel: parallel.c rng.h
	$(CC) $(CFLAGS) -fopenmp -o $@ $< -lm
parallel_mpi: parallel_mpi.c rng.h
	$(MPICC) $(CFLAGS) -o $@ $< -lm
clean:
	rm -f serial parallel parallel_mpi
