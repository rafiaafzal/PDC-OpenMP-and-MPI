CXX    = g++
MPICXX = mpicxx
FLAGS  = -O3 -march=native -std=gnu++17
all: serial parallel pingpong
serial: serial.cpp keygen.h
	$(CXX) $(FLAGS) -o $@ $<
parallel: parallel.cpp keygen.h
	$(MPICXX) $(FLAGS) -o $@ $<
pingpong: pingpong.cpp
	$(MPICXX) $(FLAGS) -o $@ $<
clean:
	rm -f serial parallel pingpong
