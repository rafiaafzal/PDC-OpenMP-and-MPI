# Top-level: build everything.  `make` builds all; `make clean` cleans all.
DIRS = cs1_heat_mpi cs2_kmer_openmp cs3_montecarlo cs4_pagerank_openmp cs5_samplesort_mpi
all:
	@for d in $(DIRS); do $(MAKE) -C $$d || exit 1; done
clean:
	@for d in $(DIRS); do $(MAKE) -C $$d clean; done
