// CS5 distributed sample sort with MPI.
// usage: mpirun -np P ./parallel N dist OVERSAMPLE [verbose]
//   N          total number of 32-bit keys (each rank generates its own N/P share)
//   dist       0 uniform, 1 skewed (Zipf-like), 2 many duplicates
//   OVERSAMPLE each rank contributes OVERSAMPLE*(P-1) regular samples (1 = textbook p-1)
//   verbose    1 = print the final element count of every rank
// Phases: (1) local sort  (2) sampling: pick samples, Gather, choose splitters, Bcast, partition
//         (3) Alltoall (counts) + Alltoallv (data)  (4) local p-way merge of received runs.
#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <numeric>
#include "keygen.h"

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, p;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    uint64_t N = argc > 1 ? strtoull(argv[1], 0, 10) : 100000000ULL;
    int dist = argc > 2 ? atoi(argv[2]) : 0;
    int over = argc > 3 ? atoi(argv[3]) : 1;
    int verbose = argc > 4 ? atoi(argv[4]) : 0;
    if (over < 1) over = 1;

    std::vector<uint32_t> a;
    gen_keys(a, my_count(N, p, rank), rank, dist);
    uint64_t sum_in_local = 0; for (uint32_t x : a) sum_in_local += x;
    uint64_t nloc = a.size();

    MPI_Barrier(MPI_COMM_WORLD);
    double T0 = wall();

    /* (1) local sort */
    std::sort(a.begin(), a.end());
    double T1 = wall();

    /* (2) sampling, splitters, partition */
    std::vector<uint32_t> splitters(p > 1 ? p - 1 : 0);
    std::vector<int> scount(p, 0), sdisp(p, 0);
    if (p > 1) {
        int ns = over * (p - 1);
        std::vector<uint32_t> samples(ns);
        for (int i = 0; i < ns; i++) {
            uint64_t idx = nloc ? ((uint64_t)(i + 1) * nloc) / (ns + 1) : 0;
            samples[i] = nloc ? a[idx] : 0;
        }
        std::vector<uint32_t> all;
        if (rank == 0) all.resize((size_t)ns * p);
        MPI_Gather(samples.data(), ns, MPI_UNSIGNED, all.data(), ns, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            std::sort(all.begin(), all.end());
            for (int i = 0; i < p - 1; i++) splitters[i] = all[((size_t)(i + 1) * all.size()) / p];
        }
        MPI_Bcast(splitters.data(), p - 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        size_t prev = 0;
        for (int i = 0; i < p; i++) {
            size_t end = (i < p - 1) ? (size_t)(std::upper_bound(a.begin(), a.end(), splitters[i]) - a.begin())
                                     : a.size();
            if (end < prev) end = prev;
            scount[i] = (int)(end - prev);
            prev = end;
        }
        for (int i = 1; i < p; i++) sdisp[i] = sdisp[i - 1] + scount[i - 1];
    } else {
        scount[0] = (int)a.size();
    }
    double T2 = wall();

    /* (3) redistribution */
    std::vector<int> rcount(p), rdisp(p, 0);
    MPI_Alltoall(scount.data(), 1, MPI_INT, rcount.data(), 1, MPI_INT, MPI_COMM_WORLD);
    for (int i = 1; i < p; i++) rdisp[i] = rdisp[i - 1] + rcount[i - 1];
    size_t nrecv = (size_t)rdisp[p - 1] + rcount[p - 1];
    std::vector<uint32_t> b(nrecv);
    double T2b = wall();
    MPI_Alltoallv(a.data(), scount.data(), sdisp.data(), MPI_UNSIGNED,
                  b.data(), rcount.data(), rdisp.data(), MPI_UNSIGNED, MPI_COMM_WORLD);
    double T3 = wall();
    std::vector<uint32_t>().swap(a);

    /* (4) merge the p sorted runs (pairwise tree of inplace_merge) */
    std::vector<size_t> bounds(rdisp.begin(), rdisp.end());
    bounds.push_back(nrecv);
    while (bounds.size() > 2) {
        std::vector<size_t> nb;
        for (size_t i = 0; i + 1 < bounds.size(); i += 2) {
            if (i + 2 < bounds.size()) {
                std::inplace_merge(b.begin() + bounds[i], b.begin() + bounds[i + 1], b.begin() + bounds[i + 2]);
                nb.push_back(bounds[i]);
            } else nb.push_back(bounds[i]);
        }
        nb.push_back(nrecv);
        bounds.swap(nb);
    }
    double T4 = wall();

    /* verification: local order, boundary order, count and checksum preserved */
    int ok = std::is_sorted(b.begin(), b.end()) ? 1 : 0;
    uint32_t mymax = b.empty() ? 0 : b.back(), mymin = b.empty() ? 0xFFFFFFFFu : b.front();
    uint32_t leftmax = 0;
    int has_left = 0;
    uint32_t mm[3] = {mymin, mymax, b.empty() ? 0u : 1u};
    std::vector<uint32_t> allmm(3 * p);
    MPI_Allgather(mm, 3, MPI_UNSIGNED, allmm.data(), 3, MPI_UNSIGNED, MPI_COMM_WORLD);
    for (int r = 0; r < rank; r++) if (allmm[3 * r + 2]) { leftmax = std::max(leftmax, allmm[3 * r + 1]); has_left = 1; }
    if (has_left && !b.empty() && leftmax > mymin) ok = 0;
    uint64_t sum_out_local = 0; for (uint32_t x : b) sum_out_local += x;
    uint64_t cnt_loc = b.size(), cnt, sin_, sout;
    int ok_all;
    MPI_Allreduce(&cnt_loc, &cnt, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&sum_in_local, &sin_, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&sum_out_local, &sout, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&ok, &ok_all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    int verified = ok_all && cnt == N && sin_ == sout;

    double ph[5] = {T1 - T0, T2 - T1, T2b - T2, T3 - T2b, T4 - T3}, phmax[5], tot = T4 - T0, totmax;
    /* ph: sort, sampling+partition, count exchange, alltoallv, merge */
    MPI_Reduce(ph, phmax, 5, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&tot, &totmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    std::vector<uint64_t> counts(p);
    MPI_Gather(&cnt_loc, 1, MPI_UINT64_T, counts.data(), 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        uint64_t mx = *std::max_element(counts.begin(), counts.end());
        uint64_t mn = *std::min_element(counts.begin(), counts.end());
        double mean = (double)N / p;
        printf("RESULT cs5 variant=samplesort N=%llu dist=%d p=%d oversample=%d time=%.6f t_sort=%.6f t_sample=%.6f "
               "t_counts=%.6f t_alltoallv=%.6f t_merge=%.6f max_per_rank=%llu min_per_rank=%llu max_over_mean=%.3f "
               "checksum=%llu %s\n",
               (unsigned long long)N, dist, p, over, totmax, phmax[0], phmax[1], phmax[2], phmax[3], phmax[4],
               (unsigned long long)mx, (unsigned long long)mn, mx / mean, (unsigned long long)sin_,
               verified ? "VERIFIED" : "VERIFY_FAIL");
        if (verbose) {
            printf("ELEMENTS_PER_RANK");
            for (int i = 0; i < p; i++) printf(" %llu", (unsigned long long)counts[i]);
            printf("\n");
        }
    }
    MPI_Finalize();
    return verified ? 0 : 2;
}
