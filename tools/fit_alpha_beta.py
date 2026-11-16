#!/usr/bin/env python3
"""Least-squares fit of T(m) = alpha + beta*m from the pingpong CSV, plus Alltoallv prediction.
usage: python3 fit_alpha_beta.py pingpong.csv [N p]
If N and p are given, predicts the Alltoallv time for N 4-byte keys on p ranks, assuming each rank
sends (p-1)/p of its N/p keys in p-1 messages:  T ~ (p-1)*alpha + beta*(N/p)*4*(p-1)/p   (no contention)."""
import sys, csv
pts = [(float(r["bytes"]), float(r["one_way_seconds"])) for r in csv.DictReader(open(sys.argv[1]))]
# weight small and large messages evenly by fitting on relative error (weights 1/t)
sw = sum(1/t**2 for m, t in pts); sx = sum(m/t**2 for m, t in pts); sy = sum(1/t for m, t in pts)
sxx = sum(m*m/t**2 for m, t in pts); sxy = sum(m/t for m, t in pts)
beta = (sw*sxy - sx*sy) / (sw*sxx - sx*sx); alpha = (sy - beta*sx) / sw
print(f"alpha = {alpha*1e6:.2f} us   beta = {beta*1e9:.4f} ns/byte   (bandwidth ~ {1/beta/1e9:.2f} GB/s)")
if len(sys.argv) >= 4:
    N, p = int(sys.argv[2]), int(sys.argv[3])
    print(f"predicted Alltoallv time: {(p-1)*alpha + beta*(N/p)*4*(p-1)/p:.4f} s")
