#!/usr/bin/env python3
"""CS4: log-log degree distribution from the CSV written by  ./serial -d results/deg.csv
usage: python3 plot_degree.py results/deg.csv"""
import sys, csv
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
rows = list(csv.DictReader(open(sys.argv[1])))
d = [int(r["degree"]) for r in rows if int(r["degree"]) > 0]
plt.loglog(d, [int(r["count_in"]) for r in rows if int(r["degree"]) > 0], ".", label="in-degree")
plt.loglog(d, [int(r["count_out"]) for r in rows if int(r["degree"]) > 0], ".", label="out-degree")
plt.xlabel("degree"); plt.ylabel("number of vertices"); plt.legend(); plt.grid(alpha=.3)
plt.title("Degree distribution (log-log)"); plt.tight_layout(); plt.savefig(sys.argv[1][:-4] + ".png", dpi=150)
