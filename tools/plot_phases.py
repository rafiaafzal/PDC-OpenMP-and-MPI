#!/usr/bin/env python3
"""CS5: stacked phase-time bar chart from RESULT lines.
usage: python3 plot_phases.py results/cs5_phases.log [dist]   (default dist=0)"""
import sys, re
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
dist = sys.argv[2] if len(sys.argv) > 2 else "0"
data = {}
for line in open(sys.argv[1]):
    d = dict(re.findall(r"(\w+)=(\S+)", line))
    if d.get("variant") == "samplesort" and d.get("dist") == dist and d.get("oversample", "1") == "1":
        data.setdefault(int(d["p"]), []).append(d)
ps = sorted(data); keys = ["t_sort", "t_sample", "t_counts", "t_alltoallv", "t_merge"]
vals = {k: [sorted(float(x[k]) for x in data[p])[len(data[p]) // 2] for p in ps] for k in keys}
bot = [0] * len(ps)
for k in keys:
    plt.bar([str(p) for p in ps], vals[k], bottom=bot, label=k); bot = [a + b for a, b in zip(bot, vals[k])]
plt.xlabel("MPI ranks p"); plt.ylabel("time (s), max over ranks, median of runs"); plt.legend()
plt.title(f"CS5 phase breakdown (dist={dist})"); plt.tight_layout(); plt.savefig("results/cs5_phases_dist%s.png" % dist, dpi=150)
print("| p | " + " | ".join(keys) + " |\n|" + "---|" * (len(keys) + 1))
for i, p in enumerate(ps): print(f"| {p} | " + " | ".join(f"{vals[k][i]:.4f}" for k in keys) + " |")
