#!/usr/bin/env python3
"""Turn a benchmark CSV (group,label,p,run,time) into speedup/efficiency tables + graphs.

usage: python3 analyze.py results/cs1_strong.csv            # strong scaling (needs label "serial", p=1)
       python3 analyze.py results/cs1_weak.csv --weak        # weak scaling: E(p)=T(1)/T(p) per label
Writes <csv>_table.md and one PNG per group (speedup+efficiency).  Medians over runs are used.
"""
import csv, sys, statistics, collections
path = sys.argv[1]; weak = "--weak" in sys.argv
rows = collections.defaultdict(list)
for r in csv.DictReader(open(path)):
    rows[(r["group"], r["label"], int(r["p"]))].append(float(r["time"]))
med = {k: statistics.median(v) for k, v in rows.items()}
runs = {k: len(v) for k, v in rows.items()}
groups = sorted({k[0] for k in med})
md = []
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
except Exception:
    plt = None; print("matplotlib missing: tables only")
for g in groups:
    labels = sorted({k[1] for k in med if k[0] == g and k[1] != "serial"})
    ts = med.get((g, "serial", 1))
    md.append(f"\n### {g}" + (f"  (serial baseline T = {ts:.4f} s)" if ts and not weak else "") + "\n")
    md.append("| variant | p | median time (s) | runs | " + ("efficiency E=T1/Tp" if weak else "speedup S | efficiency E") + " |")
    md.append("|---|---|---|---|" + ("---|" if weak else "---|---|"))
    if plt: fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    for lab in labels:
        ps = sorted(k[2] for k in med if k[0] == g and k[1] == lab)
        base = med.get((g, lab, 1)) if weak else ts
        if base is None: continue
        S, E = [], []
        for p in ps:
            t = med[(g, lab, p)]
            if weak: e = base / t; E.append(e); md.append(f"| {lab} | {p} | {t:.4f} | {runs[(g,lab,p)]} | {e:.3f} |")
            else:
                s = base / t; S.append(s); E.append(s / p)
                md.append(f"| {lab} | {p} | {t:.4f} | {runs[(g,lab,p)]} | {s:.2f} | {s/p:.3f} |")
        if plt:
            if not weak: ax[0].plot(ps, S, "o-", label=lab)
            else: ax[0].plot(ps, [med[(g, lab, p)] for p in ps], "o-", label=lab)
            ax[1].plot(ps, E, "o-", label=lab)
    if plt:
        pmax = max(k[2] for k in med if k[0] == g)
        if not weak:
            ax[0].plot([1, pmax], [1, pmax], "k--", label="ideal")
            ax[0].set_ylabel("speedup S(p) = T_serial / T_p")
        else:
            ax[0].set_ylabel("time (s)")
        ax[1].axhline(0.7, color="gray", ls=":"); ax[1].set_ylabel("efficiency")
        for a in ax: a.set_xlabel("p (threads or MPI ranks)"); a.set_xscale("log", base=2); a.grid(alpha=.3); a.legend(fontsize=7)
        nr = max(runs[k] for k in runs if k[0] == g)
        fig.suptitle(f"{g}  (median of {nr} runs per point)")
        fig.tight_layout(); fig.savefig(f"{path[:-4]}_{g}.png", dpi=150); plt.close(fig)
open(path[:-4] + "_table.md", "w").write("\n".join(md) + "\n")
print("\n".join(md))
