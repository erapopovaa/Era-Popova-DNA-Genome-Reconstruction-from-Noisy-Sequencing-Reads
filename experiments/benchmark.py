"""Eksperimentet. Përdorimi (nga root-i i projektit):
    python -m experiments.benchmark ablation
    python -m experiments.benchmark scale | error | k | coverage
"""
import csv
import random
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from src import config as C
from src.io_utils import read_fasta
from src.simulator import mutate_genome, simulate_reads
from src.pipeline import run_pipeline

SEED = 42
REF_PATH = "data/reference/phix174.fasta"
COLS = ["param", "reads_total", "candidate_comparisons", "brute_force_comparisons",
        "overlaps_identified", "longest_contig", "reconstruction_length",
        "identity", "genome_fraction", "edit_distance", "pct_reads_incorporated",
        "tp", "fp", "fn", "time_seconds", "peak_memory_mb"]


def phix_data(coverage=C.COVERAGE, error=C.ERROR_RATE):
    """phiX174 real + 16 mutacione të njohura + reads të simuluara."""
    ref = list(read_fasta(REF_PATH).values())[0]
    mutated, truth = mutate_genome(ref, n_snp=10, n_ins=3, n_del=3, seed=SEED)
    reads = simulate_reads(mutated, coverage=coverage, error_rate=error, seed=SEED)
    return ref, mutated, truth, reads


def run(reads, ref, mutated, truth, sim_error=C.ERROR_RATE, **kw):
    max_err = max(C.MAX_ERR, 4 * sim_error)  # overlap-i midis dy reads ka gabimet e të dyve
    report, _, _ = run_pipeline(reads, ref, mutated, truth, max_err=max_err, **kw)
    return report


def show(rep):
    print("  ->", {k: rep.get(k) for k in
                   ("reads_total", "longest_contig", "identity", "tp", "fp", "fn",
                    "time_seconds", "error")})


def save(name, rows):
    path = "results/benchmark_%s.csv" % name
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(COLS)
        for r in rows:
            w.writerow([r.get(c, "") for c in COLS])
    print("\nU ruajt:", path)
    print(" | ".join(COLS))
    for r in rows:
        print(" | ".join(str(r.get(c, "")) for c in COLS))


def plot_params(name, rows, xlabel):
    labels = [str(r["param"]) for r in rows]
    x = range(len(rows))
    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    ax[0].plot(x, [100 * (r.get("identity") or 0) for r in rows], "o-")
    ax[0].set_ylabel("Identiteti (%)")
    ax[1].plot(x, [r.get("tp", 0) for r in rows], "o-", label="TP")
    ax[1].plot(x, [r.get("fp", 0) for r in rows], "s-", label="FP")
    ax[1].plot(x, [r.get("fn", 0) for r in rows], "^-", label="FN")
    ax[1].set_ylabel("Mutacione"); ax[1].legend()
    ax[2].plot(x, [r.get("candidate_comparisons", 0) for r in rows], "o-")
    ax[2].set_ylabel("Kandidatë të krahasuar")
    for a in ax:
        a.set_xticks(list(x)); a.set_xticklabels(labels)
        a.set_xlabel(xlabel); a.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("results/benchmark_%s.png" % name, dpi=150)
    print("U ruajt: results/benchmark_%s.png" % name)


def sweep(name, params, xlabel, make):
    """make(param) -> (raporti)."""
    rows = []
    for p in params:
        print("\n%s = %s" % (xlabel, p))
        rep = make(p)
        rep["param"] = p
        show(rep)
        rows.append(rep)
    save(name, rows)
    plot_params(name, rows, xlabel)


def exp_scale():
    rows = []
    for n in [1000, 2000, 5000, 10000, 20000]:
        g_len = int(n * 125 / 30)
        rng = random.Random(SEED)
        genome = "".join(rng.choice("ACGT") for _ in range(g_len))
        mutated, truth = mutate_genome(genome, n_snp=max(10, g_len // 500),
                                       n_ins=max(3, g_len // 1500),
                                       n_del=max(3, g_len // 1500), seed=SEED)
        reads = simulate_reads(mutated, coverage=30, seed=SEED)
        print("\nN afërsisht %d | gjenomi %d bp | reads %d" % (n, g_len, len(reads)))
        rep = run(reads, genome, mutated, truth)
        rep["param"] = len(reads)
        show(rep)
        rows.append(rep)
    save("scale", rows)

    ns = [r["reads_total"] for r in rows]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    ax[0].plot(ns, [r["time_seconds"] for r in rows], "o-")
    ax[0].set_ylabel("Koha (s)")
    ax[1].plot(ns, [r.get("peak_memory_mb", 0) for r in rows], "o-")
    ax[1].set_ylabel("Memoria maksimale (MB)")
    ax[2].loglog(ns, [r["brute_force_comparisons"] for r in rows], "o--", label="Brute force N(N-1)/2")
    ax[2].loglog(ns, [r["candidate_comparisons"] for r in rows], "o-", label="Me k-mer index")
    ax[2].set_ylabel("Krahasime çiftesh"); ax[2].legend()
    for a in ax:
        a.set_xlabel("Numri i reads (N)"); a.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("results/benchmark_scale.png", dpi=150)
    print("U ruajt: results/benchmark_scale.png")


def exp_error():
    def make(e):
        ref, mutated, truth, reads = phix_data(error=e)
        return run(reads, ref, mutated, truth, sim_error=e)
    sweep("error", [0.001, 0.01, 0.02, 0.05], "Gabimi i sekuencimit", make)


def exp_k():
    ref, mutated, truth, reads = phix_data()
    sweep("k", [11, 15, 21, 31], "k",
          lambda k: run(reads, ref, mutated, truth, k=k))


def exp_coverage():
    def make(c):
        ref, mutated, truth, reads = phix_data(coverage=c)
        return run(reads, ref, mutated, truth)
    sweep("coverage", [10, 20, 30, 50], "Coverage (x)", make)


def exp_ablation():
    ref, mutated, truth, reads = phix_data()
    configs = [("OLC i plotë", {}),
               ("pa transitive reduction", {"reduce_graph": False}),
               ("pa quality scores", {"use_quality": False}),
               ("greedy (baseline)", {"method": "greedy"})]
    rows = []
    for name, kw in configs:
        print("\n" + name)
        rep = run(reads, ref, mutated, truth, **kw)
        rep["param"] = name
        show(rep)
        rows.append(rep)
    save("ablation", rows)


EXPERIMENTS = {"scale": exp_scale, "error": exp_error, "k": exp_k,
               "coverage": exp_coverage, "ablation": exp_ablation}

if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in EXPERIMENTS:
        print("Përdorimi: python -m experiments.benchmark", "|".join(EXPERIMENTS))
        sys.exit(1)
    EXPERIMENTS[sys.argv[1]]()