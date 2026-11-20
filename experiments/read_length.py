"""Eksperimenti me gjatësinë e reads: a e kapërcejnë repeat-in reads më të gjata?
    python -m experiments.read_length
"""
import csv

from src.simulator import simulate_reads
from experiments.benchmark import run
from experiments.repeats import make_genome, mean

SEEDS = [42, 43, 44]
METHODS = ["olc", "olc_best", "greedy"]
# (repeat_len, min_len, max_len)
CONFIGS = [
    (300, 100, 150),
    (300, 500, 700),
    (1000, 100, 150),
    (1000, 1200, 1500),
]


def main():
    rows = []
    for r_len, lo, hi in CONFIGS:
        for method in METHODS:
            reps = []
            for seed in SEEDS:
                genome = make_genome(r_len, seed)
                reads = simulate_reads(genome, coverage=30, min_len=lo, max_len=hi,
                                       error_rate=0.01, seed=seed)
                rep = run(reads, genome, None, None, method=method, eval_contigs=True)
                reps.append(rep)
                print("R=%-5d reads=%d-%d %-8s seed=%d -> contigs=%s longest=%s "
                      "covered=%s bad=%s error=%s"
                      % (r_len, lo, hi, method, seed, rep.get("contigs"),
                         rep.get("longest_contig"), rep.get("covered_fraction"),
                         rep.get("contigs_bad"), rep.get("error")))
            rows.append({
                "repeat_len": r_len,
                "read_len": "%d-%d" % (lo, hi),
                "method": method,
                "contigs": round(mean([x.get("contigs") for x in reps]), 1),
                "longest_contig": round(mean([x.get("longest_contig") for x in reps])),
                "genome_fraction": round(mean([x.get("genome_fraction") for x in reps]), 4),
                "edit_distance": round(mean([x.get("edit_distance") for x in reps]), 1),
                "time_seconds": round(mean([x.get("time_seconds") for x in reps]), 1),
                "covered_fraction": round(mean([x.get("covered_fraction") for x in reps]), 4),
                "n50": round(mean([x.get("n50") for x in reps])),
                "na50": round(mean([x.get("na50") for x in reps])),
                "contigs_bad": round(mean([x.get("contigs_bad") for x in reps]), 1),
            })

    cols = ["repeat_len", "read_len", "method", "contigs", "longest_contig",
            "genome_fraction", "edit_distance", "time_seconds",
            "covered_fraction", "n50", "na50", "contigs_bad"]
    with open("results/benchmark_readlen.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in rows:
            w.writerow([r[c] for c in cols])

    print("\nMesataret mbi %d seed-a" % len(SEEDS))
    print(" | ".join(cols))
    for r in rows:
        print(" | ".join(str(r[c]) for c in cols))
    print("\nU ruajt: results/benchmark_readlen.csv")


if __name__ == "__main__":
    main()