"""Eksperimenti me repeats: OLC kundrejt greedy në gjenom me segment të kopjuar 3 herë.
    python -m experiments.repeats
"""
import csv
import random

from src.simulator import simulate_reads
from experiments.benchmark import run

GENOME_LEN = 10000
REPEAT_LENS = [0, 80, 300, 1000]
SEEDS = [42, 43, 44]
METHODS = ["olc", "olc_best", "greedy"]


def make_genome(repeat_len, seed):
    """Gjenom i rastësishëm; segmenti 2000..2000+R kopjohet te 5000 dhe 8000."""
    rng = random.Random(seed)
    g = [rng.choice("ACGT") for _ in range(GENOME_LEN)]
    if repeat_len > 0:
        seg = g[2000:2000 + repeat_len]
        for pos in (5000, 8000):
            g[pos:pos + repeat_len] = seg
    return "".join(g)


def mean(xs):
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else 0.0


def main():
    rows = []
    for r_len in REPEAT_LENS:
        for method in METHODS:
            reps = []
            for seed in SEEDS:
                genome = make_genome(r_len, seed)
                reads = simulate_reads(genome, coverage=30, error_rate=0.01, seed=seed)
                rep = run(reads, genome, None, None, method=method)
                reps.append(rep)
                print("R=%-5d %-7s seed=%d -> contigs=%s longest=%s fraction=%s error=%s"
                      % (r_len, method, seed, rep.get("contigs"), rep.get("longest_contig"),
                         rep.get("genome_fraction"), rep.get("error")))
            rows.append({
                "repeat_len": r_len,
                "method": method,
                "contigs": round(mean([x.get("contigs") for x in reps]), 1),
                "longest_contig": round(mean([x.get("longest_contig") for x in reps])),
                "genome_fraction": round(mean([x.get("genome_fraction") for x in reps]), 4),
                "edit_distance": round(mean([x.get("edit_distance") for x in reps]), 1),
                "time_seconds": round(mean([x.get("time_seconds") for x in reps]), 1),
            })

    cols = ["repeat_len", "method", "contigs", "longest_contig", "genome_fraction",
            "edit_distance", "time_seconds"]
    with open("results/benchmark_repeats.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in rows:
            w.writerow([r[c] for c in cols])

    print("\nMesataret mbi %d seed-a (genome_fraction = contig-u me i gjate / gjenomi)" % len(SEEDS))
    print(" | ".join(cols))
    for r in rows:
        print(" | ".join(str(r[c]) for c in cols))
    print("\nU ruajt: results/benchmark_repeats.csv")


if __name__ == "__main__":
    main()