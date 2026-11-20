"""Eksperimenti me 10 seed-a për dy konfigurimet e zhurmshme (reads të shkurtra, R=300 dhe R=1000).
Jep mesataren dhe devijimin standard të metrikave të mbulimit.
    python -m experiments.repeats_seeds
"""
import csv
import statistics

from src.simulator import simulate_reads
from experiments.benchmark import run
from experiments.repeats import make_genome

SEEDS = list(range(42, 52))
METHODS = ["olc", "olc_best", "greedy"]
REPEAT_LENS = [300, 1000]
METRICS = ["covered_fraction", "contigs_bad", "n50", "na50", "contigs"]


def stats(xs):
    xs = [x for x in xs if x is not None]
    if not xs:
        return 0.0, 0.0
    sd = statistics.stdev(xs) if len(xs) > 1 else 0.0
    return statistics.mean(xs), sd


def main():
    rows = []
    for r_len in REPEAT_LENS:
        for method in METHODS:
            reps = []
            for seed in SEEDS:
                genome = make_genome(r_len, seed)
                reads = simulate_reads(genome, coverage=30, error_rate=0.01, seed=seed)
                rep = run(reads, genome, None, None, method=method, eval_contigs=True)
                reps.append(rep)
                print("R=%-5d %-8s seed=%d -> contigs=%s covered=%s bad=%s n50=%s error=%s"
                      % (r_len, method, seed, rep.get("contigs"),
                         rep.get("covered_fraction"), rep.get("contigs_bad"),
                         rep.get("n50"), rep.get("error")))
            row = {"repeat_len": r_len, "method": method, "n_seeds": len(SEEDS)}
            for m in METRICS:
                mu, sd = stats([x.get(m) for x in reps])
                row[m + "_mean"] = round(mu, 4)
                row[m + "_sd"] = round(sd, 4)
            # sa seed-a e rindërtuan gjenomin pothuajse plotësisht (mbulim i saktë >= 0.95)
            row["seeds_covered_ge_095"] = sum(
                1 for x in reps if (x.get("covered_fraction") or 0) >= 0.95)
            rows.append(row)

    cols = ["repeat_len", "method", "n_seeds"]
    for m in METRICS:
        cols += [m + "_mean", m + "_sd"]
    cols.append("seeds_covered_ge_095")

    with open("results/benchmark_repeats_seeds.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in rows:
            w.writerow([r[c] for c in cols])

    print("\nMesatarja dhe devijimi standard mbi %d seed-a (reads 100-150 bp)" % len(SEEDS))
    print(" | ".join(cols))
    for r in rows:
        print(" | ".join(str(r[c]) for c in cols))
    print("\nU ruajt: results/benchmark_repeats_seeds.csv")


if __name__ == "__main__":
    main()