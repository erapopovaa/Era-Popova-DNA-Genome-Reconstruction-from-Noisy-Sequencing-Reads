"""Shkallëzimi: ekzekuton një madhësi të vetme (p.sh. 50000 reads).
    python -m experiments.scale_big 50000
"""
import csv
import os
import random
import sys

from src.simulator import mutate_genome, simulate_reads
from experiments.benchmark import run

MEAN_LEN = 125
COVERAGE = 30


def main():
    n_reads = int(sys.argv[1]) if len(sys.argv) > 1 else 50000
    glen = int(n_reads * MEAN_LEN / COVERAGE)
    rng = random.Random(42)
    ref = "".join(rng.choice("ACGT") for _ in range(glen))

    total = max(16, glen // 300)
    n_snp = total * 10 // 16
    n_ins = total * 3 // 16
    n_del = total - n_snp - n_ins
    mutated, truth = mutate_genome(ref, n_snp=n_snp, n_ins=n_ins, n_del=n_del, seed=42)
    reads = simulate_reads(mutated, coverage=COVERAGE, error_rate=0.01, seed=42)
    print("genomi=%d bp, reads=%d, variante te vertetat=%d" % (glen, len(reads), len(truth)))

    rep = run(reads, ref, mutated, truth, method="olc")
    for k, v in rep.items():
        print("  %s: %s" % (k, v))

    path = "results/benchmark_scale_big.csv"
    new = not os.path.exists(path)
    with open(path, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["n_reads", "genome_bp"] + list(rep.keys()))
        w.writerow([len(reads), glen] + list(rep.values()))
    print("U shtua te", path)


if __name__ == "__main__":
    main()