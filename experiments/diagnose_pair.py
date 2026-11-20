"""Diagnostikë: pse mungon overlap-i read_13 <-> read_735 (seed 43, R=0).
    python -m experiments.diagnose_pair
"""
from src import config as C
from src.io_utils import reverse_complement
from src.kmer_index import build_index, find_candidates
from src.overlap import compute_overlaps
from src.graph import find_contained
from src.simulator import simulate_reads
from experiments.repeats import make_genome

SEED = 43
A, B = 13, 735


def edit_distance(s, t):
    prev = list(range(len(t) + 1))
    for i in range(1, len(s) + 1):
        cur = [i] + [0] * len(t)
        for j in range(1, len(t) + 1):
            cost = 0 if s[i - 1] == t[j - 1] else 1
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + cost)
        prev = cur
    return prev[-1]


def forward_seq(read):
    """Kthen sekuencën në orientimin e gjenomit, sipas strand-it në emër."""
    name, seq = read[0], read[1]
    return reverse_complement(seq) if name.endswith("_-") else seq


def main():
    genome = make_genome(0, SEED)
    reads = simulate_reads(genome, coverage=30, error_rate=0.01, seed=SEED)
    print("Emrat:", reads[A][0], "|", reads[B][0])

    index, _ = build_index(reads, k=C.K, max_occ=C.MAX_OCC)
    cands, _ = find_candidates(reads, index, k=C.K, min_shared=C.MIN_SHARED,
                               min_overlap=C.MIN_OVERLAP, band=C.BAND)
    overlaps, _ = compute_overlaps(reads, cands, min_overlap=C.MIN_OVERLAP,
                                   max_err=C.MAX_ERR, band=C.ALIGN_BAND)
    contained = find_contained(reads, overlaps)
    print("read_%d i përmbajtur? %s | read_%d i përmbajtur? %s"
          % (A, A in contained, B, B in contained))

    pair = [o for o in overlaps if {o["a"], o["b"]} == {A, B}]
    print("Overlap-et e gjetura mes dy reads:", pair)

    for r in (B, A):
        print("\nOverlap-et e gjetura që përfshijnë read_%d (%s):" % (r, reads[r][0]))
        for o in overlaps:
            if r in (o["a"], o["b"]):
                other = o["b"] if o["a"] == r else o["a"]
                print("   me %s | diag=%s length=%s score=%s same=%s"
                      % (reads[other][0], o["diag"], o["length"], o["score"], o["same"]))

    # gabimet reale në rajonin e përbashkët 8888..9009
    sA, sB = 8879, 8888
    fa = forward_seq(reads[A])[sB - sA:]            # nga 8888
    fb = forward_seq(reads[B])
    L = min(len(fa), len(fb), 9009 - 8888)
    fa, fb = fa[:L], fb[:L]
    ref = genome[8888:8888 + L]
    print("\nRajoni i përbashkët (~%d bp):" % L)
    print("  edit(read_13, gjenomi) =", edit_distance(fa, ref))
    print("  edit(read_735, gjenomi) =", edit_distance(fb, ref))
    d = edit_distance(fa, fb)
    print("  edit(read_13, read_735) = %d  -> raporti %.3f (MAX_ERR=%s)" % (d, d / L, C.MAX_ERR))


if __name__ == "__main__":
    main()