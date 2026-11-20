"""Diagnostikë: pse olc_best humbet ~1000 bp te seed 43, R=0.
    python -m experiments.diagnose_best
"""
import copy
import re

from src import config as C
from src.io_utils import reverse_complement
from src.kmer_index import build_index, find_candidates
from src.overlap import compute_overlaps
from src.graph import find_contained, build_graph, transitive_reduction
from src.assemble import simplify, build_contigs, resolve_branches
from src.simulator import simulate_reads
from experiments.repeats import make_genome

SEED = 43
PROBE = 25


def true_pos(name):
    m = re.search(r"pos(\d+)_([+-])", name)
    return (int(m.group(1)), m.group(2)) if m else (None, None)


def locate(seq, genome):
    """Gjen (fillim, fund, strand) të contig-ut në gjenom, me probe 25-mer."""
    for strand, s in (("+", seq), ("-", reverse_complement(seq))):
        for off in range(0, len(s) - PROBE, PROBE):
            hit = genome.find(s[off:off + PROBE])
            if hit != -1:
                start = max(hit - off, 0)
                return start, min(start + len(s), len(genome)), strand
    return None


def build(reads, graph, active):
    g = copy.deepcopy(graph)
    simplify(g, list(active), C.TIP_READS, C.BUBBLE_LEN, C.SIMPLIFY_ROUNDS)
    return build_contigs(reads, g, list(active))


def main():
    genome = make_genome(0, SEED)
    reads = simulate_reads(genome, coverage=30, error_rate=0.01, seed=SEED)
    lens = [len(r[1]) for r in reads]

    index, _ = build_index(reads, k=C.K, max_occ=C.MAX_OCC)
    cands, _ = find_candidates(reads, index, k=C.K, min_shared=C.MIN_SHARED,
                               min_overlap=C.MIN_OVERLAP, band=C.BAND)
    overlaps, _ = compute_overlaps(reads, cands, min_overlap=C.MIN_OVERLAP,
                                   max_err=C.MAX_ERR, band=C.ALIGN_BAND)
    contained = find_contained(reads, overlaps)
    active = [r for r in range(len(reads)) if r not in contained]
    graph = build_graph(reads, overlaps, contained)
    transitive_reduction(graph, lens, fuzz=C.FUZZ)

    g_before = copy.deepcopy(graph)
    g_after = copy.deepcopy(graph)
    n_removed = resolve_branches(g_after)
    print("Skaje para: %d | hequr nga resolve_branches: %d"
          % (sum(len(o) for o in g_before.values()), n_removed))

    for label, g in (("olc", g_before), ("olc_best", g_after)):
        contigs = build(reads, g, active)
        print("\n== %s: %d contig-e ==" % (label, len(contigs)))
        for i, c in enumerate(contigs):
            print("  contig %d: gjatësia %d, vendndodhja në gjenom: %s"
                  % (i, len(c[1]), locate(c[1], genome)))

    # pikat e thyerjes: skajet e contig-eve të olc_best (jo skajet e gjenomit)
    breaks = []
    for c in build(reads, g_after, active):
        loc = locate(c[1], genome)
        if loc:
            for b in (loc[0], loc[1]):
                if 50 < b < len(genome) - 50:
                    breaks.append(b)
    print("\nPika thyerjeje (koordinata të gjenomit):", sorted(set(breaks)))

    def name(n):
        r, s = n >> 1, n & 1
        p, st = true_pos(reads[r][0])
        return "%s[len=%d, pos=%s%s, orient=%d]" % (reads[r][0], len(reads[r][1]), p, st, s)

    # 1) simetria e binjakëve pas resolve_branches
    asym = [(u, v) for u, outs in g_after.items() for v in outs
            if (v ^ 1) not in g_after or (u ^ 1) not in g_after[v ^ 1]]
    print("\nSkaje pa binjak pas resolve_branches: %d" % len(asym))

    # 2) çdo skaj i hequr (u->v që ishte para dhe s'është pas)
    removed = [(u, v, e) for u, outs in g_before.items()
               for v, e in outs.items() if v not in g_after.get(u, {})]
    print("Skaje të hequra (duke numëruar edhe binjakët): %d" % len(removed))

    for u, v, e in removed:
        print("\nHEQUR: %s\n    -> %s | (gjatësia, score)=%s" % (name(u), name(v), e))
        print("   daljet e u para:")
        for w, ew in g_before[u].items():
            print("      -> %s | %s" % (name(w), ew))
        print("   hyrjet e v para:")
        for x, outs in g_before.items():
            if v in outs:
                print("      %s -> | %s" % (name(x), outs[v]))
        print("   daljet e u pas: %s" % [name(w) for w in g_after.get(u, {})])


if __name__ == "__main__":
    main()