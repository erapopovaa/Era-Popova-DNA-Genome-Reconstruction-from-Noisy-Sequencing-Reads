"""Diagnostikim: pse contig-u fragmentohet.
    python -m experiments.diagnose phix30
    python -m experiments.diagnose phix50
    python -m experiments.diagnose scale5000
"""
import random
import sys
from collections import defaultdict

from src import config as C
from src.io_utils import reverse_complement
from src.simulator import mutate_genome, simulate_reads
from src.kmer_index import build_index, find_candidates
from src.overlap import compute_overlaps
from src.graph import find_contained, build_graph, transitive_reduction, graph_stats
from src.assemble import simplify, build_contigs
from src.evaluate import _read_info
from experiments.benchmark import phix_data, SEED


def get_data(name):
    if name == "phix30":
        ref, mutated, truth, reads = phix_data(coverage=30)
    elif name == "phix50":
        ref, mutated, truth, reads = phix_data(coverage=50)
    elif name.startswith("scale"):
        n = int(name[5:])
        g_len = int(n * 125 / 30)
        rng = random.Random(SEED)
        genome = "".join(rng.choice("ACGT") for _ in range(g_len))
        mutated, truth = mutate_genome(genome, n_snp=max(10, g_len // 500),
                                       n_ins=max(3, g_len // 1500),
                                       n_del=max(3, g_len // 1500), seed=SEED)
        reads = simulate_reads(mutated, coverage=30, seed=SEED)
    else:
        raise SystemExit("Emër i panjohur: " + name)
    return mutated, reads


def locate(seq, genome, k=25):
    """Kthen (strand, fillimi, fundi) i contig-ut në gjenom, ose None."""
    for strand, s in (("+", seq), ("-", reverse_complement(seq))):
        if len(s) < k:
            continue
        for i in range(0, len(s) - k + 1, 10):
            idx = genome.find(s[i:i + k])
            if idx != -1:
                return strand, idx - i, idx - i + len(s)
    return None


def true_coverage(reads, glen):
    """Coverage i vërtetë për çdo pozicion, nga pozicionet në emrat e reads."""
    diff = [0] * (glen + 2)
    for start, length, _ in _read_info(reads).values():
        diff[start] += 1
        diff[min(glen, start + length)] -= 1
    cov, cur = [], 0
    for i in range(glen):
        cur += diff[i]
        cov.append(cur)
    return cov


def edges_of(g):
    return {(u, v) for u, outs in g.items() for v in outs}


def main(name):
    mutated, reads = get_data(name)
    glen = len(mutated)
    lens = [len(r[1]) for r in reads]
    print("Dataset:", name, "| gjenomi:", glen, "| reads:", len(reads))

    cov = true_coverage(reads, glen)
    print("\n--- Mbulimi i vërtetë (nga simulatori) ---")
    print("  mesatar: %.1f | minimal: %d" % (sum(cov) / glen, min(cov)))
    print("  pozicione me coverage 0:", sum(1 for c in cov if c == 0))
    print("  pozicione me coverage < 3:", sum(1 for c in cov if c < 3))

    max_err = max(C.MAX_ERR, 4 * C.ERROR_RATE)
    index, _ = build_index(reads, k=C.K, max_occ=C.MAX_OCC)
    cands, cst = find_candidates(reads, index, k=C.K, min_shared=C.MIN_SHARED,
                                 min_overlap=C.MIN_OVERLAP, band=C.BAND)
    overlaps, _ = compute_overlaps(reads, cands, min_overlap=C.MIN_OVERLAP,
                                   max_err=max_err, band=C.ALIGN_BAND)
    contained = find_contained(reads, overlaps)
    active = [r for r in range(len(reads)) if r not in contained]
    graph = build_graph(reads, overlaps, contained)
    snap0 = edges_of(graph)
    transitive_reduction(graph, lens, fuzz=C.FUZZ)
    snap1 = edges_of(graph)

    def indegrees():
        ind = defaultdict(int)
        for outs in graph.values():
            for v in outs:
                ind[v] += 1
        return ind

    def degree_report(label):
        ind = indegrees()
        nodes = [2 * r + s for r in active for s in (0, 1)]
        out_multi = sum(1 for u in nodes if len(graph.get(u, {})) > 1)
        in_multi = sum(1 for u in nodes if ind[u] > 1)
        dead = sum(1 for u in nodes if len(graph.get(u, {})) == 0)
        print("  %s: nyje=%d | dalje>1: %d | hyrje>1: %d | pa dalje: %d" %
              (label, len(nodes), out_multi, in_multi, dead))

    print("\n--- Grafi ---")
    print("  reads aktive (jo të përmbajtura):", len(active), "nga", len(reads))
    degree_report("pas transitive reduction")
    simplify(graph, active, C.TIP_READS, C.BUBBLE_LEN, C.SIMPLIFY_ROUNDS)
    snap2 = edges_of(graph)
    degree_report("pas tips/bubbles       ")
    print("  komponente:", graph_stats(graph, len(reads), contained)["connected_components"])

    contigs = build_contigs(reads, graph, active)
    ind = indegrees()
    lens_c = [len(s) for _, s in contigs]
    print("\n--- Contigs ---")
    print("  gjithsej:", len(contigs), "| me gjatësi >= 300:", sum(1 for l in lens_c if l >= 300),
          "| gjatësia totale:", sum(lens_c), "| gjenomi:", glen)

    rows = []
    for path, seq in contigs:
        if len(seq) < 300:
            continue
        loc = locate(seq, mutated)
        rows.append((loc[1] if loc else -1, loc[2] if loc else -1,
                     loc[0] if loc else "?", len(seq), path))
    rows.sort()
    print("\n  #   gjatësia strand   fillimi     fundi  hyrje(fillim) dalje(fund)")
    for i, (s, e, st, L, path) in enumerate(rows[:40]):
        print("  %-3d %-9d %-6s %-10d %-9d %-13d %d" %
              (i, L, st, s, e, ind[path[0]], len(graph.get(path[-1], {}))))

    print("\n--- Boshllëqet midis contigs të njëpasnjëshme ---")
    mapped = [r for r in rows if r[0] >= 0]
    for a, b in zip(mapped, mapped[1:]):
        gap = b[0] - a[1]
        lo, hi = max(0, min(a[1], b[0]) - 100), min(glen, max(a[1], b[0]) + 100)
        window = cov[lo:hi] or [0]
        print("  contig fund %-7d -> fillim %-7d | boshllëk: %-6d | coverage i vërtetë aty: "
              "min=%d mesatar=%.1f" % (a[1], b[0], gap, min(window), sum(window) / len(window)))
    if len(mapped) < 2:
        print("  (një contig i vetëm i gjatë)")

    print("\n--- Nyjet me degëzim të mbetur ---")
    info = _read_info(reads)

    def describe(x):
        r, s = x >> 1, x & 1
        if r not in info:
            return "read %d (dublikatë)" % r
        st, ln, strand = info[r]
        fwd = (strand == "+") != (s == 1)
        return "read %d @%d-%d %s" % (r, st, st + ln, "fwd" if fwd else "rev")

    def run_length(x):
        n, seen = 1, {x >> 1}
        while len(graph.get(x, {})) == 1:
            y = next(iter(graph[x]))
            if ind[y] != 1 or (y >> 1) in seen:
                break
            n += 1
            seen.add(y >> 1)
            x = y
        return n, len(graph.get(x, {}))

    shown = 0
    for u in sorted(graph):
        outs = graph[u]
        if len(outs) > 1 and shown < 6:
            shown += 1
            print("  nyja", describe(u), "| hyrje:", ind[u], "| dalje:", len(outs))
            for v, (ov, sc) in outs.items():
                n, d = run_length(v)
                print("     ->", describe(v),
                      "| overlap=%d score=%.1f | rruga pa degëzim: %d reads, pastaj dalje=%d"
                      % (ov, sc, n, d))

    print("\n--- Ku humbet skaji ---")
    print("\n--- Pse refuzohen 273-1421 dhe 1067-1239 ---")
    from src.overlap import align_overlap
    cand_map = {}
    for a, b, same, diag, votes, ovl in cands:
        cand_map[frozenset((a, b))] = (a, b, same, diag, votes, ovl)
    for a, b in [(273, 1421), (1067, 1239)]:
        c = cand_map.get(frozenset((a, b)))
        if c is None:
            print("  reads %d-%d | NUK është kandidat (indeksi nuk e gjeti)" % (a, b))
            continue
        ca, cb, same, diag, votes, ovl = c
        sa = reads[ca][1]
        sb = reads[cb][1] if same else reverse_complement(reads[cb][1])
        for w in (8, 16):
            ed, length, dp = align_overlap(sa, sb, diag, w)
            err = (ed / length) if ed is not None and length else None
            print("  reads %d-%d | kandidat: votat=%d diag=%d | band=%d -> ed=%s gjatësia=%s err=%s"
                  % (a, b, votes, diag, w, ed, length, None if err is None else round(err, 3)))
    for a, b in [(273, 1421), (1067, 1239), (971, 1421), (198, 1239)]:
        if a >= len(reads) or b >= len(reads):
            continue
        ov = [o for o in overlaps if {o["a"], o["b"]} == {a, b}]
        ovs = ("length=%d err=%.3f" % (ov[0]["length"], ov[0]["err"])) if ov \
            else "NUK ka overlap të pranuar"

        def has(snap):
            return any({u >> 1, v >> 1} == {a, b} for u, v in snap)

        print("  reads %d-%d | %s | të përmbajtura: %s/%s | skaj: para TR=%s, pas TR=%s, "
              "pas thjeshtimit=%s"
              % (a, b, ovs, a in contained, b in contained, has(snap0), has(snap1), has(snap2)))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Përdorimi: python -m experiments.diagnose phix30|phix50|scale5000")
    main(sys.argv[1])