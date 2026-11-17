import time

from src.config import (K, MAX_OCC, MIN_SHARED, BAND, MIN_OVERLAP, MAX_ERR,
                        ALIGN_BAND, FUZZ, TIP_READS, BUBBLE_LEN, SIMPLIFY_ROUNDS)
from src.io_utils import read_fastq, read_fasta, write_fasta
from src.kmer_index import build_index, find_candidates
from src.overlap import compute_overlaps
from src.graph import (find_contained, build_graph, count_edges,
                       transitive_reduction, graph_stats)
from src.assemble import simplify, build_contigs, contig_stats
from src.greedy import greedy_contigs
from src.evaluate import contig_vs_reference

reads = list(read_fastq("data/reads/sim_phix174.fastq"))
ref = list(read_fasta("data/reference/phix174_mutated.fasta").values())[0]
lens = [len(r[1]) for r in reads]
print("Reads:", len(reads), "| gjenomi i referencës (i mutuar):", len(ref))

index, _ = build_index(reads, k=K, max_occ=MAX_OCC)
candidates, _ = find_candidates(reads, index, k=K, min_shared=MIN_SHARED,
                                min_overlap=MIN_OVERLAP, band=BAND)
overlaps, _ = compute_overlaps(reads, candidates, min_overlap=MIN_OVERLAP,
                               max_err=MAX_ERR, band=ALIGN_BAND)
contained = find_contained(reads, overlaps)
active = [r for r in range(len(reads)) if r not in contained]
graph = build_graph(reads, overlaps, contained)
greedy_graph = {u: dict(o) for u, o in graph.items()}   # kopje para reduktimit

transitive_reduction(graph, lens, fuzz=FUZZ)
print("\nGrafi pas transitive reduction: skaje =", count_edges(graph))

t0 = time.time()
sstats = simplify(graph, active, TIP_READS, BUBBLE_LEN, SIMPLIFY_ROUNDS)
olc = build_contigs(reads, graph, active)
t1 = time.time()
gre = greedy_contigs(reads, greedy_graph, active)
t2 = time.time()

print("\n--- Thjeshtimi i grafit ---")
print("  tips të hequra:", sstats["tips_removed"], "| bubbles të hequra:", sstats["bubbles_popped"])
for key, val in graph_stats(graph, len(reads), contained).items():
    print(" ", key, ":", round(val, 3) if isinstance(val, float) else val)


def report(name, contigs, seconds):
    st = contig_stats(contigs)
    used = st["reads_in_contigs"] + len(contained)
    print("\n--- %s ---" % name)
    print("  koha:", round(seconds, 3), "s")
    print("  contigs:", st["contigs"], "| gjatesia totale:", st["total_length"],
          "| me i gjati:", st["longest"], "| N50:", st["n50"])
    print("  5 me te gjatet:", [len(s) for _, s in contigs[:5]])
    print("  reads te perdorura: %d (%.1f%%)" % (used, 100 * used / len(reads)))
    q = contig_vs_reference(contigs[0][1], ref)
    if q is None:
        print("  contig-u me i gjate: nuk u vendos dot ne referenca")
    else:
        print("  contig-u me i gjate: edit distance =", q["edit_distance"],
              "| identiteti =", round(q["identity"], 4),
              "| strand:", q["strand"], "| pozicioni:", q["offset"])


report("OLC (string graph + tips/bubbles)", olc, t1 - t0)
report("Greedy (baseline)", gre, t2 - t1)

write_fasta("results/contigs_olc.fasta",
            {"contig_%d" % (i + 1): s for i, (_, s) in enumerate(olc[:20])})
write_fasta("results/contigs_greedy.fasta",
            {"contig_%d" % (i + 1): s for i, (_, s) in enumerate(gre[:20])})