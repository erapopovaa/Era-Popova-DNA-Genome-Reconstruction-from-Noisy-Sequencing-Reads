import time

from src.config import (K, MAX_OCC, MIN_SHARED, BAND, MIN_OVERLAP,
                        MAX_ERR, ALIGN_BAND, FUZZ)
from src.io_utils import read_fastq
from src.kmer_index import build_index, find_candidates
from src.overlap import compute_overlaps
from src.graph import (find_contained, build_graph, count_edges,
                       transitive_reduction, graph_stats)
from src.evaluate import graph_edge_accuracy

reads = list(read_fastq("data/reads/sim_phix174.fastq"))
lens = [len(r[1]) for r in reads]
print("Reads:", len(reads))

t0 = time.time()
index, _ = build_index(reads, k=K, max_occ=MAX_OCC)
candidates, _ = find_candidates(reads, index, k=K, min_shared=MIN_SHARED,
                                min_overlap=MIN_OVERLAP, band=BAND)
overlaps, _ = compute_overlaps(reads, candidates, min_overlap=MIN_OVERLAP,
                               max_err=MAX_ERR, band=ALIGN_BAND)
t1 = time.time()
print("Overlap-e të pranuara:", len(overlaps), "| koha:", round(t1 - t0, 2), "s")

contained = find_contained(reads, overlaps)
print("\n--- Containment ---")
print("  reads të përmbajtura (hiqen nga grafi):", len(contained),
      "(%.1f%%)" % (100 * len(contained) / len(reads)))

graph = build_graph(reads, overlaps, contained)
edges_before = count_edges(graph)
acc_before = graph_edge_accuracy(reads, graph)

t2 = time.time()
removed = transitive_reduction(graph, lens, fuzz=FUZZ)
t3 = time.time()
edges_after = count_edges(graph)
acc_after = graph_edge_accuracy(reads, graph)

print("\n--- Transitive reduction ---")
print("  skaje para:", edges_before, "| pas:", edges_after,
      "| hequr:", removed, "(%.1f%%)" % (100 * removed / max(1, edges_before)))
print("  saktësia e skajeve para:", round(acc_before, 4), "| pas:", round(acc_after, 4))
print("  koha e reduktimit:", round(t3 - t2, 3), "s")

print("\n--- Grafi final ---")
for key, val in graph_stats(graph, len(reads), contained).items():
    print(" ", key, ":", round(val, 3) if isinstance(val, float) else val)