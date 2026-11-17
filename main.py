import json
import os
import time

import psutil

from src.config import (K, MAX_OCC, MIN_SHARED, BAND, MIN_OVERLAP, MAX_ERR, ALIGN_BAND,
                        FUZZ, TIP_READS, BUBBLE_LEN, SIMPLIFY_ROUNDS, POLISH_ROUNDS,
                        POLISH_BAND, MIN_ALT, MIN_ALT_FRAC, SEQ_ERROR, PVALUE_MAX)
from src.io_utils import read_fastq, read_fasta, write_fasta
from src.kmer_index import build_index, find_candidates
from src.overlap import compute_overlaps
from src.graph import find_contained, build_graph, transitive_reduction
from src.assemble import simplify, build_contigs, contig_stats
from src.consensus import (polish, orient_contig, call_variants, minor_alleles,
                           edit_vs_reference)
from src.evaluate import read_truth, compare_variants,merge_adjacent_indels

T_START = time.time()

reads = list(read_fastq("data/reads/sim_phix174.fastq"))
orig = list(read_fasta("data/reference/phix174.fasta").values())[0]
mutated = list(read_fasta("data/reference/phix174_mutated.fasta").values())[0]
truth = read_truth("data/reference/phix174_truth.tsv")
lens = [len(r[1]) for r in reads]
print("Reads:", len(reads))

# 1-3. Indeksi, kandidatët, overlap-et
index, _ = build_index(reads, k=K, max_occ=MAX_OCC)
candidates, cstats = find_candidates(reads, index, k=K, min_shared=MIN_SHARED,
                                     min_overlap=MIN_OVERLAP, band=BAND)
overlaps, ostats = compute_overlaps(reads, candidates, min_overlap=MIN_OVERLAP,
                                    max_err=MAX_ERR, band=ALIGN_BAND)
print("Kandidatë:", cstats["candidate_pairs"], "| overlap-e:", len(overlaps))

# 4. Grafi dhe rindërtimi i contig-ut
contained = find_contained(reads, overlaps)
active = [r for r in range(len(reads)) if r not in contained]
graph = build_graph(reads, overlaps, contained)
transitive_reduction(graph, lens, fuzz=FUZZ)
simplify(graph, active, TIP_READS, BUBBLE_LEN, SIMPLIFY_ROUNDS)
contigs = build_contigs(reads, graph, active)
contig = contigs[0][1]
t_asm = time.time()
print("Contig-u kryesor:", len(contig), "baza | koha deri këtu: %.1f s" % (t_asm - T_START))

# 5. Konsensusi
oriented, strand, off = orient_contig(contig, orig)
ed_before = edit_vs_reference(oriented, mutated)
print("\n--- Para konsensusit ---")
print("  strand:", strand, "| edit distance kundrejt gjenomit të mutuar:", ed_before)

t0 = time.time()
cons, pile, pst = polish(reads, oriented, rounds=POLISH_ROUNDS, k=K, w=POLISH_BAND)
t1 = time.time()
ed_after = edit_vs_reference(cons, mutated)
print("\n--- Pas konsensusit ---")
print("  gjatësia:", len(cons), "| koha: %.1f s" % (t1 - t0))
print("  edit distance kundrejt gjenomit të mutuar:", ed_after,
      "| identiteti: %.5f" % (1 - ed_after / len(cons)))
print("  reads të alinjuara: %d / %d (%.1f%%)" %
      (pst["reads_aligned"], pst["reads_total"], pst["pct_aligned"]))
print("  coverage mesatar: %.1f | gabim mesatar i reads: %.4f" %
      (pst["mean_coverage"], pst["avg_read_error"]))

# 6. Variantet
called, filtered, ed_ref = call_variants(cons, orig, pile, offset=max(0, off),
                                         min_alt=MIN_ALT, min_frac=MIN_ALT_FRAC,
                                         eps=SEQ_ERROR, max_p=PVALUE_MAX)
called = merge_adjacent_indels(called)
minors = minor_alleles(cons, pile, MIN_ALT, MIN_ALT_FRAC, SEQ_ERROR, PVALUE_MAX)
res = compare_variants(called, truth)

print("\n--- Variantet ---")
print("  edit distance konsensus -> referencë origjinale:", ed_ref)
print("  të pranuara:", len(called), "| të filtruara (dëshmi e dobët):", len(filtered),
      "| aleleshë minore:", len(minors))
print("  POS    TYPE  REF  ALT  SUPPORT  COVERAGE  PVALUE")
for v in called:
    print("  %-6d %-5s %-4s %-4s %-8d %-9d %.1e" % v)
print("\n  TP: %d | FP: %d | FN: %d  (nga %d mutacione të vërteta)" %
      (res["tp"], res["fp"], res["fn"], len(truth)))
if res["missed"]:
    print("  humbur:", res["missed"])
if res["false_positives"]:
    print("  false positives:", [v[:4] for v in res["false_positives"]])

# Skedarët dalës
os.makedirs("results", exist_ok=True)
write_fasta("results/reconstructed.fasta", {"reconstructed_phix174": cons})
with open("results/variants.tsv", "w") as f:
    f.write("POS\tTYPE\tREF\tALT\tSUPPORT\tCOVERAGE\tPVALUE\n")
    for v in called:
        f.write("%d\t%s\t%s\t%s\t%d\t%d\t%.2e\n" % v)

mem = psutil.Process(os.getpid()).memory_info()
total_time = time.time() - T_START
report = {
    "reads_total": len(reads),
    "reads_aligned_to_reconstruction": pst["reads_aligned"],
    "pct_reads_incorporated": round(pst["pct_aligned"], 2),
    "brute_force_comparisons": cstats["brute_force_pairs"],
    "candidate_comparisons": cstats["candidate_pairs"],
    "overlaps_identified": len(overlaps),
    "avg_overlap_quality": round(ostats["avg_quality"], 4),
    "reconstruction_length": len(cons),
    "edit_distance_vs_true_genome": ed_after,
    "identity_vs_true_genome": round(1 - ed_after / len(cons), 5),
    "variants_called": len(called),
    "true_positives": res["tp"], "false_positives": res["fp"], "false_negatives": res["fn"],
    "time_seconds": round(total_time, 2),
    "peak_memory_mb": round(getattr(mem, "peak_wset", mem.rss) / 1e6, 1),
}
with open("results/report.json", "w") as f:
    json.dump(report, f, indent=2)
print("\n--- Raporti (results/report.json) ---")
for key, val in report.items():
    print(" ", key, ":", val)