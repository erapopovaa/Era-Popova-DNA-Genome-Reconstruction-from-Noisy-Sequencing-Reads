import time

from src.config import K, MAX_OCC, MIN_SHARED, BAND, MIN_OVERLAP, MAX_ERR, ALIGN_BAND
from src.io_utils import read_fastq
from src.kmer_index import build_index, find_candidates
from src.overlap import compute_overlaps
from src.evaluate import candidate_quality, overlap_precision

reads = list(read_fastq("data/reads/sim_phix174.fastq"))
print("Reads:", len(reads))

t0 = time.time()
index, istats = build_index(reads, k=K, max_occ=MAX_OCC)
candidates, cstats = find_candidates(reads, index, k=K, min_shared=MIN_SHARED,
                                     min_overlap=MIN_OVERLAP, band=BAND)
t1 = time.time()
print("\n--- Kandidatët ---")
print("  brute force:", cstats["brute_force_pairs"], "| kandidatë:", cstats["candidate_pairs"])
q = candidate_quality(reads, candidates, min_overlap=MIN_OVERLAP)
print("  precision:", round(q["precision"], 4), "| recall:", round(q["recall"], 4))
print("  koha:", round(t1 - t0, 2), "s")

overlaps, ostats = compute_overlaps(reads, candidates, min_overlap=MIN_OVERLAP,
                                    max_err=MAX_ERR, band=ALIGN_BAND)
t2 = time.time()
print("\n--- Overlap-et e verifikuara ---")
for key, val in ostats.items():
    print(" ", key, ":", round(val, 4) if isinstance(val, float) else val)
print("  precision i overlap-eve:", round(overlap_precision(reads, overlaps, MIN_OVERLAP), 4))
print("  koha:", round(t2 - t1, 2), "s")

kinds = {}
for o in overlaps:
    kinds[o["kind"]] = kinds.get(o["kind"], 0) + 1
print("  llojet:", kinds)