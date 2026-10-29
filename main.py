import time

from src.config import K, MAX_OCC, MIN_SHARED, BAND, MIN_OVERLAP
from src.io_utils import read_fastq
from src.kmer_index import build_index, find_candidates
from src.evaluate import candidate_quality

reads = list(read_fastq("data/reads/sim_phix174.fastq"))
print("Reads:", len(reads))

t0 = time.time()
index, istats = build_index(reads, k=K, max_occ=MAX_OCC)
t1 = time.time()
candidates, cstats = find_candidates(reads, index, k=K, min_shared=MIN_SHARED,
                                     min_overlap=MIN_OVERLAP, band=BAND)
t2 = time.time()

print("\n--- Indeksi ---")
for key, val in istats.items():
    print(" ", key, ":", val)
print("  koha:", round(t1 - t0, 2), "s")

print("\n--- Kandidatët ---")
for key, val in cstats.items():
    print(" ", key, ":", val)
print("  koha:", round(t2 - t1, 2), "s")
print("  zvogelimi kundrejt brute force: %.1fx" %
      (cstats["brute_force_pairs"] / max(1, cstats["candidate_pairs"])))

print("\n--- Cilësia (kundrejt të vërtetës nga simulatori) ---")
q = candidate_quality(reads, candidates, min_overlap=MIN_OVERLAP)
for key, val in q.items():
    print(" ", key, ":", round(val, 4) if isinstance(val, float) else val)