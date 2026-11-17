import json

from src.io_utils import read_fastq, read_fasta, write_fasta
from src.pipeline import run_pipeline

READS_FILE = "data/reads/real_phix.fastq.gz"
N_READS = 1200      # për ~30x me reads 150 bp; rrite nëse reads janë më të shkurtra
MIN_LEN = 50

ref = list(read_fasta("data/reference/phix174.fasta").values())[0]
reads = [r for r in read_fastq(READS_FILE, max_reads=N_READS) if len(r[1]) >= MIN_LEN]

lens = [len(r[1]) for r in reads]
print("Reads:", len(reads), "| gjatësia min/mesatare/max:",
      min(lens), round(sum(lens) / len(lens)), max(lens))
print("Coverage i pritur: %.1fx" % (sum(lens) / len(ref)))

report, cons, called = run_pipeline(reads, ref)

print("\n--- Raporti (reads reale) ---")
for key, val in report.items():
    print(" ", key, ":", val)

print("\n--- Variantet kundrejt referencës ---")
print("  POS    TYPE  REF  ALT  SUPPORT  COVERAGE  PVALUE")
for v in called:
    print("  %-6d %-5s %-4s %-4s %-8d %-9d %.1e" % v)

if cons:
    write_fasta("results/reconstructed_real.fasta", {"reconstructed_real_phix": cons})
with open("results/report_real.json", "w") as f:
    json.dump(report, f, indent=2)