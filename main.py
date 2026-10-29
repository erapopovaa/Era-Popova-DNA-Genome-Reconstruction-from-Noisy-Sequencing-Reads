from src.io_utils import read_fasta, write_fasta, read_fastq
from src.simulator import mutate_genome, simulate_reads, write_fastq, write_truth

ref = list(read_fasta("data/reference/phix174.fasta").values())[0]

mutated, truth = mutate_genome(ref, n_snp=10, n_ins=3, n_del=3, seed=42)
print("Gjatesia referenca:", len(ref), "| e mutuar:", len(mutated))
print("Mutacione:", len(truth))
for v in truth:
    print("  ", v)

reads = simulate_reads(mutated, coverage=30, min_len=100, max_len=150,
                       error_rate=0.01, seed=42)
print("Numri i reads:", len(reads))
print("Shembull:", reads[0][0], reads[0][1][:50])

write_fasta("data/reference/phix174_mutated.fasta", {"phix174_mutated": mutated})
write_fastq("data/reads/sim_phix174.fastq", reads)
write_truth("data/reference/phix174_truth.tsv", truth)

# kontroll: lexo FASTQ-in mbrapsht
back = list(read_fastq("data/reads/sim_phix174.fastq"))
print("Lexuar nga FASTQ:", len(back), "reads")