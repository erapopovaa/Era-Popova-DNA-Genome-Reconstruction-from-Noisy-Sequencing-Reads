from src.io_utils import read_fasta, reverse_complement

records = read_fasta("data/reference/phix174.fasta")

for name, seq in records.items():
    print("Emri:", name)
    print("Gjatesia:", len(seq))
    print("Bazat:", {b: seq.count(b) for b in "ACGT"})
    print("Fillimi:", seq[:60])
    print("RC e dyfishte e barabarte me origjinalin:",
          reverse_complement(reverse_complement(seq)) == seq)