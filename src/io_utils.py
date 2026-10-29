"""Funksione për lexim/shkrim të FASTA dhe FASTQ, dhe reverse-complement."""

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq):
    """Kthen reverse-complement të një sekuence ADN-je."""
    return seq.translate(_COMPLEMENT)[::-1]


def read_fasta(path):
    """Lexon një FASTA dhe kthen fjalor {emri: sekuenca}."""
    records = {}
    name = None
    chunks = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        records[name] = "".join(chunks).upper()
    return records


def write_fasta(path, records, width=70):
    """Shkruan fjalor {emri: sekuenca} në FASTA."""
    with open(path, "w") as f:
        for name, seq in records.items():
            f.write(">" + name + "\n")
            for i in range(0, len(seq), width):
                f.write(seq[i:i + width] + "\n")


def read_fastq(path, max_reads=None):
    """Gjenerator: jep (emri, sekuenca, quality) për çdo read.
    quality është listë numrash Phred (Q-score)."""
    count = 0
    with open(path, "r") as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip().upper()
            f.readline()  # rreshti "+"
            qual = f.readline().strip()
            yield header[1:].strip().split()[0], seq, [ord(c) - 33 for c in qual]
            count += 1
            if max_reads is not None and count >= max_reads:
                break