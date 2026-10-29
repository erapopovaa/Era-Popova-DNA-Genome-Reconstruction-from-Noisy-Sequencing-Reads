"""Simulator: fut mutacione në një gjenom real dhe gjeneron reads me gabime."""
import random

from src.io_utils import reverse_complement


def mutate_genome(genome, n_snp=10, n_ins=3, n_del=3, seed=42):
    """Fut mutacione të njohura. Kthen (gjenomi_i_mutuar, lista_e_variantave).
    Çdo variant është (pos_1based, lloji, ref, alt) në koordinatat e referencës.
    Për INS, pos është baza para futjes."""
    rng = random.Random(seed)
    n_total = n_snp + n_ins + n_del
    slot = 100  # mutacionet ndahen në "kutia" 100bp që të mos mbivendosen
    n_slots = len(genome) // slot - 1
    slots = sorted(rng.sample(range(1, n_slots), n_total))
    types = ["SNP"] * n_snp + ["INS"] * n_ins + ["DEL"] * n_del
    rng.shuffle(types)

    variants = []  # (pos0, lloji, ref, alt)
    for s, t in zip(slots, types):
        pos = s * slot + rng.randint(10, 80)
        if t == "SNP":
            ref = genome[pos]
            alt = rng.choice([b for b in "ACGT" if b != ref])
        elif t == "INS":
            ref = "-"
            alt = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 3)))
        else:
            ref = genome[pos:pos + rng.randint(1, 3)]
            alt = "-"
        variants.append((pos, t, ref, alt))

    g = list(genome)
    for pos, t, ref, alt in reversed(variants):  # nga fundi, që pozicionet mos zhvendosen
        if t == "SNP":
            g[pos] = alt
        elif t == "INS":
            g[pos + 1:pos + 1] = list(alt)
        else:
            del g[pos:pos + len(ref)]

    truth = [(pos + 1, t, ref, alt) for pos, t, ref, alt in variants]
    return "".join(g), truth


def _add_errors(seq, rng, error_rate):
    """Shton gabime: 90% zëvendësim, 5% insertion, 5% deletion.
    Kthen (sekuenca, quality-t). Bazat me gabim marrin quality të ulët."""
    out, quals = [], []
    for b in seq:
        if rng.random() < error_rate:
            r = rng.random()
            if r < 0.90:
                out.append(rng.choice([x for x in "ACGT" if x != b]))
                quals.append(rng.randint(5, 15))
            elif r < 0.95:
                out.append(b)
                out.append(rng.choice("ACGT"))
                quals.append(rng.randint(5, 15))
                quals.append(rng.randint(5, 15))
            # else: deletion, baza anashkalohet
        else:
            out.append(b)
            quals.append(rng.randint(30, 40))
    return "".join(out), quals


def simulate_reads(genome, coverage=30, min_len=100, max_len=150,
                   error_rate=0.01, rc_fraction=0.5, dup_fraction=0.02, seed=42):
    """Gjeneron reads. Kthen listë (emri, sekuenca, quality_list)."""
    rng = random.Random(seed)
    mean_len = (min_len + max_len) / 2
    n_reads = int(coverage * len(genome) / mean_len)
    reads = []
    for i in range(n_reads):
        if reads and rng.random() < dup_fraction:
            _, seq, qual = rng.choice(reads)  # dublikatë
            reads.append(("read_%d_dup" % i, seq, list(qual)))
            continue
        length = rng.randint(min_len, max_len)
        start = rng.randint(0, len(genome) - length)
        frag = genome[start:start + length]
        seq, qual = _add_errors(frag, rng, error_rate)
        strand = "+"
        if rng.random() < rc_fraction:
            seq = reverse_complement(seq)
            qual = qual[::-1]
            strand = "-"
        reads.append(("read_%d_pos%d_%s" % (i, start, strand), seq, qual))
    return reads


def write_fastq(path, reads):
    """Shkruan reads në FASTQ."""
    with open(path, "w") as f:
        for name, seq, qual in reads:
            f.write("@" + name + "\n" + seq + "\n+\n")
            f.write("".join(chr(q + 33) for q in qual) + "\n")


def write_truth(path, truth):
    """Shkruan mutacionet e njohura në TSV."""
    with open(path, "w") as f:
        f.write("POS\tTYPE\tREF\tALT\n")
        for pos, t, ref, alt in truth:
            f.write("%d\t%s\t%s\t%s\n" % (pos, t, ref, alt))