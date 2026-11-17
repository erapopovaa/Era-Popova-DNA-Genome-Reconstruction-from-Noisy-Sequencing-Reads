import random

from src.consensus import banded_align, build_pileup, call_variants, polish
from src.io_utils import reverse_complement


def rand_seq(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


ref = rand_seq(200, 5)
x = ref[10:60]

# 1. Zëvendësim: kosto 1, pa indels
y = list(x)
y[20] = "A" if y[20] != "A" else "C"
cost, ops, start = banded_align("".join(y), ref, 10, 5)
assert cost == 1 and start == 10 and all(o[0] == "M" for o in ops), (cost, start)

# 2. Deletion në x (ref ka një bazë më shumë)
cost, ops, _ = banded_align(x[:25] + x[26:], ref, 10, 5)
assert cost == 1 and [o[0] for o in ops].count("D") == 1, (cost, ops)

# 3. Insertion në x
cost, ops, _ = banded_align(x[:25] + "A" + x[25:], ref, 10, 5)
assert cost == 1 and [o[0] for o in ops].count("I") == 1, (cost, ops)


def make_reads(g, length=80, step=10):
    reads = []
    for n, s in enumerate(range(0, len(g) - length + 1, step)):
        seq, q = g[s:s + length], [30] * length
        if n % 2:
            seq, q = reverse_complement(seq), q[::-1]
        reads.append(("r%d" % n, seq, q))
    return reads


# 4. Variant calling: SNP, INS, DEL
g = rand_seq(300, 11)
sub = "A" if g[100] != "A" else "C"
cons = g[:100] + sub + g[101:200] + g[201:250] + "T" + g[250:]
pile, st = build_pileup(make_reads(cons), cons)
assert st["reads_aligned"] == st["reads_total"], st
called, filtered, cost = call_variants(cons, g, pile, 0)
assert sorted(v[1] for v in called) == ["DEL", "INS", "SNP"], called
snp = [v for v in called if v[1] == "SNP"][0]
assert snp[0] == 101 and snp[3] == sub and snp[4] >= 3, snp

# 5. Korrigjimi: contig me një SNP dhe një deletion rikthehet në gjenomin e vërtetë
bad = g[:150] + ("A" if g[150] != "A" else "C") + g[151:220] + g[221:]
fixed, _, st = polish(make_reads(g), bad, rounds=2)
assert fixed == g, (len(fixed), len(g))

print("Të gjitha testet kaluan.")