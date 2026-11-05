import random

from src.overlap import align_overlap, compute_overlaps


def rand_seq(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


genome = rand_seq(300, 1)
a = genome[0:100]
b = genome[60:160]          # b fillon në pozicionin 60 të a, overlap = 40

# 1. Overlap i saktë
ed, length, dp = align_overlap(a, b, 60, 8)
assert (ed, length) == (0, 40), (ed, length)

# 2. Një zëvendësim
b2 = list(b)
b2[10] = "A" if b2[10] != "A" else "C"
ed, length, dp = align_overlap(a, "".join(b2), 60, 8)
assert ed == 1 and length == 40, (ed, length)

# 3. Një deletion në b (kërkon DP)
b3 = b[:10] + b[11:]
ed, length, dp = align_overlap(a, b3, 60, 8)
assert dp is True and ed == 1, (ed, length, dp)

# 4. Reads të palidhura duhet të refuzohen
other = rand_seq(100, 99)
reads = [("a", a, []), ("o", other, [])]
ovs, st = compute_overlaps(reads, [(0, 1, True, 60, 5, 40)], min_overlap=30, max_err=0.05)
assert len(ovs) == 0, ovs

# 5. Reverse-complement: çifti pranohet kur same=False
from src.io_utils import reverse_complement
reads = [("a", a, []), ("b", reverse_complement(b), [])]
ovs, st = compute_overlaps(reads, [(0, 1, False, 60, 5, 40)], min_overlap=30, max_err=0.05)
assert len(ovs) == 1 and ovs[0]["kind"] == "a_to_b", ovs

print("Të gjitha testet kaluan.")