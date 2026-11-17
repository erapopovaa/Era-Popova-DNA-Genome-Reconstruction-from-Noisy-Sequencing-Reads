"""Pileup, konsensus (korrigjim gabimesh) dhe variant calling."""
import math
from collections import Counter, defaultdict

from src.io_utils import reverse_complement

INF = 10 ** 9
IDX = {"A": 0, "C": 1, "G": 2, "T": 3}
BASES = "ACGT"


def banded_align(x, ref, d, w):
    """Alinjim me traceback në brez rreth diagonales d. x konsumohet plotësisht,
    fillimi dhe fundi në ref janë të lirë.
    Kthen (kosto, ops, fillimi_në_ref) ose None. ops janë tuple:
      ('M', i, c)  x[i] alinjohet me ref[c] (match ose mismatch)
      ('I', i, c)  x[i] është insertion (ref[:c] është konsumuar)
      ('D', i, c)  ref[c] është deletion (x[:i] është konsumuar)"""
    n, m = len(x), len(ref)
    width = 2 * w + 1
    prev = [INF] * width
    for k in range(width):
        c = d + k - w
        if 0 <= c <= m:
            prev[k] = 0
    tb = []
    for i in range(1, n + 1):
        cur = [INF] * width
        row = [-1] * width
        xi = x[i - 1]
        base = d + i - w
        for k in range(width):
            c = base + k
            if c < 0 or c > m:
                continue
            v = INF
            t = -1
            if c >= 1 and prev[k] < INF:
                v = prev[k] + (0 if xi == ref[c - 1] else 1)
                t = 0
            if k >= 1 and cur[k - 1] < INF and cur[k - 1] + 1 < v:
                v = cur[k - 1] + 1
                t = 2
            if k + 1 < width and prev[k + 1] < INF and prev[k + 1] + 1 < v:
                v = prev[k + 1] + 1
                t = 1
            cur[k] = v
            row[k] = t
        tb.append(row)
        prev = cur

    best_k = min(range(width), key=lambda q: prev[q])
    if prev[best_k] >= INF:
        return None
    cost = prev[best_k]
    k, i, ops = best_k, n, []
    while i > 0:
        c = d + i - w + k
        t = tb[i - 1][k]
        if t == 0:
            ops.append(("M", i - 1, c - 1))
            i -= 1
        elif t == 2:
            ops.append(("D", i, c - 1))
            k -= 1
        elif t == 1:
            ops.append(("I", i - 1, c))
            i -= 1
            k += 1
        else:
            return None
    ops.reverse()
    return cost, ops, d - w + k


def orient_contig(contig, ref, k=25):
    """Kthen (contig në strand-in e ref, strand, offset i përafërt në ref)."""
    for strand, seq in (("+", contig), ("-", reverse_complement(contig))):
        for s in range(0, max(1, min(len(seq) - k, 600)), 20):
            idx = ref.find(seq[s:s + k])
            if idx != -1:
                return seq, strand, idx - s
    return contig, "?", 0


def edit_vs_reference(seq, ref, w=60, k=25):
    """Edit distance e seq kundrejt ref (skajet e ref të lira). None nëse s'gjendet."""
    for s in range(0, max(1, min(len(seq) - k, 600)), 20):
        idx = ref.find(seq[s:s + k])
        if idx != -1:
            res = banded_align(seq, ref, max(0, idx - s), w)
            if res:
                return res[0]
    return None


def _kmer_positions(seq, k):
    d = defaultdict(list)
    for i in range(len(seq) - k + 1):
        d[seq[i:i + k]].append(i)
    return d


def _locate(seq, kmap, k, step=2):
    votes = Counter()
    for i in range(0, len(seq) - k + 1, step):
        pos = kmap.get(seq[i:i + k])
        if pos and len(pos) <= 5:
            for p in pos:
                votes[p - i] += 1
    best, best_v = None, 0
    for dd in votes:
        s = sum(votes.get(dd + t, 0) for t in range(-3, 4))
        if s > best_v:
            best, best_v = dd, s
    return best, best_v


def _apply(ops, x, xq, counts, raw, ins, use_quality):
    j, nops = 0, len(ops)
    while j < nops:
        op, i, c = ops[j]
        if op == "M":
            b = IDX.get(x[i])
            if b is not None:
                counts[c][b] += max(1, xq[i]) if use_quality else 1
                raw[c][b] += 1
            j += 1
        elif op == "D":
            wgt = max(1, xq[min(i, len(xq) - 1)]) if use_quality else 1
            counts[c][4] += wgt
            raw[c][4] += 1
            j += 1
        else:
            start = j
            while j < nops and ops[j][0] == "I":
                j += 1
            if c == 0:
                continue
            idxs = [ops[t][1] for t in range(start, j)]
            s = "".join(x[t] for t in idxs)
            if "N" in s:
                continue
            wgt = sum(xq[t] for t in idxs) / len(idxs) if use_quality else 1
            if ins[c - 1] is None:
                ins[c - 1] = {}
            e = ins[c - 1].setdefault(s, [0.0, 0])
            e[0] += max(1, wgt)
            e[1] += 1


def build_pileup(reads, contig, k=15, w=10, use_quality=True, min_votes=3):
    """Alinjon çdo read me contig-un dhe ndërton pileup-in.
    Kthen ((counts, raw, ins), statistikat)."""
    L = len(contig)
    kmap = _kmer_positions(contig, k)
    counts = [[0.0] * 5 for _ in range(L)]
    raw = [[0] * 5 for _ in range(L)]
    ins = [None] * L
    aligned = total_cost = total_len = 0
    for _, seq, qual in reads:
        if len(qual) != len(seq):
            qual = [30] * len(seq)
        best = None
        for strand in (0, 1):
            s = reverse_complement(seq) if strand else seq
            d, v = _locate(s, kmap, k)
            if d is not None and (best is None or v > best[0]):
                best = (v, s, strand, d)
        if best is None or best[0] < min_votes:
            continue
        _, s, strand, d = best
        q = qual[::-1] if strand else qual
        lo, hi = max(0, -d), min(len(s), L - d)
        if hi - lo < 30:
            continue
        x, xq = s[lo:hi], q[lo:hi]
        res = banded_align(x, contig, d + lo, w)
        if res is None:
            continue
        cost, ops, _ = res
        aligned += 1
        total_cost += cost
        total_len += len(x)
        _apply(ops, x, xq, counts, raw, ins, use_quality)
    stats = {
        "reads_total": len(reads),
        "reads_aligned": aligned,
        "pct_aligned": 100.0 * aligned / len(reads) if reads else 0.0,
        "avg_read_error": total_cost / total_len if total_len else 0.0,
        "mean_coverage": sum(sum(r) for r in raw) / L if L else 0.0,
    }
    return (counts, raw, ins), stats


def call_consensus(contig, pile):
    """Vota e shumicës (e peshuar) për çdo pozicion, plus insertions të mbështetura."""
    counts, _, ins = pile
    out = []
    for c in range(len(contig)):
        col = counts[c]
        tot = sum(col)
        if tot == 0:
            out.append(contig[c])
        else:
            b = max(range(5), key=lambda t: col[t])
            if b < 4:
                out.append(BASES[b])
        if ins[c] and tot > 0:
            s, e = max(ins[c].items(), key=lambda kv: kv[1][0])
            if e[0] > 0.5 * tot:
                out.append(s)
    return "".join(out)


def polish(reads, contig, rounds=2, k=15, w=10, use_quality=True):
    """Korrigjon contig-un me `rounds` raunde, pastaj ndërton pileup-in final
    kundrejt konsensusit. Kthen (konsensusi, pileup, statistikat)."""
    for _ in range(rounds):
        pile, _ = build_pileup(reads, contig, k, w, use_quality)
        contig = call_consensus(contig, pile)
    pile, stats = build_pileup(reads, contig, k, w, use_quality)
    return contig, pile, stats


def binom_tail(k, n, p):
    """P(X >= k) për X ~ Binomial(n, p)."""
    if n <= 0 or k <= 0:
        return 1.0
    if k > n:
        return 0.0
    return min(1.0, sum(math.comb(n, x) * p ** x * (1 - p) ** (n - x)
                        for x in range(k, n + 1)))


def _judge(kind, support, cov, eps, min_alt, min_frac, max_p):
    p_err = eps / 3 if kind == "SNP" else eps * 0.05
    p = binom_tail(support, cov, p_err)
    frac = support / cov if cov else 0.0
    return (support >= min_alt and frac >= min_frac and p <= max_p), p


def call_variants(cons, ref, pile, offset=0, w=60, min_alt=3, min_frac=0.25,
                  eps=0.01, max_p=1e-6):
    """Alinjon konsensusin me referencën dhe raporton ndryshimet.
    Kthen (të pranuarat, të filtruarat, edit distance). Çdo variant është
    (pos_1based, lloji, ref, alt, support, coverage, pvalue)."""
    _, raw, ins = pile
    res = banded_align(cons, ref, offset, w)
    if res is None:
        return [], [], None
    cost, ops, _ = res
    accepted, filtered = [], []

    def emit(pos, kind, r, a, support, cov):
        ok, p = _judge(kind, support, cov, eps, min_alt, min_frac, max_p)
        (accepted if ok else filtered).append((pos, kind, r, a, support, cov, p))

    j, n = 0, len(ops)
    while j < n:
        op, i, c = ops[j]
        if op == "M":
            if cons[i] != ref[c]:
                b = IDX.get(cons[i])
                if b is not None:
                    emit(c + 1, "SNP", ref[c], cons[i], raw[i][b], sum(raw[i]))
            j += 1
        elif op == "I":
            start = j
            while j < n and ops[j][0] == "I":
                j += 1
            ii = [ops[t][1] for t in range(start, j)]
            alt = "".join(cons[t] for t in ii)
            sup = min(raw[t][IDX[cons[t]]] for t in ii if cons[t] in IDX)
            emit(c, "INS", "-", alt, sup, sum(raw[ii[0]]))
        else:
            start = j
            while j < n and ops[j][0] == "D":
                j += 1
            c_first, c_last = ops[start][2], ops[j - 1][2]
            i2 = min(ops[start][1], len(cons) - 1)
            ins_reads = 0
            if i2 >= 1 and ins[i2 - 1]:
                ins_reads = sum(e[1] for e in ins[i2 - 1].values())
            sup = max(0, sum(raw[i2][:4]) - ins_reads)
            emit(c_first + 1, "DEL", ref[c_first:c_last + 1], "-", sup, sum(raw[i2]))
    accepted.sort()
    filtered.sort()
    return accepted, filtered, cost


def minor_alleles(cons, pile, min_alt=3, min_frac=0.25, eps=0.01, max_p=1e-6):
    """Pozicione ku një alel jo-konsensus ka mbështetje mbi nivelin e gabimit."""
    raw = pile[1]
    out = []
    for c, base in enumerate(cons):
        r = raw[c]
        tot = sum(r)
        if tot == 0:
            continue
        cb = IDX.get(base)
        for b in range(4):
            if b == cb:
                continue
            cnt = r[b]
            if cnt >= min_alt and cnt / tot >= min_frac \
                    and binom_tail(cnt, tot, eps / 3) <= max_p:
                out.append((c + 1, BASES[b], cnt, tot))
    return out