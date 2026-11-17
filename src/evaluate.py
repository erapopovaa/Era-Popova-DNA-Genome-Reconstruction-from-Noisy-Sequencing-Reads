"""Vlerësim: sa të mira janë kandidatët krahasuar me overlap-et e vërteta (nga simulatori)."""
import re

_NAME = re.compile(r"_pos(\d+)_([+-])$")


def _read_info(reads):
    """{read_id: (start, gjatësia, strand)} për reads me pozicion të njohur (jo dublikatat)."""
    info = {}
    for rid, (name, seq, _) in enumerate(reads):
        m = _NAME.search(name)
        if m:
            info[rid] = (int(m.group(1)), len(seq), m.group(2))
    return info


def _true_pairs(info, min_ovl):
    ids = sorted(info, key=lambda r: info[r][0])
    pairs = set()
    for x, a in enumerate(ids):
        sa, la, _ = info[a]
        for b in ids[x + 1:]:
            sb, lb, _ = info[b]
            if sb > sa + la - min_ovl:
                break
            if min(sa + la, sb + lb) - sb >= min_ovl:
                pairs.add((min(a, b), max(a, b)))
    return pairs


def candidate_quality(reads, candidates, min_overlap, recall_overlap=60):
    """Llogarit precision, recall dhe saktësinë e orientimit të kandidatëve.
    Recall matet për overlap-e të vërteta me gjatësi >= recall_overlap."""
    info = _read_info(reads)
    cand = {(a, b): same for a, b, same, *_ in candidates if a in info and b in info}

    true_loose = _true_pairs(info, min_overlap)
    true_strict = _true_pairs(info, recall_overlap)

    correct = sum(1 for p in cand if p in true_loose)
    found = sum(1 for p in true_strict if p in cand)
    orient_ok = sum(1 for (a, b), same in cand.items()
                    if (p := (a, b)) in true_loose and same == (info[a][2] == info[b][2]))
    return {
        "precision": correct / len(cand) if cand else 0.0,
        "recall": found / len(true_strict) if true_strict else 0.0,
        "true_overlaps_strict": len(true_strict),
        "orientation_accuracy": orient_ok / correct if correct else 0.0,
    }



def overlap_precision(reads, overlaps, min_overlap):
    """Sa nga overlap-et e pranuara janë të vërteta (sipas pozicioneve nga simulatori)."""
    info = _read_info(reads)
    true_loose = _true_pairs(info, min_overlap)
    checked = [(o["a"], o["b"]) for o in overlaps if o["a"] in info and o["b"] in info]
    ok = sum(1 for p in checked if p in true_loose)
    return ok / len(checked) if checked else 0.0



def graph_edge_accuracy(reads, graph):
    """Sa nga skajet e grafit janë gjeometrikisht të sakta sipas pozicioneve nga simulatori:
    të dy nyjet në të njëjtin orientim, drejtimi përgjatë gjenomit i duhur, dhe reads mbivendosen."""
    info = _read_info(reads)
    total = ok = 0
    for u, outs in graph.items():
        ru, su = u >> 1, u & 1
        if ru not in info:
            continue
        s_u, l_u, st_u = info[ru]
        fwd_u = (st_u == "+") != (su == 1)
        for v in outs:
            rv, sv = v >> 1, v & 1
            if rv not in info:
                continue
            s_v, l_v, st_v = info[rv]
            fwd_v = (st_v == "+") != (sv == 1)
            total += 1
            if fwd_u != fwd_v:
                continue
            direction_ok = s_v > s_u if fwd_u else s_v < s_u
            overlap_ok = s_v < s_u + l_u and s_u < s_v + l_v
            if direction_ok and overlap_ok:
                ok += 1
    return ok / total if total else 0.0



def _global_edit_distance(x, y, w=60):
    """Edit distance global në brez |i-j| <= w. None nëse gjatësitë ndryshojnë më shumë se w."""
    INF = 10 ** 9
    n, m = len(x), len(y)
    if abs(n - m) > w:
        return None
    width = 2 * w + 1
    prev = [INF] * width
    for k in range(width):
        j = k - w
        if 0 <= j <= m:
            prev[k] = j
    for i in range(1, n + 1):
        cur = [INF] * width
        xi = x[i - 1]
        for k in range(width):
            j = i + k - w
            if j < 0 or j > m:
                continue
            v = INF
            if j >= 1:
                v = prev[k] + (0 if xi == y[j - 1] else 1)
            if k + 1 < width and prev[k + 1] + 1 < v:
                v = prev[k + 1] + 1
            if k >= 1 and cur[k - 1] + 1 < v:
                v = cur[k - 1] + 1
            cur[k] = v
        prev = cur
    k = m - n + w
    return prev[k] if 0 <= k < width and prev[k] < INF else None


def contig_vs_reference(contig, ref, k=25, band=60):
    """Vendos contig-un në referencë (provon të dy strand-et), pastaj llogarit edit distance.
    Kthen fjalor ose None nëse nuk gjendet pozicioni."""
    from src.io_utils import reverse_complement
    best = None
    for strand, seq in (("+", contig), ("-", reverse_complement(contig))):
        off = None
        for s in range(0, max(1, min(len(seq) - k, 600)), 20):
            idx = ref.find(seq[s:s + k])
            if idx != -1:
                off = idx - s
                break
        if off is None:
            continue
        off = max(0, off)
        seg = ref[off:off + len(seq)]
        ed = _global_edit_distance(seq, seg, band)
        if ed is None:
            continue
        if best is None or ed < best["edit_distance"]:
            best = {"edit_distance": ed, "strand": strand, "offset": off,
                    "contig_length": len(seq),
                    "identity": 1 - ed / max(len(seq), len(seg))}
    return best