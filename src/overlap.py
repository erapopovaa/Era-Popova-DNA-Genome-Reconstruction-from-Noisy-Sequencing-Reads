"""Overlap i përafërt: Hamming i shpejtë + banded edit distance."""
from src.io_utils import reverse_complement

INF = 10 ** 9


def _segments(a, b, diag):
    """Pjesët e a dhe b që fillojnë bashkë në fillim të overlap-it."""
    return a[max(0, diag):], b[max(0, -diag):]


def banded_edit_distance(x, y, w):
    """Edit distance me fillim të fiksuar dhe fund të lirë, në brez |i-j| <= w.
    Overlap-i mbaron kur x ose y konsumohet plotësisht.
    Kthen (edit_distance, gjatësia_e_overlap-it) ose (None, 0) nëse s'ka rrugë."""
    n, m = len(x), len(y)
    width = 2 * w + 1
    N = min(n, m + w)

    prev = [INF] * width
    for k in range(width):
        j = k - w
        if 0 <= j <= m:
            prev[k] = j

    best = None  # (ed, -gjatësia)

    def consider(ed, i, j):
        nonlocal best
        if ed >= INF:
            return
        cand = (ed, -(i + j) // 2)
        if best is None or cand < best:
            best = cand

    for i in range(1, N + 1):
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
            if j == m or i == n:
                consider(v, i, j)
        prev = cur

    if best is None:
        return None, 0
    return best[0], -best[1]


def align_overlap(a, b, diag, w):
    """Alinjon overlap-in e a dhe b (b tashmë i orientuar, fillon te diag në a).
    Kthen (edit_distance, gjatësia, përdorur_dp)."""
    x, y = _segments(a, b, diag)
    L = min(len(x), len(y))
    if L <= 0:
        return None, 0, False
    ham = sum(1 for i in range(L) if x[i] != y[i])
    if ham <= 0.025 * L:
        return ham, L, False
    ed, ovl = banded_edit_distance(x, y, w)
    return ed, ovl, True


def overlap_kind(la, lb, diag):
    """Lloji i overlap-it sipas gjeometrisë."""
    if diag >= 0:
        return "b_in_a" if diag + lb <= la else "a_to_b"
    return "a_in_b" if diag + lb >= la else "b_to_a"


def compute_overlaps(reads, candidates, min_overlap=40, max_err=0.05, band=8):
    """Verifikon çdo kandidat. Kthen (overlaps, statistikat).
    Çdo overlap është fjalor me: a, b, same, diag, kind, length, ed, err, score."""
    overlaps = []
    dp_used = 0
    for a, b, same, diag, votes, _ in candidates:
        sa = reads[a][1]
        sb = reads[b][1]
        if not same:
            sb = reverse_complement(sb)
        ed, length, used_dp = align_overlap(sa, sb, diag, band)
        if used_dp:
            dp_used += 1
        if ed is None or length < min_overlap:
            continue
        err = ed / length
        if err > max_err:
            continue
        overlaps.append({
            "a": a, "b": b, "same": same, "diag": diag,
            "kind": overlap_kind(len(sa), len(sb), diag),
            "length": length, "ed": ed, "err": err,
            "score": length * (1 - err),
        })

    n = len(overlaps)
    stats = {
        "candidates_checked": len(candidates),
        "overlaps_accepted": n,
        "dp_alignments": dp_used,
        "hamming_fast_path": len(candidates) - dp_used,
        "avg_overlap_length": sum(o["length"] for o in overlaps) / n if n else 0.0,
        "avg_error_rate": sum(o["err"] for o in overlaps) / n if n else 0.0,
        "avg_quality": sum(o["score"] / o["length"] for o in overlaps) / n if n else 0.0,
    }
    return overlaps, stats