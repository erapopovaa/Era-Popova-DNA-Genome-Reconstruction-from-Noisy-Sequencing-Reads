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