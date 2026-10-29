"""K-mer inverted index me k-mers kanonikë dhe gjenerim kandidatësh sipas diagonales."""
from collections import defaultdict

from src.io_utils import reverse_complement


def canonical(kmer):
    """Kthen (k-mer kanonik, True nëse kmer është vetë kanoniku)."""
    rc = reverse_complement(kmer)
    if kmer <= rc:
        return kmer, True
    return rc, False


def build_index(reads, k=15, max_occ=100):
    """Ndërton indeksin k-mer -> [(read_id, pozicioni, forward_flag)].
    reads: listë (emri, sekuenca, quality).
    K-mers me më shumë se max_occ shfaqje (repeats) hidhen.
    Kthen (indeksi, statistikat)."""
    index = defaultdict(list)
    total_kmers = 0
    for rid, (_, seq, _) in enumerate(reads):
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            if "N" in kmer:
                continue
            c, fwd = canonical(kmer)
            if c == kmer and fwd and reverse_complement(kmer) == kmer:
                continue  # k-mer palindromik: orientimi i paqartë
            index[c].append((rid, i, fwd))
            total_kmers += 1

    distinct_before = len(index)
    filtered = {c: occ for c, occ in index.items() if len(occ) <= max_occ}
    stats = {
        "k": k,
        "total_kmer_occurrences": total_kmers,
        "distinct_kmers": distinct_before,
        "kept_kmers": len(filtered),
        "dropped_repeat_kmers": distinct_before - len(filtered),
    }
    return filtered, stats


def _best_cluster(diags, band):
    """Gjen grupin më të madh të diagonaleve brenda një brezi me gjerësi band.
    Kthen (numri i votave, diagonalja mediane e grupit)."""
    diags.sort()
    best_cnt, best_diag, j = 0, None, 0
    for i in range(len(diags)):
        while diags[i] - diags[j] > band:
            j += 1
        cnt = i - j + 1
        if cnt > best_cnt:
            best_cnt = cnt
            best_diag = diags[(i + j) // 2]
    return best_cnt, best_diag


def find_candidates(reads, index, k=15, min_shared=3, min_overlap=40, band=5):
    """Gjen çiftet kandidate të reads që mund të mbivendosen.
    Kthen (kandidatët, statistikat). Çdo kandidat është një tuple:
    (a, b, same_strand, diagonal, votat, overlap_i_vlerësuar)
    ku diagonal = pozicioni ku fillon b (ose revcomp(b) nëse same_strand=False)
    brenda a, dhe a < b."""
    lens = [len(r[1]) for r in reads]
    hits = defaultdict(list)  # (a, b, same) -> [diagonale]

    for occ in index.values():
        n = len(occ)
        for x in range(n):
            ra, pa, fa = occ[x]
            for y in range(x + 1, n):
                rb, pb, fb = occ[y]
                if ra == rb:
                    continue
                if ra < rb:
                    a, b, pos_a, pos_b, fl_a, fl_b = ra, rb, pa, pb, fa, fb
                else:
                    a, b, pos_a, pos_b, fl_a, fl_b = rb, ra, pb, pa, fb, fa
                same = (fl_a == fl_b)
                if same:
                    diag = pos_a - pos_b
                else:
                    diag = pos_a - (lens[b] - k - pos_b)
                hits[(a, b, same)].append(diag)

    candidates = []
    for (a, b, same), diags in hits.items():
        votes, diag = _best_cluster(diags, band)
        if votes < min_shared:
            continue
        ovl = min(lens[a], diag + lens[b]) - max(0, diag)
        if ovl >= min_overlap:
            candidates.append((a, b, same, diag, votes, ovl))

    n = len(reads)
    stats = {
        "reads": n,
        "brute_force_pairs": n * (n - 1) // 2,
        "pairs_sharing_a_kmer": len(hits),
        "candidate_pairs": len(candidates),
    }
    return candidates, stats