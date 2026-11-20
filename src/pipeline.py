"""Pipeline i plotë si funksion i ripërdorshëm (për main, reads reale dhe benchmark)."""
import os
import time

import psutil

from src import config as C
from src.kmer_index import build_index, find_candidates
from src.overlap import compute_overlaps
from src.graph import find_contained, build_graph, transitive_reduction
from src.assemble import simplify, build_contigs, resolve_branches
from src.greedy import greedy_contigs
from src.consensus import polish, orient_contig, call_variants, edit_vs_reference
from src.evaluate import compare_variants, merge_adjacent_indels


def run_pipeline(reads, ref, mutated=None, truth=None, k=C.K,
                 min_shared=C.MIN_SHARED, min_overlap=C.MIN_OVERLAP,
                 max_err=C.MAX_ERR, fuzz=C.FUZZ, reduce_graph=True,
                 use_quality=True, method="olc"):
    """Ekzekuton pipeline-in. ref = referenca origjinale; mutated = gjenomi i vërtetë
    (vetëm në simulim); truth = lista e mutacioneve të vërteta (vetëm në simulim).
    Kthen (report, konsensusi, variantet)."""
    t0 = time.time()
    lens = [len(r[1]) for r in reads]

    index, _ = build_index(reads, k=k, max_occ=C.MAX_OCC)
    candidates, cst = find_candidates(reads, index, k=k, min_shared=min_shared,
                                      min_overlap=min_overlap, band=C.BAND)
    overlaps, ost = compute_overlaps(reads, candidates, min_overlap=min_overlap,
                                     max_err=max_err, band=C.ALIGN_BAND)
    contained = find_contained(reads, overlaps)
    active = [r for r in range(len(reads)) if r not in contained]
    graph = build_graph(reads, overlaps, contained)

    if method == "greedy":
        contigs = greedy_contigs(reads, graph, active)
    else:
        if reduce_graph:
            transitive_reduction(graph, lens, fuzz=fuzz)
        if method == "olc_best":
            resolve_branches(graph)
        simplify(graph, active, C.TIP_READS, C.BUBBLE_LEN, C.SIMPLIFY_ROUNDS)
        contigs = build_contigs(reads, graph, active)

    report = {
        "method": method,
        "reads_total": len(reads),
        "brute_force_comparisons": cst["brute_force_pairs"],
        "candidate_comparisons": cst["candidate_pairs"],
        "overlaps_identified": len(overlaps),
        "avg_overlap_quality": round(ost["avg_quality"], 4),
        "contigs": len(contigs),
        "longest_contig": len(contigs[0][1]) if contigs else 0,
    }
    if not contigs:
        report["error"] = "nuk u formua asnjë contig"
        return report, "", []

    oriented, strand, off = orient_contig(contigs[0][1], ref)
    report["contig_strand"] = strand
    if strand == "?":
        report["error"] = "contig-u nuk u vendos dot në referencë"
        report["time_seconds"] = round(time.time() - t0, 2)
        return report, oriented, []

    cons, pile, pst = polish(reads, oriented, rounds=C.POLISH_ROUNDS, k=k,
                             w=C.POLISH_BAND, use_quality=use_quality)
    target = mutated if mutated else ref
    ed = edit_vs_reference(cons, target)
    called, filtered, _ = call_variants(cons, ref, pile, offset=max(0, off),
                                        min_alt=C.MIN_ALT, min_frac=C.MIN_ALT_FRAC,
                                        eps=C.SEQ_ERROR, max_p=C.PVALUE_MAX)
    called = merge_adjacent_indels(called)

    report.update({
        "reconstruction_length": len(cons),
        "edit_distance": ed,
        "identity": round(1 - ed / len(cons), 5) if ed is not None else None,
                "genome_fraction": round(len(cons) / len(target), 4),
        "reads_aligned": pst["reads_aligned"],
        "pct_reads_incorporated": round(pst["pct_aligned"], 2),
        "mean_coverage": round(pst["mean_coverage"], 1),
        "variants_called": len(called),
    })
    if truth is not None:
        res = compare_variants(called, truth)
        report.update({"tp": res["tp"], "fp": res["fp"], "fn": res["fn"]})

    mem = psutil.Process(os.getpid()).memory_info()
    report["time_seconds"] = round(time.time() - t0, 2)
    report["peak_memory_mb"] = round(getattr(mem, "peak_wset", mem.rss) / 1e6, 1)
    return report, cons, called