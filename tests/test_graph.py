import random

from src.graph import find_contained, build_graph, transitive_reduction, count_edges


def rand_seq(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


g = rand_seq(400, 7)
reads = [
    ("r0", g[0:100], []),
    ("r1", g[30:130], []),
    ("r2", g[60:160], []),
    ("r3", g[40:90], []),    # i përmbajtur në r1 dhe r0
]


def ov(a, b, diag, length):
    return {"a": a, "b": b, "same": True, "diag": diag, "length": length,
            "score": float(length)}


overlaps = [
    ov(0, 1, 30, 70),
    ov(0, 2, 60, 40),
    ov(1, 2, 30, 70),
    ov(0, 3, 40, 50),
    ov(1, 3, 10, 50),
]

contained = find_contained(reads, overlaps)
assert contained == {3}, contained

graph = build_graph(reads, overlaps, contained)
assert count_edges(graph) == 6, count_edges(graph)   # 3 skaje + 3 binjakë

lens = [len(r[1]) for r in reads]
removed = transitive_reduction(graph, lens, fuzz=10)
assert removed == 2, removed                          # 0->2 dhe binjaku i tij

edges = {(u, v) for u, outs in graph.items() for v in outs}
assert edges == {(0, 2), (2, 4), (5, 3), (3, 1)}, edges

print("Të gjitha testet kaluan.")