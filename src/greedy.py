"""Baseline greedy: skajet me peshë më të madhe pranohen nëse nuk krijojnë degëzim ose cikël."""
from src.assemble import path_to_sequence


def greedy_contigs(reads, graph, active):
    """graph: grafi para transitive reduction. Kthen listë (rruga, sekuenca)."""
    edges = []
    for u, outs in graph.items():
        for v, (_, score) in outs.items():
            if (u, v) <= (v ^ 1, u ^ 1):  # çdo çift binjakësh vetëm një herë
                edges.append((score, u, v))
    edges.sort(reverse=True)

    parent = {r: r for r in active}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    nxt, prv = {}, {}
    for score, u, v in edges:
        if (u >> 1) not in parent or (v >> 1) not in parent or (u >> 1) == (v >> 1):
            continue
        if u in nxt or v in prv or (v ^ 1) in nxt or (u ^ 1) in prv:
            continue
        ru, rv = find(u >> 1), find(v >> 1)
        if ru == rv:
            continue  # do të krijonte cikël
        nxt[u], prv[v] = v, u
        nxt[v ^ 1], prv[u ^ 1] = u ^ 1, v ^ 1
        parent[ru] = rv

    visited, contigs = set(), []
    for r in sorted(active):
        if r in visited:
            continue
        u = 2 * r
        while u in prv:
            u = prv[u]
        path = [u]
        visited.add(u >> 1)
        while path[-1] in nxt:
            v = nxt[path[-1]]
            path.append(v)
            visited.add(v >> 1)
        contigs.append((path, path_to_sequence(reads, graph, path)))
    contigs.sort(key=lambda c: len(c[1]), reverse=True)
    return contigs