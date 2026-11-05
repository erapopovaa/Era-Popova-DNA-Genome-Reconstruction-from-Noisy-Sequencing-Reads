"""Grafi i overlap-it (string graph): containment, nyje me dy orientime, transitive reduction."""
from collections import defaultdict


def classify(la, lb, diag):
    """Lloji i overlap-it sipas gjeometrisë (b i orientuar, fillon te diag brenda a)."""
    b_in_a = diag >= 0 and diag + lb <= la
    a_in_b = diag <= 0 and diag + lb >= la
    if b_in_a and a_in_b:
        return "equal"
    if b_in_a:
        return "b_in_a"
    if a_in_b:
        return "a_in_b"
    return "a_to_b" if diag > 0 else "b_to_a"


def find_contained(reads, overlaps):
    """Kthen bashkësinë e reads që janë të përmbajtura në një read tjetër."""
    lens = [len(r[1]) for r in reads]
    contained = set()
    for o in overlaps:
        a, b = o["a"], o["b"]
        kind = classify(lens[a], lens[b], o["diag"])
        if kind == "equal":
            contained.add(max(a, b))
        elif kind == "b_in_a":
            contained.add(b)
        elif kind == "a_in_b":
            contained.add(a)
    return contained


def _add_edge(graph, u, v, o):
    val = (o["length"], o["score"])
    old = graph[u].get(v)
    if old is None or val[1] > old[1]:
        graph[u][v] = val


def build_graph(reads, overlaps, contained):
    """Ndërton grafin. Nyja = 2*read + strand (0 si është, 1 reverse-complement).
    graph[u][v] = (gjatësia_e_overlap, score). Çdo skaj ka binjakun e tij (v^1 -> u^1)."""
    lens = [len(r[1]) for r in reads]
    graph = defaultdict(dict)
    for o in overlaps:
        a, b = o["a"], o["b"]
        if a in contained or b in contained:
            continue
        kind = classify(lens[a], lens[b], o["diag"])
        sb = 0 if o["same"] else 1
        if kind == "a_to_b":
            u, v = 2 * a, 2 * b + sb
        elif kind == "b_to_a":
            u, v = 2 * b + sb, 2 * a
        else:
            continue
        _add_edge(graph, u, v, o)
        _add_edge(graph, v ^ 1, u ^ 1, o)
    return graph


def count_edges(graph):
    return sum(len(outs) for outs in graph.values())


def transitive_reduction(graph, lens, fuzz=10):
    """Myers: heq u->w nëse ekziston u->v->w me të njëjtën zhvendosje (brenda fuzz).
    lens = gjatësitë e reads. Kthen numrin e skajeve të hequra."""
    to_remove = set()
    for u, outs in graph.items():
        if len(outs) < 2:
            continue
        lu = lens[u >> 1]
        offs = sorted((lu - e[0], v) for v, e in outs.items())
        for i in range(1, len(offs)):
            off_w, w = offs[i]
            for j in range(i):
                off_v, v = offs[j]
                if off_v >= off_w:
                    continue
                e = graph.get(v, {}).get(w)
                if e is None:
                    continue
                off_vw = lens[v >> 1] - e[0]
                if abs(off_v + off_vw - off_w) <= fuzz:
                    to_remove.add((u, w))
                    to_remove.add((w ^ 1, u ^ 1))
                    break
    removed = 0
    for u, w in to_remove:
        if graph[u].pop(w, None) is not None:
            removed += 1
    return removed


def graph_stats(graph, n_reads, contained):
    """Statistika të grafit (pas reduktimit)."""
    active = [r for r in range(n_reads) if r not in contained]
    out_deg = [len(graph.get(2 * r + s, {})) for r in active for s in (0, 1)]

    parent = {r: r for r in active}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for u, outs in graph.items():
        for v in outs:
            ru, rv = find(u >> 1), find(v >> 1)
            if ru != rv:
                parent[ru] = rv
    components = len({find(r) for r in active})

    return {
        "reads_in_graph": len(active),
        "nodes": 2 * len(active),
        "edges": count_edges(graph),
        "avg_out_degree": sum(out_deg) / len(out_deg) if out_deg else 0.0,
        "max_out_degree": max(out_deg) if out_deg else 0,
        "branching_nodes": sum(1 for d in out_deg if d > 1),
        "dead_end_nodes": sum(1 for d in out_deg if d == 0),
        "connected_components": components,
    }