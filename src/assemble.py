"""Unitigs, heqja e tips, bubble popping dhe ndërtimi i contigs nga grafi."""
from collections import defaultdict

from src.io_utils import reverse_complement


def _indeg(graph):
    ind = defaultdict(int)
    for outs in graph.values():
        for v in outs:
            ind[v] += 1
    return ind


def _preds(graph):
    preds = defaultdict(list)
    for u, outs in graph.items():
        for v in outs:
            preds[v].append(u)
    return preds


def _remove_edge(graph, u, v, ind=None):
    """Heq skajin u->v dhe binjakun e tij (v^1 -> u^1)."""
    if graph.get(u, {}).pop(v, None) is not None and ind is not None:
        ind[v] -= 1
    if graph.get(v ^ 1, {}).pop(u ^ 1, None) is not None and ind is not None:
        ind[u ^ 1] -= 1


def _walk(graph, ind, u, visited):
    """Ndjek nyjen u përpara për sa kohë rruga është pa degëzim."""
    path = [u]
    visited.add(u >> 1)
    while True:
        outs = graph.get(u, {})
        if len(outs) != 1:
            break
        v = next(iter(outs))
        if ind[v] != 1 or (v >> 1) in visited:
            break
        path.append(v)
        visited.add(v >> 1)
        u = v
    return path


def unitigs(graph, active):
    """Kthen listën e unitigs (secili është listë nyjesh). Çdo read del në një unitig të vetëm."""
    ind = _indeg(graph)
    preds = _preds(graph)
    nodes = sorted(2 * r + s for r in active for s in (0, 1))
    visited = set()
    result = []
    for u in nodes:
        if (u >> 1) in visited:
            continue
        p = preds.get(u, [])
        if len(p) == 1 and len(graph.get(p[0], {})) == 1:
            continue  # nuk është fillim zinxhiri
        result.append(_walk(graph, ind, u, visited))
    for u in nodes:  # cikle të pastra (të rralla)
        if (u >> 1) not in visited:
            result.append(_walk(graph, ind, u, visited))
    return result


def remove_tips(graph, active, max_tip):
    """Heq degët e shkurtra që mbarojnë pa dalje (ose fillojnë pa hyrje) dhe lidhen me një degëzim.
    Kthen numrin e tips të hequra."""
    ind = _indeg(graph)
    preds = _preds(graph)
    paths = [p for p in unitigs(graph, active) if len(p) <= max_tip]
    paths.sort(key=len)  # më të shkurtrat së pari
    removed = 0
    for path in paths:
        head, tail = path[0], path[-1]
        live_preds = [p for p in preds.get(head, []) if head in graph.get(p, {})]
        tip = False
        if not graph.get(tail):
            # fund pa dalje, i varur te një nyje me degëzim
            tip = any(len(graph[p]) > 1 for p in live_preds)
        elif not live_preds:
            # fillim pa hyrje, që bashkohet me një nyje me më shumë se një hyrje
            tip = any(ind[v] > 1 for v in graph[tail])
        if not tip:
            continue
        for n in path:
            for m in (n, n ^ 1):
                for v in list(graph.get(m, {})):
                    _remove_edge(graph, m, v, ind)
        removed += 1
    return removed


def pop_bubbles(graph, max_len):
    """Për çdo nyje me degëzim: rrugët paralele që bashkohen në të njëjtën nyje
    krahasohen dhe mbahet më e mira (më shumë reads, pastaj peshë më e madhe).
    Kthen numrin e rrugëve të hequra."""
    ind = _indeg(graph)
    popped = 0
    for u in list(graph.keys()):
        outs = graph.get(u, {})
        if len(outs) < 2:
            continue
        groups = defaultdict(list)
        for v in list(outs):
            interior, edges, cur = [], [(u, v)], v
            while (ind[cur] == 1 and len(graph.get(cur, {})) == 1
                   and len(interior) < max_len):
                interior.append(cur)
                nxt = next(iter(graph[cur]))
                edges.append((cur, nxt))
                cur = nxt
            score = sum(graph[a][b][1] for a, b in edges)
            groups[cur].append((len(interior), score, edges))
        for w, plist in groups.items():
            if len(plist) < 2:
                continue
            plist.sort(key=lambda t: (t[0], t[1]), reverse=True)
            for _, _, edges in plist[1:]:
                for a, b in edges:
                    _remove_edge(graph, a, b, ind)
                popped += 1
    return popped


def simplify(graph, active, max_tip=4, max_bubble=6, rounds=5):
    tips = bubbles = 0
    for _ in range(rounds):
        t = remove_tips(graph, active, max_tip)
        b = pop_bubbles(graph, max_bubble)
        tips += t
        bubbles += b
        if t == 0 and b == 0:
            break
    return {"tips_removed": tips, "bubbles_popped": bubbles}


def path_to_sequence(reads, graph, path):
    """Ndërton sekuencën e një rruge nyjesh: read-i i parë + pjesët pa overlap."""
    def oriented(u):
        s = reads[u >> 1][1]
        return reverse_complement(s) if u & 1 else s

    parts = [oriented(path[0])]
    for u, v in zip(path, path[1:]):
        ov = graph[u][v][0]
        parts.append(oriented(v)[ov:])
    return "".join(parts)


def build_contigs(reads, graph, active):
    """Kthen listë (rruga, sekuenca), nga më i gjati te më i shkurti."""
    contigs = [(p, path_to_sequence(reads, graph, p)) for p in unitigs(graph, active)]
    contigs.sort(key=lambda c: len(c[1]), reverse=True)
    return contigs


def contig_stats(contigs):
    lens = sorted((len(s) for _, s in contigs), reverse=True)
    total = sum(lens)
    n50, acc = 0, 0
    for l in lens:
        acc += l
        if acc * 2 >= total:
            n50 = l
            break
    return {
        "contigs": len(lens),
        "total_length": total,
        "longest": lens[0] if lens else 0,
        "n50": n50,
        "reads_in_contigs": sum(len(p) for p, _ in contigs),
    }