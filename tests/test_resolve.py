from collections import defaultdict

from src.assemble import resolve_branches


def add(graph, u, v, ov, score):
    graph[u][v] = (ov, score)
    graph[v ^ 1][u ^ 1] = (ov, score)


graph = defaultdict(dict)
add(graph, 0, 2, 40, 100.0)
add(graph, 0, 4, 40, 50.0)
resolve_branches(graph)
edges = {(u, v) for u, outs in graph.items() for v in outs}
assert edges == {(0, 2), (3, 1)}, edges

print("Të gjitha testet kaluan.")

def test_resolve_two_pass_keeps_valid_edge():
    # u->v (10), w->v (20), w->x (30): u->v duhet të mbetet
    graph = {
        0: {2: (0, 10)},              # u=0 -> v=2
        1: {2: (0, 20), 3: (0, 30)},  # w=1 -> v=2, w=1 -> x=3
        2: {},
        3: {},
    }
    resolve_branches(graph)
    assert 2 in graph[0], "u->v u hoq gabimisht"
    assert 3 in graph[1] and 2 not in graph[1]