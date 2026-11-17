from collections import defaultdict

from src.assemble import unitigs, remove_tips, pop_bubbles, path_to_sequence
from src.io_utils import reverse_complement
import random


def add(graph, u, v, ov, score):
    graph[u][v] = (ov, score)
    graph[v ^ 1][u ^ 1] = (ov, score)


rng = random.Random(3)
g = "".join(rng.choice("ACGT") for _ in range(300))
reads = [("r0", g[0:100], []), ("r1", g[60:160], []), ("r2", g[120:220], [])]

# 1. Unitig + sekuenca (strand përpara dhe mbrapsht)
graph = defaultdict(dict)
add(graph, 0, 2, 40, 100.0)
add(graph, 2, 4, 40, 100.0)
assert unitigs(graph, [0, 1, 2]) == [[0, 2, 4]], unitigs(graph, [0, 1, 2])
assert path_to_sequence(reads, graph, [0, 2, 4]) == g[0:220]
assert path_to_sequence(reads, graph, [5, 3, 1]) == reverse_complement(g[0:220])

# 2. Bubble: rruga me peshë më të madhe mbetet
graph = defaultdict(dict)
add(graph, 0, 2, 40, 100.0)
add(graph, 2, 6, 40, 100.0)
add(graph, 0, 4, 40, 50.0)
add(graph, 4, 6, 40, 50.0)
assert pop_bubbles(graph, 5) == 1
edges = {(u, v) for u, outs in graph.items() for v in outs}
assert edges == {(0, 2), (2, 6), (7, 3), (3, 1)}, edges

# 3. Tip: dega e shkurtër pa dalje hiqet
graph = defaultdict(dict)
add(graph, 0, 2, 40, 100.0)
add(graph, 2, 4, 40, 100.0)
add(graph, 4, 6, 40, 100.0)
add(graph, 2, 8, 40, 10.0)
assert remove_tips(graph, range(5), 1) == 1
assert 8 not in graph[2] and 4 in graph[2]

print("Të gjitha testet kaluan.")