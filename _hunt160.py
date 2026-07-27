"""Targeted hunt for a counterexample to conjecture 160 in the C4-free case.

The tension is explicit. Ls >= Delta always, and the right-hand side is
max l + max T <= Delta + max T, so a counterexample needs

    Ls - Delta  <  max_w T(w),

i.e. a graph whose maximum-leaf spanning tree is barely better than the star at
a maximum-degree vertex, while some *other* vertex carries many triangles. Those
pull against each other: triangles add structure, and structure raises Ls.
"""
import random
import networkx as nx
import wowii as W

ml = lambda g: max(W.local_independence(g, v) for v in g)
mT = W.max_triangles_at_vertex
Ls = W.max_leaves_via_cds
D = lambda g: max(d for _, d in g.degree())


def margin(g):
    """Positive means a counterexample."""
    if not nx.is_connected(g) or g.number_of_nodes() < 3:
        return -99
    if W.has_c4_subgraph(g):          # c_C4 = 0, the trivial case
        return -99
    return ml(g) + mT(g) - Ls(g)


def repair(g, rng):
    """Delete edges until C4-free."""
    guard = 0
    while W.has_c4_subgraph(g) and g.number_of_edges() and guard < 300:
        cyc = None
        for c in nx.simple_cycles(g, length_bound=4):
            if len(c) == 4:
                cyc = c
                break
        if cyc is None:
            break
        u, v = cyc[0], cyc[1]
        if g.has_edge(u, v):
            g.remove_edge(u, v)
        guard += 1
    return g


def hunt(n, iters=800, restarts=14, seed0=0):
    best = (-99, None)
    for s in range(seed0, seed0 + restarts):
        rng = random.Random(s * 3571 + n)
        g = nx.gnp_random_graph(n, rng.choice([0.18, 0.25, 0.32]),
                                seed=rng.randrange(10 ** 6))
        g = repair(g, rng)
        if not nx.is_connected(g):
            continue
        cur = margin(g)
        for _ in range(iters):
            u, v = rng.sample(range(n), 2)
            h = g.copy()
            if h.has_edge(u, v):
                h.remove_edge(u, v)
            else:
                h.add_edge(u, v)
                if W.has_c4_subgraph(h):
                    continue
            if not nx.is_connected(h):
                continue
            m = margin(h)
            if m >= cur:
                g, cur = h, m
            if cur > 0:
                return g, cur
        if cur > best[0]:
            best = (cur, g)
    return best[1], best[0]


if __name__ == "__main__":
    print("hunting conjecture 160, C4-free case\n")
    print("   a counterexample needs  Ls - Delta < max T\n")
    for n in range(9, 17):
        g, m = hunt(n)
        tag = "   *** COUNTEREXAMPLE ***" if m > 0 else ""
        print(f"   n={n:>2}  best margin {m:>3}{tag}", flush=True)
        if m > 0:
            print(f"      Ls={Ls(g)} maxl={ml(g)} maxT={mT(g)} Delta={D(g)}")
            print(f"      edges {sorted(g.edges())}")
            print(f"      graph6 "
                  f"{nx.to_graph6_bytes(g, header=False).strip().decode()}")
            break
