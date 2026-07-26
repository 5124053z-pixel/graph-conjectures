"""Independent search for a counterexample to WOWII 291.

DeLaVina's page records one found by Zyad Tamimi on July 23, 2026: 12 vertices,
gamma_t = 4, freqMinTriangles = 1, Havel-Hakimi zero step k = 2. This searches
for one without using that graph, which checks the whole pipeline against a
conjecture now known to be false.
"""
import random
import networkx as nx
import wowii as W


def score(g):
    """How badly 291 fails: gamma_t - (k + freq). Positive means counterexample."""
    if not nx.is_connected(g) or g.number_of_nodes() <= 2:
        return -99
    try:
        return (W.total_domination_number(g)
                - W.havel_hakimi_zero_step(g) - W.freq_min_triangles(g))
    except Exception:
        return -99


def search(n, iters=4000, seed=0):
    rng = random.Random(seed)
    g = nx.gnp_random_graph(n, 0.35, seed=seed)
    while not nx.is_connected(g):
        g = nx.gnp_random_graph(n, 0.35, seed=rng.randrange(10 ** 6))
    best = score(g)
    for _ in range(iters):
        u, v = rng.sample(range(n), 2)
        h = g.copy()
        h.remove_edge(u, v) if h.has_edge(u, v) else h.add_edge(u, v)
        s = score(h)
        if s >= best:
            g, best = h, s
        if best > 0:
            return g, best
    return g, best


if __name__ == "__main__":
    print("independent search for a counterexample to WOWII 291\n")
    for n in (11, 12, 13):
        hit = None
        for seed in range(60):
            g, s = search(n, seed=seed)
            if s > 0:
                hit = (g, s, seed)
                break
        if hit:
            g, s, seed = hit
            print(f"n={n}: FOUND (seed {seed})")
            print(f"   gamma_t={W.total_domination_number(g)} "
                  f"k={W.havel_hakimi_zero_step(g)} "
                  f"freq={W.freq_min_triangles(g)}  slack={-s}")
            print(f"   m={g.number_of_edges()} "
                  f"degrees={sorted([d for _, d in g.degree()], reverse=True)}")
            print(f"   graph6 {nx.to_graph6_bytes(g, header=False).strip().decode()}")
            print(f"   conj_291 -> {W.conj_291(g)}")
            print(f"   edges {sorted(g.edges())}")
        else:
            print(f"n={n}: none found")
