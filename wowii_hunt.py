"""Counterexample hunt above the exhaustive range, for every conjecture still open.

The exhaustive range here is 8 vertices. Conjecture 200 looked settled over all
12,112 connected graphs on at most 8 -- every applicable one held, and classical
Hamiltonicity theory covered all 671 of them -- and it is **false**, with an
11-vertex counterexample. Conjecture 291 is false at 12.

So "settled in range" is worth exactly as much as the range, and the range is
too small. This searches 9..14 vertices by hill-climbing on the violation
margin, which is how the 291 counterexample was reproduced here in seconds.

Each conjecture is turned into a margin `lhs - rhs` (or `rhs - lhs`) that is
positive exactly on a counterexample, and the search climbs it by toggling one
edge at a time, restarting from random graphs. The margin is only a heuristic --
it can be flat or misleading -- so a run finding nothing is weak evidence, while
a run finding something is a proof.
"""
from __future__ import annotations

import math
import random
import sys

import networkx as nx

import wowii as W


def _girth(g):
    return 0 if nx.is_forest(g) else nx.girth(g)


def _eccB(g):
    e = nx.eccentricity(g)
    hi = max(e.values())
    return W.ecc_set(g, [v for v in e if e[v] == hi], include_S=True)


def _maxl(g):
    return max(W.local_independence(g, v) for v in g)


def _margin_145(g):
    """Conjecture 145 needs lMin(complement) > 0, not merely a nonempty
    complement -- omitting that guard produced 202 false positives on the
    995-graph corpus, which is why every margin is checked against the corpus
    before any hunt is run."""
    gc = nx.complement(g)
    if gc.number_of_edges() == 0:
        return None
    lmin = min(W.local_independence(gc, v) for v in gc)
    if lmin <= 0:
        return None
    return 2 * _eccB(g) - W.largest_induced_tree(g) * lmin


# margin(g) > 0  <=>  g is a counterexample.  None means the hypothesis fails.
MARGINS = {
    2:   lambda g: 2 * (W.avg_l(g) - 1) - W.max_leaves_via_cds(g),
    19:  lambda g: math.floor(sum(nx.eccentricity(g).values()) / g.number_of_nodes()
                              + _maxl(g)) - W.largest_induced_bipartite(g),
    40:  lambda g: math.ceil((W.path_cover_number(g)
                              + W.largest_induced_bipartite(g) + 1) / 2)
                   - W.largest_induced_forest(g),
    59:  lambda g: math.ceil(math.sqrt(W.residue(g)
                                       * W.largest_induced_bipartite(g)))
                   - W.largest_induced_forest(g),
    61:  lambda g: W.residue(g) + math.ceil(nx.diameter(g) / 3)
                   - W.largest_induced_forest(g),
    100: lambda g: (W.indep_number(g)
                    - math.ceil((_maxl(g)
                                 + 0.5 * W.degree_l2_norm(nx.complement(g))) / 2)
                    if nx.is_connected(nx.complement(g)) else None),
    133: lambda g: nx.radius(g) + math.floor(W.avg_l(g))
                   ** (0 if W.has_c4_subgraph(g) else 1) - W.largest_induced_path(g),
    141: lambda g: _girth(g) / 2 - 1 + _maxl(g) - W.largest_induced_tree(g),
    142: lambda g: (2 / 3) * _girth(g) + _eccB(g) - W.largest_induced_tree(g),
    144: lambda g: _girth(g) - 1 + W.ecc_set(g, nx.center(g))
                   - W.largest_induced_tree(g),
    145: lambda g: _margin_145(g),
}

HAM = {
    194: lambda g: W.indep_number(g) <= 1 + W.avg_l(g),
    "198a": lambda g: (W.largest_induced_bipartite(g)
                       <= 2 + sum(nx.eccentricity(g).values()) / g.number_of_nodes()),
    217: lambda g: W.max_leaves_via_cds(g) <= 4 * (1 if W.residue(g) == 2 else 0) + 2,
}


def margin(num, g):
    """Positive means counterexample; None means the hypothesis does not apply."""
    if num in HAM:
        try:
            if not HAM[num](g):
                return None
        except Exception:
            return None
        return 1 if not W.has_hamiltonian_path(g) else -1
    try:
        return MARGINS[num](g)
    except Exception:
        return None


def hunt(num, n, iters=1500, restarts=25, seed0=0, verbose=False):
    best_seen = None
    for seed in range(seed0, seed0 + restarts):
        rng = random.Random(seed * 7919 + n)
        p = rng.choice([0.2, 0.3, 0.4, 0.5, 0.6])
        g = nx.gnp_random_graph(n, p, seed=rng.randrange(10 ** 6))
        if not nx.is_connected(g):
            continue
        cur = margin(num, g)
        cur = -99 if cur is None else cur
        for _ in range(iters):
            u, v = rng.sample(range(n), 2)
            h = g.copy()
            h.remove_edge(u, v) if h.has_edge(u, v) else h.add_edge(u, v)
            if not nx.is_connected(h):
                continue
            m = margin(num, h)
            m = -99 if m is None else m
            if m >= cur:
                g, cur = h, m
            if cur > 0:
                return g, cur
        if best_seen is None or cur > best_seen[1]:
            best_seen = (g, cur)
    return best_seen if best_seen else (None, -99)


def main():
    nums = list(MARGINS) + list(HAM)
    lo, hi = 9, 13
    if len(sys.argv) > 1:
        nums = [int(x) if x.isdigit() else x for x in sys.argv[1].split(",")]
    if len(sys.argv) > 2:
        lo, hi = (int(x) for x in sys.argv[2].split("-"))
    print(f"counterexample hunt, {lo}..{hi} vertices\n")
    print(f"   {'conj':>6} {'n':>3} {'best margin':>12}   verdict")
    for num in nums:
        for n in range(lo, hi + 1):
            g, m = hunt(num, n)
            if m > 0:
                print(f"   {str(num):>6} {n:>3} {m:>12}   *** COUNTEREXAMPLE ***",
                      flush=True)
                print(f"          graph6 "
                      f"{nx.to_graph6_bytes(g, header=False).strip().decode()}")
                print(f"          edges {sorted(g.edges())}")
                break
            print(f"   {str(num):>6} {n:>3} {m:>12}", flush=True)


if __name__ == "__main__":
    main()
