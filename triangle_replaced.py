"""The frontier family of the cubic zero-forcing conjecture, and a scan of it.

Conjecture: connected cubic G != K_4  =>  Z(G) <= alpha(G) + 1.

Reading the exhaustive scan rather than just its verdict: every graph attaining
equality has the *smallest possible* independence number for its order. For
connected cubic graphs other than K_4 that minimum is alpha = n/3, and equality
forces the vertex set to split into n/3 triangles. A cubic graph whose vertices
partition into triangles has, in each triangle, two internal neighbours and one
external, so the triangles are joined by a perfect matching -- i.e. G is
obtained from a cubic multigraph H on k = n/3 vertices by replacing every vertex
with a triangle.

    G = truncate(H),   alpha(G) = k,   n = 3k.

So the conjecture, restricted to the only place equality ever happens, becomes

    Z(truncate(H))  <=  k + 1     for every cubic multigraph H on k vertices,

and the search space collapses from all cubic graphs on n vertices to all cubic
multigraphs on n/3 of them. That is a cube-root of the problem: the exhaustive
scan reaches n = 14 in minutes, whereas H on 14 vertices covers G on 42.

This does not prove anything about graphs off the frontier, and a counterexample
could in principle have alpha > n/3. But equality never happens there in the
exhaustive range, so this is where to look first.
"""
from __future__ import annotations

import argparse
import time

import networkx as nx

from cubic_graphs import generate_multigraphs
from zero_forcing import independence_number, has_forcing_set_of_size, \
    zero_forcing_number


def truncate(edges, k):
    """Replace every vertex of the cubic multigraph H by a triangle.

    Vertex v of H becomes the triangle (v,0), (v,1), (v,2); the three edge-ends
    at v are attached to its three corners, one each."""
    slots = {v: [] for v in range(k)}
    for i, (u, v) in enumerate(edges):
        slots[u].append(i)
        slots[v].append(i)          # a loop lands in the same list twice
    corner = {}
    g = nx.Graph()
    for v in range(k):
        for c in range(3):
            g.add_node((v, c))
        g.add_edges_from([((v, 0), (v, 1)), ((v, 1), (v, 2)), ((v, 2), (v, 0))])
        for c, e in enumerate(slots[v]):
            corner.setdefault(e, []).append((v, c))
    for e, ends in corner.items():
        if len(ends) == 2:
            g.add_edge(ends[0], ends[1])
    return nx.convert_node_labels_to_integers(g)


def scan(kmax, verbose=True):
    results = []
    for k, level in generate_multigraphs(kmax, verbose=False).items():
        t0 = time.time()
        n = 3 * k
        worst, tight, bad = 0, 0, []
        for edges in level:
            g = truncate(edges, k)
            if g.number_of_nodes() != n or not nx.is_connected(g):
                continue
            if set(dict(g.degree()).values()) != {3}:
                continue
            a = independence_number(g)
            ok, _ = has_forcing_set_of_size(g, a + 1)
            if not ok:
                z = zero_forcing_number(g, lo=a + 2)
                bad.append((n, a, z, sorted(g.edges())))
            elif not has_forcing_set_of_size(g, a, tries=1200)[0]:
                tight += 1
            worst = max(worst, a)
        results.append((k, n, len(level), tight, bad))
        if verbose:
            print(f"  H on k = {k:>2} vertices ({len(level):>7,} multigraphs) "
                  f"-> G on n = {n:>3}   Z = alpha+1 exactly: {tight:>6,}   "
                  f"counterexamples: {len(bad)}   [{time.time() - t0:.1f}s]")
        for b in bad:
            print(f"\n*** COUNTEREXAMPLE  n={b[0]}  alpha={b[1]}  Z={b[2]}")
            print(f"    edges: {b[3]}\n")
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-k", type=int, default=10,
                    help="largest H; the graphs tested have n = 3k vertices")
    args = ap.parse_args()
    print("conjecture on its frontier: Z(truncate(H)) <= |V(H)| + 1")
    print("H ranges over all connected cubic multigraphs\n")
    scan(args.k)


if __name__ == "__main__":
    main()
