"""Exhaustive check of a graph-invariant inequality over ALL trees on n vertices.

For small n the number of non-isomorphic trees is modest (317,955 at n=19), so
an inequality conjectured for all connected graphs can be settled *definitively*
on the tree case rather than searched heuristically. When the extremal examples
are trees -- as they are for the eigenvalue/matching conjecture below -- this is
both faster and more conclusive than a stochastic search.

Used first as a validation of the tooling: the answer for wagner-2.1 at n=19 is
known independently (Wagner reports a 19-vertex counterexample), so if this
script disagrees, the tooling is wrong rather than the literature.

Usage:
    python exhaustive_trees.py -n 19
"""
from __future__ import annotations

import argparse
import math
import time

import numpy as np
import networkx as nx


def tree_matching_number(g):
    """Maximum matching of a tree, greedily from the leaves. O(n).

    Standard: repeatedly match a leaf to its neighbour and delete both. This is
    optimal on trees (any maximum matching can be rearranged to include such an
    edge), and far faster than the general blossom algorithm.
    """
    adj = {v: set(nbrs) for v, nbrs in g.adjacency()}
    alive = set(adj)
    size = 0
    leaves = [v for v in alive if len(adj[v]) <= 1]
    while leaves:
        v = leaves.pop()
        if v not in alive:
            continue
        if not adj[v]:                      # isolated: cannot be matched
            alive.discard(v)
            continue
        u = next(iter(adj[v]))
        size += 1
        for w in (v, u):                    # delete both endpoints
            for x in adj[w]:
                if x in alive:
                    adj[x].discard(w)
            alive.discard(w)
            adj[w] = set()
        for w in list(alive):
            if len(adj[w]) <= 1:
                leaves.append(w)
    return size


def largest_eigenvalue(g, n):
    a = nx.to_numpy_array(g, nodelist=range(n))
    return float(np.linalg.eigvalsh(a)[-1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=19)
    ap.add_argument("--report-top", type=int, default=5,
                    help="how many smallest-slack trees to print")
    args = ap.parse_args()
    n = args.n

    threshold = math.sqrt(n - 1) + 1.0
    print(f"conjecture: connected G on n>=3 vertices has lambda_1 + mu >= sqrt(n-1)+1")
    print(f"n = {n}, threshold = {threshold:.9f}")
    print("enumerating all non-isomorphic trees ...")

    t0 = time.time()
    best = []          # (slack, edges)
    count = 0
    counterexamples = 0

    for g in nx.nonisomorphic_trees(n):
        count += 1
        lam = largest_eigenvalue(g, n)
        mu = tree_matching_number(g)
        slack = lam + mu - threshold
        if slack < 0:
            counterexamples += 1
        best.append((slack, lam, mu, tuple(sorted(g.edges()))))
        if len(best) > 4000:                      # keep memory bounded
            best.sort(key=lambda t: t[0])
            del best[args.report_top * 4:]

    best.sort(key=lambda t: t[0])
    elapsed = time.time() - t0

    print(f"trees checked      : {count:,}   ({elapsed:.1f}s)")
    print(f"counterexamples    : {counterexamples}")
    print()
    print(f"{args.report_top} smallest slacks:")
    for slack, lam, mu, edges in best[:args.report_top]:
        mark = "  <-- COUNTEREXAMPLE" if slack < 0 else ""
        print(f"  slack={slack:+.9f}  lambda_1={lam:.6f}  mu={mu}{mark}")
        print(f"      edges: {list(edges)}")


if __name__ == "__main__":
    main()
